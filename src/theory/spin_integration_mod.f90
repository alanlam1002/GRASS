module spin_integration_mod
  use precision_mod, only: wp
  use iso_fortran_env, only: int32
  use para_mod, only: wp, SDIV, MDIV, LMAX, s_gp, mu, sin_theta, s_pwr, &
                      rho, gama, alpha, ww, omg, sphi, &
                      energy, pressure, enthalpy, velocity_sq, &
                      P_2n, P1_2n_1, l_uni, KAPPA, C, G, MSUN, MB, pi, KSCALE, &
                      r_ratio, r_e, enthalpy_min, s_e, SMAX, &
                      Omega_c, Omega_e, sphi_m, &
                      B_coup, mphi_r, mass, mass_0, ang_mom, &
                      A_diff, lambda1, lambda2, solver_type, output, timing, eos_file, inc_p, inc_lambda
  use cheb_mod, only: cheb_fit_stats, cheb_std_base, cheb_get_deriv_point
  use eos_mod, only: n0_at_e
  use exporter_mod, only: initial_data_for_sacra_aei
  use toolkit_mod, only: bessel_even_tables
  use spin_derivatives_mod, only: deriv_s_sub, deriv_m_sub
  use spin_workspace_mod
  use lapack_interfaces_mod, only: dgemm, dpbsv
  implicit none
  private
  public :: get_all_targets, update_alpha_potential, output_helper, write_restart_file

contains

  pure logical function is_massive_scalar()
    is_massive_scalar = mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp
  end function

  pure logical function is_spherical()
    is_spherical = abs(r_ratio - 1.0_wp) < epsilon(r_ratio)
  end function

  elemental real(wp) function pade_exp(x) result(y)
    real(wp), intent(in) :: x
    real(wp) :: t, t2
    t = x * 0.00390625_wp
    t2 = t * t
    y = (12.0_wp + 6.0_wp * t + t2) / (12.0_wp - 6.0_wp * t + t2)
    y = y*y; y = y*y; y = y*y; y = y*y
    y = y*y; y = y*y; y = y*y; y = y*y
  end function

  subroutine precompute(r_e_new, root_mphi_re)
    real(wp), intent(in) :: r_e_new, root_mphi_re
    integer :: s, m
    real(wp) :: r2, t0, t1
    real(wp) :: dt_bessel, dt_ds, dt_dm, dt_d2s, dt_d2m, dt_exp
    integer, parameter :: NT = 5
    integer, save :: nc = 0
    real(wp), save :: st(6) = 0.0_wp

    if (timing) then
      nc = nc + 1
      if (nc == 1) st = 0.0_wp
      call cpu_time(t0)
    end if

    if (is_massive_scalar()) then
      mr_cache(:) = root_mphi_re * s_gp(:) / (1.0_wp - s_gp(:))
      do s = 1, SDIV
        call bessel_even_tables(mr_cache(s), LMAX, besseli_cache(:,s), besselk_cache(:,s))
      end do
    end if
    if (timing) then; call cpu_time(t1); dt_bessel = t1 - t0; call cpu_time(t0); end if

    call deriv_s_sub(gama, dg_s_cache)
    call deriv_s_sub(rho,  dr_s_cache)
    call deriv_s_sub(ww,   dww_s_cache)
    call deriv_s_sub(sphi, ds_s_cache)
    if (timing) then; call cpu_time(t1); dt_ds = t1 - t0; call cpu_time(t0); end if

    call deriv_m_sub(gama, dg_m_cache)
    call deriv_m_sub(rho,  dr_m_cache)
    call deriv_m_sub(ww,   dww_m_cache)
    call deriv_m_sub(sphi, ds_m_cache)
    if (timing) then; call cpu_time(t1); dt_dm = t1 - t0; call cpu_time(t0); end if

    call deriv_s_sub(dg_s_cache, d2g_ss_cache)
    do m = 1, MDIV
      d2g_ss_cache(:,m) = s1_geom * d2g_ss_cache(:,m) + (1.0_wp - 2.0_wp * s_gp) * dg_s_cache(:,m)
    end do
    if (timing) then; call cpu_time(t1); dt_d2s = t1 - t0; call cpu_time(t0); end if

    call deriv_m_sub(dg_m_cache, d2g_mm_cache)
    do m = 1, MDIV
      d2g_mm_cache(:,m) = m1_geom(m) * d2g_mm_cache(:,m) - 2.0_wp * mu(m) * dg_m_cache(:,m)
    end do
    if (timing) then; call cpu_time(t1); dt_d2m = t1 - t0; call cpu_time(t0); end if

    r2 = r_e_new * r_e_new
    do m = 1, MDIV
      do s = 1, SDIV
        e_gsm_cache(s,m)      = pade_exp(0.5_wp * gama(s,m))
        e_rsm_cache(s,m)      = pade_exp(-rho(s,m))
        e2alpha_r2_cache(s,m) = pade_exp(2.0_wp * alpha(s,m)) * r2
        Acoup4_cache(s,m)     = pade_exp(-sphi(s,m)**2 * B_coup)
      end do
    end do
    if (timing) then; call cpu_time(t1); dt_exp = t1 - t0; end if

    if (timing) then
      st = st + [dt_bessel, dt_ds, dt_dm, dt_d2s, dt_d2m, dt_exp]
      if (nc >= NT) then
        write(*,'(A,I0,A)') 'precompute avg over ', NT, ':'
        write(*,'(A,1X,ES12.5)') '  bessel',       st(1)/NT
        write(*,'(A,1X,ES12.5)') '  first_deriv_s', st(2)/NT
        write(*,'(A,1X,ES12.5)') '  first_deriv_m', st(3)/NT
        write(*,'(A,1X,ES12.5)') '  second_deriv_s',st(4)/NT
        write(*,'(A,1X,ES12.5)') '  second_deriv_m',st(5)/NT
        write(*,'(A,1X,ES12.5)') '  exp_tables',    st(6)/NT
        write(*,'(A,1X,ES12.5)') '  total', sum(st)/NT
      end if
    end if
  end subroutine precompute

  subroutine build_source_terms(r_e_new)
    real(wp), intent(in) :: r_e_new
    integer :: m
    real(wp) :: mum, m1
    real(wp), dimension(SDIV) :: esm, psm, vphi, scal_p
    real(wp), dimension(SDIV) :: vel_fac, matter_sum, matter_trace
    real(wp), dimension(SDIV) :: e2alpha_s2, src_common, e_rsm2
    real(wp), dimension(SDIV) :: dg_s_sc, dg_m_sc, rho_bracket
    real(wp), dimension(SDIV) :: omega_matter, omega_bracket, sphi_src
    real(wp), dimension(SDIV) :: vsq, one_plus_vsq

    do m = 1, MDIV
      mum = mu(m)
      m1 = m1_geom(m)

      associate( gs => dg_s_cache(:,m), gm => dg_m_cache(:,m), &
                 rs => dr_s_cache(:,m), rm => dr_m_cache(:,m), &
                 wws => dww_s_cache(:,m), wwm => dww_m_cache(:,m), &
                 ss => ds_s_cache(:,m), sm => ds_m_cache(:,m), &
                 gss => d2g_ss_cache(:,m), gmm => d2g_mm_cache(:,m), &
                 egsm => e_gsm_cache(:,m), ersm => e_rsm_cache(:,m), &
                 e2ar2 => e2alpha_r2_cache(:,m), Ac4 => Acoup4_cache(:,m), &
                 ww_col => ww(:,m), omg_col => omg(:,m), &
                 rho_col => rho(:,m), gama_col => gama(:,m), sphi_col => sphi(:,m) )

      esm = energy(:,m) * Ac4
      psm = pressure(:,m) * Ac4
      vphi = sphi_col**2 * mphi_r * 0.5_wp * e2ar2
      scal_p = sphi_col * egsm
      vsq = velocity_sq(:,m)
      one_plus_vsq = 1.0_wp + vsq
      vel_fac = 1.0_wp / (1.0_wp - vsq)
      matter_sum = esm + psm
      matter_trace = esm - 3.0_wp * psm
      e_rsm2 = ersm**2
      e2alpha_s2 = e2ar2 * s2_geom
      src_common = 16.0_wp * pi * e2alpha_s2 * psm - 4.0_wp * vphi * s2_geom
      dg_s_sc = s1_geom * gs
      dg_m_sc = m1 * gm

      rho_bracket = src_common &
        - dg_s_sc * (0.5_wp * dg_s_sc + 1.0_wp) &
        - gm * (0.5_wp * dg_m_sc - mum)

      S_metric_rho(:,m) = egsm * ( &
          8.0_wp * pi * e2alpha_s2 * matter_sum * one_plus_vsq * vel_fac &
        + s2_geom * m1 * e_rsm2 * (s1_sq_geom * wws**2 + m1 * wwm**2) &
        + dg_s_sc - mum * gm &
        + rho_col * 0.5_wp * rho_bracket )

      S_metric_gama(:,m) = egsm * (src_common &
        + gama_col * 0.5_wp * (src_common - 0.5_wp * dg_s_sc**2 &
        - 0.5_wp * dg_m_sc * gm) )

      if (is_spherical()) then
        S_metric_omega(:,m) = 0.0_wp
      else
        omega_matter = (one_plus_vsq * esm + 2.0_wp * vsq * psm) * vel_fac
        omega_bracket = -8.0_wp * pi * e2alpha_s2 * omega_matter &
          - s1_geom * (2.0_wp * rs + 0.5_wp * gs) &
          + mum * (2.0_wp * rm + 0.5_wp * gm) &
          + 0.25_wp * s1_sq_geom * (4.0_wp * rs**2 - gs**2) &
          + 0.25_wp * m1 * (4.0_wp * rm**2 - gm**2) &
          - m1 * e_rsm2 * (sgp4_geom * wws**2 + s2_geom * m1 * wwm**2) &
          - 2.0_wp * vphi * s2_geom

        S_metric_omega(:,m) = egsm * ersm * ( &
          -16.0_wp * pi * e2alpha_s2 * (omg_col - ww_col) * matter_sum * vel_fac &
          + ww_col * omega_bracket )
      endif

      if (is_massive_scalar()) then
        sphi_src = -2.0_wp * pi * B_coup * matter_trace + mphi_r
        S_metric_sphi(:,m) = -r_e_new**2 * s2_geom * scal_p * mphi_r &
          + e2alpha_s2 * scal_p * sphi_src &
          + scal_p * (s1_one_minus_s_geom * gs + s1_sq_geom * (0.5_wp * gss + 0.25_wp * gs**2) &
          + m1 * (0.5_wp * gmm + 0.25_wp * gm**2) - mum * gm)
      else
        sphi_src = -2.0_wp * pi * B_coup * matter_trace
        S_metric_sphi(:,m) = -s1_sq_geom * gs * ss - dg_m_sc * sm &
          + sphi_col * sphi_src * e2alpha_s2
      endif
      end associate
    end do
  end subroutine build_source_terms

  subroutine project_integrate(root_mphi_re)
    real(wp), intent(in) :: root_mphi_re
    integer, parameter :: TILE = 8
    integer :: ncols, c1, c2, ic, n
    real(wp) :: f2n(SDIV), lp(SDIV), rs(SDIV), mscale(SDIV)

    ncols = min(TILE, max(1, LMAX+1))

    D2_metric_gama(:,1) = 0.0_wp
    D2_metric_omega(:,1) = 0.0_wp
    if (is_massive_scalar()) mscale = radial_quad_weights * wfac_cache * root_mphi_re

    f2n = 1.0_wp
    do c1 = 1, LMAX+1, ncols
      c2 = min(LMAX+1, c1 + ncols - 1)
      call dgemm('N', 'N', SDIV, c2-c1+1, MDIV, 1.0_wp, S_metric_rho, SDIV, &
                 weighted_even_basis(:,c1:c2), MDIV, 0.0_wp, proj_work(:,1:c2-c1+1), SDIV)
      do ic = 1, c2 - c1 + 1
        n = c1 + ic - 2
        call green_rho(n, f2n, proj_work(:,ic), radial_quad_weights, D2_metric_rho(:,c1+ic-1), lp, rs)
        f2n(2:) = f2n(2:) * rad_ratio_g(2:)
      end do
    end do

    if (is_massive_scalar()) then
      do c1 = 1, LMAX+1, ncols
        c2 = min(LMAX+1, c1 + ncols - 1)
        call dgemm('N', 'N', SDIV, c2-c1+1, MDIV, 1.0_wp, S_metric_sphi, SDIV, &
                   weighted_even_basis(:,c1:c2), MDIV, 0.0_wp, proj_work(:,1:c2-c1+1), SDIV)
        do ic = 1, c2 - c1 + 1
          call green_bessel(c1+ic-1, proj_work(:,ic), mscale, D2_metric_sphi(:,c1+ic-1), lp, rs)
        end do
      end do
    else
      f2n = 1.0_wp
      do c1 = 1, LMAX+1, ncols
        c2 = min(LMAX+1, c1 + ncols - 1)
        call dgemm('N', 'N', SDIV, c2-c1+1, MDIV, 1.0_wp, S_metric_sphi, SDIV, &
                   weighted_even_basis(:,c1:c2), MDIV, 0.0_wp, proj_work(:,1:c2-c1+1), SDIV)
        do ic = 1, c2 - c1 + 1
          n = c1 + ic - 2
          call green_rho(n, f2n, proj_work(:,ic), radial_quad_weights, D2_metric_sphi(:,c1+ic-1), lp, rs)
          f2n(2:) = f2n(2:) * rad_ratio_g(2:)
        end do
      end do
    end if

    if (LMAX > 0) then
      f2n(1) = 1.0_wp
      f2n(2:) = rad_ratio_g(2:)
      do c1 = 1, LMAX, ncols
        c2 = min(LMAX, c1 + ncols - 1)
        call dgemm('N', 'N', SDIV, c2-c1+1, MDIV, 1.0_wp, S_metric_gama, SDIV, &
                   weighted_gama_basis(:,c1:c2), MDIV, 0.0_wp, proj_work(:,1:c2-c1+1), SDIV)
        do ic = 1, c2 - c1 + 1
          n = c1 + ic - 1
          call green_multipole(n, f2n, proj_work(:,ic), radial_quad_weights, &
                               rad_inv_s1, .false., D2_metric_gama(:,n+1), lp, rs)
          f2n(2:) = f2n(2:) * rad_ratio_g(2:)
        end do
      end do

      f2n(1) = 1.0_wp
      f2n(2:) = rad_ratio_g(2:)
      do c1 = 1, LMAX, ncols
        c2 = min(LMAX, c1 + ncols - 1)
        call dgemm('N', 'N', SDIV, c2-c1+1, MDIV, 1.0_wp, S_metric_omega, SDIV, &
                   weighted_omega_basis(:,c1:c2), MDIV, 0.0_wp, proj_work(:,1:c2-c1+1), SDIV)
        do ic = 1, c2 - c1 + 1
          n = c1 + ic - 1
          call green_multipole(n, f2n, proj_work(:,ic), radial_quad_weights, &
                               rad_left_rho, .true., D2_metric_omega(:,n+1), lp, rs)
          f2n(2:) = f2n(2:) * rad_ratio_g(2:)
        end do
      end do
    end if

  contains

    subroutine green_rho(n_phys, f2n, src, wgt, out, lp, rs)
      integer, intent(in) :: n_phys
      real(wp), intent(in) :: f2n(SDIV), src(SDIV), wgt(SDIV)
      real(wp), intent(out) :: out(SDIV)
      real(wp), intent(inout) :: lp(SDIV), rs(SDIV)
      integer :: k
      real(wp) :: wv

      lp(1) = 0.0_wp
      if (n_phys == 0) then
        do k = 2, SDIV
          wv = wgt(k) * src(k)
          lp(k) = lp(k-1) + rad_left_rho(k) * wv
        end do
        rs(SDIV) = rad_inv_s1(SDIV) * wgt(SDIV) * src(SDIV)
        do k = SDIV-1, 2, -1
          rs(k) = rs(k+1) + rad_inv_s1(k) * wgt(k) * src(k)
        end do
        rs(1) = rs(2)
        out(1) = rs(1)
        do k = 2, SDIV
          out(k) = rad_ratio_s(k) * lp(k-1) + rs(k)
        end do
      else
        do k = 2, SDIV
          wv = wgt(k) * src(k)
          lp(k) = lp(k-1) + rad_left_rho(k) * wv / f2n(k)
        end do
        rs(SDIV) = rad_inv_s1(SDIV) * f2n(SDIV) * wgt(SDIV) * src(SDIV)
        do k = SDIV-1, 2, -1
          rs(k) = rs(k+1) + rad_inv_s1(k) * f2n(k) * wgt(k) * src(k)
        end do
        rs(1) = rs(2)
        out(1) = 0.0_wp
        do k = 2, SDIV
          out(k) = f2n(k) * rad_ratio_s(k) * lp(k-1) + rs(k) / f2n(k)
        end do
      end if
    end subroutine green_rho

    subroutine green_multipole(n_phys, f2n, src, wgt, lp_coeff, use_ratio_s, out, lp, rs)
      integer, intent(in) :: n_phys
      real(wp), intent(in) :: f2n(SDIV), src(SDIV), wgt(SDIV), lp_coeff(SDIV)
      logical, intent(in) :: use_ratio_s
      real(wp), intent(out) :: out(SDIV)
      real(wp), intent(inout) :: lp(SDIV), rs(SDIV)
      integer :: k
      real(wp) :: boundary_sum, wv

      lp(1) = 0.0_wp
      if (n_phys == 1) then
        boundary_sum = 0.0_wp
        do k = 2, SDIV
          wv = wgt(k) * src(k)
          lp(k) = lp(k-1) + lp_coeff(k) * wv / f2n(k)
          boundary_sum = boundary_sum + rad_inv_s1(k) * wv
        end do
      else
        do k = 2, SDIV
          wv = wgt(k) * src(k)
          lp(k) = lp(k-1) + lp_coeff(k) * wv / f2n(k)
        end do
      end if
      rs(SDIV) = rad_right_gama(SDIV) * f2n(SDIV) * wgt(SDIV) * src(SDIV)
      do k = SDIV-1, 2, -1
        rs(k) = rs(k+1) + rad_right_gama(k) * f2n(k) * wgt(k) * src(k)
      end do
      rs(1) = rs(2)

      out(1) = merge(boundary_sum, 0.0_wp, n_phys == 1)
      do k = 2, SDIV
        if (use_ratio_s) then
          out(k) = f2n(k) * rad_ratio_s(k) * lp(k-1) + rad_ratio_g(k) * rs(k) / f2n(k)
        else
          out(k) = real(s_pwr, wp) * f2n(k) * lp(k-1) + rad_ratio_g(k) * rs(k) / f2n(k)
        end if
      end do
    end subroutine green_multipole

    subroutine green_bessel(n_idx, src, wgt, out, lp, rs)
      integer, intent(in) :: n_idx
      real(wp), intent(in) :: src(SDIV), wgt(SDIV)
      real(wp), intent(out) :: out(SDIV)
      real(wp), intent(inout) :: lp(SDIV), rs(SDIV)
      integer :: s
      real(wp) :: wv

      lp(1) = wgt(1) * src(1) * besseli_cache(n_idx,1)
      do s = 2, SDIV
        wv = wgt(s) * src(s)
        lp(s) = lp(s-1) + wv * besseli_cache(n_idx,s)
      end do

      rs(SDIV) = wgt(SDIV) * src(SDIV) * besselk_cache(n_idx,SDIV)
      do s = SDIV-1, 1, -1
        wv = wgt(s) * src(s)
        rs(s) = rs(s+1) + wv * besselk_cache(n_idx,s)
      end do

      out(1) = besseli_cache(n_idx,1) * rs(1)
      do s = 2, SDIV
        out(s) = besselk_cache(n_idx,s) * lp(s-1) + besseli_cache(n_idx,s) * rs(s)
      end do
    end subroutine green_bessel

  end subroutine project_integrate

  subroutine reconstruct()
    call dgemm('N','T', SDIV, MDIV, LMAX+1, 1.0_wp, D2_metric_rho, SDIV, P_2n, MDIV, &
               0.0_wp, target_rho, SDIV)

    if (is_massive_scalar()) then
      call dgemm('N','T', SDIV, MDIV, LMAX+1, -1.0_wp, D2_metric_sphi, SDIV, &
                 recon_even_massive_basis, MDIV, 0.0_wp, target_sphi, SDIV)
    else
      call dgemm('N','T', SDIV, MDIV, LMAX+1, -1.0_wp, D2_metric_sphi, SDIV, P_2n, MDIV, &
                 0.0_wp, target_sphi, SDIV)
    endif

    if (is_spherical()) then
      call dgemm('N','T', SDIV, MDIV, 1, -(2.0_wp/pi), D2_metric_gama(:,2:2), SDIV, &
                 recon_gama_basis(:,1:1), MDIV, 0.0_wp, target_gama, SDIV)
      target_rho = -target_rho / e_gsm_cache
      target_gama = target_gama / e_gsm_cache
      if (is_massive_scalar()) then
        target_sphi = target_sphi / e_gsm_cache
      end if
      target_ww = 0.0_wp
      return
    endif

    call dgemm('N','T', SDIV, MDIV, LMAX, 1.0_wp, D2_metric_omega(:,2:LMAX+1), SDIV, &
               recon_omega_basis, MDIV, 0.0_wp, target_ww, SDIV)
    call dgemm('N','T', SDIV, MDIV, LMAX, -(2.0_wp/pi), D2_metric_gama(:,2:LMAX+1), SDIV, &
               recon_gama_basis, MDIV, 0.0_wp, target_gama, SDIV)

    target_rho = -target_rho / e_gsm_cache
    target_gama = target_gama / e_gsm_cache
    target_ww = target_ww / (e_gsm_cache * e_rsm_cache)
    if (is_massive_scalar()) then
      target_sphi = target_sphi / e_gsm_cache
    end if
  end subroutine reconstruct

  subroutine get_all_targets(r_e_new, root_mphi_re, &
                            out_target_rho, out_target_gama, out_target_ww, out_target_sphi)
    real(wp), intent(in) :: r_e_new, root_mphi_re
    real(wp), intent(out) :: out_target_rho(SDIV,MDIV), out_target_gama(SDIV,MDIV)
    real(wp), intent(out) :: out_target_ww(SDIV,MDIV), out_target_sphi(SDIV,MDIV)
    real(wp) :: t0, t1, dt(4)
    integer, parameter :: NT = 5
    integer, save :: nc = 0
    real(wp), save :: st(4) = 0.0_wp

    if (timing) then
      nc = nc + 1
      if (nc == 1) st = 0.0_wp
      dt = 0.0_wp
      call cpu_time(t0)
    end if

    call precompute(r_e_new, root_mphi_re)
    if (timing) then; call cpu_time(t1); dt(1) = t1 - t0; call cpu_time(t0); end if

    call build_source_terms(r_e_new)
    if (timing) then; call cpu_time(t1); dt(2) = t1 - t0; call cpu_time(t0); end if

    call project_integrate(root_mphi_re)
    if (timing) then; call cpu_time(t1); dt(3) = t1 - t0; call cpu_time(t0); end if

    call reconstruct()
    if (timing) then; call cpu_time(t1); dt(4) = t1 - t0; end if

    out_target_rho  = target_rho
    out_target_gama = target_gama
    out_target_ww   = target_ww
    out_target_sphi = target_sphi

    if (timing) then
      st = st + dt
      if (nc >= NT) then
        write(*,'(A,I0,A)') 'get_all_targets avg over ', NT, ':'
        write(*,'(A,1X,ES12.5)') '  precompute',          st(1)/NT
        write(*,'(A,1X,ES12.5)') '  build_source_terms',  st(2)/NT
        write(*,'(A,1X,ES12.5)') '  project_integrate',   st(3)/NT
        write(*,'(A,1X,ES12.5)') '  reconstruct',         st(4)/NT
        write(*,'(A,1X,ES12.5)') '  total',               sum(st)/NT
      end if
    end if
  end subroutine get_all_targets

  subroutine update_alpha_potential(r_e_new)
    real(wp), intent(in) :: r_e_new
    integer :: m, s
    real(wp) :: m1, mu_m
    real(wp) :: t0, t1, dt(4)
    integer, save :: nc = 0
    real(wp), save :: st(4) = 0.0_wp
    real(wp), dimension(SDIV,MDIV) :: da_dm, d_gama_sm
    real(wp) :: adj_const(SDIV)
    real(wp) :: gs, gm, rs, rm, ss, sm, gsm, wws, wwm, gss, gmm, e_cache
    real(wp) :: sg, s1, sg_ratio, sg4, numer_m, one_plus_s1dgs, inv_denom
    real(wp) :: temp1, temp3, temp4, temp5, temp6, temp7, temp8, temp9

    if (timing) then
      nc = nc + 1
      if (nc == 1) st = 0.0_wp
      dt = 0.0_wp
      call cpu_time(t0)
    end if

    alpha(:,:) = 0.0_wp
    if (is_spherical()) return

    da_dm(1,:) = 0.0_wp
    call deriv_m_sub(dg_s_cache, d_gama_sm)
    if (timing) then; call cpu_time(t1); dt(1) = t1 - t0; call cpu_time(t0); end if

    do m = 1, MDIV
      mu_m = mu(m)
      m1   = m1_geom(m)
      da_dm(1,m) = 0.0_wp
      do s = 2, SDIV
        gs = dg_s_cache(s,m)
        gm = dg_m_cache(s,m)
        rs = dr_s_cache(s,m)
        rm = dr_m_cache(s,m)
        ss = ds_s_cache(s,m)
        sm = ds_m_cache(s,m)
        gsm = d_gama_sm(s,m)
        wws = dww_s_cache(s,m)
        wwm = dww_m_cache(s,m)
        gss = d2g_ss_cache(s,m)
        gmm = d2g_mm_cache(s,m)
        e_cache = e_rsm_cache(s,m)
        sg = s_gp(s)
        s1 = s1_geom(s)
        sg_ratio = sgp_term_2d_cache(s,1)
        sg4 = sgp4_geom(s)

        numer_m = -mu_m + m1 * gm
        one_plus_s1dgs = 1.0_wp + s1 * gs
        inv_denom = 1.0_wp / (m1 * one_plus_s1dgs**2 + numer_m**2)

        temp1 = 2.0_wp * sg**2 * sg_ratio * m1 * wws * wwm * one_plus_s1dgs &
          - (sg4 * wws**2 - (sg * wwm * sg_ratio)**2 * m1) * numer_m
        temp3 = s1 * gss + (s1 * gs)**2
        temp4 = gm * numer_m
        temp5 = ((s1 * (rs + gs))**2 - m1 * (rm + gm)**2) * numer_m
        temp6 = s1 * m1 * ((rs + gs) * (rm + gm) / 2.0_wp + gsm + gs * gm) * one_plus_s1dgs
        temp7 = s1 * mu_m * gs * one_plus_s1dgs
        temp8 = m1 * e_cache * e_cache
        temp9 = -inv_denom * numer_m * ((s1 * ss)**2 - m1 * sm**2) &
          - 2.0_wp * m1 * s1 * one_plus_s1dgs * sm * ss

        da_dm(s,m) = -(rm + gm) / 2.0_wp &
          - inv_denom * ((temp3 - gmm - temp4) * numer_m / 2.0_wp &
          + temp5 / 4.0_wp - temp6 + temp7 + temp8 * temp1 / 4.0_wp) + temp9
      end do
    end do
    if (timing) then; call cpu_time(t1); dt(2) = t1 - t0; call cpu_time(t0); end if

    do m = 1, MDIV-1
      alpha(:,m+1) = alpha(:,m) + (mu(m+1) - mu(m)) * ( da_dm(:,m+1) + da_dm(:,m) ) * 0.5_wp
    enddo
    if (timing) then; call cpu_time(t1); dt(3) = t1 - t0; call cpu_time(t0); end if

    alpha(SDIV,:) = 0.0_wp
    adj_const = alpha(:,MDIV) - ( gama(:,MDIV) - rho(:,MDIV) )/2.0_wp
    alpha = alpha - spread(adj_const, dim=2, ncopies=MDIV)
    if (timing) then; call cpu_time(t1); dt(4) = t1 - t0; end if

    if (timing) then
      st = st + dt
      if (nc >= 4) then
        write(*,'(A,I0,A,5(1X,ES12.5))') 'update_alpha_potential call ', nc, ':', dt, sum(dt)
        write(*,'(A,I0,A)') 'update_alpha_potential running avg over ', nc, ':'
        write(*,'(A,1X,ES12.5)') '  deriv_m_sub',       st(1)/nc
        write(*,'(A,1X,ES12.5)') '  column_loop',       st(2)/nc
        write(*,'(A,1X,ES12.5)') '  integrate_mu',      st(3)/nc
        write(*,'(A,1X,ES12.5)') '  adjust_boundary',   st(4)/nc
        write(*,'(A,1X,ES12.5)') '  total',             sum(st)/nc
      end if
    end if
    if (any(alpha .ge. 300.0)) then
      write(*,*) "Error: Alpha fails in at least one row."
      error stop "alpha fails"
    end if
  end subroutine update_alpha_potential

  function cheb_r_coeff(D2_col, n_pow) result(coeff)
    real(wp), intent(in) :: D2_col(SDIV)
    integer,  intent(in) :: n_pow
    real(wp) :: coeff
    integer :: s0, n_ext, n_deg, deriv_order
    real(wp) :: R0
    real(wp), allocatable :: coeffs(:)
    type(cheb_fit_stats) :: stats
    integer, parameter :: orders = 31

    s0 = min(SDIV, count(s_gp < s_e) + 1)
    n_ext = SDIV - s0 + 1
    n_deg = min(n_ext - 1, orders)
    deriv_order = s_pwr * n_pow
    R0 = r_e * sqrt(KAPPA)
    allocate(coeffs(0:n_deg))

    call cheb_std_base(n_deg, n_ext, s_gp(s0:SDIV), D2_col(s0:SDIV), s_e, SMAX, coeffs, stats)
    if (.false.) then
      write(*,'(A,I0,A,I0,A,I0)') ' cheb_r_coeff: n_pow=', n_pow, ', deriv_order=', deriv_order, &
                                  ', n_deg=', stats%n_deg
      write(*,'(A,1X,ES12.5)') ' fitting error=', stats%rms_rel
      write(*,'(A,1X,ES12.5)') '         L_inf=', stats%linf_rel
      write(*,'(A,1X,ES12.5)') 'weight at tail=', stats%coeff_tail_l2_ratio
    end if
    ! With r = R0 * (s / (1-s))^s_pwr, the exterior tail obeys
    ! D2(s) = a_n * R0^{-n} * ((1-s)/s)^{s_pwr*n} + higher orders.
    ! Recover a_n from the matching endpoint derivative of the fitted polynomial.
    coeff = real((-1)**deriv_order, wp) * R0**n_pow * &
            cheb_get_deriv_point(n_deg, coeffs, s_e, SMAX, 1.0_wp, deriv_order) / &
            real(factorial_int(deriv_order), wp)
    deallocate(coeffs)
  end function cheb_r_coeff

  pure integer function factorial_int(n) result(val)
    integer, intent(in) :: n
    integer :: k
    val = 1
    do k = 2, n
      val = val * k
    end do
  end function factorial_int

  subroutine write_moment_tail(D2_rho, D2_omega, D2_gama)
    use para_mod, only:  M2, S3, M4, S5, M6
    real(wp), intent(in) :: D2_rho(SDIV,LMAX+1), D2_omega(SDIV,LMAX+1), D2_gama(SDIV,LMAX+1)
    character(64) :: tail_fmt
    integer :: s, unit, ios, sig_digits, field_width
    real(wp) :: r_inf, nu_monopole(SDIV)
    nu_monopole = -0.5_wp * D2_rho(:,1) - (1.0_wp / pi) * D2_gama(:,2)

    r_inf = r_e * sqrt(KAPPA) * (s_gp(SDIV - 1) / ( 1._wp - s_gp(SDIV - 1) ))**s_pwr
    M2 = - D2_metric_rho  (SDIV-1,1+1 ) / 2._wp * r_inf**3 * ( C**2 / G / Mass )**3
    S3 = - D2_metric_omega(SDIV-1,2+1 ) / 2._wp * r_inf**5 * ( C**2 / G / Mass )**4 / sqrt(KAPPA)
    M4 =   D2_metric_rho  (SDIV-1,2+1 ) / 2._wp * r_inf**5 * ( C**2 / G / Mass )**5
    
    ! alternative method 
    !M2 = - cheb_r_coeff(D2_rho  (:,1+1), 3) / 2.0_wp * ( C**2 / G / Mass )**3
    !S3 = - cheb_r_coeff(D2_omega(:,2+1), 5) / 2.0_wp * ( C**2 / G / Mass )**4 / sqrt(KAPPA)
    !M4 =   cheb_r_coeff(D2_rho  (:,2+1), 5) / 2.0_wp * ( C**2 / G / Mass )**5
    !S5 = - cheb_r_coeff(D2_omega(:,3+1), 7) / 2.0_wp * ( C**2 / G / Mass )**6 / sqrt(KAPPA)
    !M6 = - cheb_r_coeff(D2_rho  (:,3+1), 7) / 2.0_wp * ( C**2 / G / Mass )**7

    open(newunit=unit, file="Cont/moment_tail.dat", status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim("Cont/moment_tail.dat")
      return
    end if
    sig_digits = precision(1.0_wp) - 1
    field_width = sig_digits + 10
    write(tail_fmt,'("(30es",i0,".",i0,"e3)")') field_width, sig_digits
    do s = 1, SDIV-1
      write(unit,tail_fmt) s_gp(s), D2_rho(s,:), D2_omega(s,:)
    end do
    close(unit)
  end subroutine write_moment_tail

  subroutine output_helper(D2_rho, D2_omega, D2_gama)
    use para_mod, only: run_task, MRbuild, donut
    use donu_mod, only: donutization_number
    real(wp), intent(in) :: D2_rho(SDIV,LMAX+1), D2_omega(SDIV,LMAX+1), D2_gama(SDIV,LMAX+1)
    real(wp) :: radial_geom(SDIV), volume_density(SDIV,MDIV), baryon_dens(SDIV,MDIV)
    real(wp) :: rho_0, t0, t1
    character(512) :: fname
    character(len=*), parameter :: restart_binary_path = "./Res/res.rst"

    if (.not. output .or. run_task == MRbuild) return

    call cpu_time(t0); write(*,*) " "
    write(*,"(A)",advance='no') " Off-loading data ..."

    call write_moment_tail(D2_rho, D2_omega, D2_gama)

    radial_geom = radial_quad_weights * real(s_pwr, wp) * (s_gp / (1.0_wp - s_gp))**(3*s_pwr - 1) / (1.0_wp - s_gp)**2
    volume_density = exp(2.0_wp * alpha + 0.5_wp * (gama - rho))
    block
      integer :: s_, m_
      real(wp) :: n0_val, vel_safe
      do s_ = 1, SDIV
        do m_ = 1, MDIV
          n0_val = n0_at_e(energy(s_, m_))
          vel_safe = min(max(velocity_sq(s_, m_), 0.0_wp), 1.0_wp - 1.0e-12_wp)
          if (n0_val > 0.0_wp) then
            baryon_dens(s_, m_) = n0_val * MB * KSCALE * C**2 &
                                * exp(-0.75_wp * sphi(s_, m_)**2 * B_coup) &
                                / sqrt(1.0_wp - vel_safe)
          else
            baryon_dens(s_, m_) = 0.0_wp
          end if
        end do
      end do
    end block
    donut = donutization_number(sphi * sqrt(B_coup), volume_density, baryon_dens, &
                                radial_geom, angular_quad_weights)

    rho_0 = n0_at_e(energy(1,1)) * MB
    write(fname,"(A, A, A, f0.2, A, f0.3, A, es0.2e2, A,es0.2e2, A, es0.3e2, A, f0.3)") &
      "./Cont/", trim(eos_file), &
      "_J",    ang_mom, &
      "_Mb",   mass_0/MSUN, &
      "_B",    B_coup, &
      "_mphi", sqrt(mphi_r*1.e10_wp/KAPPA)*l_uni, &
      "_rhoc", rho_0, &
      "_sphim",sphi_m
    fname = adjustl(fname)
    select case(trim(solver_type))
    case("uniform")
    case("const_j");  write(fname,"(A,A,f5.2)") trim(fname), "_A",  A_diff
    case("uryu");     write(fname,"(A,A,f5.2,A,f5.2)") trim(fname), "_L1", lambda1, "_L2", lambda2
    case("uryu_inc"); write(fname,"(A,A,f5.2,A,f5.3)") trim(fname), "_Pinc", inc_p, "_Lam", inc_lambda
    case default;     stop "Unknown solver type"
    end select
    fname = trim(fname)//".dat"

    !call initial_data_for_sacra_aei(trim(fname))
    call write_restart_file(restart_binary_path)

    call cpu_time(t1); write(*,"(A, f12.6, A)") "   took ", t1-t0, " [s]"
  end subroutine output_helper

  subroutine write_restart_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios
    integer(int32) :: header_ints(6)
    real(wp) :: header_meta(5)
    real(wp), allocatable :: restart_data(:,:,:)
    character(len=*), parameter :: restart_magic = "GRASSRST01"
    integer(int32), parameter :: restart_format_version = 1_int32
    integer(int32), parameter :: restart_field_count = 10_int32

    header_ints = [restart_format_version, int(storage_size(1.0_wp), int32), restart_field_count, &
                   int(SDIV, int32), int(MDIV, int32), int(s_pwr, int32)]
    header_meta = [r_e, energy(1,1), r_ratio, Omega_e, Omega_c]

    allocate(restart_data(restart_field_count, SDIV, MDIV), source=0.0_wp)
    restart_data(1,:,:)  = alpha
    restart_data(2,:,:)  = gama
    restart_data(3,:,:)  = rho
    restart_data(4,:,:)  = ww
    restart_data(5,:,:)  = pressure
    restart_data(6,:,:)  = energy
    restart_data(7,:,:)  = enthalpy
    restart_data(8,:,:)  = velocity_sq
    restart_data(9,:,:)  = omg
    restart_data(10,:,:) = sphi

    open(newunit=unit, file=filename, status='replace', action='write', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_restart_file: failed to open file ", trim(filename)
      deallocate(restart_data)
      return
    end if

    write(unit, iostat=ios) restart_magic
    if (ios == 0) write(unit, iostat=ios) header_ints
    if (ios == 0) write(unit, iostat=ios) header_meta
    if (ios == 0) write(unit, iostat=ios) s_gp
    if (ios == 0) write(unit, iostat=ios) mu
    if (ios == 0) write(unit, iostat=ios) restart_data
    if (ios /= 0) write(*,*) "write_restart_file: failed while writing ", trim(filename)

    close(unit)
    deallocate(restart_data)
  end subroutine write_restart_file

end module spin_integration_mod
