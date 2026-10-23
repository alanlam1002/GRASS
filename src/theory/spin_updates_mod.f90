module spin_updates_mod
  use precision_mod, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: update_equatorial_radius, update_angular_velocity, update_eos_and_velocity, reset_uryu_peak_cache

  ! Persistent state for narrow peak search in uryu_rotation
  real(wp), allocatable, save :: omg_mu_0_saved(:)
  integer, save :: s_peak_prev = 0

contains

  subroutine reset_uryu_peak_cache()
    s_peak_prev = 0
  end subroutine reset_uryu_peak_cache

  subroutine update_equatorial_radius(r_e_old, r_e_new, dif, sphi_pole_h, gama_pole_h, rho_pole_h, &
                                      gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h, &
                                      sphi_center_h, gama_center_h, rho_center_h)
    use para_mod, only : wp, SDIV, MDIV, s_gp, s_e, s_pwr, r_ratio, &
                         sphi, gama, rho, ww, has_scalar, B_coup, h_center, enthalpy_min
    use toolkit_mod, only : interp
    real(wp), intent(in)    :: r_e_old
    real(wp), intent(out)   :: r_e_new, dif
    real(wp), intent(out)   :: sphi_pole_h, gama_pole_h, rho_pole_h
    real(wp), intent(out)   :: gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h
    real(wp), intent(out)   :: sphi_center_h, gama_center_h, rho_center_h
    real(wp) :: r_e_new_sq, grgr
    real(wp), save :: s_p_cached = -1.0_wp
    real(wp), save :: r_ratio_prev = -1.0_wp

    if (r_ratio /= r_ratio_prev) then
      r_ratio_prev = r_ratio
      s_p_cached = r_ratio**(1.0_wp/real(s_pwr, wp)) / (1.0_wp + r_ratio**(1.0_wp/real(s_pwr, wp)))
    end if

    call interp(s_gp, sphi(:,MDIV), SDIV, s_p_cached, sphi_pole_h   )
    call interp(s_gp, gama(:,MDIV), SDIV, s_p_cached, gama_pole_h   )
    call interp(s_gp, rho(:,MDIV),  SDIV, s_p_cached, rho_pole_h    )
    call interp(s_gp, gama(:,1),    SDIV, s_e,        gama_equator_h)
    call interp(s_gp, rho(:,1),     SDIV, s_e,        rho_equator_h )
    call interp(s_gp, ww(:,1),      SDIV, s_e,        ww_equator_h  )
    call interp(s_gp, sphi(:,1),    SDIV, s_e,        sphi_equator_h)
    sphi_center_h = sphi(1,1)
    gama_center_h = gama(1,1)
    rho_center_h  = rho(1,1)

    grgr = gama_pole_h + rho_pole_h - gama_center_h - rho_center_h
    if (has_scalar) then
      grgr = grgr + B_coup / 2.0_wp * ( sphi_center_h**2 - sphi_pole_h**2 )
    end if

    r_e_new_sq = ( 2.0_wp * ( h_center - enthalpy_min ) ) / grgr
    r_e_new = sqrt( r_e_new_sq )
    dif = abs(r_e_old - r_e_new) / r_e_new

    if (r_e_new / r_e_old > 2 .or. ieee_is_nan(r_e_new) ) then
      write(*,"(a,es22.14)") " r_e_old        :", r_e_old
      write(*,"(a,es22.14)") " r_e_new        :", r_e_new
      write(*,"(a,es22.14)") " grgr           :", grgr
      write(*,"(a,es22.14)") " r_e_new_sq     :", r_e_new_sq
      write(*,"(a,es22.14)") " gama_pole_h    :", gama_pole_h
      write(*,"(a,es22.14)") " rho_pole_h     :", rho_pole_h
      write(*,"(a,es22.14)") " gama_center_h  :", gama_center_h
      write(*,"(a,es22.14)") " rho_center_h   :", rho_center_h
      write(*,"(a,es22.14)") " gama_equator_h :", gama_equator_h
      write(*,"(a,es22.14)") " rho_equator_h  :", rho_equator_h
      write(*,"(a,es22.14)") " sphi_pole_h    :", sphi_pole_h
      write(*,"(a,es22.14)") " sphi_center_h  :", sphi_center_h
      write(*,"(a,es22.14)") " h_center       :", h_center
      write(*,"(a,es22.14)") " enthalpy_min   :", enthalpy_min
      write(*,"(a,es22.14)") " s_e            :", s_e
      write(*,"(a,es22.14)") " s_p_cached     :", s_p_cached
      write(*,"(a,es22.14)") " r_ratio        :", r_ratio
      if (has_scalar) write(*,"(a,es22.14)") " B_coup         :", B_coup
      error stop 'r_e cannot be found.'
    endif
  end subroutine update_equatorial_radius

  subroutine update_angular_velocity(r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, &
                                    sphi_pole_h, sphi_equator_h, ww_equator_h)
    use para_mod, only : wp, SDIV, MDIV, s_gp, mu, rho, ww, omg, r_ratio, &
                         Omega_c, Omega_e, Omg, F_j, has_scalar, B_coup, A_diff, solver_type, &
                         lambda1, lambda2, uryu_p, uryu_q, uryu_fmax_mode, enthalpy, enthalpy_min, &
                         Fmax_h, F_equator_h, timing
    use rotation_law_mod, only: diff_rotation_const_j, rotation_law_const_j, &
                                diff_rotation_uryu, rotation_law_uryu, &
                                cache_uryu_ab, ctx_exp2re2rho, &
                                diff_rotation_uryu_inc, uryu_inc_omega, &
                                uryu_fmax_ratio, uryu_y_star, diff_rotation_uryu_an, uryu_omega
    use brent_mod, only : find_omega_e, zbrent_rot
    real(wp), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h
    real(wp), intent(in) :: sphi_pole_h, sphi_equator_h, ww_equator_h
    real(wp) :: metric_diff, term_in_Omega_h
    real(wp), parameter :: TOLERANCE_SPHERICAL = 1.0e-3_wp
    real(wp), parameter :: GUESS_FACTOR = 0.8_wp
    real(wp), parameter :: TOLERANCE_ROOT = 1.e-5_wp
    real(wp), parameter :: TOLERANCE_FMAX = 1.e-7_wp
    real(wp), parameter :: FMAX_INITIAL = 1.0e-2_wp
    real(wp), parameter :: FMAX_STEP = 1.e-2_wp
    integer :: s, m

    if (abs(r_ratio - 1.0_wp) < TOLERANCE_SPHERICAL) then
      Omega_c = 0.0_wp; Omega_e = 0.0_wp; Omg = 0.0_wp
      return
    end if

    associate(re2_val => r_e_new**2)
      metric_diff = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h &
                  + B_coup / 2.0_wp * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.0_wp - exp( re2_val * metric_diff )
      if (term_in_Omega_h >= 0.0_wp) then
        Omega_e = ww_equator_h + exp(re2_val * rho_equator_h) * sqrt(term_in_Omega_h)
      else
        write(*,"('Solving for axis ratio: ', f12.5)") r_ratio
        write(*,"(10A15)") "gama_pole", "rho_pole", "gama_equator", "rho_equator", "sphi_pole", "sphi_equator"
        write(*,"(10es15.3)") gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, sphi_pole_h, sphi_equator_h
        error stop "Omega can't be found; Line 99 of spin helper"
      endif

      select case(trim(solver_type))
      case("uniform")
        Omega_c = Omega_e; Omg = Omega_e
      case("const_j")
        call const_j_rotation(re2_val)
      case("uryu")
        call uryu_rotation(re2_val)
      case("uryu_inc")
        call uryu_inc_rotation(re2_val)
      case default
        error stop "Unknown solver type"
      end select
    end associate
  contains
    subroutine uryu_inc_rotation(re2_val)
      use para_mod, only: inc_p, inc_lambda, inc_A2, s_e
      real(wp), intent(in) :: re2_val
      real(wp) :: guess, F_e
      logical :: found

      guess = Omega_e * GUESS_FACTOR
      call find_omega_e(guess, r_e_new, rho_equator_h, gama_equator_h, &
                      ww_equator_h, rho_pole_h, gama_pole_h, TOLERANCE_ROOT, Omega_e, diff_rotation_uryu_inc)
      F_e = (Omega_e - ww_equator_h) / ( exp(2.0_wp * re2_val * rho_equator_h) - (Omega_e - ww_equator_h)**2 )
      if (F_e <= 0.0_wp) error stop "uryu_inc_rotation: F_e <= 0 at the equator"
      Omega_c = Omega_e / inc_lambda
      inc_A2 = F_e / ( Omega_e * (inc_lambda - 1.0_wp)**inc_p )

      Omg(1,:) = Omega_c
      Omg(1:2*SDIV/3,MDIV) = Omega_c
      F_j(1,:) = 0.0_wp
      F_j(1:2*SDIV/3,MDIV) = 0.0_wp
      do s = 2, 2*SDIV/3
        do m = 1, MDIV-1
          associate(sg => s_gp(s), mum_loc => mu(m), &
                    rsm_loc => rho(s,m), wwsm_loc => ww(s,m), o => omg(s,m))
            call uryu_inc_omega(r_e_new, rsm_loc, wwsm_loc, sg, mum_loc, Omega_c, inc_A2, inc_p, &
                                Omg(s-1,m) / Omega_c - 1.0_wp, o, found)
            if (.not. found) then
              ! Beyond the equatorial surface the increasing law can need Omega above
              ! the light cylinder; no matter is there, so keep the neighbour's value.
              if (sg <= s_e) error stop "uryu_inc_rotation: no Omega root inside the star"
              o = Omg(s-1,m)
            end if
            ctx_exp2re2rho = exp(2.0_wp * re2_val * rsm_loc)
            F_j(s,m) = (o - wwsm_loc) * sg**2 * (1.0_wp - mum_loc**2) &
                  / (ctx_exp2re2rho * (1.0_wp - sg)**2 - (o - wwsm_loc)**2 * sg**2 * (1.0_wp - mum_loc**2))
          end associate
        end do
      end do
    end subroutine uryu_inc_rotation
    subroutine const_j_rotation(re2_val)
      real(wp), intent(in) :: re2_val
      real(wp) :: guess, term_omega_diff

      guess = Omega_e * GUESS_FACTOR
      call find_omega_e(guess, r_e_new, rho_equator_h, gama_equator_h, &
                      ww_equator_h, rho_pole_h, gama_pole_h, TOLERANCE_ROOT, Omega_e, diff_rotation_const_j)

      term_omega_diff = abs(Omega_e - ww_equator_h)
      term_in_Omega_h = term_omega_diff * exp(-2.0_wp * re2_val * rho_equator_h)
      Omega_c = Omega_e + term_in_Omega_h / (1.0_wp - term_in_Omega_h * term_omega_diff) / A_diff**2

      Omg(1,:) = Omega_c
      Omg(1:2*SDIV/3,MDIV) = Omega_c
      do s = 2, 2*SDIV/3
        do m = 1, MDIV-1
          associate(sg => s_gp(s), mum_loc => mu(m), &
                    rsm_loc => rho(s,m), wwsm_loc => ww(s,m), o => omg(s,m))
            ctx_exp2re2rho = exp(2.0_wp * re2_val * rsm_loc)
            call zbrent_rot(Omg(s-1,m) * GUESS_FACTOR, r_e_new, rsm_loc, wwsm_loc, sg, mum_loc, TOLERANCE_ROOT, o, rotation_law_const_j)
            F_j(s,m) = (o - wwsm_loc) * sg**2 * (1.0_wp - mum_loc**2) &
                  / (ctx_exp2re2rho * (1.0_wp - sg)**2 - (o - wwsm_loc)**2 * sg**2 * (1.0_wp - mum_loc**2))
          end associate
        end do
      end do
    end subroutine const_j_rotation
    subroutine uryu_rotation(re2_val)
      real(wp), intent(in) :: re2_val
      real(wp) :: diff_Fmax, guess, Fa, omg_max_h, exp_term_eq, mum, e2rr
      integer :: s_lo, s_hi, s_peak
      integer, parameter :: DELTA_PEAK = 3

      if (trim(uryu_fmax_mode) == "analytic") then
        call uryu_rotation_analytic(re2_val); return
      end if
      if (.not. allocated(omg_mu_0_saved)) allocate(omg_mu_0_saved(SDIV), source=0.0_wp)
      exp_term_eq = exp(2.0_wp * re2_val * rho_equator_h)

      diff_Fmax = 1.0_wp
      Fmax_h = FMAX_INITIAL
      do while(abs(diff_Fmax) > TOLERANCE_FMAX)
        guess = Omega_e
        call find_omega_e(guess, r_e_new, rho_equator_h, gama_equator_h, ww_equator_h, &
                        rho_pole_h, gama_pole_h, TOLERANCE_ROOT, Fa, diff_rotation_uryu)
        Omega_e = Fa
        F_equator_h = (Omega_e - ww_equator_h) / ( exp_term_eq - (Omega_e - ww_equator_h)**2 )
        if ( F_equator_h < 0.0_wp ) stop "negative F_equator_h; L120 in uryu"

        Omega_c = Omega_e / lambda2
        call cache_uryu_ab()
        mum = 0.0_wp
        if (s_peak_prev == 0) then
          omg_mu_0_saved(1) = Omega_c
          omg_max_h = Omega_c
          s_peak = 1
          do s = 2, (SDIV-1)/2
            associate(o => omg_mu_0_saved(s), o_prev => omg_mu_0_saved(s-1), sg => s_gp(s), &
                      rsm_loc => rho(s,1), wwsm_loc => ww(s,1))
              guess = o_prev
              ctx_exp2re2rho = exp(2.0_wp * re2_val * rsm_loc)
              call zbrent_rot(guess, r_e_new, rsm_loc, wwsm_loc, sg, mum, TOLERANCE_ROOT, o, rotation_law_uryu)
              if (o > omg_max_h) then
                omg_max_h = o
                s_peak = s
              else
                exit
              end if
            end associate
          end do
          s_peak_prev = s_peak
        else
          s_lo = max(2, s_peak_prev - DELTA_PEAK)
          s_hi = min((SDIV-1)/2, s_peak_prev + DELTA_PEAK)
          omg_max_h = Omega_c
          s_peak = 1
          do s = s_lo, s_hi
            associate(o => omg_mu_0_saved(s), sg => s_gp(s), &
                      rsm_loc => rho(s,1), wwsm_loc => ww(s,1))
              guess = o
              ctx_exp2re2rho = exp(2.0_wp * re2_val * rsm_loc)
              call zbrent_rot(guess, r_e_new, rsm_loc, wwsm_loc, sg, mum, TOLERANCE_ROOT, o, rotation_law_uryu)
              if (o > omg_max_h) then
                omg_max_h = o
                s_peak = s
              end if
            end associate
          end do
          s_peak_prev = s_peak
        end if
        diff_Fmax = ( lambda1 - omg_max_h / Omega_c )
        Fmax_h = Fmax_h - diff_Fmax * FMAX_STEP
      end do

      Omg(1,:) = Omega_c
      Omg(1:3*SDIV/4,MDIV) = Omega_c
      do s = 2, SDIV*3/4
        e2rr = exp(2.0_wp * re2_val * rho(s,1))
        do m = 1, MDIV-1
          associate(sg => s_gp(s), mum_loc => mu(m), &
                    rsm_loc => rho(s,m), wwsm_loc => ww(s,m), o => omg(s,m))
            guess = Omg(s-1,m)
            ctx_exp2re2rho = exp(2.0_wp * re2_val * rsm_loc)
            call zbrent_rot(guess, r_e_new, rsm_loc, wwsm_loc, sg, mum_loc, TOLERANCE_ROOT, o, rotation_law_uryu)
            F_j(s,m) = (o - wwsm_loc) * sg**2 * (1.0_wp - mum_loc**2) &
                  / (ctx_exp2re2rho * (1.0_wp - sg)**2 - (o - wwsm_loc)**2 * sg**2 * (1.0_wp - mum_loc**2))
          end associate
        end do
      end do
    end subroutine uryu_rotation

    subroutine uryu_rotation_analytic(re2_val)
      ! Law (8) with Fmax = y* F_e (rotation_law_mod::uryu_fmax_ratio): one equator root, no Fmax iteration, then the
      ! bulk Omega root per point by rotation_law_mod::uryu_omega (bracketed inside the physical interval).
      real(wp), intent(in) :: re2_val
      real(wp) :: exp_term_eq, guess
      logical :: found
      exp_term_eq = exp(2.0_wp * re2_val * rho_equator_h)
      uryu_y_star = uryu_fmax_ratio(lambda1, lambda2, uryu_p, uryu_q)
      guess = Omega_e
      call find_omega_e(guess, r_e_new, rho_equator_h, gama_equator_h, ww_equator_h, &
                        rho_pole_h, gama_pole_h, TOLERANCE_ROOT, Omega_e, diff_rotation_uryu_an)
      F_equator_h = (Omega_e - ww_equator_h) / ( exp_term_eq - (Omega_e - ww_equator_h)**2 )
      if ( F_equator_h <= 0.0_wp ) error stop "uryu_rotation_analytic: F_e <= 0 at the equator"
      Fmax_h = uryu_y_star * F_equator_h
      Omega_c = Omega_e / lambda2
      call cache_uryu_ab()
      Omg(1,:) = Omega_c
      Omg(1:3*SDIV/4,MDIV) = Omega_c
      do s = 2, SDIV*3/4
        do m = 1, MDIV-1
          associate(sg => s_gp(s), mum_loc => mu(m), &
                    rsm_loc => rho(s,m), wwsm_loc => ww(s,m), o => omg(s,m))
            ctx_exp2re2rho = exp(2.0_wp * re2_val * rsm_loc)
            call uryu_omega(r_e_new, rsm_loc, wwsm_loc, sg, mum_loc, Omg(s-1,m), o, found)
            if (.not. found) then
              if (enthalpy(s,m) > enthalpy_min) error stop "uryu_rotation_analytic: no Omega root inside the star"
              o = Omg(s-1,m)
            end if
            F_j(s,m) = (o - wwsm_loc) * sg**2 * (1.0_wp - mum_loc**2) &
                  / (ctx_exp2re2rho * (1.0_wp - sg)**2 - (o - wwsm_loc)**2 * sg**2 * (1.0_wp - mum_loc**2))
          end associate
        end do
      end do
    end subroutine uryu_rotation_analytic
  end subroutine update_angular_velocity

  subroutine update_eos_and_velocity(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, &
                                     sgp_term_2d_cache_arg, sin_theta_2d_cache_arg, sgp_2d_cache_arg)
    use para_mod, only : wp, SDIV, MDIV, r_ratio, h_center, &
                         rho, gama, alpha, sphi, velocity_sq, enthalpy, pressure, energy, &
                         Omg, ww, Omega_c, omg, F_j, &
                         s_gp, has_scalar, B_coup, A_diff, solver_type, enthalpy_min, s_e, timing, &
                         inc_p, inc_A2
    use rotation_law_mod, only: intF, uryu_inc_G
    use eos_mod, only : pe_at_h
    real(wp), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(wp), intent(in) :: sgp_term_2d_cache_arg(:,:), sin_theta_2d_cache_arg(:,:), sgp_2d_cache_arg(:,:)
    integer :: s, m, eos_idx_hint, s_active_max
    real(wp) :: re2
    real(wp) :: t0, t1, dt_hydro, dt_rotlaw, dt_eos, dt_rescale
    integer, parameter :: timing_calls = 5
    integer, save :: eos_update_call_count = 0
    real(wp), save :: sum_dt_hydro = 0.0_wp, sum_dt_rotlaw = 0.0_wp
    real(wp), save :: sum_dt_eos = 0.0_wp, sum_dt_rescale = 0.0_wp

    if (timing) then
      if (eos_update_call_count == 0) then
        sum_dt_hydro = 0.0_wp
        sum_dt_rotlaw = 0.0_wp
        sum_dt_eos = 0.0_wp
        sum_dt_rescale = 0.0_wp
      end if
      eos_update_call_count = eos_update_call_count + 1
      dt_hydro = 0.0_wp
      dt_rotlaw = 0.0_wp
      dt_eos = 0.0_wp
      dt_rescale = 0.0_wp
      call cpu_time(t0)
    end if

    re2 = r_e_new**2

    if (abs(r_ratio - 1.0_wp) < epsilon(r_ratio)) then
      velocity_sq = 0.0_wp
      enthalpy = enthalpy_min + 0.5_wp * re2 * &
               ( gama_pole_h + rho_pole_h - gama - rho + ( sphi**2 - sphi_pole_h**2 ) * B_coup / 2.0_wp )
    else
      velocity_sq = ((Omg - ww) * sgp_term_2d_cache_arg * sin_theta_2d_cache_arg * exp(-rho * re2))**2
      where (velocity_sq > 1.0_wp) velocity_sq = 0.0_wp
      enthalpy = enthalpy_min + 0.5_wp * ( &
            re2 * ( gama_pole_h + rho_pole_h - gama - rho + ( sphi**2 - sphi_pole_h**2 ) * B_coup / 2.0_wp ) &
            - log( max(1.e-30_wp, 1.0_wp-velocity_sq) )  )
    end if

    if (timing) then
      call cpu_time(t1); dt_hydro = t1 - t0; call cpu_time(t0)
    end if

    s_active_max = count(s_gp <= s_e)

    if ( trim(solver_type) == "const_j" ) then
      enthalpy = enthalpy + 0.5_wp * A_diff**2 * (Omg - Omega_c)**2
    elseif ( trim(solver_type) == "uryu" ) then
      if (s_active_max < SDIV) F_j(s_active_max+1:SDIV,:) = 0.0_wp
      do m = 1, MDIV
        do s = 1, s_active_max
          enthalpy(s,m) = enthalpy(s,m) - intF(omg(s,m), F_j(s,m))
        end do
      end do
    elseif ( trim(solver_type) == "uryu_inc" ) then
      do m = 1, MDIV
        do s = 1, s_active_max
          enthalpy(s,m) = enthalpy(s,m) - uryu_inc_G(omg(s,m) / Omega_c - 1.0_wp, inc_A2, Omega_c, inc_p)
        end do
      end do
    endif

    if (timing) then
      call cpu_time(t1); dt_rotlaw = t1 - t0; call cpu_time(t0)
    end if

    pressure = 0.0_wp
    energy = 0.0_wp
    if (s_active_max < SDIV) enthalpy(s_active_max+1:SDIV,:) = enthalpy_min

    do m = 1, MDIV
      eos_idx_hint = 1
      do s = 1, s_active_max
        associate( h => enthalpy(s,m), p => pressure(s,m), e => energy(s,m) )
        if (h > enthalpy_min) then
          call pe_at_h(h, p, e, eos_idx_hint)  ! Hint carries forward within meridian m
        else
          h = enthalpy_min
          eos_idx_hint = 1  ! Reset on invalid point to keep hint accurate
        end if
        end associate
      end do
    end do

    if (timing) then
      call cpu_time(t1); dt_eos = t1 - t0; call cpu_time(t0)
    end if

    rho   = rho   * re2
    gama  = gama  * re2
    alpha = alpha * re2
    sphi  = sphi  * r_e_new

    if (timing) then
      call cpu_time(t1); dt_rescale = t1 - t0
      sum_dt_hydro = sum_dt_hydro + dt_hydro
      sum_dt_rotlaw = sum_dt_rotlaw + dt_rotlaw
      sum_dt_eos = sum_dt_eos + dt_eos
      sum_dt_rescale = sum_dt_rescale + dt_rescale
      if (eos_update_call_count >= timing_calls) then
        write(*,'(A,I0,A)') 'update_eos_and_velocity avg over ', timing_calls, ':'
        write(*,'(A,1X,ES12.5)') '  hydro', sum_dt_hydro / timing_calls
        write(*,'(A,1X,ES12.5)') '  rotation_law', sum_dt_rotlaw / timing_calls
        write(*,'(A,1X,ES12.5)') '  eos_loop', sum_dt_eos / timing_calls
        write(*,'(A,1X,ES12.5)') '  rescale', sum_dt_rescale / timing_calls
        write(*,'(A,1X,ES12.5)') '  total', &
          (sum_dt_hydro + sum_dt_rotlaw + sum_dt_eos + sum_dt_rescale) / timing_calls
      end if
    end if
  end subroutine update_eos_and_velocity

end module spin_updates_mod
