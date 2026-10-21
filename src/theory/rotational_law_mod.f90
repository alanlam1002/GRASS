module rotation_law_mod
  use precision_mod, only: wp
  implicit none
  public
  real(wp), save :: ctx_rho_equator_h, ctx_gama_equator_h, ctx_ww_equator_h
  real(wp), save :: ctx_rho_pole_h, ctx_gama_pole_h
  real(wp), save :: cached_aa = 0._wp, cached_bb = 0._wp
  real(wp), save :: ctx_exp2re2rho = 1._wp
  private :: AA_h, BB_h
contains

  ! AA_h's defining exponent 1/(uryu_p+uryu_q) is fractional, so its
  ! bracketed base -- (F_e*Fmax_h)**uryu_p * tmp1/tmp2, tmp1 = lambda2*
  ! F_e**uryu_q - lambda1*Fmax_h**uryu_q, tmp2 = (lambda1-1)*F_e**uryu_p -
  ! (lambda2-1)*Fmax_h**uryu_p -- must stay non-negative or AA_h silently
  ! returns NaN (Fortran's real**noninteger is undefined for a negative
  ! base). For lambda1>1>lambda2>0 (the defaults), tmp2>0 whenever F_e >
  ! -(1-lambda2)/(lambda1-1)*Fmax_h (comfortably true once F_e is
  ! positive), but tmp1 is negative for the whole band
  ! 0 <= F_e < (lambda1/lambda2)**(1/uryu_q)*Fmax_h -- i.e. NaN hits for
  ! any F_e near zero, not just F_e<0. This band widens as Fmax_h grows
  ! during the outer Fmax-shooting loop. F_e<0 is also already treated as
  ! outside this rotation law's physical domain elsewhere (the "negative
  ! F_equator_h" guard in spin_updates_mod.f90's uryu_rotation), so
  ! clamping F_e comfortably above the whole unsafe band -- rather than
  ! just above zero -- keeps AA_h/BB_h evaluable everywhere they're called
  ! with an F_e that could land in this band, holding the nearest physical
  ! boundary value flat instead. Shared by diff_rotation_uryu's search (any
  ! trial x Brent should never accept as the final root anyway) and
  ! cache_uryu_ab (the raw, unclamped F_equator_h global can land here too
  ! once find_omega_e's own root is found, which previously left
  ! cached_aa/cached_bb NaN and poisoned every downstream rotation_law_uryu
  ! call in the peak search).
  pure real(wp) function safe_F_e(F_e, Fmax_h) result(F_e_safe)
    use para_mod, only: wp, lambda1, lambda2, uryu_q
    real(wp), intent(in) :: F_e, Fmax_h
    real(wp), parameter :: MARGIN = 1.0_wp + 1.0e-6_wp
    real(wp) :: F_e_min
    F_e_min = MARGIN * (lambda1/lambda2)**(1.0_wp/real(uryu_q, wp)) * Fmax_h
    F_e_safe = max(F_e, F_e_min)
  end function safe_F_e

  subroutine cache_uryu_ab()
    use para_mod, only: wp, F_equator_h, Fmax_h
    real(wp) :: F_e_safe
    F_e_safe = safe_F_e(F_equator_h, Fmax_h)
    cached_aa = AA_h(F_e_safe, Fmax_h)
    cached_bb = BB_h(F_e_safe, Fmax_h)
  end subroutine cache_uryu_ab

  subroutine diff_rotation_const_j(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    use para_mod, only: wp, A_diff
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out):: fx
    real(wp) :: LHS
    LHS =  re**2 * (g_h + rho_h - g_p - rho_p) + log(1._wp - ((x-w_h)*exp(-re**2*rho_h))**2 )
    fx = LHS * (1._wp - ((x-w_h)*exp(-re**2*rho_h))**2 )**2 &
        - ( (x-w_h) * exp(-2._wp*re**2*rho_h) / A_diff )**2
  end subroutine diff_rotation_const_j

  subroutine rotation_law_const_j(x, fx, re, rho_p, ww_p, sgp, mugp)
    use para_mod, only: wp, A_diff, Omega_c
    real(wp), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(wp), intent(out):: fx
    real(wp) :: sin2m, inv_e2
    sin2m = 1._wp - mugp * mugp
    inv_e2 = 1._wp / ctx_exp2re2rho
    fx = (Omega_c - x) * ( (1._wp-sgp)**2 - ( sgp*(x-ww_p) )**2 * inv_e2 * sin2m ) &
        - ( (x-ww_p) * sgp * sgp * sin2m * inv_e2 ) / A_diff**2
  end subroutine rotation_law_const_j

  subroutine diff_rotation_uryu(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    use para_mod, only: wp, lambda2, Fmax_h
    implicit none
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out):: fx
    real(wp) :: F_e, ocre, RHS, aa, bb, log_arg
    real(wp), parameter :: LOG_ARG_FLOOR = 1.0e-6_wp

    ocre = x / lambda2
    F_e  = (x - w_h) / ( exp(2._wp*re**2*rho_h) - (x-w_h)**2 )
    log_arg = 1._wp - ((x-w_h)*exp(-re**2*rho_h))**2

    ! See safe_F_e's docstring: keeps AA_h/BB_h clear of their near-F_e=0
    ! NaN band, which find_omega_e's bracket search (starting right at the
    ! F_e~0 equator-seed regime for a near-spherical restart) would
    ! otherwise walk straight into.
    F_e = safe_F_e(F_e, Fmax_h)

    aa = AA_h(F_e, Fmax_h)
    bb = BB_h(F_e, Fmax_h)

    RHS = F_e * x - aa * ocre / 4._wp * &
          (2._wp * aa / bb * atan(F_e**2/aa**2) &
          - sqrt(2._wp) * ( atan(1._wp - F_e*sqrt(2._wp)/aa) - atan(1._wp + F_e*sqrt(2._wp)/aa) ) &
          + sqrt(2._wp) * ATANH( aa * F_e * sqrt(2._wp) / (F_e**2 + aa**2) ) )

    fx = re**2 * (g_h + rho_h - g_p - rho_p) + log( max(log_arg, LOG_ARG_FLOOR) ) &
        + 2._wp * RHS
  end subroutine diff_rotation_uryu

  subroutine rotation_law_uryu(x, fx, re, rho_p, ww_p, sgp, mugp)
    use para_mod, only: wp, uryu_p, uryu_q, Omega_c
    implicit none
    real(wp), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(wp), intent(out):: fx
    real(wp) :: Fj, sin2m

    sin2m = 1._wp - mugp**2
    Fj = (x-ww_p) * sgp**2 * sin2m / ( ctx_exp2re2rho * (1._wp-sgp)**2 - (x-ww_p)**2 * sgp**2 * sin2m )
    fx = x / Omega_c * ( 1._wp + (Fj / cached_aa)**(uryu_p+uryu_q) ) - ( 1._wp + (Fj / cached_bb)**uryu_p )
  end subroutine rotation_law_uryu

  elemental real(wp) function intF(x, F_at_x)
    use para_mod, only: wp, Omega_c
    implicit none
    real(wp), intent(in) :: x, F_at_x

    if ( abs(x) < epsilon(x) .and. abs(F_at_x) < epsilon(F_at_x) ) then
      intF = 0._wp
    else
      intF = F_at_x * x - cached_aa * omega_c / 4._wp * &
          (2._wp * cached_aa / cached_bb * atan(F_at_x**2/cached_aa**2) &
          - sqrt(2._wp) * ( atan(1._wp-F_at_x*sqrt(2._wp)/cached_aa) - atan(1._wp+F_at_x*sqrt(2._wp)/cached_aa) ) &
          + sqrt(2._wp) * ATANH(cached_aa*F_at_x*sqrt(2._wp)/(F_at_x**2 + cached_aa**2)) )
    endif
  end function intF

  pure elemental real(wp) function AA_h(F_e, F_m)
    use para_mod, only: wp, lambda1, lambda2, uryu_p, uryu_q
    implicit none
    real(wp), intent(in) :: F_e, F_m
    real(wp) :: tmp1, tmp2

    tmp1 = lambda2 * F_e**uryu_q &
        - lambda1 * F_m**uryu_q
    tmp2 = (lambda1-1._wp) * F_e**uryu_p &
        - (lambda2-1._wp) * F_m**uryu_p
    AA_h = ( (F_e*F_m)**uryu_p * tmp1 / tmp2 )**(1._wp/real(uryu_p+uryu_q, wp))
  end function

  pure elemental real(wp) function BB_h(F_e, F_m)
    use para_mod, only: wp, lambda1, lambda2, uryu_p, uryu_q
    implicit none
    real(wp), intent(in) :: F_e, F_m
    real(wp) :: tmp1, tmp2

    tmp1 = lambda2 * F_e**uryu_q &
        - lambda1 * F_m**uryu_q
    tmp2 = lambda2 * (lambda1-1._wp) * F_e**(uryu_p+uryu_q) &
        - lambda1 * (lambda2-1._wp) * F_m**(uryu_p+uryu_q)
    BB_h = F_e * F_m * ( tmp1 / tmp2 )**(1._wp/real(uryu_p, wp))
  end function

  subroutine set_shoot_context(rho_equator_h, gama_equator_h, ww_equator_h, rho_pole_h, gama_pole_h)
    real(wp), intent(in) :: rho_equator_h, gama_equator_h, ww_equator_h, rho_pole_h, gama_pole_h
    ctx_rho_equator_h  = rho_equator_h
    ctx_gama_equator_h = gama_equator_h
    ctx_ww_equator_h   = ww_equator_h
    ctx_rho_pole_h     = rho_pole_h
    ctx_gama_pole_h    = gama_pole_h
  end subroutine set_shoot_context

  subroutine shoot_Fmax(current_Fmax_h, diff_sq)
    use para_mod, only: wp, Fmax_h, Omega_e, Omega_c, r_e, rho, ww, F_equator_h, &
                        SDIV, s_gp, lambda1, lambda2
    use brent_mod, only: find_omega_e, zbrent_rot
    real(wp), intent(in)  :: current_Fmax_h
    real(wp), intent(out) :: diff_sq
    real(wp) :: guess, Fa, rsm, wwsm, sgp, mum, omg_max_h, diff_Fmax
    real(wp), dimension(SDIV) :: omg_mu_0
    integer :: imax, s
    real(wp), parameter :: TOLERANCE = 1.e-5

    Fmax_h = current_Fmax_h

    guess = Omega_e
    call find_omega_e(guess, r_e, ctx_rho_equator_h, ctx_gama_equator_h, ctx_ww_equator_h, &
                    ctx_rho_pole_h, ctx_gama_pole_h, tolerance, Fa, diff_rotation_uryu)
    Omega_e = Fa
    F_equator_h = (Omega_e - ww(1,1)) / (exp(2._wp*r_e**2*rho(1,1)) - (Omega_e-ww(1,1))**2)
    if (F_equator_h < 0._wp) then
      diff_sq = 1.e3_wp
      return
    endif

    Omega_c = Omega_e / lambda2
    call cache_uryu_ab()
    omg_mu_0(1) = Omega_c
    mum = 0._wp
    do s = 2, SDIV*2/3
      rsm = rho(s,1)
      wwsm= ww (s,1)
      sgp = s_gp(s)
      guess  = omg_mu_0(s-1)
      ctx_exp2re2rho = exp(2._wp * r_e**2 * rsm)
      call zbrent_rot( guess, r_e, rsm, wwsm, sgp, mum, 1.e-5_wp, omg_mu_0(s), rotation_law_uryu)
    enddo
    imax = maxloc(omg_mu_0(1:SDIV*2/3), 1)
    omg_max_h = omg_mu_0(imax)
    diff_Fmax = (lambda1 - omg_max_h / Omega_c)
    diff_sq = diff_Fmax**2
  end subroutine shoot_Fmax
end module rotation_law_mod
