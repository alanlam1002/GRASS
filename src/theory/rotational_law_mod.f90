module rotation_law_mod
  use precision_mod, only: wp
  implicit none
  public
  real(wp), save :: ctx_rho_equator_h, ctx_gama_equator_h, ctx_ww_equator_h
  real(wp), save :: ctx_rho_pole_h, ctx_gama_pole_h
  real(wp), save :: cached_aa = 0._wp, cached_bb = 0._wp
  real(wp), save :: ctx_exp2re2rho = 1._wp
  real(wp), save :: uryu_y_star = 0._wp        ! Fmax/F_e of law (8) in uryu_fmax_mode = "analytic" (uryu_fmax_ratio)
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

  ! ------------------------------------------------------------------------
  ! Uryu law (8) with an analytic Fmax (uryu_fmax_mode = "analytic").
  ! Omega/Omega_c = (1 + (F/B)^p) / (1 + (F/A)^(p+q)), with A and B fixed by
  ! Omega(Fmax) = lambda1 Omega_c and Omega(F_e) = lambda2 Omega_c (AA_h, BB_h).
  ! The relaxation in uryu_rotation iterates Fmax until the equatorial peak
  ! equals lambda1. Since Omega(Fmax) = lambda1 Omega_c by construction, the
  ! peak is >= lambda1 for every Fmax, with equality only where Fmax is the
  ! law's own maximum: a tangent (double) root, which the fixed-step update
  ! can miss. Imposing dOmega/dF = 0 at Fmax on the closed-form A, B instead
  ! gives p/B^p = lambda1 (p+q) Fmax^q / A^(p+q), i.e. with y = Fmax/F_e
  !   h(y) = p lambda2 (lambda1-1) - (p+q) lambda1 (lambda1-1) y^q
  !          + q lambda1 (lambda2-1) y^(p+q) = 0,
  ! independent of F_e. A^(p+q) > 0 needs y < y_max = (lambda2/lambda1)^(1/q);
  ! h(0) > 0, h(y_max) = q lambda2 ((lambda2-1) y_max^p - (lambda1-1)) < 0,
  ! and h is decreasing for lambda2 < 1: one root, so Fmax = y* F_e with no
  ! iteration, and F_e > (lambda1/lambda2)^(1/q) Fmax holds automatically
  ! (safe_F_e is never active).
  ! ------------------------------------------------------------------------
  pure real(wp) function uryu_h(y, lam1, lam2, p, q) result(h)
    real(wp), intent(in) :: y, lam1, lam2
    integer, intent(in)  :: p, q
    h = real(p, wp) * lam2 * (lam1 - 1._wp) - real(p+q, wp) * lam1 * (lam1 - 1._wp) * y**q &
        + real(q, wp) * lam1 * (lam2 - 1._wp) * y**(p+q)
  end function uryu_h

  real(wp) function uryu_fmax_ratio(lam1, lam2, p, q) result(y)
    ! y* = Fmax/F_e: bisection on (0, y_max) to machine precision
    real(wp), intent(in) :: lam1, lam2
    integer, intent(in)  :: p, q
    real(wp) :: a, b, ha, hm
    integer :: it
    if (lam1 <= 1._wp .or. lam2 <= 0._wp .or. lam2 >= 1._wp) &
      error stop "uryu_fmax_ratio: analytic Fmax needs lambda1 > 1 and 0 < lambda2 < 1"
    a = 0._wp; b = (lam2 / lam1)**(1._wp / real(q, wp))
    ha = uryu_h(a, lam1, lam2, p, q)
    if (.not. (ha > 0._wp .and. uryu_h(b, lam1, lam2, p, q) < 0._wp)) &
      error stop "uryu_fmax_ratio: h(y) does not change sign on (0, y_max)"
    do it = 1, 200
      y = 0.5_wp * (a + b)
      hm = uryu_h(y, lam1, lam2, p, q)
      if ((hm > 0._wp) .eqv. (ha > 0._wp)) then
        a = y; ha = hm
      else
        b = y
      end if
      if (b - a <= 4._wp * epsilon(1._wp) * b) exit
    end do
    y = 0.5_wp * (a + b)
  end function uryu_fmax_ratio

  subroutine diff_rotation_uryu_an(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    ! Equator equation of law (8) with Fmax = uryu_y_star * F_e (analytic mode); otherwise as diff_rotation_uryu,
    ! without the safe_F_e clamp (never needed here). F_e <= 0 (x at or below the frame-dragging rate) is outside the
    ! law's domain and returns NaN, which brent_core steps back from.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use para_mod, only: wp, lambda2, Fmax_h
    implicit none
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out):: fx
    real(wp) :: F_e, ocre, RHS, aa, bb, log_arg, Fm
    F_e  = (x - w_h) / ( exp(2._wp*re**2*rho_h) - (x-w_h)**2 )
    log_arg = 1._wp - ((x-w_h)*exp(-re**2*rho_h))**2
    if (F_e <= 0._wp .or. log_arg <= 0._wp) then
      fx = ieee_value(fx, ieee_quiet_nan); return
    end if
    ocre = x / lambda2
    Fm = uryu_y_star * F_e
    Fmax_h = Fm
    aa = AA_h(F_e, Fm)
    bb = BB_h(F_e, Fm)
    RHS = F_e * x - aa * ocre / 4._wp * &
          (2._wp * aa / bb * atan(F_e**2/aa**2) &
          - sqrt(2._wp) * ( atan(1._wp - F_e*sqrt(2._wp)/aa) - atan(1._wp + F_e*sqrt(2._wp)/aa) ) &
          + sqrt(2._wp) * ATANH( aa * F_e * sqrt(2._wp) / (F_e**2 + aa**2) ) )
    fx = re**2 * (g_h + rho_h - g_p - rho_p) + log(log_arg) + 2._wp * RHS
  end subroutine diff_rotation_uryu_an

  subroutine uryu_omega(re, rho_p, ww_p, sgp, mugp, x_guess, omega, found)
    ! Bulk Omega root of law (8) at one grid point, restricted to the physical interval w < Omega < Omega_lc (above the
    ! local frame-dragging rate, below the light cylinder, where F_j runs from 0 to +infinity):
    !   f(x) = x/Omega_c (1 + (F/A)^(p+q)) - (1 + (F/B)^p),
    ! with f -> +infinity at the light cylinder. Bracket outward from the guess (geometric steps in the distance from
    ! the guess, clipped to the interval), else a linear scan of the interval; of the sign changes, the one nearest
    ! the guess is taken (continuity with the neighbour). Refinement: Illinois regula falsi with bisection every 8
    ! iterations. found = .false. if the interval has no sign change.
    use para_mod, only: wp, uryu_p, uryu_q, Omega_c
    real(wp), intent(in)  :: re, rho_p, ww_p, sgp, mugp, x_guess
    real(wp), intent(out) :: omega
    logical, intent(out)  :: found
    integer, parameter :: NSCAN = 4000, NGROW = 200
    real(wp) :: vfac, e2, xlo, xhi, g, dlt, a, b, fa, fb, c, fc, xs(0:NSCAN), fs(0:NSCAN)
    integer :: k, it, side, best
    real(wp) :: dbest
    found = .false.; omega = x_guess
    vfac = sgp**2 * (1._wp - mugp**2)
    if (vfac <= 0._wp) then
      omega = Omega_c; found = .true.; return
    end if
    e2 = ctx_exp2re2rho
    xlo = ww_p
    xhi = ww_p + sqrt(e2) * (1._wp - sgp) / sqrt(vfac)
    dlt = (xhi - xlo) * 1.e-12_wp
    xlo = xlo + dlt; xhi = xhi - dlt
    g = min(max(x_guess, xlo), xhi)
    ! outward bracketing from the guess on both sides
    do k = 1, NGROW
      dlt = (xhi - xlo) * 1.e-6_wp * 1.5_wp**(k-1)
      a = max(xlo, g - dlt); b = min(xhi, g + dlt)
      fa = f(a); fb = f(b)
      if (fa * fb <= 0._wp .and. fa == fa .and. fb == fb) then
        ! a sign change in [a, b]: locate the one nearest g by checking each half
        if (f(g) * fa <= 0._wp) then
          b = g; fb = f(g)
        else if (f(g) * fb <= 0._wp) then
          a = g; fa = f(g)
        end if
        found = .true.; exit
      end if
      if (a <= xlo .and. b >= xhi) exit
    end do
    if (.not. found) then
      do k = 0, NSCAN
        xs(k) = xlo + (xhi - xlo) * real(k, wp) / real(NSCAN, wp); fs(k) = f(xs(k))
      end do
      best = -1; dbest = huge(1._wp)
      do k = 1, NSCAN
        if (fs(k-1) * fs(k) <= 0._wp .and. fs(k-1) == fs(k-1) .and. fs(k) == fs(k)) then
          if (abs(0.5_wp * (xs(k-1) + xs(k)) - g) < dbest) then
            dbest = abs(0.5_wp * (xs(k-1) + xs(k)) - g); best = k
          end if
        end if
      end do
      if (best < 0) return
      a = xs(best-1); fa = fs(best-1); b = xs(best); fb = fs(best); found = .true.
    end if
    if (fa == 0._wp) then; omega = a; return; end if
    if (fb == 0._wp) then; omega = b; return; end if
    side = 0
    do it = 1, 200
      if (mod(it, 8) == 0) then
        c = 0.5_wp * (a + b)
      else
        c = (a * fb - b * fa) / (fb - fa)
      end if
      fc = f(c)
      if (fc == 0._wp .or. abs(b - a) <= 1.e-14_wp * abs(c)) exit
      if (fc * fb < 0._wp) then
        a = b; fa = fb; b = c; fb = fc; side = 0
      else
        b = c; fb = fc
        if (side == -1) fa = 0.5_wp * fa
        side = -1
      end if
    end do
    omega = c
  contains
    real(wp) function f(x)
      real(wp), intent(in) :: x
      real(wp) :: Fj
      Fj = (x - ww_p) * vfac / (e2 * (1._wp - sgp)**2 - (x - ww_p)**2 * vfac)
      f = x / Omega_c * (1._wp + (Fj / cached_aa)**(uryu_p + uryu_q)) - (1._wp + (Fj / cached_bb)**uryu_p)
    end function f
  end subroutine uryu_omega
  ! Increasing Uryu law (Uryu et al. 2017 eq. 5, plus sign): j = A^2 Omega x^p,
  ! x = Omega/Omega_c - 1 >= 0, parameterized by inc_lambda = Omega_e/Omega_c and
  ! inc_p > 0 (real). G = int_{Omega_c}^{Omega} j dOmega' is the log-enthalpy term.

  pure real(wp) function uryu_inc_c(lam, p) result(c)
    real(wp), intent(in) :: lam, p
    c = ( (lam - 1._wp) / (p + 1._wp) + (lam - 1._wp)**2 / (p + 2._wp) ) / lam**2
  end function uryu_inc_c

  pure elemental real(wp) function uryu_inc_G(x, A2, omc, p) result(g)
    real(wp), intent(in) :: x, A2, omc, p
    real(wp) :: xx
    xx = max(x, 0._wp)
    g = A2 * omc**2 * ( xx**(p + 1._wp) / (p + 1._wp) + xx**(p + 2._wp) / (p + 2._wp) )
  end function uryu_inc_G

  subroutine diff_rotation_uryu_inc(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    use para_mod, only: wp, inc_p, inc_lambda
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out):: fx
    real(wp) :: F_e
    F_e = (x - w_h) / ( exp(2._wp*re**2*rho_h) - (x - w_h)**2 )
    fx = re**2 * (g_h + rho_h - g_p - rho_p) + log( 1._wp - ((x - w_h)*exp(-re**2*rho_h))**2 ) &
        + 2._wp * F_e * x * uryu_inc_c(inc_lambda, inc_p)
  end subroutine diff_rotation_uryu_inc

  ! Bulk root of F_j(Omega) - A^2 Omega x^p = 0 at one (s, mu). F_j grows without
  ! bound at the light cylinder, so the function returns to positive values there:
  ! the physical root is the first sign change above Omega_c. It is bracketed by
  ! marching up geometrically from far below the guess (no Newton: x^p has
  ! infinite slope at 0 for p < 1) and refined by Illinois regula falsi with
  ! bisection safeguards.
  subroutine uryu_inc_omega(re, rho_p, ww_p, sgp, mugp, omc, A2, p, x_guess, omega, found)
    real(wp), intent(in)  :: re, rho_p, ww_p, sgp, mugp, omc, A2, p, x_guess
    real(wp), intent(out) :: omega
    logical, intent(out), optional :: found
    real(wp) :: e2, sin2, vfac, x_up, xa, xb, xc, fa, fb, fc
    integer :: it, side
    real(wp), parameter :: XTOL = 1.e-14_wp
    integer, parameter :: NSCAN = 20000
    if (present(found)) found = .true.
    sin2 = 1._wp - mugp**2
    vfac = sgp**2 * sin2
    if (vfac <= tiny(1._wp)) then
      omega = omc; return
    end if
    e2 = exp(2._wp * re**2 * rho_p)
    x_up = ( ww_p + sqrt(e2) * (1._wp - sgp) / sqrt(vfac) ) / omc - 1._wp
    xa = 0._wp; fa = f(xa)
    if (fa < 0._wp) error stop "uryu_inc_omega: F_j < 0 at Omega_c (Omega_c below the frame-dragging rate)"
    if (fa == 0._wp .or. x_up <= 0._wp) then
      omega = omc; return
    end if
    ! Start well below the guess: a guess beyond the second (light-cylinder) root
    ! also has f > 0, and marching up from there would miss the physical root.
    xb = min(max(1.e-3_wp * x_guess, 1.e-12_wp), 0.5_wp * x_up)
    do it = 1, 2000
      fb = f(xb)
      if (fb <= 0._wp) exit
      xa = xb; fa = fb
      xb = min(1.5_wp * xb, 0.5_wp * (xb + x_up))
      if (x_up - xb < XTOL * (1._wp + x_up)) exit
    end do
    if (fb > 0._wp) then
      ! The geometric march can step over a narrow negative interval when the two
      ! roots are close (fast rotation): scan [0, x_up) linearly for the first one.
      xa = 0._wp; fa = f(xa)
      do it = 1, NSCAN
        xb = x_up * real(it, wp) / real(NSCAN + 1, wp)
        fb = f(xb)
        if (fb <= 0._wp) exit
        xa = xb; fa = fb
      end do
      if (fb > 0._wp) then
        if (present(found)) then
          found = .false.; omega = omc * (1._wp + x_guess); return
        end if
        error stop "uryu_inc_omega: no root below the light cylinder"
      end if
    end if
    side = 0
    do it = 1, 300
      if (fb == 0._wp .or. abs(xb - xa) <= XTOL * (1._wp + abs(xb))) exit
      xc = (xa * fb - xb * fa) / (fb - fa)
      if (.not. (xc > min(xa, xb) .and. xc < max(xa, xb)) .or. mod(it, 8) == 0) xc = 0.5_wp * (xa + xb)
      fc = f(xc)
      if (fc * fb < 0._wp) then
        xa = xb; fa = fb
        side = 0
      else
        if (side == 1) fa = 0.5_wp * fa
        side = 1
      end if
      xb = xc; fb = fc
    end do
    omega = omc * (1._wp + xb)
  contains
    real(wp) function f(x)
      real(wp), intent(in) :: x
      real(wp) :: o, d
      o = omc * (1._wp + x)
      d = e2 * (1._wp - sgp)**2 - (o - ww_p)**2 * vfac
      f = (o - ww_p) * vfac / d - A2 * o * x**p
    end function f
  end subroutine uryu_inc_omega

end module rotation_law_mod
