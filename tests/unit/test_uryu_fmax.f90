program test_uryu_fmax
  ! Uryu law (8) with the analytic Fmax (uryu_fmax_mode = "analytic"):
  !  (1) y* = uryu_fmax_ratio(lambda1, lambda2, p, q) makes Fmax = y* F_e the law's maximum: with A, B from the
  !      value conditions (the closed forms of rotation_law_mod's AA_h/BB_h), Omega(Fmax) = lambda1 Omega_c,
  !      Omega(F_e) = lambda2 Omega_c and dOmega/dF(Fmax) = 0 (central difference), for several (lambda1, lambda2)
  !      and F_e (y* does not depend on F_e), and Omega(F) <= lambda1 Omega_c on a dense F grid;
  !  (2) uryu_omega against a brute-force scan of the same f(x) in a flat metric (ww = 0, exp(2 re^2 rho) = 1): the
  !      returned root satisfies |f| small and is the sign change nearest the guess; the axis returns Omega_c.
  use precision_mod, only: wp
  use test_utils
  use para_mod, only: Omega_c, uryu_p, uryu_q
  use rotation_law_mod, only: uryu_fmax_ratio, uryu_h, uryu_omega, cached_aa, cached_bb, ctx_exp2re2rho
  implicit none
  real(wp), parameter :: lam(2, 5) = reshape([1.5_wp, 0.3_wp, 1.5_wp, 0.45_wp, 2.0_wp, 0.5_wp, 2.2_wp, 0.6_wp, &
                                              1.8_wp, 0.4_wp], [2, 5])
  real(wp), parameter :: Fes(3) = [0.05_wp, 0.25_wp, 1.3_wp], gfac(3) = [1.0_wp, 1.6_wp, 0.6_wp]
  real(wp) :: l1, l2, y, Fe, Fm, A, B, d, dO, Om_max, Fk, omega, sgp, mu, guess, xlo, xhi, fbest, x, fx, fprev, dist
  integer :: i, j, k, nbest
  logical :: found
  character(len=80) :: lab

  uryu_p = 1; uryu_q = 3; Omega_c = 0.7_wp
  do i = 1, size(lam, 2)
    l1 = lam(1, i); l2 = lam(2, i)
    y = uryu_fmax_ratio(l1, l2, 1, 3)
    write(lab, '(a,f4.2,a,f4.2)') 'h(y*) = 0, lambda ', l1, '/', l2
    call assert_near(trim(lab), 0.0_wp, uryu_h(y, l1, l2, 1, 3), 1.e-13_wp)
    call assert_true('y* inside (0, y_max)', y > 0.0_wp .and. y < (l2 / l1)**(1.0_wp / 3.0_wp))
    do j = 1, size(Fes)
      Fe = Fes(j); Fm = y * Fe
      call ab(Fe, Fm, l1, l2, A, B)
      write(lab, '(a,f4.2,a,f4.2,a,f4.2)') 'Omega(Fmax) = lambda1, ', l1, '/', l2, ' F_e ', Fe
      call assert_near(trim(lab), l1, law(Fm, A, B), 1.e-12_wp)
      call assert_near('Omega(F_e) = lambda2', l2, law(Fe, A, B), 1.e-12_wp)
      d = 1.e-5_wp * Fm
      dO = (law(Fm + d, A, B) - law(Fm - d, A, B)) / (2.0_wp * d)
      call assert_near('dOmega/dF = 0 at Fmax (relative)', 0.0_wp, dO * Fm, 1.e-8_wp)
      Om_max = 0.0_wp
      do k = 1, 20000
        Fk = 3.0_wp * Fe * real(k, wp) / 20000.0_wp; Om_max = max(Om_max, law(Fk, A, B))
      end do
      call assert_true('Omega(F) <= lambda1 Omega_c everywhere', Om_max <= l1 * (1.0_wp + 1.e-12_wp))
    end do
  end do

  ! (2) bulk root in a flat metric
  l1 = 2.0_wp; l2 = 0.5_wp; y = uryu_fmax_ratio(l1, l2, 1, 3); Fe = 0.8_wp
  call ab(Fe, y * Fe, l1, l2, A, B); cached_aa = A; cached_bb = B; ctx_exp2re2rho = 1.0_wp
  do k = 1, 6
    sgp = 0.05_wp * k; mu = 0.3_wp
    xlo = 0.0_wp; xhi = (1.0_wp - sgp) / (sgp * sqrt(1.0_wp - mu**2))
    do j = 1, 3
      guess = Omega_c * gfac(j)
      call uryu_omega(0.0_wp, 0.0_wp, 0.0_wp, sgp, mu, guess, omega, found)
      write(lab, '(a,f4.2,a,f4.2)') 'uryu_omega found, s ', sgp, ' guess/Omega_c ', guess / Omega_c
      call assert_true(trim(lab), found)
      call assert_near('f(root) = 0', 0.0_wp, fflat(omega, sgp, mu), 1.e-10_wp)
      ! brute force: the sign change nearest the guess
      nbest = -1; dist = huge(1.0_wp); fprev = fflat(xlo + 1.e-12_wp * xhi, sgp, mu)
      do i = 1, 200000
        x = xlo + (xhi - xlo) * real(i, wp) / 200001.0_wp; fx = fflat(x, sgp, mu)
        if (fx * fprev <= 0.0_wp .and. abs(x - guess) < dist) then
          dist = abs(x - guess); fbest = x; nbest = i
        end if
        fprev = fx
      end do
      call assert_true('brute force finds a root', nbest > 0)
      call assert_near('root nearest the guess', fbest, omega, 1.e-5_wp * xhi)
    end do
  end do
  call uryu_omega(0.0_wp, 0.0_wp, 0.0_wp, 0.3_wp, 1.0_wp, Omega_c, omega, found)
  call assert_near('axis: Omega = Omega_c', Omega_c, omega, 1.e-15_wp)

  call test_summary('test_uryu_fmax')

contains
  subroutine ab(Fe, Fm, l1, l2, A, B)          ! closed forms of AA_h, BB_h (p = 1, q = 3)
    real(wp), intent(in)  :: Fe, Fm, l1, l2
    real(wp), intent(out) :: A, B
    A = ((Fe * Fm) * (l2 * Fe**3 - l1 * Fm**3) / ((l1 - 1.0_wp) * Fe - (l2 - 1.0_wp) * Fm))**0.25_wp
    B = Fe * Fm * (l2 * Fe**3 - l1 * Fm**3) / (l2 * (l1 - 1.0_wp) * Fe**4 - l1 * (l2 - 1.0_wp) * Fm**4)
  end subroutine ab
  real(wp) function law(F, A, B)               ! Omega/Omega_c of law (8), p = 1, q = 3
    real(wp), intent(in) :: F, A, B
    law = (1.0_wp + F / B) / (1.0_wp + (F / A)**4)
  end function law
  real(wp) function fflat(x, sgp, mu)          ! uryu_omega's f(x) with ww = 0, exp(2 re^2 rho) = 1
    real(wp), intent(in) :: x, sgp, mu
    real(wp) :: vfac, Fj
    vfac = sgp**2 * (1.0_wp - mu**2)
    Fj = x * vfac / ((1.0_wp - sgp)**2 - x**2 * vfac)
    fflat = x / Omega_c * (1.0_wp + (Fj / cached_aa)**4) - (1.0_wp + Fj / cached_bb)
  end function fflat
end program test_uryu_fmax
