program test_uryu_inc
  use precision_mod, only: wp
  use test_utils
  use rotation_law_mod, only: uryu_inc_G, uryu_inc_c, uryu_inc_omega
  implicit none

  integer, parameter :: NGL = 40, M_SUB = 10
  real(wp) :: xg(NGL), wg(NGL)
  real(wp), parameter :: ps(4) = [0.5_wp, 1.0_wp, 1.3_wp, 2.0_wp]
  real(wp), parameter :: xs(3) = [0.05_wp, 0.3_wp, 1.0_wp]
  real(wp), parameter :: A2 = 2.7_wp, omc = 0.31_wp
  real(wp) :: p, x, quad, G, lam, F_e, Oe, xstar, vp2, s, mu, omega, guess, sinth
  integer :: i, j, k
  character(len=64) :: lab

  call gauss_legendre(NGL, xg, wg)

  ! G(x) = A^2 Omega_c^2 int_0^x (1+u) u^p du against quadrature. With u = x y^m the
  ! integrand is x^(p+1) m y^(m p + m - 1) (1 + x y^m): a polynomial in y for m p
  ! integer, which NGL-point Gauss-Legendre integrates exactly up to round-off.
  do i = 1, size(ps)
    p = ps(i)
    do j = 1, size(xs)
      x = xs(j)
      quad = 0.0_wp
      do k = 1, NGL
        associate(y => 0.5_wp * (xg(k) + 1.0_wp))
          quad = quad + 0.5_wp * wg(k) * x**(p + 1.0_wp) * real(M_SUB, wp) &
               * y**(real(M_SUB, wp) * p + real(M_SUB, wp) - 1.0_wp) * (1.0_wp + x * y**M_SUB)
        end associate
      end do
      quad = A2 * omc**2 * quad
      G = uryu_inc_G(x, A2, omc, p)
      write(lab, '(A,F4.2,A,F4.2)') "G vs quadrature, p=", p, " x=", x
      call assert_rel_near(trim(lab), quad, G, 1.0e-12_wp)
    end do
  end do
  call assert_near("G(0) = 0", 0.0_wp, uryu_inc_G(0.0_wp, A2, omc, 0.5_wp), 1.0e-300_wp)

  ! Equator identity: with A^2 = F_e/(Omega_e (lambda-1)^p) and Omega_c = Omega_e/lambda,
  ! G(lambda - 1) = F_e Omega_e c(lambda, p).
  F_e = 0.37_wp; Oe = 0.42_wp
  do i = 1, size(ps)
    p = ps(i)
    do j = 1, 3
      lam = 1.0_wp + 0.25_wp * real(j, wp)
      G = uryu_inc_G(lam - 1.0_wp, F_e / (Oe * (lam - 1.0_wp)**p), Oe / lam, p)
      write(lab, '(A,F4.2,A,F4.2)') "equator identity, p=", p, " lambda=", lam
      call assert_rel_near(trim(lab), F_e * Oe * uryu_inc_c(lam, p), G, 1.0e-13_wp)
    end do
  end do

  ! Bulk root, analytic case: flat metric (rho = w = 0), F_j = Omega vp^2/(1 - Omega^2 vp^2)
  ! with vp = s sin(theta)/(1-s). For a chosen root x*, vp^2 = A^2 x*^p/(1 + A^2 x*^p Omega*^2);
  ! the (s, mu) built from it must return x*.
  do i = 1, size(ps)
    p = ps(i)
    do j = 1, size(xs)
      xstar = xs(j)
      do k = 1, 2
        mu = merge(0.0_wp, 0.6_wp, k == 1)
        sinth = sqrt(1.0_wp - mu**2)
        vp2 = A2 * xstar**p / (1.0_wp + A2 * xstar**p * (omc * (1.0_wp + xstar))**2)
        s = sqrt(vp2) / (sinth + sqrt(vp2))
        guess = merge(1.0e-3_wp * xstar, 3.0_wp * xstar, j == 2)
        call uryu_inc_omega(1.0_wp, 0.0_wp, 0.0_wp, s, mu, omc, A2, p, guess, omega)
        ! x* is the physical root only if it is the first sign change, f'(x*) < 0;
        ! otherwise the solver must return the first root, below x*.
        if (fflat(xstar + 1.0e-7_wp, vp2, p) - fflat(xstar - 1.0e-7_wp, vp2, p) < 0.0_wp) then
          write(lab, '(A,F4.2,A,F4.2,A,F3.1)') "bulk root, p=", p, " x*=", xstar, " mu=", mu
          call assert_rel_near(trim(lab), xstar, omega / omc - 1.0_wp, 1.0e-10_wp)
        else
          write(lab, '(A,F4.2,A,F4.2,A,F3.1)') "first root below a second-root x*, p=", p, " x*=", xstar, " mu=", mu
          call assert_true(trim(lab), omega / omc - 1.0_wp < xstar .and. &
                           abs(fflat(omega / omc - 1.0_wp, vp2, p)) < 1.0e-10_wp * A2 * omc)
        end if
      end do
    end do
  end do
  call uryu_inc_omega(1.0_wp, 0.0_wp, 0.0_wp, 0.3_wp, 1.0_wp, omc, A2, 0.5_wp, 0.1_wp, omega)
  call assert_near("bulk root on the axis = Omega_c", omc, omega, 1.0e-15_wp)

  call test_summary("uryu_inc")

contains

  real(wp) function fflat(x, vp2, p)
    real(wp), intent(in) :: x, vp2, p
    real(wp) :: o
    o = omc * (1.0_wp + x)
    fflat = o * vp2 / (1.0_wp - o**2 * vp2) - A2 * o * x**p
  end function fflat

  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(wp), intent(out) :: x(n), w(n)
    real(wp), parameter :: pi = 3.14159265358979323846_wp
    real(wp) :: z, z1, p1, p2, p3, pp
    integer :: i, jj, it
    do i = 1, (n + 1) / 2
      z = cos(pi * (real(i, wp) - 0.25_wp) / (real(n, wp) + 0.5_wp))
      do it = 1, 100
        p1 = 1.0_wp; p2 = 0.0_wp
        do jj = 1, n
          p3 = p2; p2 = p1
          p1 = ((2.0_wp * jj - 1.0_wp) * z * p2 - (jj - 1.0_wp) * p3) / real(jj, wp)
        end do
        pp = n * (z * p1 - p2) / (z * z - 1.0_wp)
        z1 = z; z = z1 - p1 / pp
        if (abs(z - z1) < 1.0e-15_wp) exit
      end do
      x(i) = -z; x(n + 1 - i) = z
      w(i) = 2.0_wp / ((1.0_wp - z * z) * pp * pp); w(n + 1 - i) = w(i)
    end do
  end subroutine gauss_legendre

end program test_uryu_inc
