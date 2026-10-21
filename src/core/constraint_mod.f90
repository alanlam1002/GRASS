module constrain_mod
  use precision_mod, only: wp
  use para_mod, only: SDIV, MDIV, DS, DM, r_e, s_gp, mu, rho, gama, alpha, ww, &
                      energy, pressure, sphi, mphi_r, B_coup, pi, has_scalar, velocity_sq, &
                      angular_collocation, COLLOCATION_UNI, w_mu, output, enthalpy
  use ope_eq_mod, only: laplacian_operator, gradient_vector
  implicit none
  private

  public :: hamiltonian

  type(laplacian_operator) :: op   ! allocated once on first hamiltonian call

contains

  subroutine hamiltonian(hamL2)
    real(wp), intent(out) :: hamL2
    real(wp), dimension(SDIV,MDIV) :: ham
    real(wp), dimension(SDIV,MDIV) :: logPsi4, acoup4
    real(wp), dimension(SDIV,MDIV) :: ricci, ricci_lap, ricci_scal
    real(wp), dimension(SDIV,MDIV) :: psi4, gutt, gurr, twist, KK, rhoH, dphidphi, Vphi
    real(wp), dimension(SDIV) :: r_phys
    real(wp), dimension(SDIV,MDIV) :: r2_2d, m1_2d
    real(wp), dimension(MDIV) :: mu_weights
    integer :: unit, ios, s, m
    type(gradient_vector) :: grad_logPsi4, grad_ww, grad_sphi, grad_alp

    r_phys = r_e * s_gp / max(1.e-30_wp, 1.0_wp - s_gp)
    r2_2d  = spread(r_phys**2, dim=2, ncopies=MDIV)
    m1_2d  = spread(1.0_wp - mu**2, dim=1, ncopies=SDIV)
    psi4   = exp(gama - rho)
    logPsi4= gama - rho
    acoup4 = exp(-sphi**2 * B_coup)
    gurr   = exp(-2.0_wp * alpha)
    twist  = psi4 * r2_2d * m1_2d * ww
    gutt   = -exp(gama + rho)
    rhoH   = (energy + pressure) / (1.0_wp - velocity_sq) - pressure
    
    call op%init()
    grad_logPsi4 = op%grad(logPsi4)
    grad_ww      = op%grad(ww)
    grad_sphi    = op%grad(sphi)
    grad_alp     = op%grad(alpha)

    ricci_lap = op%lap2(alpha) + 0.5_wp * op%laplacian(logPsi4) !op%laplacian(alpha + 0.5_wp * logPsi4)
    ricci_scal = 0.25_wp * op%scal(grad_logPsi4, grad_logPsi4) + 0.5_wp * op%divr(grad_logPsi4%r)
    ricci = -2.0_wp * gurr * (ricci_lap + ricci_scal )

    if (has_scalar .and. abs(B_coup) > 1.e-30_wp) then
      dphidphi = op%scal(grad_sphi, grad_sphi) * gurr
      Vphi     = 0.5_wp * mphi_r * sphi**2 !/ B_coup
    else
      dphidphi = 0.0_wp
      Vphi     = 0.0_wp
    end if

    !KK = 0.5_wp * gurr / psi4 / max(lapsesq, 1.e-30_wp) * op%scal(grad_ww, grad_ww) &
    !   / max(r2_2d, 1.e-30_wp) / max(m1_2d, 1.e-30_wp)
    KK = 0.5_wp * gurr * exp(-2.0_wp * rho) * r2_2d * m1_2d * op%scal(grad_ww, grad_ww)
    ham = ricci + KK - 16.0_wp * pi * rhoH * acoup4 - 2.0_wp * (dphidphi + 2.0_wp * Vphi)

    if (angular_collocation /= COLLOCATION_UNI) then
      mu_weights = w_mu
    else
      mu_weights = DM
      if (MDIV > 1) then
        mu_weights(1) = 0.5_wp * DM
        mu_weights(MDIV) = 0.5_wp * DM
      end if
    end if

    hamL2 = 4.0_wp * pi * r_e**3 * DS * &
            sum(ham**2 * spread(s_gp**2 / max(1.e-30_wp, 1.0_wp - s_gp)**4, dim=2, ncopies=MDIV) * &
                spread(mu_weights, dim=1, ncopies=SDIV))
    !write(*,'(16es8.1)') ricci(:,1)
    !write(*,*) " "
    !write(*,'(16es8.1)') 16.0_wp * pi * rhoH(:,1) * acoup4(:,1)
    !write(*,*) " "
    !write(*,'(16es8.1)') Vphi(:,1)

    if (.true.) then
      open(newunit=unit, file="./Cont/hamiltonian.dat", status="replace", action="write", iostat=ios)
      if (ios == 0) then
        write(unit,"(2(i0,2X))") SDIV, MDIV
        do s = 1, SDIV
          do m = 1, MDIV
            write(unit, "(13es24.15)") s_gp(s), mu(m), r_phys(s), ham(s,m), enthalpy(s,m), &
                                        ww(s,m), velocity_sq(s,m), rho(s,m), gama(s,m), alpha(s,m), &
                                        rhoH(s,m), ricci(s,m), KK(s,m)
          end do
        end do
        close(unit)
      else
        write(*,*) "hamiltonian: failed to open ./Cont/hamiltonian.dat"
      end if
    end if
  end subroutine hamiltonian

end module constrain_mod
