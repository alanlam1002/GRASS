program grass
  use iso_fortran_env, only: error_unit
  use eos_mod, only: loadEos
  use grid_mod, only: make_grid, GridTrig
  use shoot_mod, only: shoot_v2
  use toolkit_mod, only: debug_mod_bessel
  use para_mod, only: wp, initialize_theory, angular_collocation, run_task, &
                      has_scalar, mphi_goal, scalarton, B_goal, &
                      solver_type, A_diff, lambda1, lambda2, inc_p, inc_lambda, &
                      s_gp, SDIV, MDIV, s_pwr, e_surface, C, KSCALE
  use constrain_mod, only: hamiltonian
  use starting_model_mod, only: initialize_starting_model
  use MRcurve_mod, only: MRcurve
  implicit none
  real(wp) :: hamL2
  character(len=10) :: mu_grid_name

  call initialize_theory()
  
  print *, " "
  print *, " "
  if (has_scalar) then
    write(*,"(A32,es12.6,A11,es12.6,A2)") &
    "Theory: Scalar-Tensor  ( mphi = ", mphi_goal*scalarton, "eV,  B = ", B_goal, " )"
  else
    write(*,*) "Theory: General Relativity"
  end if
  write(*,*) "Relaxation scheme: Anderson"
  block
    character(len=10), parameter :: mu_names(3) = ["Uniform   ", "Legendre  ", "Chebyshev "]
  write(*,*) "          mu-grid: ", trim(mu_names(angular_collocation))
  end block

  call loadEos

  write(*,'(1X,A)') repeat('-', 36)
  select case (trim(adjustl(solver_type)))
  case ("const_j")
    write(*,'(1X,A,F5.3,A)') "Rotation Law: Constant-J   (A^-1 = ", 1._wp / A_diff, ")"
  case ("uryu")
    write(*,'(1X,A,F5.3,A,F5.3,A)') "Rotation Law: Uryu   (lambda1 = ", lambda1, "  lambda2 = ", lambda2, ")"
  case ("uryu_inc")
    write(*,'(1X,A,F6.3,A,F6.3,A)') "Rotation Law: Uryu increasing   (p = ", inc_p, "  lambda = ", inc_lambda, ")"
  case default
    write(*,'(1X,A)') "Rotation Law: Uniform"
  end select
  write(*,'(1X,A)') repeat('-', 36)

  call make_grid
  call GridTrig

  write(*,fmt="(A18,es18.9)",Advance='NO') " Outer boundary:", (s_gp(SDIV-1) / ( 1._wp - s_gp(SDIV-1) ))**s_pwr
  write(*,"(A18,i0,A,i0)") "Resolution: ", SDIV, " x ", MDIV
  write(*,"(A18,es18.9)") "Surface eps:", e_surface / (C*C*KSCALE)
  
  !call debug_mod_bessel; stop "debug_mod_bessel output written"
  
  select case (run_task)
  case(1)
    call shoot_v2
  case(2)
    call MRcurve
  case(3)
    call initialize_starting_model()
  case default
    write(error_unit,*) "Unknown task."; stop 3221
  end select

  call hamiltonian(hamL2)
  write(*,"(A18,es27.16)") "Ham L2:", hamL2
  write(*,*) " "
  write(*,*) "Completed!"

end program grass
