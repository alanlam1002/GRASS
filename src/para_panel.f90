module para_mod
  use precision_mod, only: wp
  implicit none
  ! -- Theory selection ------------------------------------------------------
  integer, parameter :: THEORY_GR = 0
  integer, parameter :: THEORY_ST = 1
  integer :: active_theory = THEORY_GR

  ! -- Model-family selection -------------------------------------------------
  integer, parameter :: MODEL_NS = 1
  integer, parameter :: MODEL_BH_TOROID = 2
  integer :: model_family = MODEL_NS

  ! hybrid / anderson
  ! -- Running option --------------------------------------------------------
  integer, parameter :: MODE_REGRID  = 1, MODE_DEFAULT = 2

  integer :: run_mode = MODE_REGRID  ! isolate mass vs shape at default e_center

  ! -- Running mode ----------------------------------------------------------
  integer, parameter :: shoot = 1, MRbuild = 2, OneModel = 3
  integer :: run_task = OneModel

  ! -- Rotation configuration ------------------------------------------------
  ! uniform / const_j / uryu
  character(len=20) :: solver_type = "uryu"
  integer, parameter :: COLLOCATION_UNI = 1, COLLOCATION_LEG = 2, COLLOCATION_CHEB = 3
  integer :: angular_collocation = COLLOCATION_LEG

  ! -- Solver state ----------------------------------------------------------
  logical :: output = .false.
  logical :: timing = .false.
  integer, parameter :: SHOOT_FIX1_HC = 1, SHOOT_FIX1_RP = 2, SHOOT_2D = 3
  integer :: shooting = SHOOT_FIX1_HC
  character(len=20) :: FIX1 = "Mb_goal"
  character(len=20) :: FIX2 = "chi_goal"

  ! -- Resolutions -----------------------------------------------------------
  integer, parameter :: res  = 400
  integer, parameter :: s_pwr = 1
  integer :: SDIV = 2 * res + 1
  integer :: MDIV = 41

  ! -- Target quantities -----------------------------------------------------
  character(len=128) :: eos_file = "DD2_hot_equal_pinned"
  real(wp) :: M_goal   = 1.2_wp
  real(wp) :: Mb_goal  = 2._wp
  real(wp) :: J_goal   = 1.6_wp
  real(wp) :: chi_goal = 0.1_wp
  real(wp) :: omc_goal = 30.0_wp

  real(wp) :: B_goal   = 3.5e1_wp
  real(wp) :: mphi_goal = 0.2_wp

  ! -- Rotation-law parameters (KEH, Uryu enabled) --------------------------
  real(wp) :: A_diff  = 0.5_wp
  real(wp) :: lambda1 = 1.28_wp
  real(wp) :: lambda2 = 0.72_wp
  integer :: uryu_p  = 1
  integer :: uryu_q  = 3

  real(wp) :: parA = 1.0_wp
  real(wp) :: parB = 1.0_wp

  real(wp) :: cofA = 30.0_wp
  real(wp) :: cofB = 0.0_wp
  real(wp) :: cofp = 0.8_wp
  real(wp) :: cofq = 0.9_wp

  ! -- Equation of state -----------------------------------------------------
  logical :: phase_transition = .false.
  integer :: num_tab = 0
  integer :: n_PT = 0
  integer, allocatable :: p_at_PT(:)
  real(wp), allocatable :: log_p(:), log_e(:), log_h(:), log_n0(:)
  real(wp) :: p_center = 0.0_wp
  real(wp) :: h_center = 0.0_wp
  real(wp) :: e_center = 0.0_wp
  real(wp) :: e_center_default = .8e15_wp
  real(wp) :: enthalpy_min = 0.0_wp

  ! -- Grid configuration ----------------------------------------------------
  integer, parameter :: LMAX = 10
  real(wp) :: SMAX  = 1.0_wp - (1.e-1_wp)**(9.0_wp / real(s_pwr, wp))
  integer, parameter :: RDIV = 1800

  real(wp) :: DS = 0.0_wp
  real(wp) :: DM = 0.0_wp
  real(wp), parameter :: S_E = 0.5_wp

  real(wp), allocatable :: s_gp(:), mu(:), sin_theta(:)
  real(wp), allocatable :: D_mu(:,:), D_mu_t(:,:), D2_mu(:,:), w_mu(:)

  ! -- Disk helper quantities (kept for compatibility) ----------------------
  logical :: disk_present = .false.
  real(wp) :: edge_in = 800.0_wp
  real(wp) :: s_inner = 0.5_wp
  real(wp) :: j_disk  = 4.5_wp
  real(wp) :: p_max_disk = 0.0_wp
  real(wp) :: h_max_disk = 0.0_wp
  integer :: i_isco_p = 1
  integer :: i_isco_m = 1

  ! Fluid
  real(wp), allocatable :: pressure(:,:), enthalpy(:,:), velocity_sq(:,:), &
                          energy(:,:), omg(:,:), F_j(:,:)
  real(wp), allocatable :: v_plus(:), v_minus(:), V_rr_p(:), V_rr_m(:), sound_speed(:)

  ! Metric
  real(wp), allocatable :: gama(:,:), rho(:,:), ww(:,:), alpha(:,:), sphi(:,:)

  ! Scalar field
  logical :: has_scalar = .false.
  real(wp) :: B_coup = 0.0_wp
  real(wp) :: mphi_r = 0.0_wp
  real(wp) :: sphi_c     = 0.0_wp
  real(wp) :: sphi_m     = 0.0_wp
  real(wp) :: r_sphi_max = 0.0_wp  ! physical equatorial radius at max(sphi)
  real(wp) :: donut      = 0.0_wp

  real(wp) :: B_burn_init = 16.0_wp
  real(wp) :: mphi_burn_seed = 0.1_wp
  integer, parameter :: scalar_burn_max_iter = 200
  real(wp), parameter :: MPHI_BURN_THRESHOLD = 0.1_wp

  ! Bulk properties
  real(wp) :: Omega_c = 0.0_wp
  real(wp) :: Omega_e = 0.0_wp
  real(wp) :: Omega_K = 0.0_wp
  real(wp) :: r_e     = 1.0_wp
  real(wp) :: r_ratio = 0.85_wp  ! lambda strengthening at r_ratio=0.85
  real(wp) :: r_circ  = 0.0_wp
  real(wp) :: ang_mom = 0.0_wp
  real(wp) :: mass    = 0.0_wp
  real(wp) :: mass_0  = 0.0_wp
  real(wp) :: mass_s  = 0.0_wp  ! Baryon mass in scalar-active regions (sphi > 1e-3)
  real(wp) :: mass_p  = 0.0_wp
  real(wp) :: chi     = 0.0_wp
  real(wp) :: T_kin   = 0.0_wp
  real(wp) :: I_inertia = 0.0_wp
  real(wp) :: Love2   = 0.0_wp
  real(wp) :: Fmax_h  = 0.0_wp
  real(wp) :: F_equator_h = 0.0_wp
  ! Multipole information
  real(wp) :: M2 = 0.0_wp
  real(wp) :: M4 = 0.0_wp
  real(wp) :: M6 = 0.0_wp
  real(wp) :: S3 = 0.0_wp
  real(wp) :: S5 = 0.0_wp

  integer :: n_of_relaxation_steps = 0

  ! Green's functions (Legendre weights — precomputed on grid)
  real(wp), allocatable :: P_2n(:,:), P1_2n_1(:,:), sin_2n_1_theta(:,:)

  ! Timing
  real(wp) :: start = 0.0_wp
  real(wp) :: finish = 0.0_wp

  ! -- Physical constants ----------------------------------------------------
  real(wp), parameter :: C    = 2.99792458e10_wp
  real(wp), parameter :: G    = 6.67408e-8_wp
  real(wp), parameter :: MSUN = 1.98847e33_wp
  real(wp), parameter :: MB   = 1.6749286e-24_wp
  real(wp), parameter :: PI = acos(-1.0_wp)

  real(wp), parameter :: ACCURACY  = 1.e-5_wp
  real(wp), parameter :: TOV_RMIN  = 1.e-15_wp
  real(wp), parameter :: KAPPA     = 1.e-15_wp * C**2 / G
  real(wp), parameter :: KSCALE    = KAPPA * G / C**4
  ! DD2_hot_equal(_pinned) is a hot merger-remnant profile, not a cold NS
  ! crust EOS -- it doesn't extend down anywhere near the standard
  ! iron-density vacuum surface (7.8 g/cm^3, 1.01e8 dyn/cm^2) these two
  ! defaults assume. sphere_mod's TOV integrator has no iteration cap, so
  ! with the standard defaults it runs (near-)forever trying to reach a
  ! pressure the table was never built to cover, extrapolating off the
  ! table's first cell the whole way. Re-pointed at this table's own
  ! low-density floor instead (row 5 of eos/DD2_hot_equal_pinned.dat,
  ! e=2.404e3 g/cm^3, p=7.253e21 dyn/cm^2 -- a few rows in from the table's
  ! edge, not the very first row, so this is reached by direct
  ! interpolation, not by extrapolating past the table).
  ! Test: does MPA1 also degrade under rotation if ITS surface is raised to
  ! match DD2_hot_equal_pinned's absolute cutoff level? Row 130 of MPA1.dat
  ! (e=1.9967e3, p=1.0345e18), directly tabulated -- no extrapolation.
  ! Now within reach of the standard near-vacuum convention (~7.8 g/cm^3)
  ! thanks to the polytropic low-density extension -- row 52 of the
  ! extended DD2_hot_equal_pinned.dat, directly tabulated, no extrapolation
  ! needed. Gamma for this extension is now FIXED at 1.58425 (the Read
  ! et al. 2009 / MPA1-floor literature low-density-crust value) rather
  ! than matched to DD2's own softer (Gamma~0.95) edge slope -- the old,
  ! sub-linear Gamma<1 tail gave low-mass stars an unboundedly-growing
  ! radius (Newtonian Lane-Emden n<0 regime; several low-e_center TOV
  ! solves failed to converge at all). P_SURFACE moved accordingly: the
  ! stiffer tail drops pressure much faster with density, so P at the same
  ! e is ~27x smaller than under the old Gamma.
  real(wp), parameter :: E_SURFACE = 8.95e0_wp * C**2 * KSCALE
  real(wp), parameter :: P_SURFACE = 1.307030e18_wp * KSCALE
  real(wp), parameter :: RHO_UNI   = 7.4259154861063358e-19_wp
  real(wp), parameter :: PRS_UNI = rho_uni/(C*1.e5_wp)**2
  real(wp), parameter :: F_UNI     = 2.029739818539300e5_wp
  real(wp), parameter :: HBAR      = 6.582119569e-16_wp
  real(wp), parameter :: L_UNI     = 1.4769994423016508_wp
  real(wp), parameter :: N_SAT     = 2.7e14_wp
  real(wp), parameter :: SCALARTON = HBAR * C / l_uni / 1.e5_wp
contains

  pure function to_lower_str(str) result(out)
    character(*), intent(in) :: str
    character(len(str)) :: out
    integer :: i, code

    out = str
    do i = 1, len(str)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end function to_lower_str

  subroutine initialize_theory()
    select case (active_theory)
    case (THEORY_GR)
      call apply_gr_defaults()
    case (THEORY_ST)
      call apply_st_defaults()
    case default
      error stop "initialize_theory: unknown theory mode"
    end select

    solver_type = trim(to_lower_str(adjustl(solver_type)))

    DS   = SMAX / (real(SDIV, wp) - 1.0_wp)
    DM   = 1.0_wp  / (real(MDIV, wp) - 1.0_wp)
    s_inner = edge_in**(1.0_wp / real(s_pwr, wp)) / (edge_in**(1.0_wp / real(s_pwr, wp)) + 1.0_wp)

    call apply_runtime_overrides()

    call allocate_fields()
  end subroutine initialize_theory

  subroutine apply_runtime_overrides()
    character(len=64) :: buf
    integer :: stat, length, ios
    real(wp) :: rval

    call get_environment_variable("GRASS_TIMING", buf, length, stat)
    if (stat == 0 .and. length > 0) then
      select case (trim(adjustl(to_lower_str(buf))))
      case ("1", "true", "yes", "on", "t")
        timing = .true.
      case ("0", "false", "no", "off", "f")
        timing = .false.
      end select
    end if

    call get_environment_variable("GRASS_E_CENTER", buf, length, stat)
    if (stat == 0 .and. length > 0) then
      read(buf, *, iostat=ios) rval
      if (ios == 0) e_center_default = rval
    end if
  end subroutine apply_runtime_overrides

  subroutine apply_gr_defaults()
    has_scalar = .false.
    B_goal  = 0.0_wp
    mphi_goal = 0.0_wp
    B_coup = 0.0_wp
    mphi_r = 0.0_wp
    sphi_c = 0.0_wp
    sphi_m = 0.0_wp
  end subroutine apply_gr_defaults

  subroutine apply_st_defaults()
    has_scalar = .true.
  end subroutine apply_st_defaults

  subroutine allocate_fields()
    call deallocate_fields()

    allocate(s_gp(SDIV), source=0.0_wp)
    allocate(mu(MDIV), source=0.0_wp)
    allocate(sin_theta(MDIV), source=0.0_wp)
    allocate(D_mu(MDIV, MDIV), source=0.0_wp)
    allocate(D_mu_t(MDIV, MDIV), source=0.0_wp)
    allocate(D2_mu(MDIV, MDIV), source=0.0_wp)
    allocate(w_mu(MDIV), source=0.0_wp)

    allocate(pressure(SDIV,MDIV), source=0.0_wp)
    allocate(enthalpy(SDIV,MDIV), source=0.0_wp)
    allocate(velocity_sq(SDIV,MDIV), source=0.0_wp)
    allocate(energy(SDIV,MDIV), source=0.0_wp)
    allocate(omg(SDIV,MDIV), source=0.0_wp)
    allocate(F_j(SDIV,MDIV), source=0.0_wp)
    allocate(v_plus(SDIV), source=0.0_wp)
    allocate(v_minus(SDIV), source=0.0_wp)
    allocate(V_rr_p(SDIV), source=0.0_wp)
    allocate(V_rr_m(SDIV), source=0.0_wp)
    allocate(sound_speed(SDIV), source=0.0_wp)

    allocate(gama(SDIV,MDIV), source=0.0_wp)
    allocate(rho(SDIV,MDIV), source=0.0_wp)
    allocate(ww(SDIV,MDIV), source=0.0_wp)
    allocate(alpha(SDIV,MDIV), source=0.0_wp)
    allocate(sphi(SDIV,MDIV), source=0.0_wp)

    allocate(P_2n(MDIV, LMAX+1), source=0.0_wp)
    allocate(P1_2n_1(MDIV, LMAX+1), source=0.0_wp)
    allocate(sin_2n_1_theta(MDIV, LMAX), source=0.0_wp)
  end subroutine allocate_fields

  subroutine deallocate_fields()
    if (allocated(s_gp))             deallocate(s_gp)
    if (allocated(mu))               deallocate(mu)
    if (allocated(sin_theta))        deallocate(sin_theta)
    if (allocated(D_mu))             deallocate(D_mu)
    if (allocated(D_mu_t))           deallocate(D_mu_t)
    if (allocated(D2_mu))            deallocate(D2_mu)
    if (allocated(w_mu))             deallocate(w_mu)
    if (allocated(pressure))         deallocate(pressure)
    if (allocated(enthalpy))         deallocate(enthalpy)
    if (allocated(velocity_sq))      deallocate(velocity_sq)
    if (allocated(energy))           deallocate(energy)
    if (allocated(omg))              deallocate(omg)
    if (allocated(F_j))              deallocate(F_j)
    if (allocated(v_plus))           deallocate(v_plus)
    if (allocated(v_minus))          deallocate(v_minus)
    if (allocated(V_rr_p))           deallocate(V_rr_p)
    if (allocated(V_rr_m))           deallocate(V_rr_m)
    if (allocated(gama))             deallocate(gama)
    if (allocated(rho))              deallocate(rho)
    if (allocated(ww))               deallocate(ww)
    if (allocated(alpha))            deallocate(alpha)
    if (allocated(sphi))             deallocate(sphi)
    if (allocated(P_2n))             deallocate(P_2n)
    if (allocated(P1_2n_1))          deallocate(P1_2n_1)
    if (allocated(sin_2n_1_theta))   deallocate(sin_2n_1_theta)
  end subroutine deallocate_fields

end module para_mod
