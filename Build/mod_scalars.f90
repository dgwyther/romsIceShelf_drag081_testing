      MODULE mod_scalars
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        USE mod_param
!
        implicit none
!
!-----------------------------------------------------------------------
! Multiple grid structure.
!-----------------------------------------------------------------------
!
!    Fstate        Logical switches to control computations of the
!                    Forcing Singular Vectors or Stochastic Optimals.
!    Lstate        Logical switches to control computations of the
!                    model state.
!    Cs_r          Set of S-curves used to stretch the vertical grid
!                    that follows the bathymetry at vertical RHO-points.
!    Cs_w          Set of S-curves used to stretch the vertical grid
!                    that follows the bathymetry at vertical W-points.
!    sc_r          S-coordinate independent variable, [-1 < sc < 0] at
!                    vertical RHO-points.
!    sc_w          S-coordinate independent variable, [-1 < sc < 0] at
!                    vertical W-points.
!
        TYPE T_SCALARS
          logical, pointer :: Fstate(:)
          logical, pointer :: Lstate(:)
          real(r8), pointer :: Cs_r(:)
          real(r8), pointer :: Cs_w(:)
          real(r8), pointer :: sc_r(:)
          real(r8), pointer :: sc_w(:)
        END TYPE T_SCALARS
        TYPE (T_SCALARS), allocatable :: SCALARS(:)
!
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!
        integer :: itemp              ! Potential temperature
        integer :: isalt              ! Salinity
!
!-----------------------------------------------------------------------
!  Time stepping indices, variables, and clocks.
!-----------------------------------------------------------------------
!
!    indx1         2D timestep rolling counter.
!    iic           Timestep counter for 3D primitive equations.
!    iif           Timestep counter for 2D primitive equations.
!    ndtfast       Number of barotropic timesteps between each
!                    baroclinic timestep.
!    nfast         Number of barotropic timesteps needed to compute
!                    time-averaged barotropic variables centered at
!                    time level n+1.
!    dt            Size baroclinic timestep (s).
!    dtfast        Size barotropic timestep (s).
!    run_time      Total run time for all nested grids (s).
!    tdays         Model time clock (days).
!    time          Model time clock (s).
!    time_code     Model time clock (string, Day HH:MM:SS)
!    AVGtime       Model time clock for averages output (s).
!    DIAtime       Model time clock for diagnostics output (s).
!    IMPtime       Impulse forcing time (s) to process.
!    ObsTime       Observation time (s) to process.
!    FrcTime       Adjoint or tangent linear Impulse forcing time (s).
!    r_date        Model initialization reference date (vector):
!                    r_date(1) => reference date (YYYYMMDD.dd).
!                    r_date(2) => year.
!                    r_date(3) => year day.
!                    r_date(4) => month.
!                    r_date(5) => day.
!                    r_date(6) => hour.
!                    r_date(7) => minute.
!                    r_date(8) => second.
!    dstart        Time stamp assigned to model initialization (usually
!                    a Calendar day, like modified Julian Day).
!    tide_start    Reference time for tidal forcing (days).
!    time_ref      Reference time for "units" attribute (YYYYMMDD.dd).
!    r_text        Model initialization reference text (string).
!
        logical, allocatable :: PerfectRST(:)
        logical, allocatable :: PREDICTOR_2D_STEP(:)
        integer, allocatable :: indx1(:)
        integer, allocatable :: iic(:)
        integer, allocatable :: iif(:)
        integer, allocatable :: ndtfast(:)
        integer, allocatable :: nfast(:)
        real(r8), allocatable :: tdays(:)                ! days
        real(r8), allocatable :: time(:)                 ! seconds
        real(r8), allocatable :: dt(:)                   ! seconds
        real(r8), allocatable :: dtfast(:)               ! seconds
        real(r8), allocatable :: TimeEnd(:)              ! seconds
        real(r8), allocatable :: AVGtime(:)              ! seconds
        real(r8), allocatable :: DIAtime(:)              ! seconds
        real(r8), allocatable :: IMPtime(:)              ! seconds
        real(r8), allocatable :: ObsTime(:)              ! seconds
        real(r8), allocatable :: FrcTime(:)              ! seconds
        real(r8), dimension(8) :: r_date
        real(r8) :: dstart = 0.0_r8                      ! days
        real(r8) :: run_time = 0.0_r8                    ! seconds
        real(r8) :: tide_start = 0.0_r8                  ! days
        real(r8) :: time_ref = 0.0_r8                    ! YYYYMMDD.dd
        character (len=14), allocatable :: time_code(:)  ! DD HH:MM:SS
        character (len=19) :: r_text
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  Fields.  The power-law shape filters are given by:
!
!     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi
!
!  Possible settings of parameters to yield the second-order accuracy:
!
!     Falpha  Fbeta      Fgamma
!     ------------------------------
!      2.0     1.0    0.1181  0.169     The problem here is setting
!      2.0     2.0    0.1576  0.234     Fgamma. Its value here is
!      2.0     3.0    0.1772  0.266     understood as the MAXIMUM
!      2.0     4.0    0.1892  0.284     allowed. It is computed using
!      2.0     5.0    0.1976  0.296     a Newton iteration scheme.
!      2.0     6.0    0.2039  0.304
!      2.0     8.0    0.2129  0.314
!
!  NOTE: Theoretical values of Fgamma presented in the table above are
!  derived assuming "exact" barotropic mode stepping. Consequently, it
!  does not account for effects caused by Forward-Euler (FE) startup
!  of the barotropic mode at every 3D time step.  As the result, the
!  code may become unstable if the theoretical value of Fgamma is used
!  when mode splitting ratio "ndtfast" is small, thus yielding non-
!  negligible start up effects.  To compensate this, the accepted
!  value of Fgamma is reduced relatively to theoretical one, depending
!  on splitting ratio "ndtfast".  This measure is empirical. It is
!  shown to work with setting of "ndtfast" as low as 15, which is
!  more robust that the Hamming Window the squared cosine weights
!  options in "set_weights".
!
        real(r8) :: Falpha = 2.0_r8
        real(r8) :: Fbeta  = 4.0_r8
        real(r8) :: Fgamma = 0.284_r8
!
!  Total number timesteps in current run. In 3D configurations, "ntimes"
!  is the total of baroclinic timesteps. In 2D configuration, "ntimes"
!  is the total of barotropic timesteps.
!
        integer, allocatable :: ntimes(:)
!
!  Step counter for current execution time-window.
!
        integer, allocatable :: step_counter(:)
!
!  Number of time interval divisions for Stochastic Optimals.  It must
!  a multiple of "ntimes".
!
        integer :: Nintervals = 1
!
!  Starting, current, and ending ensemble run parameters.
!
        integer :: ERstr = 1                    ! Starting value
        integer :: ERend = 1                    ! Ending value
        integer :: Ninner = 1                   ! number of inner loops
        integer :: Nouter = 1                   ! number of outer loops
        integer :: Nrun = 1                     ! Current counter
        integer :: inner = 0                    ! inner loop counter
        integer :: outer = 0                    ! outer loop counter
!
!  First, starting, and ending timestepping parameters
!
        integer, allocatable :: ntfirst(:)      ! Forward-Euler step
        integer, allocatable :: ntstart(:)      ! Start step
        integer, allocatable :: ntend(:)        ! End step
!
!  Adjoint model or tangent linear model impulse forcing time record
!  counter and number of records available.
!
        integer, allocatable :: FrcRec(:)
        integer, allocatable :: NrecFrc(:)
!
!-----------------------------------------------------------------------
!  Control switches.
!-----------------------------------------------------------------------
!
!  Switch to proccess nudging coefficients.
!
        logical, allocatable :: NudgingCoeff(:)
!
!  Switch to proccess input boundary data.
!
        logical, allocatable :: ObcData(:)
!
!  These switches are designed to control computational options within
!  nested and/or multiple connected grids.  They are .TRUE. by default.
!  They can turned off for a particular grind in input scripts.
!
        logical, allocatable :: Lbiology(:)
        logical, allocatable :: Lfloats(:)
        logical, allocatable :: Lsediment(:)
        logical, allocatable :: Lstations(:)
!
!-----------------------------------------------------------------------
!  Physical constants.
!-----------------------------------------------------------------------
!
!    Cp            Specific heat for seawater (Joules/Kg/degC).
!    Csolar        Solar irradiantion constant (W/m2).
!    Eradius       Earth equatorial radius (m).
!    Infinity      Value resulting when dividing by zero.
!    StefBo        Stefan-Boltzmann constant (W/m2/K4).
!    emmiss        Infrared emmissivity.
!    g             Acceleration due to gravity (m/s2).
!    gorho0        gravity divided by mean density anomaly.
!    rhow          fresh water density (kg/m3).
!    vonKar        von Karman constant.
!
        real(r8) :: Cp = 3985.0_r8              ! Joules/kg/degC
        real(r8) :: Csolar = 1353.0_r8          ! 1360-1380 W/m2
        real(r8) :: Infinity                    ! Infinity = 1.0/0.0
        real(r8) :: Eradius = 6371315.0_r8      ! m
        real(r8) :: StefBo = 5.67E-8_r8         ! Watts/m2/K4
        real(r8) :: emmiss = 0.97_r8            ! non_dimensional
        real(r8) :: rhow = 1000.0_r8            ! kg/m3
        real(r8) :: g = 9.81_r8                 ! m/s2
        real(r8) :: gorho0                      ! m4/s2/kg
        real(r8) :: vonKar = 0.41_r8            ! non-dimensional
!
!-----------------------------------------------------------------------
!  Various model parameters.  Some of these parameters are overwritten
!  with the values provided from model standard input script.
!-----------------------------------------------------------------------
!
!  Switch for spherical grid (lon,lat) configurations.
!
        logical :: spherical = .FALSE.
!
!  Switch to compute the grid stiffness.
!
        logical :: Lstiffness = .TRUE.
!
!  Composite grid a refined grids switches. They are .FALSE. by default.
!
        logical, allocatable :: CompositeGrid(:,:)
        logical, allocatable :: RefinedGrid(:)
!
!  Refinement grid scale factor from donor grid.
!
        integer, allocatable :: RefineScale(:)
!
!  Switch to extract donor grid (coarse) data at the refinement grid
!  contact point locations. The coarse data is extracted at the first
!  sub-refined time step.  Recall that the finer grid time-step is
!  smaller than the coarser grid by a factor of RefineScale(:). This
!  switch is only relevant during refinement nesting.
!
        logical, allocatable :: GetDonorData(:)
!
!  Periodic boundary swiches for distributed-memory exchanges.
!
        logical, allocatable :: EWperiodic(:)
        logical, allocatable :: NSperiodic(:)
!
!  Lateral open boundary edges volume conservation switches.
!
        logical, allocatable :: VolCons(:,:)
!
!  Switch indicating which tracer field needs to be processed for
!  climatology.
!
        logical, allocatable :: LtracerCLM(:,:)
!
!  Switch indicating which tracer field needs to be processed for
!  Sources/Sinks terms.
!
        logical, allocatable :: LtracerSrc(:,:)
!
!  Execution termination flag.
!
!    exit_flag = 0   No error
!    exit_flag = 1   Blows up
!    exit_flag = 2   Input error
!    exit_flag = 3   Output error
!    exit_flag = 4   IO error
!    exit_flag = 5   Configuration error
!    exit_flag = 6   Partition error
!    exit_flag = 7   Illegal input parameter
!    exit_flag = 8   Fatal algorithm result
!
        integer :: exit_flag = 0
        integer :: blowup = 0
        integer :: NoError = 0
!
!  Set threshold maximum speed (m/s) and density anomaly (kg/m3) to
!  test if the model is blowing-up.
!
        real(r8), allocatable :: maxspeed(:)
        real(r8), allocatable :: maxrho(:)
!
        real(r8) :: max_speed = 20.0_r8         ! m/s
        real(r8) :: max_rho = 200.0_r8          ! kg/m3
!
!  Interpolation scheme.
!
        integer, parameter :: linear = 0        ! linear interpolation
        integer, parameter :: cubic  = 1        ! cubic  interpolation
!
        integer :: InterpFlag = linear          ! interpolation flag
!
!  Shallowest and Deepest levels to apply bottom momemtum stresses as
!  a bodyforce
!
        integer, allocatable :: levsfrc(:)
        integer, allocatable :: levbfrc(:)
!
!  Vertical coordinates transform.  Currently, there are two vertical
!  transformation equations (see set_scoord.F for details):
!
!    Original transform (Vtransform=1):
!
!         z_r(x,y,s,t) = Zo_r + zeta(x,y,t) * [1.0 + Zo_r / h(x,y)]
!
!                 Zo_r = hc * [s(k) - C(k)] + C(k) * h(x,y)
!
!    New transform (Vtransform=2):
!
!         z_r(x,y,s,t) = zeta(x,y,t) + [zeta(x,y,t)+ h(x,y)] * Zo_r
!
!                 Zo_r = [hc * s(k) + C(k) * h(x,y)] / [hc + h(x,y)]
!
        integer, allocatable :: Vtransform(:)
!
!  Vertical grid stretching function flag:
!
!    Vstretcing = 1   Original function (Song and Haidvogel, 1994)
!               = 2   A. Shchepetkin (ROMS-UCLA) function
!               = 3   R. Geyer BBL function
!
        integer, allocatable :: Vstretching(:)
!
!  Vertical grid stretching parameters.
!
!    Tcline        Width (m) of surface or bottom boundary layer in
!                    which higher vertical resolution is required
!                    during stretching.
!    hc            S-coordinate critical depth, hc=MIN(hmin,Tcline).
!    theta_s       S-coordinate surface control parameter.
!    theta_b       S-coordinate bottom control parameter.
!
        real(r8), allocatable :: Tcline(:)      ! m, positive
        real(r8), allocatable :: hc(:)          ! m, positive
        real(r8), allocatable :: theta_s(:)     ! 0 < theta_s < 20
        real(r8), allocatable :: theta_b(:)     ! 0 < theta_b < 1
!
!  Bathymetry range values.
!
        real(r8), allocatable :: hmin(:)        ! m, positive
        real(r8), allocatable :: hmax(:)        ! m, positive
!
!  Length (m) of domain box in the XI- and ETA-directions.
!
        real(r8), allocatable :: xl(:)          ! m
        real(r8), allocatable :: el(:)          ! m
!
!  Minimum and Maximum longitude and latitude at RHO-points
!
        real(r8), allocatable :: LonMin(:)      ! degrees east
        real(r8), allocatable :: LonMax(:)      ! degrees east
        real(r8), allocatable :: LatMin(:)      ! degrees north
        real(r8), allocatable :: LatMax(:)      ! degrees north
!
!  Constant used in the Shchepetkin boundary conditions for 2D momentum,
!  Co = 1.0_r8/(2.0_r8+SQRT(2.0_r8)).
!
        real(r8) :: Co
!
!  Number of digits in grid size for format statements.
!
        integer, allocatable :: Idigits(:)
        integer, allocatable :: Jdigits(:)
        integer, allocatable :: Kdigits(:)
!
!  Diagnostic volume averaged variables.
!
        integer, allocatable :: first_time(:)
        real(r8) :: avgke = 0.0_r8              ! Kinetic energy
        real(r8) :: avgpe = 0.0_r8              ! Potential energy
        real(r8) :: avgkp = 0.0_r8              ! Total energy
        real(r8) :: volume = 0.0_r8             ! diagnostics volume
        real(r8) :: ad_volume                   ! adjoint volume
        real(r8) :: TotVolume                   ! Total volume
        real(r8) :: MaxVolume                   ! Minimum cell volume
        real(r8) :: MinVolume                   ! Maximum cell volume
!
!  Minimun and maximum grid spacing
!
        real(r8), allocatable :: DXmin(:)
        real(r8), allocatable :: DXmax(:)
        real(r8), allocatable :: DYmin(:)
        real(r8), allocatable :: DYmax(:)
        real(r8), allocatable :: DZmin(:)
        real(r8), allocatable :: DZmax(:)
!
!  Maximum size of a grid node (m) over the whole curvilinear grid
!  application. Used for scaling horizontal mixing by the grid size.
!
        real(r8), allocatable :: grdmax(:)
!
!  Courant Numbers due to gravity wave speed limits.
!
        real(r8) :: Cg_min                      ! Minimun barotropic
        real(r8) :: Cg_max                      ! Maximun barotropic
        real(r8) :: Cg_Cor                      ! Maximun Coriolis
!
!  Time dependent Counrant Numbers due to velocity components and
!  indices location of maximum value.
!
        integer :: max_Ci, max_Cj, max_Ck       ! maximum location
        real(r8) :: max_C                       ! maximum total
        real(r8) :: max_Cu                      ! maximum I-component
        real(r8) :: max_Cv                      ! maximum J-component
        real(r8) :: max_Cw                      ! maximum K-component
!
!  Linear equation of state parameters.
!
!    R0            Background constant density anomaly (kg/m3).
!    Tcoef         Thermal expansion coefficient (1/Celsius).
!    Scoef         Saline contraction coefficient (1/PSU).
!
        real(r8), allocatable :: R0(:)
        real(r8), allocatable :: Tcoef(:)
        real(r8), allocatable :: Scoef(:)
!
!  Background potential temperature (Celsius) and salinity (PSU) values
!  used in analytical initializations.
!
        real(r8), allocatable :: T0(:)
        real(r8), allocatable :: S0(:)
!
!  Slipperiness variable, either 1.0 (free slip) or -1.0 (no slip).
!
        real(r8), allocatable :: gamma2(:)
!
!  Weighting coefficient for the newest (implicit) time step derivatives
!  in the Crack-Nicolson implicit scheme (usually, lambda=0.5).
!
        real(r8) :: lambda = 1.0_r8             ! backward implicit
!
!  Jerlov water type to assign everywhere, range values: 1 - 5.
!
        integer, allocatable :: lmd_Jwt(:)
!
!  Grid r-factor (non-dimensional).
!
        real(r8) :: rx0 = 0.0_r8                ! Beckmann and Haidvogel
        real(r8) :: rx1 = 0.0_r8                ! Haney
!
!  Linear (m/s) and quadratic (nondimensional) bottom drag coefficients.
!
        real(r8), allocatable :: rdrg(:)
        real(r8), allocatable :: rdrg2(:)
!
!  Minimum and maximum threshold for transfer coefficient of momentum.
!
        real(r8) :: Cdb_min = 0.000001_r8
        real(r8) :: Cdb_max = 0.5_r8
!
!  Surface and bottom roughness (m)
!
        real(r8), allocatable :: Zos(:)
        real(r8), allocatable :: Zob(:)
!
!  Minimum depth for wetting and drying (m).
!
        real(r8), allocatable :: Dcrit(:)
!
!  Mean density (Kg/m3) used when the Boussinesq approximation is
!  inferred.
!
        real(r8) :: rho0 = 1025.0_r8
!
!  Background Brunt-Vaisala frequency (1/s2)
!
        real(r8) :: bvf_bak = 0.00001_r8
!
!  Vector containing USER generic parameters.
!
        integer :: Nuser
        real(r8), dimension(25) :: user(25)
!
!  Weights for the time average of 2D fields.
!
        real(r8), allocatable :: weight(:,:,:)
!
!  Constants.
!
        real(r8), parameter :: pi = 3.14159265358979323846_r8
        real(r8), parameter :: deg2rad = pi / 180.0_r8
        real(r8), parameter :: rad2deg = 180.0_r8 / pi
        real(r8), parameter :: day2sec = 86400.0_r8
        real(r8), parameter :: sec2day = 1.0_r8 / 86400.0_r8
        real(r8), parameter :: spval = 1.0E+37_r8
        real(r8), parameter :: jul_off = 2440000.0_r8
!
!  Set special check value.  Notice that a smaller value is assigned
!  to account for both NetCDF fill value and roundoff. There are
!  many Matlab scripts out there that do not inquire correctly
!  the spval from the _FillValue attribute in single/double
!  precision.
!
        real(r8), parameter :: spval_check = 1.0E+35_r8
!
!-----------------------------------------------------------------------
!  Horizontal and vertical constant mixing coefficients.
!-----------------------------------------------------------------------
!
!    Akk_bak       Background vertical mixing coefficient (m2/s) for
!                    turbulent energy.
!    Akp_bak       Background vertical mixing coefficient (m2/s) for
!                    generic statistical field "psi".
!    Akt_bak       Background vertical mixing coefficient (m2/s) for
!                    tracers.
!    Akv_bak       Background vertical mixing coefficient (m2/s) for
!                    momentum.
!    Kdiff         Isopycnal mixing thickness diffusivity (m2/s) for
!                    tracers.
!    ad_visc2      ADM lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    nl_visc2      NLM lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    tl_visc2      TLM lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    visc2         Current lateral harmonic constant mixing coefficient
!                    (m2/s) for momentum.
!    ad_visc4      ADM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    nl_visc4      NLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    tl_visc4      TLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    visc4         Current lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for momentum.
!    ad_tnu2       ADM lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    nl_tnu2       NLM lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    tl_tnu2       TLM lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    tnu2          Current lateral harmonic constant mixing coefficient
!                    (m2/s) for tracer type variables.
!    ad_tnu4       ADM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    nl_tnu4       NLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    tl_tnu4       TLM lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    tnu4          Current lateral biharmonic (squared root) constant
!                     mixing coefficient (m2 s^-1/2) for tracers.
!    tkenu2        Lateral harmonic constant mixing coefficient
!                    (m2/s) for turbulent energy.
!    tkenu4        Lateral biharmonic (squared root) constant mixing
!                    coefficient (m2 s^-1/2) for turbulent energy.
!
        real(r8), allocatable :: Akk_bak(:)          ! m2/s
        real(r8), allocatable :: Akp_bak(:)          ! m2/s
        real(r8), allocatable :: Akv_bak(:)          ! m2/s
        real(r8), allocatable :: ad_visc2(:)         ! m2/s
        real(r8), allocatable :: nl_visc2(:)         ! m2/s
        real(r8), allocatable :: tl_visc2(:)         ! m2/s
        real(r8), allocatable :: visc2(:)            ! m2/s
        real(r8), allocatable :: ad_visc4(:)         ! m2 s-1/2
        real(r8), allocatable :: nl_visc4(:)         ! m2 s-1/2
        real(r8), allocatable :: tl_visc4(:)         ! m2 s-1/2
        real(r8), allocatable :: visc4(:)            ! m2 s-1/2
        real(r8), allocatable :: tkenu2(:)           ! m2/s
        real(r8), allocatable :: tkenu4(:)           ! m2 s-1/2
        real(r8), allocatable :: Akt_bak(:,:)        ! m2/s
        real(r8), allocatable :: Kdiff(:,:)          ! m2/s
        real(r8), allocatable :: ad_tnu2(:,:)        ! m2/s
        real(r8), allocatable :: nl_tnu2(:,:)        ! m2/s
        real(r8), allocatable :: tl_tnu2(:,:)        ! m2/s
        real(r8), allocatable :: tnu2(:,:)           ! m2/s
        real(r8), allocatable :: ad_tnu4(:,:)        ! m2 s-1/2
        real(r8), allocatable :: nl_tnu4(:,:)        ! m2 s-1/2
        real(r8), allocatable :: tl_tnu4(:,:)        ! m2 s-1/2
        real(r8), allocatable :: tnu4(:,:)           ! m2 s-1/2
!
!  Horizontal diffusive relaxation coefficients (m2/s) used to smooth
!  representer tangent linear solution during Picard iterations to
!  improve stability and convergence.
!
        real(r8), allocatable :: tl_M2diff(:)        ! 2D momentum
        real(r8), allocatable :: tl_M3diff(:)        ! 3D momentum
        real(r8), allocatable :: tl_Tdiff(:,:)       ! tracers
!
!  Basic state vertical mixing coefficient scale factors for adjoint
!  based algorithms. In some applications, a smaller/larger values of
!  vertical mixing are necessary for stability.
!
        real(r8), allocatable :: ad_Akv_fac(:)       ! ADM momentum
        real(r8), allocatable :: tl_Akv_fac(:)       ! TLM momentum
        real(r8), allocatable :: ad_Akt_fac(:,:)     ! ADM tracers
        real(r8), allocatable :: tl_Akt_fac(:,:)     ! TLM tracers
!
!-----------------------------------------------------------------------
!  IO parameters.
!-----------------------------------------------------------------------
!
!  Switches to activate creation and writing of output NetCDF files.
!
        logical, allocatable :: LdefADJ(:)       ! Adjoint file
        logical, allocatable :: LdefAVG(:)       ! Average file
        logical, allocatable :: LdefDIA(:)       ! Diagnostics file
        logical, allocatable :: LdefERR(:)       ! 4DVar error file
        logical, allocatable :: LdefFLT(:)       ! Floats file
        logical, allocatable :: LdefHIS(:)       ! History file
        logical, allocatable :: LdefHSS(:)       ! Hessian file
        logical, allocatable :: LdefINI(:)       ! Initial file
        logical, allocatable :: LdefIRP(:)       ! Initial RPM file
        logical, allocatable :: LdefITL(:)       ! Initial TLM file
        logical, allocatable :: LdefLCZ(:)       ! Lanczos file
        logical, allocatable :: LdefMOD(:)       ! 4DVAR file
        logical, allocatable :: LdefRST(:)       ! Restart file
        logical, allocatable :: LdefSTA(:)       ! Stations file
        logical, allocatable :: LdefTIDE(:)      ! tide forcing file
        logical, allocatable :: LdefTLM(:)       ! Tangent linear file
        logical, allocatable :: LdefTLF(:)       ! TLM/RPM impulse file
        logical, allocatable :: LwrtADJ(:)       ! Write adjoint file
        logical, allocatable :: LwrtAVG(:)       ! Write average file
        logical, allocatable :: LwrtDIA(:)       ! Write diagnostic file
        logical, allocatable :: LwrtHIS(:)       ! Write history file
        logical, allocatable :: LwrtPER(:)       ! Write during ensemble
        logical, allocatable :: LwrtRST(:)       ! Write restart file
        logical, allocatable :: LwrtTLM(:)       ! Write tangent file
        logical, allocatable :: LwrtTLF(:)       ! Write impulse file
        logical, allocatable :: LdefNRM(:,:)     ! Norm file
        logical, allocatable :: LwrtNRM(:,:)     ! Write norm file
!
!  Switch to write out adjoint 2D state arrays instead of IO solution
!  arrays and adjoint ocean time. This is used in 4DVAR for IO
!  maniputations.
!
        logical, allocatable :: LwrtState2d(:)
        logical, allocatable :: LwrtTime(:)
!
!  Switch to write out adjoint surface forcing fields adjusted by the
!  4DVAR algorithms.
!
        logical, allocatable :: Ladjusted(:)
!
!  Switch to write application set-up information to standard output.
!
        logical, allocatable :: LwrtInfo(:)
!
!  Switch used to create new output NetCDF files. If TRUE, new output
!  files are created. If FALSE, data is appended to an existing output
!  files.  Used only for history, average and station files.
!
        logical, allocatable :: ldefout(:)       ! New output files
!
!  Number of timesteps between creation of new output files.
!
        integer, allocatable :: ndefADJ(:)       ! Adjoint file
        integer, allocatable :: ndefAVG(:)       ! Average file
        integer, allocatable :: ndefDIA(:)       ! Diagnostics file
        integer, allocatable :: ndefHIS(:)       ! History file
        integer, allocatable :: ndefTLM(:)       ! Tangent linear file
!
!  Starting timestep for accumulation of output.
!
        integer, allocatable :: ntsAVG(:)        ! Average file
        integer, allocatable :: ntsDIA(:)        ! Diagnostics file
!
!  Number of timesteps between writing of output data.
!
        integer, allocatable :: nADJ(:)          ! Adjoint file
        integer, allocatable :: nAVG(:)          ! Average file
        integer, allocatable :: nDIA(:)          ! Diagnostics file
        integer, allocatable :: nFLT(:)          ! Floats file
        integer, allocatable :: nHIS(:)          ! History file
        integer, allocatable :: nRST(:)          ! Restart file
        integer, allocatable :: nSTA(:)          ! Stations file
        integer, allocatable :: nTLM(:)          ! Tangent linear file
!
!  Number of timesteps between print of single line information to
!  standard output.
!
        integer, allocatable :: ninfo(:)
!
!  Number of timesteps between 4DVAR adjustment of open boundaries.
!  In strong constraint 4DVAR, it is possible to open bounadies at
!  other intervals in addition to initial time. These parameters are
!  used to store the appropriate number of open boundary records in
!  output history NetCDF files.
!
!    Nbrec(:) = 1 + ntimes(:) / nOBC(:)
!
!  Here, it is assumed that nOBC is a multiple of NTIMES or greater
!  than NTIMES. If nOBC > NTIMES, only one record is stored in the
!  output history NetCDF files and the adjustment is for constant
!  open boundaries with constant correction.
!
        integer, allocatable :: nOBC(:)          ! number of timesteps
        integer, allocatable :: Nbrec(:)         ! number of records
        integer, allocatable :: OBCcount(:)      ! record counter
!
!  Number of timesteps between adjustment of 4DVAR surface forcing
!  fields. In strong constraint 4DVAR, it is possible to adjust surface
!  forcing fields at other intervals in addition to initial time.
!  These parameters are used to store the appropriate number of
!  surface forcing records in output history NetCDF files.
!
!    Nfrec(:) = 1 + ntimes(:) / nSFF(:)
!
!  Here, it is assumed that nSFF is a multiple of NTIMES or greater
!  than NTIMES. If nSFF > NTIMES, only one record is stored in the
!  output history NetCDF files and the adjustment is for constant
!  forcing with constant correction.
!
        integer, allocatable :: nSFF(:)          ! number of timesteps
        integer, allocatable :: Nfrec(:)         ! number of records
        integer, allocatable :: SFcount(:)       ! record counter
!
!  Restart time record to read from disk and use as the initial
!  conditions. Use nrrec=0 for new solutions. If nrrec is negative
!  (say, nrrec=-1), the model will restart from the most recent
!  time record. That is, the initialization record is assigned
!  internally.
!
        integer, allocatable :: nrrec(:)
!
!  Switch to activate processing of input data.  This switch becomes
!  very useful when reading input data serially in parallel
!  applications.
!
        logical, allocatable :: synchro_flag(:)
!
!  Switch to inialize model with latest time record from initial
!  (restart/history) NetCDF file.
!
        logical, allocatable :: LastRec(:)
!
!  Generalized Statbility Theory (GST) parameters.
!
        logical :: LmultiGST          ! multiple eigenvector file switch
        logical :: LrstGST            ! restart switch
        integer :: MaxIterGST         ! Number of iterations
        integer :: nGST               ! check pointing interval
!
!  Switches used to recycle time records in some output file. If TRUE,
!  only the latest two time records are maintained.  If FALSE, all
!  field records are saved.
!
        logical, allocatable :: LcycleADJ(:)
        logical, allocatable :: LcycleRST(:)
        logical, allocatable :: LcycleTLM(:)
!
!-----------------------------------------------------------------------
!  Adjoint sensitivity parameters.
!-----------------------------------------------------------------------
!
!  Starting and ending vertical levels of the 3D adjoint state whose
!  sensitivity is required.
!
        integer, allocatable :: KstrS(:)           ! starting level
        integer, allocatable :: KendS(:)           ! ending level
!
!  Starting and ending day for adjoint sensitivity forcing.
!
        real(r8), allocatable :: DstrS(:)          ! starting day
        real(r8), allocatable :: DendS(:)          ! ending day
!
!-----------------------------------------------------------------------
!  Stochastic optimals parameters.
!-----------------------------------------------------------------------
!
!  Stochastic optimals forcing records counter.
!
        integer, allocatable :: SOrec(:)
!
!  Trace of stochastic optimals matrix.
!
        real(r8), allocatable :: TRnorm(:)
!
!  Stochastic optimals time decorrelation scale (days) assumed for
!  red noise processes.
!
        real(r8), allocatable :: SO_decay(:)
!
!  Stochastic optimals surface forcing standard deviation for
!  dimensionalization.
!
        real(r8), allocatable :: SO_sdev(:,:)
!
!-----------------------------------------------------------------------
!  Nudging variables for passive (outflow) and active (inflow) oepn
!  boundary conditions.
!-----------------------------------------------------------------------
!
!    iwest         West  identification index in boundary arrays.
!    isouth        South identification index in boundary arrays.
!    ieast         East  identification index in boundary arrays.
!    inorth        North identification index in boundary arrays.
!    obcfac        Factor between passive and active open boundary
!                    conditions (nondimensional and greater than one).
!                    The nudging time scales for the active conditions
!                    are obtained by multiplying the passive values by
!                    factor.
!    FSobc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards free-surface data at  inflow.
!    FSobc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards free-surface data at outflow.
!    M2obc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards 2D momentum data at  inflow.
!    M2obc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards 2D momentum data at outflow.
!    M3obc_in      Active and strong time-scale (1/sec) coefficients
!                    for nudging towards 3D momentum data at  inflow.
!    M3obc_out     Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards 3D momentum data at outflow.
!    Tobc_in       Active and strong time-scale (1/sec) coefficients
!                    for nudging towards tracer data at  inflow.
!    Tobc_out      Passive and weak  time-scale (1/sec) coefficients
!                    for nudging towards tracer data at outflow.
!
        integer, parameter :: iwest = 1
        integer, parameter :: isouth = 2
        integer, parameter :: ieast = 3
        integer, parameter :: inorth = 4
        real(r8), allocatable :: obcfac(:)
        real(r8), allocatable :: FSobc_in(:,:)
        real(r8), allocatable :: FSobc_out(:,:)
        real(r8), allocatable :: M2obc_in(:,:)
        real(r8), allocatable :: M2obc_out(:,:)
        real(r8), allocatable :: M3obc_in(:,:)
        real(r8), allocatable :: M3obc_out(:,:)
        real(r8), allocatable :: Tobc_in(:,:,:)
        real(r8), allocatable :: Tobc_out(:,:,:)
!
!  Inverse time-scales (1/s) for nudging at open boundaries and sponge
!  areas.
!
        real(r8), allocatable :: Znudg(:)          ! Free-surface
        real(r8), allocatable :: M2nudg(:)         ! 2D momentum
        real(r8), allocatable :: M3nudg(:)         ! 3D momentum
        real(r8), allocatable :: Tnudg(:,:)        ! Tracers
!
!  Variables used to impose mass flux conservation in open boundary
!  configurations.
!
        real(r8) :: bc_area = 0.0_r8
        real(r8) :: bc_flux = 0.0_r8
        real(r8) :: ubar_xs = 0.0_r8
!
!-----------------------------------------------------------------------
!  Water clarity parameters.
!-----------------------------------------------------------------------
!
!    lmd_mu1       Reciprocal of the absorption coefficient for solar
!                    wavelength band 1 as a function of the Jerlov
!                    water type.
!    lmd_mu2       Reciprocal of the absorption coefficient for solar
!                    wavelength band 2 as a function of the Jerlov
!                    water type.
!    lmd_r1        Fraction of total radiance for wavelength band 1 as
!                    a function of the Jerlov water type.
!
        real(r8), dimension(5) :: lmd_mu1 =                             &
     &            (/ 0.35_r8, 0.6_r8, 1.0_r8, 1.5_r8, 1.4_r8 /)
        real(r8), dimension(5) :: lmd_mu2 =                             &
     &            (/ 23.0_r8, 20.0_r8, 17.0_r8, 14.0_r8, 7.9_r8 /)
        real(r8), dimension(5) :: lmd_r1 =                              &
     &            (/ 0.58_r8, 0.62_r8, 0.67_r8, 0.77_r8, 0.78_r8 /)
!
!-----------------------------------------------------------------------
!  Large et al. (1994) K-profile parameterization.
!-----------------------------------------------------------------------
!
!    lmd_Ri0       Critical gradient Richardson number below which
!                    turbulent mixing occurs.
!    lmd_Rrho0     Value of double-diffusive density ratio where
!                    mixing goes to zero in salt fingering.
!    lmd_bvfcon    Brunt-Vaisala frequency (1/s2) limit for convection.
!    lmd_fdd       Scaling factor for double diffusion of temperature
!                    in salt fingering case (lmd_fdd=0.7).
!    lmd_nu        Molecular viscosity (m2/s).
!    lmd_nu0c      Maximum interior convective viscosity and diffusivity
!                    due to shear instability.
!    lmd_nu0m      Maximum interior viscosity (m2/s) due shear
!                    instability.
!    lmd_nu0s      Maximum interior diffusivity (m2/s) due shear
!                    instability.
!    lmd_nuf       Scaling factor for double diffusion in salt
!                    fingering.
!    lmd_nuwm      Interior viscosity (m2/s) due to wave breaking.
!    lmd_nuws      Interior diffusivity (m2/s) due to wave breaking.
!    lmd_sdd1      Double diffusion constant for salinity in diffusive
!                    convection case (lmd_sdd1=0.15).
!    lmd_sdd2      Double diffusion constant for salinity in diffusive
!                    convection case (lmd_sdd2=1.85).
!    lmd_sdd3      Double diffusion constant for salinity in diffusive
!                    convection case (lmd_sdd3=0.85).
!    lmd_tdd1      Double diffusion constant for temperature
!                    in diffusive convection case (lmd_tdd1=0.909).
!    lmd_tdd2      Double diffusion constant for temperature in
!                    diffusive convection case (lmd_tdd2=4.6).
!    lmd_tdd3      Double diffusion constant for temperature in
!                    diffusive convection case (lmd_tdd3=0.54).
!
        real(r8) :: lmd_Ri0 = 0.7_r8          ! non-dimensional
        real(r8) :: lmd_Rrho0 = 1.9_r8        ! m2/s
        real(r8) :: lmd_bvfcon = -2.0E-5_r8   ! 1/s2
        real(r8) :: lmd_fdd = 0.7_r8          ! non-dimensional
        real(r8) :: lmd_nu = 1.5E-6_r8        ! m2/s
        real(r8) :: lmd_nu0c = 0.01_r8        ! m2/s
        real(r8) :: lmd_nu0m = 10.0E-4_r8     ! m2/s
        real(r8) :: lmd_nu0s = 10.0E-4_r8     ! m2/s
        real(r8) :: lmd_nuf = 10.0E-4_r8      ! m2/s
        real(r8) :: lmd_nuwm = 1.0E-5_r8      ! m2/s
        real(r8) :: lmd_nuws = 1.0E-6_r8      ! m2/s
        real(r8) :: lmd_sdd1 = 0.15_r8        ! non-dimensional
        real(r8) :: lmd_sdd2 = 1.85_r8        ! non-dimensional
        real(r8) :: lmd_sdd3 = 0.85_r8        ! non-dimensional
        real(r8) :: lmd_tdd1 = 0.909_r8       ! non-dimensional
        real(r8) :: lmd_tdd2 = 4.6_r8         ! non-dimensional
        real(r8) :: lmd_tdd3 = 0.54_r8        ! non-dimensional
!
!-----------------------------------------------------------------------
!  Large et al. (1994) oceanic boundary layer parameters.
!-----------------------------------------------------------------------
!
!    lmd_Cg        Proportionality coefficient parameterizing nonlocal
!                    transport.
!    lmd_Cstar     Proportionality coefficient parameterizing nonlocal
!                    transport.
!    lmd_Cv        Ratio of interior Brunt-Vaisala frequency to that
!                    at entrainment depth "he".
!    lmd_Ric       Critical bulk Richardson number.
!    lmd_am        Coefficient of flux profile for momentum in their
!                    1/3 power law regimes.
!    lmd_as        Coefficient of flux profile for tracers in their
!                    1/3 power law regimes.
!    lmd_betaT     Ratio of entrainment flux to surface buoyancy flux.
!    lmd_cekman    Constant used in the computation of Ekman depth.
!    lmd_cmonob    Constant used in the computation of Monin-Obukhov
!                    depth.
!    lmd_cm        Coefficient of flux profile for momentum in their
!                    1/3 power law regimes.
!    lmd_cs        Coefficient of flux profile for tracers in their
!                    1/3 power law regimes.
!    lmd_epsilon   Non-dimensional extent of the surface layer.
!    lmd_zetam     Maximum stability parameter "zeta" value of the 1/3
!                    power law regime of flux profile for momentum.
!    lmd_zetas     Maximum stability parameter "zeta" value of the 1/3
!                    power law regime of flux profile for tracers.
!
        real(r8) :: lmd_Cg
        real(r8) :: lmd_Cstar = 10.0_r8
        real(r8) :: lmd_Cv = 1.25_r8
        real(r8) :: lmd_Ric = 0.3_r8
        real(r8) :: lmd_am = 1.257_r8
        real(r8) :: lmd_as = -28.86_r8
        real(r8) :: lmd_betaT = -0.2_r8
        real(r8) :: lmd_cekman = 0.7_r8
        real(r8) :: lmd_cmonob = 1.0_r8
        real(r8) :: lmd_cm = 8.36_r8
        real(r8) :: lmd_cs = 98.96_r8
        real(r8) :: lmd_epsilon = 0.1_r8
        real(r8) :: lmd_zetam = -0.2_r8
        real(r8) :: lmd_zetas = -1.0_r8
!
!-----------------------------------------------------------------------
!  Generic Length Scale parameters.
!-----------------------------------------------------------------------
!
!    gls_Gh0
!    gls_Ghcri
!    gls_Ghmin
!    gls_Kmin      Minimum value of specific turbulent kinetic energy.
!    gls_Pmin      Minimum Value of dissipation.
!    gls_cmu0      Stability coefficient (non-dimensional).
!    gls_c1        Shear production coefficient (non-dimensional).
!    gls_c2        Dissipation coefficient (non-dimensional).
!    gls_c3m       Buoyancy production coefficient (minus).
!    gls_c3p       Buoyancy production coefficient (plus).
!    gls_E2
!    gls_m         Turbulent kinetic energy exponent (non-dimensional).
!    gls_n         Turbulent length scale exponent (non-dimensional).
!    gls_p         Stability exponent (non-dimensional).
!    gls_sigk      Constant Schmidt number (non-dimensional) for
!                    turbulent kinetic energy diffusivity.
!    gls_sigp      Constant Schmidt number (non-dimensional) for
!                    turbulent generic statistical field, "psi".
!
        real(r8), allocatable :: gls_m(:)
        real(r8), allocatable :: gls_n(:)
        real(r8), allocatable :: gls_p(:)
        real(r8), allocatable :: gls_sigk(:)
        real(r8), allocatable :: gls_sigp(:)
        real(r8), allocatable :: gls_cmu0(:)
        real(r8), allocatable :: gls_cmupr(:)
        real(r8), allocatable :: gls_c1(:)
        real(r8), allocatable :: gls_c2(:)
        real(r8), allocatable :: gls_c3m(:)
        real(r8), allocatable :: gls_c3p(:)
        real(r8), allocatable :: gls_Kmin(:)
        real(r8), allocatable :: gls_Pmin(:)
!
! Constants used in the various formulation of surface flux boundary
! conditions for the GLS vertical turbulence closure in terms of
! Charnok surface roughness (CHARNOK_ALPHA), roughness from wave
! amplitude (zos_hsig_alpha), wave dissipation (SZ_ALPHA), and
! Craig and Banner wave breaking (CRGBAN_CW).
!
        real(r8), allocatable :: charnok_alpha(:)
        real(r8), allocatable :: zos_hsig_alpha(:)
        real(r8), allocatable :: sz_alpha(:)
        real(r8), allocatable :: crgban_cw(:)
!
!-----------------------------------------------------------------------
!  Tangent linear and adjoint model parameters.
!-----------------------------------------------------------------------
!
!  Tangent linear and adjoint model control switches.
!
        logical :: TLmodel = .FALSE.
        logical :: ADmodel = .FALSE.
!
      CONTAINS
!
      SUBROUTINE allocate_scalars
!
!=======================================================================
!                                                                      !
!  This routine allocates structure and  several variables in module   !
!  that depend on the number of nested grids.                          !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: ng
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate and initialize variables in module structure.
!-----------------------------------------------------------------------
!
      allocate ( SCALARS(Ngrids) )
      DO ng=1,Ngrids
        allocate ( SCALARS(ng) % Fstate(7+2*MT) )
        SCALARS(ng) % Fstate(1:7+2*MT) = .FALSE.
        allocate ( SCALARS(ng) % Lstate(5+MT) )
        SCALARS(ng) % Lstate(1:5+MT) = .FALSE.
        allocate ( SCALARS(ng) % Cs_r(N(ng)) )
        SCALARS(ng) % Cs_r(1:N(ng)) = IniVal
        allocate ( SCALARS(ng) % Cs_w(0:N(ng)) )
        SCALARS(ng) % Cs_w(0:N(ng)) = IniVal
        allocate ( SCALARS(ng) % sc_r(N(ng)) )
        SCALARS(ng) % sc_r(1:N(ng)) = IniVal
        allocate ( SCALARS(ng) % sc_w(0:N(ng)) )
        SCALARS(ng) % sc_w(0:N(ng)) = IniVal
      END DO
!
!  Allocate variables that require special treatment in shared-memory.
!  These variables are private for each thread to avoid collisions.
!
      allocate ( PREDICTOR_2D_STEP(Ngrids) )
      allocate ( indx1(Ngrids) )
      allocate ( iic(Ngrids) )
      allocate ( iif(Ngrids) )
      allocate ( FrcRec(Ngrids) )
      allocate ( SOrec(Ngrids) )
!
!-----------------------------------------------------------------------
!  Allocate variables.
!-----------------------------------------------------------------------
!
      allocate ( PerfectRST(Ngrids) )
      allocate ( ndtfast(Ngrids) )
      allocate ( nfast(Ngrids) )
      allocate ( dt(Ngrids) )
      allocate ( dtfast(Ngrids) )
      allocate ( TimeEnd(Ngrids) )
      allocate ( AVGtime(Ngrids) )
      allocate ( DIAtime(Ngrids) )
      allocate ( IMPtime(Ngrids) )
      allocate ( ObsTime(Ngrids) )
      allocate ( FrcTime(Ngrids) )
      allocate ( ntimes(Ngrids) )
      allocate ( synchro_flag(Ngrids) )
      allocate ( first_time(Ngrids) )
      allocate ( ntfirst(Ngrids) )
      allocate ( ntstart(Ngrids) )
      allocate ( ntend(Ngrids) )
      allocate ( step_counter(Ngrids) )
      allocate ( tdays(Ngrids) )
      allocate ( time(Ngrids) )
      allocate ( time_code(Ngrids) )
      allocate ( NrecFrc(Ngrids) )
      allocate ( NudgingCoeff(Ngrids) )
      allocate ( ObcData(Ngrids) )
      allocate ( Lbiology(Ngrids) )
      allocate ( Lfloats(Ngrids) )
      allocate ( Lsediment(Ngrids) )
      allocate ( Lstations(Ngrids) )
      allocate ( CompositeGrid(4,Ngrids) )
      allocate ( RefinedGrid(Ngrids) )
      allocate ( RefineScale(Ngrids) )
      allocate ( GetDonorData(Ngrids) )
      allocate ( EWperiodic(Ngrids) )
      allocate ( NSperiodic(Ngrids) )
      allocate ( VolCons(4,Ngrids) )
      allocate ( LtracerCLM(MT,Ngrids) )
      allocate ( LtracerSrc(MT,Ngrids) )
      allocate ( maxspeed(Ngrids) )
      allocate ( maxrho(Ngrids) )
      allocate ( levsfrc(Ngrids) )
      allocate ( levbfrc(Ngrids) )
      allocate ( Vtransform(Ngrids) )
      allocate ( Vstretching(Ngrids) )
      allocate ( Tcline(Ngrids) )
      allocate ( hc(Ngrids) )
      allocate ( theta_s(Ngrids) )
      allocate ( theta_b(Ngrids) )
      allocate ( hmin(Ngrids) )
      allocate ( hmax(Ngrids) )
      allocate ( xl(Ngrids) )
      allocate ( el(Ngrids) )
      allocate ( LonMin(Ngrids) )
      allocate ( LonMax(Ngrids) )
      allocate ( LatMin(Ngrids) )
      allocate ( LatMax(Ngrids) )
      allocate ( Idigits(Ngrids) )
      allocate ( Jdigits(Ngrids) )
      allocate ( Kdigits(Ngrids) )
      allocate ( DXmin(Ngrids) )
      allocate ( DXmax(Ngrids) )
      allocate ( DYmin(Ngrids) )
      allocate ( DYmax(Ngrids) )
      allocate ( DZmin(Ngrids) )
      allocate ( DZmax(Ngrids) )
      allocate ( grdmax(Ngrids) )
      allocate ( R0(Ngrids) )
      allocate ( Tcoef(Ngrids) )
      allocate ( Scoef(Ngrids) )
      allocate ( T0(Ngrids) )
      allocate ( S0(Ngrids) )
      allocate ( gamma2(Ngrids) )
      allocate ( lmd_Jwt(Ngrids) )
      allocate ( rdrg(Ngrids) )
      allocate ( rdrg2(Ngrids) )
      allocate ( Zos(Ngrids) )
      allocate ( Zob(Ngrids) )
      allocate ( Dcrit(Ngrids) )
      allocate ( weight(2,0:256,Ngrids) )
      allocate ( Akk_bak(Ngrids) )
      allocate ( Akp_bak(Ngrids) )
      allocate ( Akv_bak(Ngrids) )
      allocate ( ad_visc2(Ngrids) )
      allocate ( nl_visc2(Ngrids) )
      allocate ( tl_visc2(Ngrids) )
      allocate ( visc2(Ngrids) )
      allocate ( ad_visc4(Ngrids) )
      allocate ( nl_visc4(Ngrids) )
      allocate ( tl_visc4(Ngrids) )
      allocate ( visc4(Ngrids) )
      allocate ( tkenu2(Ngrids) )
      allocate ( tkenu4(Ngrids) )
      allocate ( Akt_bak(MT,Ngrids) )
      allocate ( Kdiff(MT,Ngrids) )
      allocate ( ad_tnu2(MT,Ngrids) )
      allocate ( nl_tnu2(MT,Ngrids) )
      allocate ( tl_tnu2(MT,Ngrids) )
      allocate ( tnu2(MT,Ngrids) )
      allocate ( ad_tnu4(MT,Ngrids) )
      allocate ( nl_tnu4(MT,Ngrids) )
      allocate ( tl_tnu4(MT,Ngrids) )
      allocate ( tnu4(MT,Ngrids) )
      allocate ( tl_M2diff(Ngrids) )
      allocate ( tl_M3diff(Ngrids) )
      allocate ( tl_Tdiff(MT,Ngrids) )
      allocate ( ad_Akv_fac(Ngrids) )
      allocate ( tl_Akv_fac(Ngrids) )
      allocate ( ad_Akt_fac(MT,Ngrids) )
      allocate ( tl_Akt_fac(MT,Ngrids) )
      allocate ( LdefADJ(Ngrids) )
      allocate ( LdefAVG(Ngrids) )
      allocate ( LdefDIA(Ngrids) )
      allocate ( LdefERR(Ngrids) )
      allocate ( LdefFLT(Ngrids) )
      allocate ( LdefHIS(Ngrids) )
      allocate ( LdefHSS(Ngrids) )
      allocate ( LdefINI(Ngrids) )
      allocate ( LdefIRP(Ngrids) )
      allocate ( LdefITL(Ngrids) )
      allocate ( LdefLCZ(Ngrids) )
      allocate ( LdefMOD(Ngrids) )
      allocate ( LdefRST(Ngrids) )
      allocate ( LdefSTA(Ngrids) )
      allocate ( LdefTIDE(Ngrids) )
      allocate ( LdefTLM(Ngrids) )
      allocate ( LdefTLF(Ngrids) )
      allocate ( LwrtADJ(Ngrids) )
      allocate ( LwrtAVG(Ngrids) )
      allocate ( LwrtDIA(Ngrids) )
      allocate ( LwrtHIS(Ngrids) )
      allocate ( LwrtPER(Ngrids) )
      allocate ( LwrtRST(Ngrids) )
      allocate ( LwrtTLM(Ngrids) )
      allocate ( LwrtTLF(Ngrids) )
      allocate ( LdefNRM(4,Ngrids) )
      allocate ( LwrtNRM(4,Ngrids) )
      allocate ( LwrtState2d(Ngrids) )
      allocate ( LwrtTime(Ngrids) )
      allocate ( Ladjusted(Ngrids) )
      allocate ( LwrtInfo(Ngrids) )
      allocate ( ldefout(Ngrids) )
      allocate ( ndefADJ(Ngrids) )
      allocate ( ndefAVG(Ngrids) )
      allocate ( ndefDIA(Ngrids) )
      allocate ( ndefHIS(Ngrids) )
      allocate ( ndefTLM(Ngrids) )
      allocate ( ntsAVG(Ngrids) )
      allocate ( ntsDIA(Ngrids) )
      allocate ( nADJ(Ngrids) )
      allocate ( nAVG(Ngrids) )
      allocate ( nDIA(Ngrids) )
      allocate ( nFLT(Ngrids) )
      allocate ( nHIS(Ngrids) )
      allocate ( nRST(Ngrids) )
      allocate ( nSTA(Ngrids) )
      allocate ( nTLM(Ngrids) )
      allocate ( ninfo(Ngrids) )
      allocate ( nOBC(Ngrids) )
      allocate ( Nbrec(Ngrids) )
      allocate ( OBCcount(Ngrids) )
      allocate ( nSFF(Ngrids) )
      allocate ( Nfrec(Ngrids) )
      allocate ( SFcount(Ngrids) )
      allocate ( nrrec(Ngrids) )
      allocate ( LastRec(Ngrids) )
      allocate ( LcycleADJ(Ngrids) )
      allocate ( LcycleRST(Ngrids) )
      allocate ( LcycleTLM(Ngrids) )
      allocate ( KstrS(Ngrids) )
      allocate ( KendS(Ngrids) )
      allocate ( DstrS(Ngrids) )
      allocate ( DendS(Ngrids) )
      allocate ( TRnorm(Ngrids) )
      allocate ( SO_decay(Ngrids) )
      allocate ( SO_sdev(7+2*MT,Ngrids) )
      allocate ( obcfac(Ngrids) )
      allocate ( FSobc_in(Ngrids,4) )
      allocate ( FSobc_out(Ngrids,4) )
      allocate ( M2obc_in(Ngrids,4) )
      allocate ( M2obc_out(Ngrids,4) )
      allocate ( M3obc_in(Ngrids,4) )
      allocate ( M3obc_out(Ngrids,4) )
      allocate ( Tobc_in(MT,Ngrids,4) )
      allocate ( Tobc_out(MT,Ngrids,4) )
      allocate ( Znudg(Ngrids) )
      allocate ( M2nudg(Ngrids) )
      allocate ( M3nudg(Ngrids) )
      allocate ( Tnudg(MT,Ngrids) )
      allocate ( gls_m(Ngrids) )
      allocate ( gls_n(Ngrids) )
      allocate ( gls_p(Ngrids) )
      allocate ( gls_sigk(Ngrids) )
      allocate ( gls_sigp(Ngrids) )
      allocate ( gls_cmu0(Ngrids) )
      allocate ( gls_cmupr(Ngrids) )
      allocate ( gls_c1(Ngrids) )
      allocate ( gls_c2(Ngrids) )
      allocate ( gls_c3m(Ngrids) )
      allocate ( gls_c3p(Ngrids) )
      allocate ( gls_Kmin(Ngrids) )
      allocate ( gls_Pmin(Ngrids) )
      allocate ( charnok_alpha(Ngrids) )
      allocate ( zos_hsig_alpha(Ngrids) )
      allocate ( sz_alpha(Ngrids) )
      allocate ( crgban_cw(Ngrids) )
      RETURN
      END SUBROUTINE allocate_scalars
      SUBROUTINE initialize_scalars
!
!=======================================================================
!                                                                      !
!  This routine initializes several variables in module for all nested !
!  grids.                                               !              !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: i, ic, j, ng, itrc
      real(r8) :: one, zero
      real(r8), parameter :: IniVal = 0.0_r8
!
!---------------------------------------------------------------------
!  Set tracer identification indices.
!---------------------------------------------------------------------
!
      itemp=1
      isalt=2
      ic=NAT
!
!-----------------------------------------------------------------------
!  Activate all computation control switches.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        LastRec(ng)=.FALSE.
        CompositeGrid(1:4,ng)=.FALSE.
        RefinedGrid(ng)=.FALSE.
        GetDonorData(ng)=.FALSE.
        Lbiology(ng)=.TRUE.
        LcycleADJ(ng)=.FALSE.
        LcycleRST(ng)=.FALSE.
        LcycleTLM(ng)=.FALSE.
        Lfloats(ng)=.TRUE.
        Lsediment(ng)=.TRUE.
        Lstations(ng)=.TRUE.
      END DO
!
!-----------------------------------------------------------------------
!  Initialize several scalar variables.
!-----------------------------------------------------------------------
!
      one=1.0_r8
      zero=0.0_r8
      Co=1.0_r8/(2.0_r8+SQRT(2.0_r8))
      gorho0=g/rho0
      DO ng=1,Ngrids
        EWperiodic(ng)=.FALSE.
        NSperiodic(ng)=.FALSE.
        NudgingCoeff(ng)=.FALSE.
        ObcData(ng)=.FALSE.
        RefineScale(ng)=0
        ntfirst(ng)=1
        ntstart(ng)=1
        ntend(ng)=0
        step_counter(ng)=0
        gamma2(ng)=-1.0_r8
        Vtransform(ng)=1
        Vstretching(ng)=1
        first_time(ng)=0
        Idigits(ng)=INT(LOG10(REAL(Lm(ng),r8)))+1
        Jdigits(ng)=INT(LOG10(REAL(Mm(ng),r8)))+1
        Kdigits(ng)=INT(LOG10(REAL(N (ng),r8)))+1
        DO itrc=1,MT
          LtracerCLM(itrc,ng)=.FALSE.
          LtracerSrc(itrc,ng)=.FALSE.
          ad_Akt_fac(itrc,ng)=1.0_r8
          tl_Akt_fac(itrc,ng)=1.0_r8
          ad_tnu2(itrc,ng)=IniVal
          nl_tnu2(itrc,ng)=IniVal
          tl_tnu2(itrc,ng)=IniVal
          tnu2(itrc,ng)=IniVal
          ad_tnu4(itrc,ng)=IniVal
          nl_tnu4(itrc,ng)=IniVal
          tl_tnu4(itrc,ng)=IniVal
          tnu4(itrc,ng)=IniVal
        END DO
        ad_Akv_fac(ng)=1.0_r8
        tl_Akv_fac(ng)=1.0_r8
        ad_visc2(ng)=IniVal
        nl_visc2(ng)=IniVal
        tl_visc2(ng)=IniVal
        visc2(ng)=IniVal
        ad_visc4(ng)=IniVal
        nl_visc4(ng)=IniVal
        tl_visc4(ng)=IniVal
        visc4(ng)=IniVal
        DO i=1,4
          VolCons(i,ng)=.FALSE.
          FSobc_in (ng,i)=IniVal
          FSobc_out(ng,i)=IniVal
          M2obc_in (ng,i)=IniVal
          M2obc_out(ng,i)=IniVal
          M3obc_in (ng,i)=IniVal
          M3obc_out(ng,i)=IniVal
        END DO
      END DO
      Tobc_in = IniVal
      Tobc_out = IniVal
!
!  Proportionality coefficient parameterizing boundary layer
!  nonlocal transport.
!
      lmd_Cg=lmd_Cstar*                                                 &
     &       vonKar*(lmd_cs*vonKar*lmd_epsilon)**(1.0_r8/3.0_r8)
!
!  Initialize several IO flags.
!
      LmultiGST=.FALSE.
      LrstGST=.FALSE.
      DO ng=1,Ngrids
        PerfectRST(ng)=.FALSE.
        Ladjusted(ng)=.FALSE.
        LdefADJ(ng)=.FALSE.
        LdefAVG(ng)=.TRUE.
        LdefDIA(ng)=.TRUE.
        LdefERR(ng)=.FALSE.
        LdefHIS(ng)=.TRUE.
        LdefINI(ng)=.FALSE.
        LdefIRP(ng)=.FALSE.
        LdefITL(ng)=.FALSE.
        LdefMOD(ng)=.FALSE.
        LdefRST(ng)=.TRUE.
        LdefSTA(ng)=.TRUE.
        LdefTLM(ng)=.FALSE.
        LwrtADJ(ng)=.FALSE.
        LwrtAVG(ng)=.FALSE.
        LwrtDIA(ng)=.FALSE.
        LwrtHIS(ng)=.FALSE.
        LwrtPER(ng)=.FALSE.
        LwrtRST(ng)=.FALSE.
        LwrtTLM(ng)=.FALSE.
        LwrtInfo(ng)=.TRUE.
        LwrtState2d(ng)=.FALSE.
        LwrtTime(ng)=.TRUE.
        ldefout(ng)=.FALSE.
        synchro_flag(ng)=.FALSE.
      END DO
      RETURN
      END SUBROUTINE initialize_scalars
      END MODULE mod_scalars
