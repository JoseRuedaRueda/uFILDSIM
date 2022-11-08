! -----------------------------------------------------------------------------
!  TRACKER MAIN MODULE
! -----------------------------------------------------------------------------
! TITLE         : SINPA
! PROJECT       : Synthetic INPA
! MODULE        : Auxiliar routines
! AFFILIATION   : University of Sevilla
!> \author Jose Rueda - Universidad de Sevilla
!> \date 21/05/2022
!> \version 2.2
!> \see https://gitlab.mpcdf.mpg.de/poyo/fosd
!
! DESCRIPTION:
!> \brief Contains auxiliary routines for SINPA
! Full-orbit is implemented via Boris leap-frog algorithm (taken from iHIBPsim)
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
! Notes:
!  - File units:
!       - 60: Temporal files: Inputs/output. Open, read(write) and close it
module sinpa_module

  implicit none
  !----------------------------------------------------------------------------
  ! PARAMETERS
  !----------------------------------------------------------------------------
  integer, parameter:: versionID1 = 4  !< ID version number, to identify ouput
  integer, parameter:: versionID2 = 0  !< ID version 2, to identify the ouput
  real (8), parameter:: pi = 3.141592653589793 !< pi
  real (8), parameter:: twopi = 6.283185307
  real (8), parameter:: amu_to_kg = 1.66054e-27 !< amu to kg
  real (8), parameter:: qe = 1.60217662e-19 !< Electron charge C
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! STRUCTURES
  !----------------------------------------------------------------------------
  ! ABSTRACT CLASS FOR THE GEOMETRY ELEMENTS
  !> \brief Class with plate geometry.
  type :: geom_element
      character (len=50) :: name   !< Name of the element
      integer :: triangleNum       !< Number of triangles
      integer :: kind              !< Kind of plate, 0=coll, 1=foil, 2=scint
      real(8), allocatable :: triangles(:,:,:)  !< Triangle data
  end type geom_element

  ! ABSTRACT CLASS FOR THE MC MARKERS
  ! >\brief Class with the marker information
  type :: marker
      !------------------------------------------------------------------------
      ! --- Particle structure
      !------------------------------------------------------------------------
      integer :: n_t    !< Number of integration steps
      real (8), dimension(:,:), allocatable:: position !< Particle position in cm
      real (8), dimension(:,:), allocatable:: velocity !< Particle velocity in cm/s
      real (8):: energy0 !< particle energy at the entrance (keV)
      real (8):: weight  !< Weight of the markers, [part/s of pinhole, for FIDASIM, in au, for mapping]
      real (8):: type    !< Marker type: identify if is active or passive, follows FIDASIM criteria
      logical :: collision = .FALSE. !< Flag for the collision
      real (8), dimension(3) :: collision_point !< Point of collision in m
      real (8), dimension(3) :: ionization_point !< Point of ionization in m
      real (8):: qm !< q over m, SI
      ! For the mapping
      real (8):: xi  !< xi angle (or pitch for FILD), radians
      real (8):: rl !< gyroradius of the particle [cm]
      real (8):: beta  !< beta angle (or Gyrophase for FILD), radians
      ! Integration
      real (8):: dt !< dt for the integration [s]
      integer::kindOfCollision ! With whom has collide
      real (8):: cosalpha_foil ! cos alpha inpinging in the INPA foil

  end type marker

  ! ABSTRACT CLASS FOR THE FIELDS
  ! >/briefg FIDASIM markers
  type :: FIDASIM
        integer(4):: shot_number !< shot number of the performed sim
        real(4) :: time !< time point of the simulation
        integer(4) :: nr_npa !< number of NPA markers
        integer(4) :: counter !< number of markers hitting the detector
        real(8), dimension(:, :), allocatable :: ipos !< birth position
        real(8), dimension(:, :), allocatable :: fpos !< pinhole position
        real(8), dimension(:, :), allocatable :: v !< marker velocity
        real(8), dimension(:), allocatable :: wght !< weight
        integer, dimension(:), allocatable :: kind !< kind, active, passive...
  end type FIDASIM

  ! type :: scatteringData
  !   ! Sigma for the model of scattering in the carbon foil
  !   integer :: na                 ! Number of angles
  !   integer :: ne                 ! Number of energies
  !   real(8):: amin     ! Angles array
  !   real(8):: da       !
  !   real(8):: emin     ! Energies array
  !   real(8):: de       !
  !   real(8), dimension(:,:), allocatable :: sigma      ! Sigma matrix
  !   ! Coefficients of the fit for the energy. See the wiki or the carbon foil
  !   ! routines in section 4
  !   real(8), dimension(:,:), allocatable :: p1_ener,p2_ener
  !   real(8):: depth
  ! end type scatteringData
  !

  ! ABSTRACT CLASS FOR THE NBI
  ! >/briefg NBI geometry
  type :: NBI_class
    character (len=20) :: name  !< NBI name
    real(8), dimension(3) ::  p0  !< Initial point of the NBI
    real(8), dimension(3) ::  u   !< Velocity directoin of the NBI
    ! real(8), dimension(3) ::  d   !< distance between points (cm)
    ! integer:: npoints  !< Number of poins along the NBI line
    ! real(8), dimension(:,:), allocatable :: points !< 3d coordinates of the points
    ! real(8), dimension(:), allocatable :: Rpoints !< R coordinate of the points
    ! real(8), dimension(:), allocatable :: Phipoints !< Phi coordinate of the points
  end type NBI_class

  ! ABSTRACT CLASS FOR THE NBI
  !> \brief Contains the data of 1D direction in a compact type.
  type:: grid_class
    real(8), allocatable::data(:)   !< Contains the points of the grid.
    real(8)::x0                     !< Left bound of the grid.
    real(8)::x1                     !< Right bound of the grid.
    integer::size                            !< Number of grid points.
    real(8)::dx                     !< Difference between consecutive points.
  end type grid_class

  ! ABSTRACT CLASS FOR THE Pinhole reference system
  !> \brief Contains the reference system in the pinhole
  type :: pinhole_system_type
    integer:: kind  !< kind of pinhole (circle, 0) (rectangle, 1)
    real(8):: d1 ! Pinhole radius, or size along u1 (in m)
    real(8):: d2 ! Pinhole radius, or size along u2 (in m)
    real(8), dimension(3):: r  !< Position of the pinhole (m)
    real(8), dimension(3):: u1  !< Vector in the pinhole plane 1
    real(8), dimension(3):: u2  !< vector in the pinhole plane 2
    real(8), dimension(3):: u3  !< normal to the pinhole

    real(8), dimension(3):: e1  !< Normal to the B field at pinhole 1
    real(8), dimension(3):: e2  !< Normal to the B field at pinhole 2
    real(8), dimension(3):: e3  !< Parallel to the B field at pinhole

    real(8):: eta  !< angle (rad) between u3 and B
  end type pinhole_system_type
  ! ---------------------------------------------------------------------------
  ! Interfaces
  ! ---------------------------------------------------------------------------
  abstract interface
    ! Pointer for the field interpolator
    subroutine interpolateField_interface(rq, zq, phiq, tq, out)
      implicit none
      real(kind=8), intent(in)::rq, zq, phiq, tq! Query points in all directions.
      real(kind=8), intent(out)::out(6) ! Br, Bz, Bphi
    end subroutine interpolateField_interface
    ! pointer for the coordinate transform
    subroutine coordinatesTransform(cart, pol)
      implicit none
      real(kind=8), intent(in)::cart(3) ! input coordinates
      real(kind=8), intent(out)::pol(3) ! Output Coordinates
    end subroutine coordinatesTransform
    ! pointer for the data storage
    subroutine savePartToStrikeMatrix(particle, k)
      import marker
      type(marker), intent(in)::particle ! marker to be saved
      integer, intent(in)::k     ! place to be storaged in the matrix
    end subroutine savePartToStrikeMatrix
    subroutine savePartToStrikeMatrixSignal(particle, k)
      import marker
      type(marker), intent(in)::particle ! marker to be saved
      integer, intent(in)::k     ! place to be storaged in the matrix
    end subroutine savePartToStrikeMatrixSignal
  end interface
  procedure(interpolateField_interface),pointer, public::getField
  procedure(coordinatesTransform),pointer, public::cart2pol
  procedure(savePartToStrikeMatrix),pointer, public::savePartToStrike
  procedure(savePartToStrikeMatrixSignal),pointer, public::savePartToStrikeSignal
  
  ! -------------------------------------------------------------------------
  ! Variables
  ! -------------------------------------------------------------------------
  ! --- Mapping
  type(NBI_class) :: nbi
  real(8), dimension(:, :), allocatable:: closestPoints    !< closer points to NBI
  real(8), dimension(:), allocatable:: dPoints    !< distance to NBI
  type(marker) :: part, backPart, tempPart
  real(8) :: dt !< dt for the time integration, in s
  real(8) :: dt1
  real(8) :: OmegaPart !< Gyrofrequency of mapping particles
  real(8), dimension(:, :), allocatable:: Strike !< Matrix to store data
  real(8), dimension(:, :), allocatable:: StrikeMap !< Matrix to store the strike map
  real(8), dimension(:, :), allocatable:: CollimatorStrikes !< Matrix to store the strike map
  real(8), dimension(:, :), allocatable:: backCollimatorStrikes !< Matrix to store the strike map
  real(8), dimension(:, :), allocatable:: WrongMarkers !< Matrix to store the strike map
  real(8), dimension(3,3):: rotation !< Rotation matrix
  real(8), dimension(3):: ps !< reference point in the scintillator
  real(8), dimension(3):: ScintNormal !< Normal to the scintillator

  ! --- Magnetic field
  real(8), dimension(:, :, :, :), allocatable, protected::Brfield    !< Radial magnetic field.
  real(8), dimension(:, :, :, :), allocatable, protected::Bzfield    !< Vertical magnetic field.
  real(8), dimension(:, :, :, :), allocatable, protected::Bphifield  !< Toroidal magnetic field.
  real(8), dimension(:, :, :, :), allocatable, protected::Erfield    !< Radial electric field [V/cm].
  real(8), dimension(:, :, :, :), allocatable, protected::Ezfield    !< Vertical electric field [V/cm].
  real(8), dimension(:, :, :, :), allocatable, protected::Ephifield  !< Toroidal electric field [V/cm].
  ! --- grids
  type(grid_class), protected:: grr  !< Grid classes to store the R grid
  type(grid_class), protected:: gzz  !< Grid classes to store the Z grid
  type(grid_class), protected:: gphi !< Grid classes to store the phi grid
  type(grid_class), protected:: gtt  !< Grid classes to store the t grid

  ! --- Geometry structures
  type(geom_element), dimension (:), allocatable:: geometry
  type(geom_element) :: scintillator, foil, collimator
  ! --- FileNames and paths
  character (len=1000) :: BFieldFile !< Magnetic field file
  character (len=1000) :: CollimatorFile !< Collimator file
  character (len=1000) :: FoilFile !< Foil file
  character (len=1000) :: ScintillatorFile !< Scintillator file
  character (len=1000) :: input_filename !< Scintillator file
  character (len=1000) :: runFolder !< root directory of SINPA
  character (len=1000) :: FIDASIMfolder = ''!< root directory of SINPA

  ! --- Counters
  integer:: cCollimator !< Number of markers impinging the collimator
  integer:: backCollimator !< Number of markers impinging the collimator
  integer:: cWrongNeutral !< Number of markers not colliding with the collimator neither carbon foil
  integer:: cWrongIons !< Number of markers not colliding with the scintillator
  integer:: cWrong !< Number of ions not colliding with anything
  integer:: cInpingBack !< Number of ions colliding in the back part of the plate
  integer:: cSelfshadowed !< Number of markers shadowed

  integer:: cScintillator !< Number of markers colliding with the scintillator, for each alpha
  integer:: cFoil !< Number of markers colliding with the foil, for each alpha
  integer:: cOrb = 0!< Number of saved orbits
  ! --- Signal
  type(FIDASIM):: F4Markers !< FIDASIM4 markers
  real(8):: normalization_resample

  ! --- Others
  integer:: ierr !< error management
  integer:: nToLaunch !< number of markers to ba launched
  integer:: irl, iXI, imc, istep, iistep, iiistep, i, j!< dummy for loops
  real(8), dimension(3) :: rPinCyl  !< position of the pingole (r,p,z)
  real(8), dimension(3) :: Bpinhole  !< B field at pinhole, vector
  real(8) :: BpinholeMod !< B field at pinhole
  real(8) :: v0  !< initial velocity for mapping markers
  real(8), dimension(3) :: r0  !< initial position for mapping markers
  real(8), dimension(3) :: vv0  !< initial velocity for mapping markers
  real(8), dimension(3) :: ran  !< Set of random numbers
  real(8), dimension(3) :: pos1  !< Set of random numbers
  real(8), dimension(:, :), allocatable:: random_numbers  !< Set of random numbers
  real(8) :: rand_number  !< Set of random numbers
  real(8) :: beta  !< random angle for the markers
  real(8) :: dMin  !< minimum distance to NBI of markers trajectories
  real(8), dimension(3) :: posMin  !< position to minimum distance to NBI of the mapping marker
  type(pinhole_system_type) :: pinhole !< pinhole definition
  integer:: dummy_shape !< to get the size of the strike object
  real(8) :: time_sign !< sign for dt (>0 forward modelling, <0 backtrace)
  real(8), dimension(3):: foilNormal
  real(8) :: incidentProjection
  ! --- FILDSIM mode
  real(8), dimension(:), allocatable::min_beta, delta_beta, marker_factor

  ! --- Timing
  real(8) :: t_initial_orbits  !< Time were orbits calculation start
  real(8) :: t_final_orbits !< Time were orbits were calculated
  ! ---------------------------------------------------------------------------
  ! NAMELIST
  ! --- Configuration file
  ! Dummy variables for the namelist
  character (len=50):: runID !< runID
  character (len=1000):: GeomFolder !< geometryID
  logical:: FILDSIMmode = .False. !< Flag to use FILDSIM mode or not
  integer:: nGeomElements !< Number of geometrical elements
  integer:: nxi !< number of pitches (R) to simulate
  integer:: nGyroradius !< number of energies (gyroradius) to simulate
  integer:: nMap = 5000 !< number 0 of markers per energy-pitch for the map
  real(8):: n1 = 0.0d0 !< n1 factor for the calculation of the number of markers
  real(8):: r1 = 0.0d0 !< r1 factor for the calculation of the number of markers
  logical:: mapping = .True.!< flag to decide if we launch mapping markers
  logical:: signal = .False.!< flag to decide if we follow FIDASIM4 markers
  logical:: resampling = .False.!< flag to decide if we resample FIDASIM markers
  integer:: nResampling = 1 !< Resample number for FIDASIM simulations
  logical:: saveOrbits = .False.!< flag to decide if we save the orbits
  real(8):: saveRatio = 0.0d0!< Ratio of orbits to be saved
  logical:: saveOrbitLongMode = .False. !< Flag to save orbits in long or short format
  logical:: verbose = .True. !< Print information of the run
  real(8):: M = 2.0d0!< Mass of the particle, in amu
  real(8):: Zin = 1.0d0!< Initial chargeof the particle, in |e|
  real(8):: Zout = 1.0d0!< Charge after ionization, in |e|
  integer:: IpBt = -1!< Sign of the magnetic field vs the current
  logical:: flag_efield_on = .false.!< include or not electric field
  logical:: save_collimator_strike_points = .false. !< Save the collimator strike points
  logical:: save_wrong_markers_position = .false.  !< Save the end position of the wrong markers
  logical:: save_scintillator_strike_points = .true. !< Save the scintillator strike points
  integer:: kindOfstrikeScintFile = 1 !< Kind of output strike file
  logical:: save_self_shadowing_collimator_strike_points = .false. !< Save the scintillator strike points
  logical:: backtrace = .false.!< flag to trace back the orbits
  logical:: restrict_mode = .false. !< flag to restrict the initial gyrophase
  integer:: FoilElossModel = 0
  integer:: ScintillatorYieldModel = 0
  integer:: FoilYieldModel = 0
  logical:: self_shadowing = .false.

  ! Namelist
  NAMELIST /config/ runID, GeomFolder, nxi, &
    nGyroradius, nMap, n1, r1, mapping, FILDSIMmode, &
    signal, resampling, nResampling, saveOrbits, saveRatio,saveOrbitLongMode, runFolder,&
    FIDASIMfolder, verbose, M, Zin, Zout, IpBt, flag_efield_on, save_collimator_strike_points,&
    backtrace,restrict_mode, FoilElossModel, ScintillatorYieldModel, FoilYieldModel, &
    save_wrong_markers_position, save_scintillator_strike_points, &
    kindOfstrikeScintFile, self_shadowing, &
    save_self_shadowing_collimator_strike_points

  ! --- Input
  integer:: nGyro = 300 !< number of points in a complete gyrocicle of the particle

  real(8) :: minAngle = -1.8 !< collimator vertical acceptance
  real(8) :: dAngle = 0.4 !< collimator vertical acceptance

  real(8), dimension(:), allocatable:: XI, XI_input
  real(8), dimension(:), allocatable:: rL !< FILDSIM gyroradius
  real (8):: maxT = 0.0000001!< maximum time to follow particles
  NAMELIST /inputParams/ nGyro, minAngle, dAngle, XI, rL, maxT

  ! --- Mapping orientation
  real(8) :: pinRadius  !< Radius of the head
  real(8), dimension(3) :: rPin  !< poisition of the pinhole (x,y,z)
  real(8), dimension(3) :: u3  !< normal in the pinhole
  real(8), dimension(3) :: u1  !< tangent to the pinhole
  real(8), dimension(3) :: u2  !< second tangent in the pinhole

  ! --- Particle and foil parameters
  real(8), dimension(:), allocatable:: FoilElossParameters  !< Foil parameters
  real(8), dimension(:), allocatable:: ScintillatorYieldParameters  !< Scintillator parameters
  real(8), dimension(:), allocatable:: FoilYieldParameters  !< Scintillator parameters
  ! Namelist
  NAMELIST /markerinteractions/ FoilElossParameters, ScintillatorYieldParameters, FoilYieldParameters

  ! --- NBI namelist
  ! Dummy variables for namelist
  real(8), dimension(3):: p0, u
  ! Namelist
  NAMELIST /nbi_namelist/ p0, u
contains

!-------------------------------------------------------------------------------
! SECTION 1: AUXILIARY FUNCTIONS
!-------------------------------------------------------------------------------
  subroutine module_v(A,q)
    !---------------------------------------------------------------------------
    ! This subroutine calculates the module of a vector using euclidean metric
    ! INPUTS:
    !       - A: Vector of N components
    ! OUTPUTS:
    !       - q: sqrt(sum of square elements of A)
    !---------------------------------------------------------------------------
    implicit none
    ! Dummy arguments
    real (8), dimension(3),intent(in):: A
    real (8),intent(out):: q
    ! Core
    q=sqrt(dot_product(A,A))
  end subroutine module_v


  subroutine vec_product(a,b,cross)
    !---------------------------------------------------------------------------
    ! This subroutine calculates the cross product of two 3-vectors. The product
    ! is computed following the rules of a right handled axis system
    ! INPUTS:
    !       - a: Vector of 3 components
    !       - b: Vector of 3 components
    ! OUTPUTS:
    !       - cross: a x b
    !---------------------------------------------------------------------------
    real (8), dimension(3), intent(in)  :: a,b
    real (8), dimension(3), intent(out) :: cross

    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end subroutine vec_product


  function cross(a,b) !*
    !Dummy varialbes
    real (8), dimension(3), intent(in) :: a,b
    !Local variables
    real (8), dimension(3):: cross
    !Core
      cross(1)=a(2)*b(3)-a(3)*b(2)
      cross(2)=a(3)*b(1)-a(1)*b(3)
      cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross


  subroutine rand_gauss(s,r)
    !---------------------------------------------------------------------------
    ! This subroutine generates a random number according to a Gaussian
    ! distribution of standard deviation equal to s. A rejection method is used
    !
    ! INPUTS:
    !   - s: The sigma of the standard deviation
    ! OUTPUT:
    !   - r: The random number
    !---------------------------------------------------------------------------
    ! Dummy variables
    real (8), intent(in) :: s
    real (8), intent(out):: r
    ! Local variables
    real (8) :: x

    15 call random_number(x)
    x = 5*(-1+2*x)
    call random_number(r)
    if (r .le. exp(-x**2/2)) then
      r = s*x
    else
      go to 15
    endif
  end subroutine rand_gauss


  subroutine linspace2(x0, x1, N, out)
    ! -----------------------------------------------------------------------
    ! GENERATES A LINEARLY SPACED GRID USING TYPE 'BASE'
    !> \brief Generates a linearly space grid in an type(base) using the
    !! interval [x0, x1] with N-points.
    !
    ! Written by Pablo Oyola for the iHIBPsim code
    ! -----------------------------------------------------------------------
    implicit none
    ! Dummy variables
    real(8), intent(in)::x0 !< Starting point of the grid.
    real(8), intent(in)::x1 !< Ending point of the grid.
    integer, intent(in):: N  !< Number of points in the grid.
    class(grid_class), intent(inout)::out !< Class for the grid data.
    ! local variables
    integer::i
    real(8)::step


    out%x0 = x0
    out%x1 = x1
    out%size = N

    allocate(out%data(N))

    if ((x1-x0) .eq. 0.0d0) then
      return
    end if

    if ((N .eq. 1) .or. (N .eq. 0)) then
      return
    end if

    step = (x1-x0)/(dble(N)-1.0d0)

    do i = 1, N
      out%data(i) = x0 + (dble(i)-1.0d0)*step
    end do

    out%dx = step

    ! Avoid problems with non properly initialized grids.
    if((out%dx .eq. 0.0d0) .or. (N .le. 1))then
      out%dx = 1.0d0
    end if
  end subroutine linspace2


  subroutine rl2v0(rLlocal, Mlocal, Z, B, v0local)
    ! -----------------------------------------------------------------------
    ! Give the velocity module for a given larmor radius
    !
    !
    ! Written by Jose Rueda Rueda
    ! -----------------------------------------------------------------------
    implicit none
    ! Dummy variables
    real(8), intent(in):: rLlocal !< Lamor radius
    real(8), intent(in):: Mlocal !< Mass in amu
    real(8), intent(in):: Z !< charge in electron charge (abs)
    real(8), intent(in):: B !< Modulus of the field
    real(8), intent(out):: v0local !< Modulus of the field

    ! Factor 100.0 is included because the rl is inserted in cm
    v0local = rLlocal * abs(Z) * B / Mlocal * qe / amu_to_kg / 100.0d0
  end subroutine rl2v0


  subroutine omega(M, Z, B, w)
    ! -------------------------------------------------------------------------
    ! Calculate the girofrequency of a particle
    ! -------------------------------------------------------------------------
    ! Dummy variables
    implicit none
    real(8), intent(in) :: M !< Mass in amu
    real(8), intent(in) :: Z !< Charge in abs electron charges
    real(8), intent(in) :: B !< Magnetic field in T
    real(8), intent(out) :: w !< frequency in s-1

    w = Z * qe * B / M / amu_to_kg
  end subroutine omega


  subroutine cart2pol_0_2pi(cart, pol)
    ! Created by Pablo oyola, phi is forced to be between 0 and 2pi
    implicit none

    real(8), intent(in)::cart(3)   !< Cartesian points (x, y, z)
    real(8), intent(out)::pol(3)   !< Polar points (r, z, phi)

    pol(1) = sqrt(cart(1)**2+cart(2)**2)
    pol(2) = cart(3)
    pol(3) = atan2(cart(2), cart(1))

    ! Here it is important to move the phi angle between [0, 2pi]
    if(pol(3) .lt. 0.0d0)then
    pol(3) = twopi + pol(3)
    end if

    if(pol(3) .gt. twopi)then
    pol(3) = pol(3)-twopi
    end if
  end subroutine cart2pol_0_2pi

  subroutine cart2pol_m_pi_pi(cart, pol)
    ! Created by Jose Rueda to eliminate ifs in Pablos versions
    ! Notice that this consider phi always between -pi and pi
    implicit none

    real(8), intent(in)::cart(3)   !< Cartesian points (x, y, z)
    real(8), intent(out)::pol(3)   !< Polar points (r, z, phi)

    pol(1) = sqrt(cart(1)**2+cart(2)**2)
    pol(2) = cart(3)
    pol(3) = atan2(cart(2), cart(1))
  end subroutine cart2pol_m_pi_pi


  subroutine pol2cart(polar, cart)
    ! Createed by Pablo Oyola
    implicit none

    real(8), intent(out)::cart(3)   !< Cartesian points (x, y, z)
    real(8), intent(in)::polar(3)   !< Polar points (r, z, phi)

    cart(1) = polar(1)*cos(polar(3))
    cart(2) = polar(1)*sin(polar(3))
    cart(3) = polar(2)
  end subroutine pol2cart


  subroutine pol2cart_cov(vpol, vcart, r0)
    ! -----------------------------------------------------------------
    ! POLAR TO CARTESIAN COORDINATES (COVARIANT VERSION)
    !> \brief Covariant transformation from polar coordinates to cartesian coordinates
    !! (i.e. velocities).
    ! Created by Pablo Oyola
    ! -----------------------------------------------------------------
    implicit none

    real(8), intent(out)::vcart(3)  !< Cartesian components (Ax, Ay, Az)
    real(8), intent(in)::vpol(3)    !< Polar components (Ar, Aphi, Az)
    real(8), intent(in)::r0(3)      !< Point where the field point is (r, z, phi)

    vcart(1) = vpol(1)*cos(r0(3)) - vpol(2)*sin(r0(3))
    vcart(2) = vpol(1)*sin(r0(3)) + vpol(2)*cos(r0(3))
    vcart(3) = vpol(3)
  end subroutine pol2cart_cov

! -----------------------------------------------------------------

  subroutine quicksort_nr(array)
    ! QUICK sort of an array
    !> \brief Quick sort of array (for 1D array)
    ! Taken from https://www.mjr19.org.uk/IT/sorts/
    ! -----------------------------------------------------------------
    ! This version maintains its own stack, to avoid needing to call
    ! itself recursively. By always pushing the larger "half" to the
    ! stack, and moving directly to calculate the smaller "half",
    ! it can guarantee that the stack needs no more than log_2(N)
    ! entries
    implicit none
    real(8), intent(inout)::array(:)
    real(8) :: temp,pivot
    integer :: i, j, left,right,low,high
    ! If your compiler lacks storage_size(), replace
    ! storage_size(i) by 64
    integer :: stack(2,storage_size(i)),stack_ptr

    low=1
    high=size(array)
    stack_ptr=1

    do
       if (high-low.lt.50) then ! use insertion sort on small arrays
          do i=low+1,high
             temp=array(i)
             do j=i-1,low,-1
                if (array(j).le.temp) exit
                array(j+1)=array(j)
             enddo
             array(j+1)=temp
          enddo
          ! now pop from stack
          if (stack_ptr.eq.1) return
          stack_ptr=stack_ptr-1
          low=stack(1,stack_ptr)
          high=stack(2,stack_ptr)
          cycle
       endif

       ! find median of three pivot
       ! and place sentinels at first and last elements
       temp=array((low+high)/2)
       array((low+high)/2)=array(low+1)
       if (temp.gt.array(high)) then
          array(low+1)=array(high)
          array(high)=temp
       else
          array(low+1)=temp
       endif
       if (array(low).gt.array(high)) then
          temp=array(low)
          array(low)=array(high)
          array(high)=temp
       endif
       if (array(low).gt.array(low+1)) then
          temp=array(low)
          array(low)=array(low+1)
          array(low+1)=temp
       endif
       pivot=array(low+1)

       left=low+2
       right=high-1
       do
          do while(array(left).lt.pivot)
             left=left+1
          enddo
          do while(array(right).gt.pivot)
             right=right-1
          enddo
          if (left.ge.right) exit
          temp=array(left)
          array(left)=array(right)
          array(right)=temp
          left=left+1
          right=right-1
       enddo
       if (left.eq.right) left=left+1
       !          call quicksort(array(1:left-1))
       !          call quicksort(array(left:))
       if (left.lt.(low+high)/2) then
          stack(1,stack_ptr)=left
          stack(2,stack_ptr)=high
          stack_ptr=stack_ptr+1
          high=left-1
       else
          stack(1,stack_ptr)=low
          stack(2,stack_ptr)=left-1
          stack_ptr=stack_ptr+1
          low=left
       endif

    enddo
  end subroutine quicksort_nr


!-------------------------------------------------------------------------------
! SECTION 2: NBI
!-------------------------------------------------------------------------------
  subroutine minimumDistanceLines(p, v, r, w, d_local, point)
    ! -------------------------------------------------------------------------
    ! Get the minimum distance between two lines
    !
    ! Jose Rueda: jrrueda@us.es
    !
    ! Note: The case of parallel lines will not be consider
    !
    ! -------------------------------------------------------------------------

    ! Dummy variables
    real(8), dimension(3), intent(in) :: p !< Point in line 1 (x, y, z)
    real(8), dimension(3), intent(in) :: v !< unit vector line 1
    real(8), dimension(3), intent(in) :: r !< Point in line 2 (x, y, z)
    real(8), dimension(3), intent(in) :: w !< unit vector line 2

    real(8), intent(out) :: d_local !< minimum distace
    real(8), dimension(3), intent(out) :: point !< Point in line 1 closer to line 2

    ! Auxiliar variables
    real(8):: alpha, beta, gamma   ! Please see SINPA documentation
    real(8), dimension(3):: d0 ! Point in line 1 closer to line 2
    ! Calculate auxiliar quantities
    alpha = sum((p - r) * w)
    beta = sum((p - r) * v)
    gamma = sum(v * w)
    d0 = (p - r) + ((beta - alpha*gamma)*v &
      + (alpha - gamma*beta)*w)/(gamma**2-1)
    ! Calculate the minimum distance
    d_local = sqrt(d0(1)**2 + d0(2)**2 + d0(3)**2)
    point = p + (beta - gamma * alpha)/(gamma**2-1)*v
  end subroutine minimumDistanceLines


!------------------------------------------------------------------------------
! SECTION 3: COEFFICIENTS FOR INTERPOLATORS
!------------------------------------------------------------------------------
  subroutine interpolate2D_coefficients(x1, x2, x1q, x2q, aaa, idx)
    ! -----------------------------------------------------------------------
    ! OBTAIN THE NEAREST-NEIGHBOURS COEFFICIENTS FOR 2D INTERPOLATION
    !> \brief Computes, for a two given type(base), the nearest neighbours and the
    !! normalized distance to all of them (normalized to the interval size in each direction).
    !> \detail These values are used for interpolation in 2D.
    ! -----------------------------------------------------------------------
    implicit none

    class(grid_class), intent(in)::x1 !< Grid class for each direction.
    class(grid_class), intent(in)::x2
    real(8), intent(in)::x1q, x2q !< Query points in each direction.
    real(8), intent(out)::aaa(4) !< Contains the interpolating coefficients.
    integer, intent(out)::idx(4)  !< Nearest points in the grid
                                  !! (N, N+1)->closer to the left and right respectively.
    integer::ia, ia1, ja, ja1
    real(8)::ar, az, ar1, az1

    ia  = max(1, min(x1%size-1, int((x1q - x1%x0)/x1%dx  + 1)))
    ia1 = ia + 1
    ja  = max(1, min(x2%size-1, int((x2q - x2%x0)/x2%dx  + 1)))
    ja1 = ja + 1

    ar1   = max(0.0d0, min(1.0d0, (x1q-x1%data(ia))/x1%dx))
    ar    = 1.0d0 - ar1
    az1   = max(0.0d0, min(1.0d0, (x2q-x2%data(ja))/x2%dx))
    az    = 1.0d0 - az1

    ! Interpolating coefficients.
    aaa(1) = ar *az
    aaa(2) = ar1*az
    aaa(3) = ar *az1
    aaa(4) = ar1*az1

    ! Interpolation indexes.
    idx(1) = ia
    idx(2) = ia1
    idx(3) = ja
    idx(4) = ja1
  end subroutine interpolate2D_coefficients


  subroutine interpolate3D_coefficients(x1, x2, x3, x1q, x2q, x3q, aaa, idx)
    ! -----------------------------------------------------------------------
    ! OBTAIN THE NEAREST-NEIGHBOURS COEFFICIENTS FOR 3D INTERPOLATION
    !> \brief Computes, for a three given type(base), the nearest neighbours and the
    !! normalized distance to all of them (normalized to the interval size in each direction).
    !> \detail These values are used for interpolation in 3D.
    ! -----------------------------------------------------------------------
   implicit none

   class(grid_class), intent(in)::x1 !< Grid class for each direction.
   class(grid_class), intent(in)::x2
   class(grid_class), intent(in)::x3
   real(8), intent(in)::x1q, x2q, x3q !< Query points in each direction.
   real(8), intent(out)::aaa(8) !< Contains the interpolating coefficients.
   integer, intent(out)::idx(6)  !< Nearest points in the grid
                                 !! (N, N+1)->closer to the left and right respectively.
   integer::ia, ia1, ja, ja1, ka, ka1
   real(8)::ax, ay, az, ax1, ay1, az1

   ia  = max(1, min(x1%size-1, int((x1q - x1%x0)/x1%dx  + 1)))
   ia1 = ia + 1
   ja  = max(1, min(x2%size-1, int((x2q - x2%x0)/x2%dx  + 1)))
   ja1 = ja + 1
   ka  = max(1, min(x3%size-1, int((x3q - x3%x0)/x3%dx  + 1)))
   ka1 = ka + 1

   ax1   = max(0.0d0, min(1.0d0, (x1q-x1%data(ia))/x1%dx))
   ax    = 1.0d0 - ax1
   ay1   = max(0.0d0, min(1.0d0, (x2q-x2%data(ja))/x2%dx))
   ay    = 1.0d0 - ay1
   az1   = max(0.0d0, min(1.0d0, (x3q-x3%data(ka))/x3%dx))
   az    = 1.0d0 - az1

   ! Interpolating coefficients.
   aaa(1) = ax  * ay  * az
   aaa(2) = ax  * ay1 * az
   aaa(3) = ax1 * ay  * az
   aaa(4) = ax1 * ay1 * az
   aaa(5) = ax  * ay  * az1
   aaa(6) = ax  * ay1 * az1
   aaa(7) = ax1 * ay  * az1
   aaa(8) = ax1 * ay1 * az1

   ! Interpolation indexes.
   idx(1) = ia
   idx(2) = ia1
   idx(3) = ja
   idx(4) = ja1
   idx(5) = ka
   idx(6) = ka1
  end subroutine interpolate3D_coefficients


!------------------------------------------------------------------------------
! SECTION 4: MAGNETIC FIELD
!------------------------------------------------------------------------------
  subroutine parseField(Bname, Ename, verbose_flag)
    !---------------------------------------------------------------------------
    ! This subroutine load the magnetic field
    !
    ! Adaptation from Pablo Oyola iHIBPsim code
    !
    ! INPUTS:
    !       - Bname: Name of the field to read
    !       - Ename: Name of the field to read
    ! OUTPUTS:
    !       - magnetic field arrays will be filled
    !---------------------------------------------------------------------------
    implicit none
    ! Dummy variables
    character(len=*), intent(in)::Bname   !< Filename with the magnetic field.
    character(len=*), intent(in)::Ename   !< Filename with the magnetic field.
    logical, intent(in) :: verbose_flag  !< flag to print the information or not
    ! Local variables
    integer::lr, lz, lphi, ntot
    integer::lr1, lz1, lphi1, ntot1 ! Auxiliar variables to compare the field sizes.
    real(8)::rmin, rmax, zmin, zmax, phimin, phimax, timemin, timemax
    integer::send_size, ii

    ! Parse de magnetic field
    if (verbose_flag) then
      print*, 'Parsing the magnetic field...'
    endif
    open(unit = 60, file=Bname, access='stream', status='old', action='read')
    read(60) lr, lz, lphi, ntot, &
             rmin, rmax, zmin, zmax, phimin, phimax, timemin, timemax

    if (verbose_flag) then
      print*, 'The header of the magnetic field was parsed:'
      print*, 'RESULTS'
      print*, 'nr = ', lr, ' nz = ', lz, ' nphi = ', lphi, ' ntot = ', ntot
      print*, 'INTEGRATION REGION:'
      print*, 'R = (',rmin, ' , ', rmax, ')'
      print*, 'z = (',zmin, ' , ', zmax, ')'
      print*, 'phi = (',phimin, ' , ', phimax, ')'
      print*, 't = (',timemin, ' , ', timemax, ')'
      print*, 'Allocating variables... '
    endif
    ! With this, we can allocate all the variables.
    allocate(Brfield(lr, lphi, lz, ntot))
    allocate(Bzfield(lr, lphi, lz, ntot))
    allocate(Bphifield(lr, lphi, lz, ntot))
    if (verbose_flag) then
      print*, 'Allocated! Reading the magnetic field'
    endif
    read(60) Brfield
    read(60) Bphifield
    read(60) Bzfield
    close(60)
    if (flag_efield_on) then
      if (verbose_flag) then
        print*, 'reading:', Ename
      endif
      open(unit = 60, file=Ename, access='stream', status='old', action='read')
      read(60) lr1, lz1, lphi1, ntot1, &
               rmin, rmax, zmin, zmax, phimin, phimax, timemin, timemax
      if ((lr1.ne.lr).or.(lz1.ne.lz).or.(lphi1.ne.lphi).or.(ntot.ne.ntot1)) THEN
        print*, 'Different grid size found in E and B Fields'
        close(60)
        stop
      endif
      if (verbose_flag) then
        print*, 'parsing electrinc field'
      endif
      allocate(Erfield(lr, lphi, lz, ntot))
      allocate(Ezfield(lr, lphi, lz, ntot))
      allocate(Ephifield(lr, lphi, lz, ntot))
      read(60) Erfield
      read(60) Ephifield
      read(60) Ezfield
      close(60)
    endif
    ! Once we read everything, we create the grids:
    call linspace2(rmin, rmax, lr, grr)
    call linspace2(zmin, zmax, lz, gzz)

    ! If the phi-direction is neglected (i.e., axisymmetric case)
    if(lphi .eq. 1)then
      allocate(gphi%data(1))
      gphi%data = 0.0d0
      gphi%x0   = 0.0d0
      gphi%x1   = 1.0d0
      gphi%size = 1
      gphi%dx = 1.0d0
    else
      call linspace2(phimin, phimax, lphi, gphi)
    end if
    if(ntot .eq. 1)then
        allocate(gtt%data(1))
        gtt%data = 0.0d0
        gtt%x0   = 0.0d0
        gtt%x1   = 1.0d0
        gtt%size = 1
        gtt%dx = 1.0d0
    else
      call linspace2(timemin, timemax, ntot, gtt)
    end if
    ! Prepare the interfaces for the field
    if(lr .eq. 1) then  ! We have a uniform
      getField => getField_0D
      cart2pol => cart2pol_m_pi_pi
    elseif(lphi.eq.1)then  ! We have an axisymmetric magnetic field
      getField => getField_2D
      cart2pol => cart2pol_m_pi_pi
    else
      getField => getField_3D ! We have an non-axisymmetric magnetic field
      if (phimax > 3.15) then
        print*, 'Warning, this is not optimal, please use -pi, pi grid'
        cart2pol => cart2pol_0_2pi
      else
        cart2pol => cart2pol_m_pi_pi
      endif
    end if
  end subroutine parseField


  subroutine unloadField
    ! Written by Pablo Oyola for the iHIBPsim code
    implicit none

    ! Releasing the memory used by the grid.
    deallocate(grr%data)
    deallocate(gzz%data)
    deallocate(gphi%data)

    ! Releasing the memory used by the fields.
    deallocate(Brfield)
    deallocate(Bzfield)
    deallocate(Bphifield)

    ! Releasing the memory used by the fields.
    deallocate(Erfield)
    deallocate(Ezfield)
    deallocate(Ephifield)


    getField => NULL()
  end subroutine unloadField


  subroutine getField_3D(rq, zq, phiq, tq, out)

    implicit none

    real(8), intent(in)::rq      !< Query point in R major to interpolate.
    real(8), intent(in)::zq      !< Query point in Z to interpolate.
    real(8), intent(in)::phiq    !< Query point in toroidal direciton.
    real(8), intent(in)::tq      !< Time query point.
    real(8), intent(out)::out(6) !< Br, Bz, Bphi, Er, Rz, Ephi
    real(8)::a000, a100, a010, a110, a001, a101, a011, a111
    ! Interpolation coefficients.
    real(8)::aaa(8)
    integer::idx(6), ia, ia1, ja, ja1, ka, ka1

    ! Temporal arrays to transform from cylindrical to Cartesian.
    real(8)::b(3), e(3), r0(3)

    ! First, we have to check if the data point is out of range.
    if((rq .gt. grr%x1) .or.  (zq .gt. gzz%x1) .or. (zq .lt. gzz%x0))then
      out = 0.0d0
      return
    end if

    call interpolate3D_coefficients(grr, gphi, gzz, rq, phiq, zq, aaa, idx)


    ia  = idx(1)
    ia1 = idx(2)
    ja  = idx(3)
    ja1 = idx(4)
    ka  = idx(5)
    ka1 = idx(6)

    a000 = aaa(1)
    a010 = aaa(2)
    a100 = aaa(3)
    a110 = aaa(4)
    a001 = aaa(5)
    a011 = aaa(6)
    a101 = aaa(7)
    a111 = aaa(8)

    b(1) =    Brfield(ia, ja, ka,  1)*a000  + Brfield(ia1, ja, ka,  1)*a100  &
           +  Brfield(ia, ja, ka1, 1)*a001  + Brfield(ia1, ja, ka1, 1)*a101  &
           +  Brfield(ia, ja1, ka,  1)*a010 + Brfield(ia1, ja1, ka,  1)*a110  &
           +  Brfield(ia, ja1, ka1, 1)*a011 + Brfield(ia1, ja1, ka1, 1)*a111

    b(3) =    Bzfield(ia, ja, ka,  1)*a000  + Bzfield(ia1, ja, ka,  1)*a100  &
           +  Bzfield(ia, ja, ka1, 1)*a001  + Bzfield(ia1, ja, ka1, 1)*a101  &
           +  Bzfield(ia, ja1, ka,  1)*a010 + Bzfield(ia1, ja1, ka,  1)*a110  &
           +  Bzfield(ia, ja1, ka1, 1)*a011 + Bzfield(ia1, ja1, ka1, 1)*a111

    b(2) =    Bphifield(ia, ja, ka,  1)*a000  + Bphifield(ia1, ja, ka,  1)*a100  &
           +  Bphifield(ia, ja, ka1, 1)*a001  + Bphifield(ia1, ja, ka1, 1)*a101  &
           +  Bphifield(ia, ja1, ka,  1)*a010 + Bphifield(ia1, ja1, ka,  1)*a110  &
           +  Bphifield(ia, ja1, ka1, 1)*a011 + Bphifield(ia1, ja1, ka1, 1)*a111

    ! Tranforming from cylindrical to Cartesian components.
    r0(3) = phiq
    call pol2cart_cov(b, out(1:3), r0)

    if(flag_efield_on)then
        e(1) =    Erfield(ia, ja, ka,  1)*a000  + Erfield(ia1, ja, ka,  1)*a100  &
               +  Erfield(ia, ja, ka1, 1)*a001  + Erfield(ia1, ja, ka1, 1)*a101  &
               +  Erfield(ia, ja1, ka,  1)*a010 + Erfield(ia1, ja1, ka,  1)*a110  &
               +  Erfield(ia, ja1, ka1, 1)*a011 + Erfield(ia1, ja1, ka1, 1)*a111

        e(3) =    Ezfield(ia, ja, ka,  1)*a000  + Ezfield(ia1, ja, ka,  1)*a100  &
               +  Ezfield(ia, ja, ka1, 1)*a001  + Ezfield(ia1, ja, ka1, 1)*a101  &
               +  Ezfield(ia, ja1, ka,  1)*a010 + Ezfield(ia1, ja1, ka,  1)*a110  &
               +  Ezfield(ia, ja1, ka1, 1)*a011 + Ezfield(ia1, ja1, ka1, 1)*a111

        e(2) =    Ephifield(ia, ja, ka,  1)*a000  + Ephifield(ia1, ja, ka,  1)*a100  &
               +  Ephifield(ia, ja, ka1, 1)*a001  + Ephifield(ia1, ja, ka1, 1)*a101  &
               +  Ephifield(ia, ja1, ka,  1)*a010 + Ephifield(ia1, ja1, ka,  1)*a110  &
               +  Ephifield(ia, ja1, ka1, 1)*a011 + Ephifield(ia1, ja1, ka1, 1)*a111
        ! Tranforming from cylindrical to Cartesian components.
        call pol2cart_cov(e, out(4:6), r0)
    else
      out(4:6) = 0.0d0
    end if
  end subroutine getField_3D


  subroutine getField_2D(rq, zq, phiq, tq, out)
    ! -----------------------------------------------------------------------
    ! GET THE 2D FIELDS INTERPOLATION.
    !> \brief Obtains Br, Bz, Bphi, Er, Ez, Ephi interpolated in a given point,
    !! (R, z). Linear interpolation is used.
    ! Written by Pablo Oyola for the iHIBPsim code
    ! -----------------------------------------------------------------------
    implicit none

    real(8), intent(in)::rq      !< Query point in R major to interpolate.
    real(8), intent(in)::zq      !< Query point in Z to interpolate.
    real(8), intent(in)::phiq    !< Query point in toroidal direciton.
    real(8), intent(in)::tq      !< Time query point.
    real(8), intent(out)::out(6) !< Br, Bz, Bphi

    ! Interpolation coefficients.
    real(8)::aaa(4)
    integer::idx(4), ia, ia1, ja, ja1
    real(8), dimension(3):: dum

    ! First, we have to check if the data point is out of range.
    if((rq .gt. grr%x1) .or.  (zq .gt. gzz%x1) .or. (zq .lt. gzz%x0) .or. (rq .lt. grr%x0))then
      out(1:6) = 0.0d0
      return
    end if

    call interpolate2D_coefficients(grr, gzz, rq, zq, aaa, idx)
    ia  = idx(1)
    ia1 = idx(2)
    ja  = idx(3)
    ja1 = idx(4)
    dum(1) =  Brfield(ia, 1, ja, 1)*aaa(1)   + Brfield(ia1, 1, ja, 1)*aaa(2) &
            + Brfield(ia, 1, ja1, 1)*aaa(3)  + Brfield(ia1, 1, ja1, 1)*aaa(4)

    dum(2) =  Bzfield(ia, 1, ja, 1)*aaa(1)   + Bzfield(ia1, 1, ja, 1)*aaa(2) &
            + Bzfield(ia, 1, ja1, 1)*aaa(3)  + Bzfield(ia1, 1, ja1, 1)*aaa(4)

    dum(3) =  Bphifield(ia, 1, ja, 1)*aaa(1)   + Bphifield(ia1, 1, ja, 1)*aaa(2) &
            + Bphifield(ia, 1, ja1, 1)*aaa(3)  + Bphifield(ia1, 1, ja1, 1)*aaa(4)
    call pol2cart_cov((/dum(1),dum(3),dum(2)/), out(1:3), (/rq, zq, phiq/))

    if(flag_efield_on)then
      dum(1) =  Erfield(ia, 1, ja, 1)*aaa(1)   + Erfield(ia1, 1, ja, 1)*aaa(2) &
              + Erfield(ia, 1, ja1, 1)*aaa(3)  + Erfield(ia1, 1, ja1, 1)*aaa(4)

      dum(2) =  Ezfield(ia, 1, ja, 1)*aaa(1)   + Ezfield(ia1, 1, ja, 1)*aaa(2) &
              + Ezfield(ia, 1, ja1, 1)*aaa(3)  + Ezfield(ia1, 1, ja1, 1)*aaa(4)

      dum(3) =  Ephifield(ia, 1, ja, 1)*aaa(1)   + Ephifield(ia1, 1, ja, 1)*aaa(2) &
              + Ephifield(ia, 1, ja1, 1)*aaa(3)  + Ephifield(ia1, 1, ja1, 1)*aaa(4)
      call pol2cart_cov((/dum(1),dum(3),dum(2) /), out(4:6), (/rq, zq, phiq/))
    else
      out(4:6) = 0.0d0
    end if
  end subroutine getField_2D

  subroutine getField_0D(rq, zq, phiq, tq, out)
    ! -----------------------------------------------------------------------
    ! GET THE FIELD of only 0D
    !> \brief Obtains Br, Bz, Bphi, interpolated in a given point,
    !! (R, z).
    !
    ! Note: zero field is suppose to be in cartesian coordinates!!!
    !
    ! Written by Jose Rueda
    ! -----------------------------------------------------------------------
    implicit none

    real(8), intent(in)::rq      !< Query point in R major to interpolate.
    real(8), intent(in)::zq      !< Query point in Z to interpolate.
    real(8), intent(in)::phiq    !< Query point in toroidal direciton.
    real(8), intent(in)::tq    !< Query point in toroidal direciton.
    real(8), intent(out)::out(6) !< Br, Bz, Bphi

    out(1) =  Brfield(1, 1, 1, 1)

    out(2) =  Bphifield(1, 1, 1, 1)

    out(3) =  Bzfield(1, 1, 1, 1)

    if(flag_efield_on)then
      out(4) =  Erfield(1, 1, 1, 1)

      out(5) =  Ephifield(1, 1, 1, 1)

      out(6) =  Ezfield(1, 1, 1, 1)
    else
      out(4:6) = 0.0d0
    end if
  end subroutine getField_0D

  subroutine Bsystem()
    ! -----------------------------------------------------------------------
    ! Prepare the pinhole reference system
    !> \brief Calculate the reference system of the pinhole
    ! Written by Jose Rueda
    ! When executed, it set in the code workspace:
    !   - rPinCyl !< Cylindrical coordinates of the pinhole
    !   - pinhole%e3,2,1 !< reference vectors in the pinhole
    !   - Bpinhole !< Magnetic field at the inhole in cartesian coordiates
    !   - BpinholeMod !< Modulus of the magnetic field at the pinhole
    ! -----------------------------------------------------------------------

    ! --- Local variables
    real(8) :: project !< Projection of the magnetic field in the u2 vector
    real(8) :: e1mod !< Modulus of the e1 vector
    real(8), dimension(6):: field ! Field at the pinhole in cylindrical

    ! Calculate the pinhole position in cylindrical coordinates
    call cart2pol(pinhole%r, rPinCyl)

    ! Get the magnetic field in the pinhole
    ! call getField(rPinCyl(1), rPinCyl(3), rPinCyl(2), 0.5d0,  field_in_cyl)
    ! call pol2cart_cov((/field_in_cyl(1),field_in_cyl(3), field_in_cyl(2) /), Bpinhole, rPinCyl)
    call getField(rPinCyl(1), rPinCyl(2), rPinCyl(3), 0.5d0,  field)
    Bpinhole = field(1:3)
    BpinholeMod = sqrt(Bpinhole(1)**2 + Bpinhole(2)**2 + Bpinhole(3)**2)
    ! calculate the normal vectors
    call vec_product(pinhole%u3,Bpinhole,pinhole%e1)
    e1mod = sqrt(pinhole%e1(1)**2 + pinhole%e1(2)**2 &
                                   + pinhole%e1(3)**2)
    ! See handwritten SINPA notes for full explanation with 3D drawings
    if (e1mod.lt.0.001) then ! B is in line with u2
      print*, Bpinhole
      print*, pinhole%u3
      print*, 'B field normal to the pinhole. Option not considered'
      print*, 'Write to jrrueda@us.es'
      stop
    endif
    ! Now calculate e2:
    pinhole%e3 = Bpinhole / BpinholeMod
    pinhole%e1 = pinhole%e1 / e1mod
    call vec_product(pinhole%e3,pinhole%e1,pinhole%e2)
    ! save the angle eta
    pinhole%eta = acos(pinhole%e3(1)*pinhole%u3(1) &
                       + pinhole%e3(2)*pinhole%u3(2) &
                       + pinhole%e3(3)*pinhole%u3(3))
  end subroutine Bsystem


!-------------------------------------------------------------------------------
! SECTION 5: COLLISION CHECK and INTEGRATOR
!-------------------------------------------------------------------------------
  subroutine triangleRay(num_triangles, triangles, o, p, flag_col, t)
    implicit none
    ! ------------------------------------------------------------------
    ! Main routine to compute the many triangles-single ray collision.
    !> \brief Calculation of the ray-triangle collision.
    !> \detail Compute the collision, if exists, of a certain array
    !! of triangles with the segment defined by positions vectors O and P.
    !! IMPORTANT: inputs and outputs must be in CARTESIAN COORDINATES.
    !
    ! Written by Pablo Oyola, addapted for SINPA by Jose Rueda
    ! ------------------------------------------------------------------
    ! Accuracy of the ray-triangle calculation in the same units 'triangles'
    ! coordinates are given.
    real(8), parameter :: eps = 1.0d-12

    ! Ray coordinates: origin (o) and ending (p) points and collision point (t).
    ! Cartesian coordinates are required
    real(8), intent(in)::o(3)    !< Coordinates of the origin of the ray.
    real(8), intent(in)::p(3)    !< Coordinates of the end of the ray.
    real(8), intent(inout)::t(3) !< If any collision is detected, this contains crossing point.
    integer, intent(in)::num_triangles  !< Number of triangles in the 'triangles' array.

    ! The 'triangle' array must be provided in a special format. Go to \see{parseGeometry}
    ! subroutine for more info.
    ! The last term of the array indicates the coordinate: 1.x | 2. y | 3.z
    real(8), intent(in)::triangles(num_triangles, 3, 3) !< Triangle array provided by \see{parseGeometry}.
    logical, intent(out)::flag_col !< Boolean determining whether there is collision with a triangle.

    ! Variables to check the collision. To see the description of this algorithm,
    ! visit \url{}.
    real(8)::edge0(3), edge1(3), v0(3) ! Vectors to describe the triangle.
    real(8)::dir(3) ! Direction of the segment.
    real(8)::cross_prod1(3), determinant, inv_det
    real(8)::svec(3), qvec(3)
    real(8)::u, v ! Barycentric coordinates.
    real(8)::tcoor    ! Line coordinate.

    ! The counters
    integer:: ii

    ! This will be common for all triangles.
    dir = p - o
    flag_col = .FALSE.

    ! check collisions
    do ii = 1,num_triangles
      ! Computing the edges of the triangle
      edge0(1) = triangles(ii, 1, 1)
      edge0(2) = triangles(ii, 1, 2)
      edge0(3) = triangles(ii, 1, 3)

      edge1(1) = triangles(ii, 2, 1)
      edge1(2) = triangles(ii, 2, 2)
      edge1(3) = triangles(ii, 2, 3)

      v0(1) = triangles(ii, 3, 1)
      v0(2) = triangles(ii, 3, 2)
      v0(3) = triangles(ii, 3, 3)

      ! First of all, we compute the determinant:
      cross_prod1(1) = dir(2)*edge1(3) - dir(3)*edge1(2)
      cross_prod1(2) = dir(3)*edge1(1) - dir(1)*edge1(3)
      cross_prod1(3) = dir(1)*edge1(2) - dir(2)*edge1(1)

      determinant = edge0(1)*cross_prod1(1) + edge0(2)*cross_prod1(2) + edge0(3)*cross_prod1(3)

      ! Let's check if the segment and the triangle are parallel:
      if(ABS(determinant) .lt. eps)then
        cycle
      end if

      inv_det = 1.0d0/determinant

      svec = o - v0

      ! First barycentric coordinate
      u = (svec(1)*cross_prod1(1) + svec(2)*cross_prod1(2) + &
           & svec(3)*cross_prod1(3))*inv_det
      if((u .lt. 0.0d0) .or. (u .gt. 1.0d0))then
        cycle
      end if

      ! Computing the second barycentric coordinate.
      qvec(1) = svec(2)*edge0(3) - svec(3)*edge0(2)
      qvec(2) = svec(3)*edge0(1) - svec(1)*edge0(3)
      qvec(3) = svec(1)*edge0(2) - svec(2)*edge0(1)

      v = (dir(1)*qvec(1) + dir(2)*qvec(2) + dir(3)*qvec(3))*inv_det

      if((v .lt. 0.0d0) .or. (u + v .gt. 1.0d0))then
        cycle
      end if

      ! Finally, we checked that there is an actual collision with the ray
      ! on which our segment is lying. However, we need to verify if the
      ! the collision is actually on the segment itself.
      tcoor = (edge1(1)*qvec(1) + edge1(2)*qvec(2) + edge1(3)*qvec(3))*inv_det
      if((tcoor .ge. 0.0d0) .and. (tcoor .le. 1.0d0))then
        ! We have a collision, let's obtain the coordinates of the collision point.
        flag_col = .TRUE.
        t = o + tcoor*dir
        exit ! Exit the loop.
      end if
    end do
  end subroutine triangleRay

  subroutine checkCollision(particle, kk)
    type(marker), intent(inout):: particle
    integer::iloop, kk
    if (particle%kindOfCollision .eq. 1) then
      ! We already collided with the foil, check just the scintillator, which is
      ! supposed to be the last element of the geometry array
      call triangleRay(geometry(nGeomElements)%triangleNum, geometry(nGeomElements)%triangles, &
                       particle%position(:, kk), particle%position(:, kk + 1), &
                       particle%collision, particle%collision_point)
      if (particle%collision) then
        particle%kindOfCollision = geometry(nGeomElements)%kind
      endif
    else   ! We need to loop over all the plates
      plates: do iloop=1, nGeomElements
        call triangleRay(geometry(iloop)%triangleNum, geometry(iloop)%triangles, &
                         particle%position(:, kk), particle%position(:, kk + 1), &
                         particle%collision, particle%collision_point)
        if (particle%collision) then
          particle%kindOfCollision = geometry(iloop)%kind
          exit plates
        endif
      enddo plates
    endif
  end subroutine checkCollision


  subroutine pushParticle(qm, r0, v0, r1, v1, dt)
    ! ----------------------------------------------------------------------------
    ! BORIS LEAP-FROG INTEGRATOR: evolve the position and velocity of the particle
    ! one step, following the Boris's Leap-frog method
    !
    ! This subroutine was created by Pablo Oyola following the article:
    ! A COMPREHENSIVE COMPARISON OF RELATIVISTIC PARTICLE INTEGRATORS: B. Ripperda
    !
    ! The necessary changes to insert this routine in INPASIM were carried on by
    ! Jose Rueda
    !
    ! This routine allow to have non homogeneous magnetic field
    !
    ! This subroutine will implement a time-step in the evolution of the particle.
    ! The relevant data will be contained within the module data, read at the
    ! beginning of the simulation.
    ! This is a Cartesian integrator!
    ! INPUTS:
    !   -qm: charge/mass ratio in SI units.
    !   -r0 = (x, y, z) initial coordinates [in cm]
    !   -v0 = (vx, vy, vz) initial velocity. [in cm/s]
    !   -E = Electric field. [in V/cm]
    !   -dt = time step of the integration
    !   - The magnetic field is already stored n memory in the structure fields
    ! OUTPUTS
    !   -r1 : final coordinates [in cm]
    !   -v1 : final velocity [in cm/s]
    ! ----------------------------------------------------------------------------

    !Dummy variables
    real (8), intent(in) :: qm  !< Charge-mass ratio
    real (8), intent(in) :: dt  !< Time step for integration.
    real (8), dimension(3), intent(in) :: r0 !< Initial position vector (x, y, z) = Cartesian system.
    real (8), dimension(3), intent(in) :: v0 !< Initial velocity vector (vx, vy, vz) = Cartesian system.
    real (8), dimension(3), intent(out) :: r1 !< Evolved position vector. Also cartesian.
    real (8), dimension(3), intent(out) :: v1 !< Evolved velocity vector. Also cartesian.
    !Local variables
    real (8), dimension(3)::r_plus_half, r_polar
    real (8), dimension(3)::vminus, vplus
    real (8):: mod_t_2
    real (8), dimension(3)::t, s, vprime, vprime2
    real (8), dimension(3):: B, B1, E
    real (8), dimension(6)::field_data ! This contains the interpolated data.
    ! Note, 't' is a really bad notation for an auxiliary vector, but in most books
    ! is used, so I will maintain the name

    ! 1. We push the particle to n -> n+1/2
    r_plus_half = r0 + 0.5d0*v0*dt
    ! 2. We now have the r_{n+1/2}, let's evaluate the fields on that position.
    ! The fields, however, are stored in cylindrical coordinates, so the position
    ! must be firstly translated into that coordinate system.
    call cart2pol(r_plus_half, r_polar)
    call getField(r_polar(1), r_polar(2), r_polar(3), 0.0d0, field_data)
    B = field_data(1:3)
    E = field_data(4:6)
    ! We have the fields in polar coordinates, but the Boris method requires the
    ! cartesian ones:
    ! call pol2cart_cov((/field_data(1),field_data(3), field_data(2) /), B, r_polar)
    ! call pol2cart_cov((/field_data(4),field_data(6), field_data(5) /), E, r_polar)
    ! 3. First half of the acceleration:
    vminus = v0 + 0.5d0*qm*dt*E

    ! 4. Calculation of the auxiliary vectors t and s
    t = 0.50d0*qm *B*dt
    mod_t_2 = t(1)*t(1) + t(2)*t(2) +t(3)*t(3)
    s = 2.0d0*t/(1+mod_t_2)

    ! 5. To get v+, we apply the rotation, step by step:
    !  5.1 vprime = vminus x t
    vprime(1) = vminus(2)*t(3) - vminus(3)*t(2)
    vprime(2) = vminus(3)*t(1) - vminus(1)*t(3)
    vprime(3) = vminus(1)*t(2) - vminus(2)*t(1)

    !  5.2 vprime2 = (vminus x t) + vminus
    vprime2 = vprime + vminus

    !  5.3 vplus = vminus + vprime2 x s
    vplus(1) = vprime2(2)*s(3) - vprime2(3)*s(2) + vminus(1)
    vplus(2) = vprime2(3)*s(1) - vprime2(1)*s(3) + vminus(2)
    vplus(3) = vprime2(1)*s(2) - vprime2(2)*s(1) + vminus(3)

    ! 6. Recover the real velocity at the point: v(n+1) = vplus + q/(2m)*E
    v1 = vplus + 0.50d0*qm*E*dt

    ! 6. The position.
    r1 = r_plus_half + 0.50d0*v1*dt
  end subroutine pushParticle


!-------------------------------------------------------------------------------
! SECTION 6: PLATE SETUP
!-------------------------------------------------------------------------------
  subroutine parseGeometry(filename, verb, g_element)
    ! -----------------------------------------------------------------------
    ! PARSE THE GEOMETRY FILES
    !> \brief Loads the wall file and parse the triangle structure to have a quick
    !! triangle-ray evaluation.
    !> \detail The structure of the file must be first the number of triangles and
    ! then a matrix of 3 columns x #triangles rows. The triangles will be recomposed
    ! by associating three consecutive rows.
    ! TRIANGLE ARRAY = (#triangles, 3, 3)
    ! 1. The first index states which triangle we have.
    ! 2. The last index states which coordinate (x, y, z)
    ! 3. The second index contains the following:
    !      (:, 1, :) ==> 3-vector of one of the edges of the triangle.
    !      (:, 2, :) ==> 3-vector of another edge of the triangle.
    !      (:, 2, :) ==> 3-vector of a point in the triangle.
    ! This choice was chosen to precompute few of the variables continously used
    ! for the triangle-ray algorithm.
    !
    ! Written originally by Pablo Oyola, adapted to INPASIM by Jose Rueda
    ! -----------------------------------------------------------------------
    implicit none
    ! ---------------------
    ! Inputs
    character(len = *), intent(in) :: filename !< File containing the triangle-structure.
    logical, intent(in):: verb !< Write to console some info.
    ! Ouputs
    type(geom_element), intent(out) :: g_element  !< class with the triangle geometry
    ! Local variables
    real(8), allocatable :: buffer(:,:)
    character(len = 80) :: line
    integer :: ii
    ! ---------------------
    ! Read the file
    if (verb) then
      print*, 'reading: ', filename
    endif
    open(unit = 60, file = filename, form = 'formatted', status='old', action='read')
    read(60, *) g_element%name
    read(60, *) line  ! Dummy descritption name
    read(60, *) line  ! Dummy description line
    read(60, *) g_element%kind
    read(60, *) g_element%triangleNum

    ! Allocate space for the triangle file:
    if(verb) PRINT*, 'Reading triangle structure. #triangle = ', g_element%triangleNum
    allocate(buffer(3,g_element%triangleNum*3))

    ! We read the file to a buffer.
    read(60,*) buffer
    close(60)
    ! ---------------------
    allocate(g_element%triangles(g_element%triangleNum, 3, 3))

    if(verb)then
      print*, 'Name: ', g_element%name
      print*, 'Kind: ', g_element%kind
      print*, 'Processing the triangles into 2-vector-and-point structure.'
    end if
    ! The buffer will contain a raw copy of the triangles. We have to convert them
    ! into the format needed by the triangle-ray algorithm.
    do ii = 1,g_element%triangleNum
      ! (:, 1, :) ==> 3-vector of one of the edges of the triangle.
      g_element%triangles(ii, 1, 1) = buffer(1, 3*(ii-1)+2) - buffer(1, 3*(ii-1)+1)
      g_element%triangles(ii, 1, 2) = buffer(2, 3*(ii-1)+2) - buffer(2, 3*(ii-1)+1)
      g_element%triangles(ii, 1, 3) = buffer(3, 3*(ii-1)+2) - buffer(3, 3*(ii-1)+1)

      ! (:, 2, :) ==> 3-vector of another edge of the triangle.
      g_element%triangles(ii, 2, 1) = buffer(1, 3*(ii-1)+3) - buffer(1, 3*(ii-1)+1)
      g_element%triangles(ii, 2, 2) = buffer(2, 3*(ii-1)+3) - buffer(2, 3*(ii-1)+1)
      g_element%triangles(ii, 2, 3) = buffer(3, 3*(ii-1)+3) - buffer(3, 3*(ii-1)+1)

      ! (:, 3, :) ==> 3-vector of a point in the triangle.
      g_element%triangles(ii, 3, 1) = buffer(1, 3*(ii-1)+1)
      g_element%triangles(ii, 3, 2) = buffer(2, 3*(ii-1)+1)
      g_element%triangles(ii, 3, 3) = buffer(3, 3*(ii-1)+1)
    end do

    if(verb) then
      print*, 'Triangle preprocessed!'
      print*, '---'
    endif
    deallocate(buffer)
  end subroutine parseGeometry

  subroutine readGeometry(GeomFolder, n, verb)
    ! -----------------------------------------------------------------------
    ! Prepare the geometry elements
    !> \brief Prepare the geometry elements
    ! Written by Jose Rueda
    ! When executed, it set in the code workspace:
    !   - geometry !< the array of geometrical elements
    !   - pinhole  !< structure with data from the pinhole
    !
    ! Note: Maximum number of files to be read is set to 9. This could be
    ! easily enanced, but, really, more than 9 files??
    ! -----------------------------------------------------------------------
    implicit none
    ! Dummy variables
    character (len=*), intent(in) :: GeomFolder  !< ID of the geometry to read
    integer, intent(in) :: n   !< Number of elements to read
    logical, intent(in):: verb !< Write to console some info.

    ! Local variables
    integer :: kk,ll  !< index for loops
    integer:: ierr  !< error storing variable
    integer:: pinholeKind !< kind of pinhole (namelist variables)
    real(8), dimension(3) :: u1, u2, u3, rPin !< variables for the namelist
    real(8) :: d1, d2 !< pinhole sizes (namelist variables)
    character (len=1000) :: dummy_string, err_str, geometry_dir !< dummy strings
    character (len=1) :: number !< number where to save the element we are reading
    type(geom_element):: dumGeometry

    NAMELIST /ExtraGeometryParams/ nGeomElements, u1, u2, u3, rPin, d1, d2, ps, rotation, pinholeKind

    ! Read namelist and configure the plates

    open (unit=60, file=trim(GeomFolder)//'/ExtraGeometryParams.txt',form='formatted',iostat=ierr)
    read(60, NML=ExtraGeometryParams, iostat=ierr)
    close(60)

    ! read the plates
    allocate(geometry(n))
    do i = 1,n
      write(number, '(I1)') i
      dummy_string = trim(GeomFolder)//'/Element'//number//'.txt'
      call parseGeometry(trim(dummy_string),verb, geometry(i))
      if (geometry(i)%kind .eq. 1) then
        call vec_product(geometry(i)%triangles(1, 1, :),  geometry(i)%triangles(1, 2, :), foilNormal)
        foilNormal = foilNormal / sqrt(sum(foilNormal**2))
      elseif (geometry(i)%kind .eq. 2) then
        call vec_product(geometry(i)%triangles(1, 1, :),  geometry(i)%triangles(1, 2, :), ScintNormal)
        ScintNormal = ScintNormal / sqrt(sum(ScintNormal**2))
      endif
    enddo
    ! Ensure the scintillator to be the last element
    if (geometry(n)%kind .ne. 2) then
      exploreGeom: do i = 1,n-1
        if (geometry(i)%kind .eq. 2) then
          dumGeometry = geometry(n)
          geometry(n) = geometry(i)
          geometry(i) = dumGeometry
          print*, 'Scintillator was not the last element'
          print*, 'last element was now placed in position ', i
          print*, 'Scintillator now placed as the last element'
        endif
      enddo exploreGeom
    endif

    ! Fill the pinhole object
    pinhole%u1 = u1
    pinhole%u2 = u2
    pinhole%u3 = u3
    pinhole%r = rPin
    pinhole%d1 = d1
    pinhole%d2 = d2
    pinhole%kind = pinholeKind
    if (verb) then
      print*, 'Pinhole position', pinhole%r
      print*, 'Pinhole u1', pinhole%u1
      print*, 'Pinhole u2', pinhole%u2
      print*, 'Pinhole u3', pinhole%u3
      print*, 'Pinhole e1', pinhole%e1
      print*, 'Pinhole e2', pinhole%e2
      print*, 'Pinhole e3', pinhole%e3
      print*, 'Pinhole d1', pinhole%d1
      print*, 'Pinhole d2', pinhole%d2
      print*, 'Scintillator reference point', ps
    endif
  end subroutine readGeometry


!------------------------------------------------------------------------------
! SECTION 7: FIDASIM compatibility
!------------------------------------------------------------------------------
  subroutine readFIDASIM4Markers(filename, verbose)
    ! -----------------------------------------------------------------------
    ! Load a FIDASIM npa file
    !> \brief Load markers from a FIDASIM simulation
    ! Written by Jose Rueda
    ! When executed, it set in the code workspace:
    !   - F4Markers !< Structure with NPA data
    !   - verbose !< flag to print some info
    ! -----------------------------------------------------------------------

    ! Dummy variables:
    character(*), intent(in) :: filename !< Name of the file with the markers
    logical, intent(in) :: verbose

    ! Local variables:
    real(4), dimension(:), allocatable:: dummy1D
    real(4), dimension(:, :), allocatable:: dummy2D

    if (verbose) then
      print*, 'Reading markers info from: ', trim(filename)
    endif
    open(60, file=trim(filename), action='read',access='stream')
    read(60) F4Markers%shot_number
    read(60) F4Markers%time
    read(60) F4Markers%nr_npa
    read(60) F4Markers%counter
    ! Read the initial position
    allocate(dummy2D(F4Markers%counter, 3))
    allocate(F4Markers%ipos(3, F4Markers%counter))
    allocate(F4Markers%fpos(3, F4Markers%counter))
    allocate(F4Markers%v(3, F4Markers%counter))
    allocate(dummy1D(F4Markers%counter))
    allocate(F4Markers%wght(F4Markers%counter))
    allocate(F4Markers%kind(F4Markers%counter))
    read(60) dummy2D
    F4Markers%ipos = transpose(dble(dummy2D)) / 100.0d0
    read(60) dummy2D
    F4Markers%fpos = transpose(dble(dummy2D)) / 100.0d0
    read(60) dummy2D
    F4Markers%v = transpose(dble(dummy2D)) / 100.0d0
    read(60) dummy1D
    F4Markers%wght = dble(dummy1D) !* pinhole%d1**2 * pi
    read(60) dummy1D
    F4Markers%kind = int(dummy1D)
    if (verbose) then
      print*, F4Markers%counter, '*', nResampling, 'markers to be followed'
    endif
  end subroutine readFIDASIM4Markers


!------------------------------------------------------------------------------
! SECTION 8: Initialization of markers
!------------------------------------------------------------------------------
 subroutine initMarker(vmod, dt_initial, xxi, beta_min, dbeta, rli)
   ! -----------------------------------------------------------------------
   ! Initialise the marker
   !> \brief Initialise the marker for the simulation
   ! Written by Jose Rueda
   ! When executed, it set in the code workspace:
   !   - part !< Structure with particle data
   !
   ! -----------------------------------------------------------------------
   implicit none
   !----
   real(8), intent(in):: vmod   !< Modulus of the velocity to be used
   real(8), intent(in):: dt_initial  !< initial dt to use in the tracking [s]
   real(8), intent(in):: xxi  !< Value of the xi variable to use
   real(8), intent(in):: beta_min !< Minimum value of the random angle [rad]
   real(8), intent(in):: dbeta !< Interval of the beta angle [rad]
   real(8), intent(in):: rli !< Gyroradius of the particle
   !---
   real(8):: distance
   ! real(8):: randnumberaux
   ! real(8):: projection
   ! Clean the marker position and trajectory
   part%position(:, :) = 0.0d0
   part%velocity(:, :) = 0.0d0
   part%collision    = .False.
   part%kindOfCollision = 9  ! Just initialise it to a dummy value
   part%dt = dt_initial
   part%qm    = Zin *qe / M /amu_to_kg
   part%xi = xxi
   part%rl = rli

   ! Init the marker position
   distance = 1000.0
   if (pinhole%kind.eq.0) then
      do while (distance.gt.pinhole%d1)
        call random_number(ran)
        part%position(:, 1) = pinhole%r + &
          pinhole%d1 * ((-1+2*ran(1))*pinhole%u1 + &
                        (-1+2*ran(2))*pinhole%u2)
        distance = sqrt(sum((part%position(:, 1) - pinhole%r)**2))
      enddo
    else
      call random_number(ran)
      part%position(:, 1) = pinhole%r + &
          0.5*pinhole%d1 * (-1+2*ran(1))*pinhole%u1 + &
          0.5*pinhole%d2 * (-1+2*ran(2))*pinhole%u2
    endif
    ! Init marker velocity
    part%beta = beta_min + ran(3) * dbeta
    if (FILDSIMmode) then
      part%velocity(:, 1) = vmod * (sqrt(1-xxi**2) * (cos(part%beta)*pinhole%e1 + &
                                                      sin(part%beta)*pinhole%e2) +&
                                    xxi * pinhole%e3)
      part%weight         = 1.0
    else   ! INPA mode
      part%velocity(:, 1) = ((cos(xxi)*pinhole%u3 + &
                              sin(xxi)*pinhole%u1) * cos(part%beta) &
                           + sin(part%beta)*pinhole%u2) * vmod

      part%weight         = abs(sin(part%beta))
    endif
 end subroutine initMarker

 subroutine calculate_betas(flag)
   ! -----------------------------------------------------------------------
   ! Calculate the beta angles to be used
   !> \brief calculate the beta intervas to be latter used
   ! Written by Jose Rueda
   ! When executed, it set in the code workspace:
   !   - min_beta !< Array with beta angles for the simulation
   !   - delta_beta !< array with deltas for the beta angle
   !
   ! For INPA should be always zero, as there is no 'safety issues' the markers
   ! enter like neutral particles, so no much to be done
   ! -----------------------------------------------------------------------
   ! Dummy variables
   logical, intent(in):: flag
   ! local variables
   integer:: jj
   real(8):: aux, aux2, dum_beta, dum_beta_max
   real(8):: epsilon=1e-6

   ! Allocate the arrays
   allocate(marker_factor(nXI))
   allocate(min_beta(nXI))
   allocate(delta_beta(nXI))

   if (flag) then  ! The 'secure mode'
     do jj=1, nXI
       if ((abs(pinhole%eta) .le. epsilon) .or. (abs(pinhole%eta-pi/2.0d0) .le. epsilon) )then
         ! Not efficient, I know, @todo
         min_beta(jj) = minAngle  ! Save the input namelist default
         delta_beta(jj) = dAngle
         marker_factor(jj) = 1
       else
         aux2 = sqrt(1-XI(jj)**2)
         if (abs(aux2) .le. epsilon) then
           min_beta(jj) = minAngle  ! Save the input namelist default
           delta_beta(jj) = dAngle
           marker_factor(jj) = 1
         else
           aux = -1.0d0*DBLE(IpBt)*XI(jj)*aux2*cos(pinhole%eta)/abs(sin(pinhole%eta))
           if (aux .ge. 1.0d0) then  ! all can enter
             min_beta(jj) = minAngle  ! Save the input namelist default
             delta_beta(jj) = dAngle
             marker_factor(jj) = 1
           else if (aux .le. -1.0d0) then ! Nothing can enter
             marker_factor(jj) = 0
             min_beta(jj) = minAngle  ! Save the input namelist default
             delta_beta(jj) = dAngle
           else ! Normal case
             dum_beta = asin(aux)
             if (dum_beta .ge. 0.0d0) then
               print*, XI_input(jj)
               print*, dum_beta
               print*, 'Option not implemented, contat jrrueda@us.es '
               stop
             else
               min_beta(jj) = max(pi + abs(dum_beta), minAngle)
               dum_beta_max = min(2.0d0*pi + dum_beta, dAngle + minAngle)
               delta_beta(jj) = dum_beta_max - min_beta(jj)
               marker_factor(jj) = 1

             endif
           endif
         endif
       endif
     enddo
   else
     min_beta(:) = minAngle  ! Save the input namelist default
     delta_beta(:) = dAngle
     marker_factor(:) = 1
   endif
 end subroutine calculate_betas


!-----------------------------------------------------------------------------
! SECTION 9: Saving
!-----------------------------------------------------------------------------
 subroutine saveOrbit(fid, marker_to_save, step_to_save, long_flag)
   ! -----------------------------------------------------------------------
   ! Save the orbit
   !> \brief Save the orbit into the binary file
   ! Written by Jose Rueda
   !
   ! INPUTS:
   !  - fid: FID of the file where to write (it should be open already)
   !  - marker to save: marker to be saved
   !  - step_to_save: Maximum step to be saved in the trajectory
   !  - long_flag: if true to so called 'long' format will be used
   ! -----------------------------------------------------------------------
   implicit none
   ! Dummy variables
   integer, intent(in):: fid
   type(marker), intent(in):: marker_to_save
   integer, intent(in):: step_to_save
   logical, intent(in):: long_flag

   write(fid) step_to_save
   write(fid) marker_to_save%rl
   write(fid) marker_to_save%xi
   write(fid) marker_to_save%kindOfCollision
   write(fid) transpose(marker_to_save%position(:, 1:step_to_save))
   if (long_flag) then
     write(fid) transpose(marker_to_save%velocity(:, 1:step_to_save))
   endif
 end subroutine saveOrbit

 subroutine writeStrikeFileHeader(filename, fid, kindOfFile, flag, nCol)
   ! -----------------------------------------------------------------------
   ! Open the strike file and prepare the strike size
   !> \brief Save the orbit into the binary file
   ! Written by Jose Rueda
   !
   ! INPUTS:
   !  - fid: FID of the file where to write (it should be open already)
   !  - marker to save: marker to be saved
   !  - step_to_save: Maximum step to be saved in the trajectory
   !  - long_flag: if true to so called 'long' format will be used
   ! -----------------------------------------------------------------------
   !Dummy variables
   character(*), intent(in) :: filename
   integer, intent(in) :: fid
   integer, intent(in) :: kindOfFile
   logical, intent(in) :: flag
   integer, intent(out) :: nCol

   if (FILDSIMmode.and.(kindOfFile.eq.1))then
    nCol = 12
    ! Set the pointers to the store functions
    savePartToStrike => savePartToStrikeFILD1
   elseif (FILDSIMmode.and.(kindOfFile.eq.2)) then
    nCol = 14
    ! Set the pointers to the store functions
    savePartToStrike => savePartToStrikeFILD2
   elseif (.not.FILDSIMmode.and.(kindOfFile.eq.1)) then
    nCol = 19
    ! Set the pointers to the store functions
    savePartToStrike => savePartToStrikeINPA1
   elseif (.not.FILDSIMmode.and.(kindOfFile.eq.101)) then
      nCol = 22
      ! Set the pointers to the store functions
      savePartToStrikeSignal => savePartToStrikeINPASignal1
   endif

   if (flag .and. (kindOfFile.lt.101)) then
    open(unit=fid, file=trim(filename), access = 'stream', action='write', status='replace')
    write(fid) versionID1, versionID2, runID, nGyroradius, rL, nxi, XI_input, transfer(FILDSIMmode, 1), nCol, kindOfFile
   elseif (flag .and. (kindOfFile.gt.100)) then 
    open(unit=fid, file=trim(filename), access = 'stream', action='write', status='replace')
    write(fid) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0, transfer(FILDSIMmode, 1), nCol, kindOfFile
   endif

  


 end subroutine writeStrikeFileHeader

 subroutine savePartToStrikeFILD1(particle, counter)
  ! Dummy variables
  type(marker), intent(in):: particle
  integer, intent(in) :: counter
  Strike(1:3, cScintillator ) = part%collision_point ! f point
  Strike(4, cScintillator) = part%weight ! weight
  Strike(5, cScintillator) = part%beta ! Beta angle
  Strike(6:8, cScintillator) = part%position(:, 1) ! i pos
  Strike(9:11, cScintillator ) = &
    MATMUL(rotation, part%collision_point - ps) ! f pos scint
  Strike(12, cScintillator) = &
    180/pi*acos(incidentProjection) ! incident angle
 end subroutine savePartToStrikeFILD1

 subroutine savePartToStrikeFILD2(particle, counter)
  ! Dummy variables
  type(marker), intent(in):: particle
  integer, intent(in) :: counter
  real(8) :: v1, v2, v3, rparticle(3), dummy34(6)
  Strike(1:3, cScintillator ) = part%collision_point ! f point
  Strike(4, cScintillator) = part%weight ! weight
  Strike(5, cScintillator) = part%beta ! Beta angle
  Strike(6:8, cScintillator) = part%position(:, 1) ! i pos
  Strike(9:11, cScintillator ) = &
    MATMUL(rotation, part%collision_point - ps) ! f pos scint
  Strike(12, cScintillator) = &
    180/pi*acos(incidentProjection) ! incident angle
  v1 = sum(pinhole%e1*part%velocity(:, istep))&
    /norm2(part%velocity(:, istep))
  v2 = sum(pinhole%e2*part%velocity(:, istep))&
    /norm2(part%velocity(:, istep))
  v3 = sum(pinhole%e3*part%velocity(:, istep))&
    /norm2(part%velocity(:, istep))
  Strike(13, cScintillator) = v1
  Strike(14, cScintillator) = v2
  call cart2pol(part%position(:, istep), rparticle)
  call getField(rparticle(1), rparticle(2), rparticle(3), 0.0d0, dummy34)
  write(114,'(6F14.6)') v1, v2, v3, dummy34(1:3)

 end subroutine savePartToStrikeFILD2

 subroutine savePartToStrikeINPA1(particle, counter)
  ! Dummy variables
  type(marker), intent(in):: particle
  integer, intent(in) :: counter
  Strike(1:3, cScintillator ) = part%collision_point ! f point
  Strike(4, cScintillator) = part%weight ! weight
  Strike(5, cScintillator) = part%beta ! Beta angle
  Strike(6:8, cScintillator) = part%position(:, 1) ! i pos
  Strike(9:11, cScintillator ) = &
    MATMUL(rotation, part%collision_point - ps) ! f pos scint
  Strike(12, cScintillator) = &
    180/pi*acos(incidentProjection) ! incident angle
  call minimumDistanceLines(part%position(:, 1), &
    part%velocity(:, 1)/v0, nbi%p0,&
    nbi%u, dMin, posMin)
  Strike(13:15,cScintillator) = posMin   ! Initial point (NBI)
  Strike(16:18,cScintillator) = part%velocity(:, 1)  ! Initial velocity
  Strike(19, cScintillator) = dMin ! Closer distance to NBI
 end subroutine savePartToStrikeINPA1

 subroutine savePartToStrikeINPASignal1(particle, counter)
  ! Dummy variables
  type(marker), intent(in):: particle
  integer, intent(in) :: counter
  Strike(1, cScintillator ) = dble(i) ! FIDASIM marker id
  Strike(2:4, cScintillator ) = part%collision_point ! f point
  Strike(5:7, cScintillator ) = &
     MATMUL(rotation, part%collision_point - ps) !f point in scintillator region
  Strike(8:10, cScintillator) = F4Markers%ipos(:, i) ! CX position
  Strike(11:13,cScintillator) = part%velocity(:, 1)   ! Initial velocity
  Strike(14,cScintillator) = part%weight ! weight
  Strike(15, cScintillator) = F4Markers%kind(i) ! Kind of signal
  Strike(16, cScintillator) = part%energy0 ! energy at entrance
  Strike(17, cScintillator) = 0.5*sum(part%velocity(:, istep)**2)*M/qe*amu_to_kg/1000.0
  Strike(18, cScintillator) = part%cosalpha_foil
  Strike(19, cScintillator) = F4Markers%wght(i) / normalization_resample
  Strike(20:22, cScintillator) = part%ionization_point 
 end subroutine savePartToStrikeINPASignal1
!-----------------------------------------------------------------------------
! SECTION 10: Foil interaction
!-----------------------------------------------------------------------------
 subroutine foilInteraction(MC_marker, normal, step)
   ! -----------------------------------------------------------------------
   ! Save the orbit
   !> \brief Apply the interaction to the marker
   ! Written by Jose Rueda
   !
   ! INPUTS:
   !  - interactionModel: integer specifiying the model which should be applied
   !  - MC_marker: markers to be modified
   !  - normal: normal vector to the carbon foil
   !  - step: step where the collision have taken place
   ! -----------------------------------------------------------------------
   ! Dummy variables
   type(marker), intent(inout):: MC_marker
   real(8), dimension(3), intent(in):: normal
   integer, intent(in):: step
   ! Local variables
   real(8):: modV  !< velocity modulus
   real(8):: Energy !< Energy of the particle, in keV
   real(8), dimension(3):: localV  !< Dummy copy of the marker velocity
   real(8):: cos_alpha  !< projection of the velocity on the normal
   real(8):: deltaE, yield

   if (FoilYieldModel.eq.1) then
     MC_marker%weight = MC_marker%weight * FoilYieldParameters(1)
   elseif (FoilYieldModel.eq.2) then
     ! estimate the ionization efficiency
     yield = FoilYieldParameters(2) - FoilYieldParameters(3) * exp(-FoilYieldParameters(4) * MC_marker%energy0 / M)
     ! Scale the weight
     MC_marker%weight = MC_marker%weight * FoilYieldParameters(1) * yield
   endif

   if (FoilElossModel.eq.1) then
     ! Scale the velocity
     localV = MC_marker%velocity(:, step+1)
     modV = sqrt(sum(localV**2))
     localV = localV/modV
     ! Get the incident angle (projection)
     cos_alpha = abs(localV(1)*normal(1) + localV(2)*normal(2) + localV(3)*normal(3))
     ! Scale the energy
     deltaE = - FoilElossParameters(1) * FoilElossParameters(2)/cos_alpha * sqrt(MC_marker%energy0)
     ! Scale the velocity
     MC_marker%velocity(:, step+1) = localV * sqrt(2 * (MC_marker%energy0 + deltaE)/(M/qe*amu_to_kg/1000.0))
     ! Save the cosalpha value
     MC_marker%cosalpha_foil = cos_alpha
   elseif (FoilElossModel.eq.2) then
     ! Scale the velocity
     localV = MC_marker%velocity(:, step+1)
     modV = sqrt(sum(localV**2))
     localV = localV/modV
     ! Get the incident angle (projection)
     cos_alpha = abs(localV(1)*normal(1) + localV(2)*normal(2) + localV(3)*normal(3))
     ! Scale the energy
     Energy = 0.5*modV**2*M/qe*amu_to_kg/1000.0
     deltaE = - FoilElossParameters(1)/cos_alpha * &
      (FoilElossParameters(2) * sqrt(Energy) + FoilElossParameters(3)*Energy)
     ! Scale the velocity
     MC_marker%velocity(:, step+1) = localV * sqrt(2 * (Energy + deltaE)/(M/qe*amu_to_kg/1000.0))
     ! Save the cosalpha value
     MC_marker%cosalpha_foil = cos_alpha
   elseif (FoilElossModel.eq.4) then
     ! Scale the velocity
     localV = MC_marker%velocity(:, step+1)
     modV = sqrt(sum(localV**2))
     localV = localV/modV
     ! Get the incident angle (projection)
     cos_alpha = abs(localV(1)*normal(1) + localV(2)*normal(2) + localV(3)*normal(3))
     ! Scale the energy
     Energy = 0.5*modV**2*M/qe*amu_to_kg/1000.0
     if (Energy-FoilElossParameters(3).gt.0.0) then
      deltaE = - FoilElossParameters(1)/cos_alpha * &
        (FoilElossParameters(2) * sqrt(Energy-FoilElossParameters(3)))
     else
      deltaE = 0.0
     endif
     ! Scale the velocity
     MC_marker%velocity(:, step+1) = localV * sqrt(2 * (Energy + deltaE)/(M/qe*amu_to_kg/1000.0))
     ! Save the cosalpha value
     MC_marker%cosalpha_foil = cos_alpha
    elseif (FoilElossModel.eq.3) then
     ! Scale the velocity
     localV = MC_marker%velocity(:, step+1)
     modV = sqrt(sum(localV**2))
     localV = localV/modV
     ! Get the incident angle (projection)
     cos_alpha = abs(localV(1)*normal(1) + localV(2)*normal(2) + localV(3)*normal(3))
     ! Scale the energy
     Energy = 0.5*modV**2*M/qe*amu_to_kg/1000.0
     deltaE = - FoilElossParameters(1)/cos_alpha * &
      (FoilElossParameters(2) * Energy + FoilElossParameters(3))
     ! Scale the velocity
     MC_marker%velocity(:, step+1) = localV * sqrt(2 * (Energy + deltaE)/(M/qe*amu_to_kg/1000.0))
     ! Save the cosalpha value
     MC_marker%cosalpha_foil = cos_alpha
   endif

 end subroutine foilInteraction


!-----------------------------------------------------------------------------
! SECTION 11: Yield modelling
!-----------------------------------------------------------------------------
 subroutine yieldScintillator(MC_marker, step)
   ! -----------------------------------------------------------------------
   ! Save the orbit
   !> \brief Apply the interaction to the marker
   ! Written by Jose Rueda
   !
   ! INPUTS:
   !  - MC_marker: markers to be modified
   !  - step: step where the collision have taken place
   ! -----------------------------------------------------------------------
   ! Dummy variables
   type(marker), intent(inout):: MC_marker
   integer, intent(in):: step
   ! Local variables
   real(8):: modV  !< velocity modulus
   real(8):: Energy !< Energy of the particle, in keV
   real(8), dimension(3):: localV  !< Dummy copy of the marker velocity
   real(8):: deltaE
  if (ScintillatorYieldModel.eq.1) then
     MC_marker%weight = MC_marker%weight * ScintillatorYieldParameters(1)
  else if (ScintillatorYieldModel.eq.2) then
     localV = MC_marker%velocity(:, step+1)
     modV = sqrt(sum(localV**2))
     localV = localV/modV
     ! Scale the weight
     Energy = 0.5*modV**2*M/qe*amu_to_kg/1000.0
     MC_marker%weight = MC_marker%weight * ScintillatorYieldParameters(1) * (Energy ** ScintillatorYieldParameters(2))
   endif
 end subroutine yieldScintillator

end module sinpa_module
