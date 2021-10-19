! -----------------------------------------------------------------------------
!  TRACKER MAIN MODULE
! -----------------------------------------------------------------------------
! TITLE         : SINPA
! PROJECT       : Synthetic INPA
! MODULE        : Auxiliar routines
! AFFILIATION   : University of Sevilla
!> \author Jose Rueda - Universidad de Sevilla
!> \date 26/07/2021
!> \version 0.0
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
  integer, parameter:: versionID1 = 0  !< ID version number, to identify ouput
  integer, parameter:: versionid2 = 0  !< ID version 2, to identify the ouput
  real (8), parameter:: pi = 3.141592653589793 !< pi
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
      real (8):: weight  !< Weight of the markers, [part/s/cm^2 of pinhole, for FIDASIM, in au, for mapping]
      real (8):: type    !< Marker type: identify if is active or passive, follows FIDASIM criteria
      logical :: collision = .FALSE. !< Flag for the collision
      real (8), dimension(3) :: collision_point !< Point of collision in cm
      real (8):: qm !< q over m, SI
      ! For the mapping
      real (8):: xi  !< xi angle (or pitch for FILD), radians
      real (8):: beta  !< beta angle (or Gyrophase for FILD), radians
      ! Integration
      real (8):: dt !< dt for the integration [s]
      integer::kindOfCollision ! With whom has collide
  end type marker

  ! ABSTRACT CLASS FOR THE FIELDS
  ! >/briefg
  type :: fields
    integer :: nr                 ! Size of the grid in the r dimension
    integer :: nz                 ! Size of the grid in the z dimension
    real(8):: rmin       ! Minimum value of the r array
    real(8):: zmin       ! Minimum value of the z array
    real(8):: dr         ! Spacing f the grid in the r direction
    real(8):: dz         ! Spacing f the grid in the r direction
    real(8), dimension(:,:,:), allocatable :: Brzt,E ! Fields
  end type fields

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

  type :: scatteringData
    ! Sigma for the model of scattering in the carbon foil
    integer :: na                 ! Number of angles
    integer :: ne                 ! Number of energies
    real(8):: amin     ! Angles array
    real(8):: da       !
    real(8):: emin     ! Energies array
    real(8):: de       !
    real(8), dimension(:,:), allocatable :: sigma      ! Sigma matrix
    ! Coefficients of the fit for the energy. See the wiki or the carbon foil
    ! routines in section 4
    real(8), dimension(:,:), allocatable :: p1_ener,p2_ener
    real(8):: depth
  end type scatteringData


  type :: sigma_optics
    ! Sigma for the model of optical resolution
    integer :: no                 ! Number of points
    real(8),dimension(3):: optical_axis ! Position of the optical axis
    real(8):: delta       !
    real(8), dimension(:), allocatable :: dist,so
  end type sigma_optics


  type :: NBI_class
    character (len=20) :: name  !< NBI name
    real(8), dimension(3) ::  p0  !< Initial point of the NBI
    real(8), dimension(3) ::  u   !< Velocity directoin of the NBI
    real(8), dimension(3) ::  d   !< distance between points (cm)
    integer:: npoints  !< Number of poins along the NBI line
    real(8), dimension(:,:), allocatable :: points !< 3d coordinates of the points
    real(8), dimension(:), allocatable :: Rpoints !< R coordinate of the points
    real(8), dimension(:), allocatable :: Phipoints !< Phi coordinate of the points
  end type NBI_class


  !> \brief Contains the data of 1D direction in a compact type.
  type:: grid_class
    real(8), allocatable::data(:)   !< Contains the points of the grid.
    real(8)::x0                     !< Left bound of the grid.
    real(8)::x1                     !< Right bound of the grid.
    integer::size                            !< Number of grid points.
    real(8)::dx                     !< Difference between consecutive points.
  end type grid_class

  !> \brief Contains the reference system in the pinhole
  type :: pinhole_system_type
    real(8), dimension(3):: rPin  !< Position of the pinhole (cm)
    integer:: kind  !< kind of pinhole (circle, 0) (rectangle, 1)
    real(8):: d1 ! Pinhole radius, or size along u1 (in cm)
    real(8):: d2 ! Pinhole radius, or size along u2 (in cm)
    real(8), dimension(3):: r  !< Position of the pinhole (cm)
    real(8), dimension(3):: u1  !< Vector in the pinhole plane 1
    real(8), dimension(3):: u2  !< vector in the pinhole plane 2
    real(8), dimension(3):: u3  !< normal to the pinhole
  end type pinhole_system_type
  ! ---------------------------------------------------------------------------
  ! Interfaces
  ! ---------------------------------------------------------------------------
  abstract interface
    subroutine interpolateField_interface(rq, zq, phiq,out)
      implicit none
      real(kind=8), intent(in)::rq, zq, phiq! Query points in all directions.
      real(kind=8), intent(out)::out(3) ! Br, Bz, Bphi
    end subroutine interpolateField_interface
  end interface
  procedure(interpolateField_interface),pointer, public::getField
  ! -------------------------------------------------------------------------
  ! Variables
  ! -------------------------------------------------------------------------
  ! --- Mapping
  type(NBI_class) :: nbi
  real(8), dimension(:, :), allocatable:: closestPoints    !< closer points to NBI
  real(8), dimension(:), allocatable:: dPoints    !< distance to NBI
  type(marker) :: part
  real(8) :: dt !< dt for the time integration, in s
  real(8) :: dtNeutral
  real(8) :: OmegaPart !< Gyrofrequency of mapping particles
  real(8), dimension(:, :), allocatable:: Strike !< Matrix to store data
  real(8), dimension(:, :), allocatable:: StrikeMap !< Matrix to store the strike map
  real(8), dimension(:, :), allocatable:: CollimatorStrikes !< Matrix to store the strike map
  real(8), dimension(3,3):: rotation !< Rotation matrix
  real(8), dimension(3):: ps !< reference point in the scintillator

  ! --- Magnetic field
  real(8), dimension(:, :, :), allocatable, protected::Brfield    !< Radial magnetic field.
  real(8), dimension(:, :, :), allocatable, protected::Bzfield    !< Vertical magnetic field.
  real(8), dimension(:, :, :), allocatable, protected::Bphifield  !< Toroidal magnetic field.
  ! --- grids
  type(grid_class), protected:: grr  !< Grid classes to store the R grid
  type(grid_class), protected:: gzz  !< Grid classes to store the Z grid
  type(grid_class), protected:: gphi !< Grid classes to store the phi grid
  ! --- Geometry structures
  type(geom_element), dimension (:), allocatable:: geometry
  type(geom_element) :: scintillator, foil, collimator
  ! --- FileNames and paths
  character (len=1000) :: BFieldFile !< Magnetic field file
  character (len=1000) :: CollimatorFile !< Collimator file
  character (len=1000) :: FoilFile !< Foil file
  character (len=1000) :: ScintillatorFile !< Scintillator file
  character (len=1000) :: input_filename !< Scintillator file
  character (len=1000) :: SINPA_dir !< root directory of SINPA
  character (len=1000) :: FIDASIMfolder !< root directory of SINPA

  ! --- Counters
  integer:: cCollimator !< Number of markers impinging the collimator
  integer:: cWrongNeutral !< Number of markers not colliding with the collimator neither carbon foil
  integer:: cWrongIons !< Number of markers not colliding with the scintillator
  integer:: cScintillator !< Number of markers colliding with the scintillator, for each alpha
  integer:: cFoil !< Number of markers colliding with the foil, for each alpha
  integer:: cOrb = 0!< Number of saved orbits
  ! --- Signal
  type(FIDASIM):: F4Markers !< FIDASIM4 markers

  ! --- Others
  integer:: ierr !< error management
  integer:: irl, ialpha, imc, istep, iistep, i!< dummy for loops
  real(8), dimension(3) :: rPinCyl  !< position of the pingole (r,p,z)
  real(8), dimension(3) :: Bpinhole  !< B field at pinhole, vector
  real(8) :: BpinholeMod !< B field at pinhole
  real(8) :: v0  !< initial velocity for mapping markers
  real(8), dimension(3) :: r0  !< initial position for mapping markers
  real(8), dimension(3) :: vv0  !< initial velocity for mapping markers
  real(8), dimension(3) :: ran  !< Set of random numbers
  real(8), dimension(:, :), allocatable:: random_numbers  !< Set of random numbers
  real(8) :: rand_number  !< Set of random numbers
  real(8) :: beta  !< random angle for the markers
  real(8) :: dMin  !< minimum distance to NBI of markers trajectories
  real(8), dimension(3) :: posMin  !< position to minimum distance to NBI of the mapping marker
  type(pinhole_system_type) :: pinhole !< pinhole definition

  ! ---------------------------------------------------------------------------
  ! NAMELIST
  ! --- Configuration file
  ! Dummy variables for the namelist
  character (len=50) :: runID !< runID
  character (len=50) :: geomID !< geometryID
  logical:: FILDSIMmode !< Flag to use FILDSIM mode or not
  integer:: nGeomElements !< Number of geometrical elements
  integer:: nxi !< number of pitches (R) to simulate
  integer:: nGyroradius !< number of energies (gyroradius) to simulate
  integer:: nMap !< number of markers per energy-pitch for the map
  logical:: mapping !< flag to decide if we launch mapping markers
  logical:: signal !< flag to decide if we follow FIDASIM4 markers
  logical:: resampling !< flag to decide if we resample FIDASIM markers
  integer:: nResampling !< number of markers per energy-pitch for the map
  logical:: saveOrbits !< flag to decide if we save the orbits
  real(8) :: saveRatio !< Ratio of orbits to be saved
  logical:: verbose !< Print information of the run
  real(8) :: M !< Mass of the particle, in amu
  real(8) :: Zin !< Initial chargeof the particle, in |e|
  real(8) :: Zout !< Charge after ionization, in |e|

  ! Namelist
  NAMELIST /config/ runID, geomID, FILDSIMmode, nGeomElements, nxi, nGyroradius, &
    nMap, mapping,&
    signal, resampling, nResampling, saveOrbits, saveRatio, SINPA_dir,&
    FIDASIMfolder, verbose, M, Zin, Zout

  ! --- Input
  integer:: nGyro !< number of points in a complete gyrocicle of the particle

  real(8) :: minAngle  !< collimator vertical acceptance
  real(8) :: dAngle  !< collimator vertical acceptance

  real(8), dimension(:), allocatable:: alphas
  real(8), dimension(:), allocatable:: rL !< FILDSIM gyroradius
  real (8):: maxT !< maximum time to follow particles
  NAMELIST /inputParams/ nGyro, minAngle, dAngle, alphas, rL, maxT

  ! --- Mapping orientation
  real(8) :: pinRadius  !< Radius of the head
  real(8), dimension(3) :: rPin  !< poisition of the pinhole (x,y,z)
  real(8), dimension(3) :: u3  !< normal in the pinhole
  real(8), dimension(3) :: u1  !< tangent to the pinhole
  real(8), dimension(3) :: u2  !< second tangent in the pinhole

  ! --- Particle and foil parameters
  logical:: energyLoss !< Apply energy loss
  real(8) :: a_SRIM !< First coefficient of SRIM energy loss model
  real(8) :: b_SRIM !< Second coefficient of SRIM energy loss model
  logical:: weightChange !< Apply energy loss
  real(8) :: a_ionization !< First coefficient of ionization model
  real(8) :: b_ionization !< Second coefficient of ionization model
  real(8) :: c_ionization !< Third coefficient of ionization model
  real(8) :: geometricTrans !< Mesh geometric transmission factor
  logical :: scattering !< apply scatering in the foil
  ! Namelist
  NAMELIST /markerinteractions/ energyLoss, a_SRIM, b_SRIM, weightChange, &
    a_ionization, b_ionization, c_ionization, geometricTrans, scattering

  ! --- NBI namelist
  ! Dummy variables for namelist
  real(8), dimension(3):: p0, u
  integer:: npoints
  real(8):: d
  ! Namelist
  NAMELIST /nbi_namelist/ p0, u, npoints, d
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


  subroutine rl2v0(rL, M, Z, B, v0)
    ! -----------------------------------------------------------------------
    ! Give the velocity module for a given larmor radius
    !
    !
    ! Written by Jose Rueda Rueda
    ! -----------------------------------------------------------------------
    implicit none
    ! Dummy variables
    real(8), intent(in):: rL !< Lamor radius
    real(8), intent(in):: M !< Mass in amu
    real(8), intent(in):: Z !< charge in electron charge (abs)
    real(8), intent(in):: B !< Modulus of the field
    real(8), intent(out):: v0 !< Modulus of the field

    v0 = rL * abs(Z) * B / M * qe / amu_to_kg
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


  subroutine cart2pol(cart, pol)
    ! Created by Pablo Oyola
    implicit none

    real(8), intent(in)::cart(3)   !< Cartesian points (x, y, z)
    real(8), intent(out)::pol(3)   !< Polar points (r, z, phi)

    pol(1) = sqrt(cart(1)**2+cart(2)**2)
    pol(2) = cart(3)
    pol(3) = atan2(cart(2), cart(1))
  end subroutine cart2pol


  subroutine pol2cart(polar, cart)
    ! Createed by Pablo Oyola
    implicit none

    real(8), intent(out)::cart(3)   !< Cartesian points (x, y, z)
    real(8), intent(in)::polar(3)   !< Polar points (r, z, phi)

    cart(1) = polar(1)*cos(polar(3))
    cart(2) = polar(1)*sin(polar(3))
    cart(3) = polar(2)
  end subroutine pol2cart

  ! -----------------------------------------------------------------
  ! POLAR TO CARTESIAN COORDINATES (COVARIANT VERSION)
  !> \brief Covariant transformation from polar coordinates to cartesian coordinates
  !! (i.e. velocities).
  ! Created by Pablo Oyola
  ! -----------------------------------------------------------------
  subroutine pol2cart_cov(vpol, vcart, r0)
    implicit none

    real(8), intent(out)::vcart(3)  !< Cartesian components (Ax, Ay, Az)
    real(8), intent(in)::vpol(3)    !< Polar components (Ar, Aphi, Az)
    real(8), intent(in)::r0(3)      !< Point where the field point is (r, z, phi)

    vcart(1) = vpol(1)*cos(r0(3)) - vpol(2)*sin(r0(3))
    vcart(2) = vpol(1)*sin(r0(3)) + vpol(2)*cos(r0(3))
    vcart(3) = vpol(3)
  end subroutine pol2cart_cov
!-------------------------------------------------------------------------------
! SECTION 2: NBI
!-------------------------------------------------------------------------------
  subroutine initNBI(filename, nbi)
    ! ------------------------------------------------------------------
    ! Initialise the NBI object
    !> \brief Initialise the NBI object
    !> \detail Load data from the NBI file and create the series of points
    ! along the NBI line to be latter used for mapping
    !
    ! Written Jose Rueda Rueda
    ! ------------------------------------------------------------------
    ! ----
    ! INPUTS:
    character (len=*), intent(in):: filename  !< file with NBI data
    ! Outputs:
    type(NBI_class), intent(out):: nbi   !< NBI object
    ! Local variables
    integer :: ii
    ! ---

    ! Open the file and read:
    open(60, file=filename, form = 'formatted', status='old', action='read')
    read(60, NML=nbi_namelist)
    close(60)
    nbi%p0 = p0
    nbi%u = u
    nbi%npoints = npoints
    nbi%d = d

    ! calculate the points along the line
    allocate(nbi%points(3, npoints))
    allocate(nbi%Rpoints(npoints))
    do ii = 1,npoints
      nbi%points(:, ii) = p0 + ii * u
      nbi%rpoints(ii) = sqrt(nbi%points(1, ii)**2 + nbi%points(2, ii)**2)
      nbi%phipoints(ii) = atan2(nbi%points(2, ii), nbi%points(1, ii))
    end do
  end subroutine initNBI

  ! subroutine getNBIpoint(nbi, axis_data, value, p)
  !   ! ToDo. Include tolerance
  !   ! INPUTS
  !   type(NBI_class), intent(in) :: nbi  !< NBI object
  !   real(8), dimension(nbi%npoints), intent(in) :: axis_data
  !   real(8), intent(in) :: value
  !   ! Ouputs
  !   real(8), dimension(3), intent(out) :: p
  !   ! Local
  !   integer :: i
  !
  !   ! Find the position of the desired value
  !   i = minloc(abs(axis_data - value))
  !   ! return the data
  !   p = nbi%points(:, i)
  ! end subroutine getNBIpoint

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
    real(8)::ar, az, aphi, ar1, az1, aphi1

    ia  = max(1, min(x1%size-1, int((x1q - x1%x0)/x1%dx  + 1)))
    ia1 = ia + 1
    ja  = max(1, min(x2%size-1, int((x2q - x2%x0)/x2%dx  + 1)))
    ja1 = ja + 1
    ka  = max(1, min(x3%size-1, int((x3q - x3%x0)/x3%dx  + 1)))
    ka1 = ka + 1

    ar1   = max(0.0d0, min(1.0d0, (x1q-x1%data(ia))/x1%dx))
    ar    = 1.0d0 - ar1
    az1   = max(0.0d0, min(1.0d0, (x2q-x2%data(ka))/x2%dx))
    az    = 1.0d0 - az1
    aphi1 = max(0.0d0, min(1.0d0, (x3q-x3%data(ja))/x3%dx))
    aphi  = 1.0d0 - aphi1

    ! Interpolating coefficients.
    aaa(1) = ar *az *aphi
    aaa(2) = ar1*az *aphi
    aaa(3) = ar *az1*aphi
    aaa(4) = ar1*az1*aphi
    aaa(5) = ar *az *aphi1
    aaa(6) = ar1*az *aphi1
    aaa(7) = ar *az1*aphi1
    aaa(8) = ar1*az1*aphi1

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
  subroutine parseField(Bname)
    !---------------------------------------------------------------------------
    ! This subroutine load the magnetic field
    !
    ! Adaptation from Pablo Oyola iHIBPsim code
    !
    ! INPUTS:
    !       - Bname: Name of the field to read
    ! OUTPUTS:
    !       - magnetic field arrays will be filled
    !---------------------------------------------------------------------------
    implicit none
    ! Dummy variables
    character(len=*), intent(in)::Bname   !< Filename with the magnetic field.
    ! Local variables
    integer::lr, lz, lphi
    integer::lr1, lz1, lphi1, ntot1 ! Auxiliar variables to compare the field sizes.
    real(8)::rmin, rmax, zmin, zmax, phimin, phimax
    integer::send_size, ii

    ! Parse de magnetic field
    print*, 'Parsing the magnetic field...'
    open(unit = 60, file=Bname, access='stream', status='old', action='read')
    read(60) lr, lz, lphi, &
             rmin, rmax, zmin, zmax, phimin, phimax
    print*, 'The header of the magnetic field was parsed:'
    print*, 'RESULTS'
    print*, 'nr = ', lr, ' nz = ', lz, ' nphi = ', lphi
    print*, 'INTEGRATION REGION:'
    print*, 'R = (',rmin, ' , ', rmax, ')'
    print*, 'z = (',zmin, ' , ', zmax, ')'
    print*, 'phi = (',phimin, ' , ', phimax, ')'
    print*, 'Allocating variables... '

    ! With this, we can allocate all the variables.
    allocate(Brfield(lr, lphi, lz))
    allocate(Bzfield(lr, lphi, lz))
    allocate(Bphifield(lr, lphi, lz))
    print*, 'Allocated! Reading the magnetic field'
    read(60) Brfield
    read(60) Bphifield
    read(60) Bzfield
    close(60)

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
    ! Prepare the interfaces for the field
    if(lr .eq. 1) then  ! We have a uniform
      getField => getField_0D
    elseif(lphi.eq.1)then  ! We have an axisymmetric magnetic field
      print*, 'We have a 2D field  '
      getField => getField_2D
    else
      getField => getField_3D ! We have an non-axisymmetric magnetic field
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

    getField => NULL()
  end subroutine unloadField

  subroutine getField_3D(rq, zq, phiq, out)
    ! -----------------------------------------------------------------------
    ! GET THE 3D FIELDS INTERPOLATION.
    !> \brief Obtains Br, Bz, Bphi interpolated in a given point,
    !! (R, phi, z). Linear interpolation is used.
    ! Written by Pablo Oyola for the iHIBPsim code
    ! -----------------------------------------------------------------------
    implicit none

    real(8), intent(in)::rq      !< Query point in R major to interpolate.
    real(8), intent(in)::zq      !< Query point in Z to interpolate.
    real(8), intent(in)::phiq    !< Query point in toroidal direciton.
    real(8), intent(out)::out(3) !< Br, Bz, Bphi

    ! Interpolation coefficients.
    real(8)::aaa(8)
    integer::idx(6), ia, ia1, ja, ja1, ka, ka1

    ! First, we have to check if the data point is out of range.
    if((rq .gt. grr%x1) .or.  (zq .gt. gzz%x1) .or. (zq .lt. gzz%x0))then
      out(1:3) = 0.0d0
      return
    end if

    call interpolate3D_coefficients(grr, gphi, gzz, rq, phiq, zq, aaa, idx)

    ia  = idx(1)
    ia1 = idx(2)
    ja  = idx(3)
    ja1 = idx(4)
    ka  = idx(5)
    ka1 = idx(6)
    out(1) =  Brfield(ia, ja,  ka )*aaa(1) + Brfield(ia1, ja,  ka )*aaa(2) + &
              Brfield(ia, ja1, ka )*aaa(3) + Brfield(ia1, ja1, ka )*aaa(4) + &
              Brfield(ia, ja,  ka1)*aaa(5) + Brfield(ia1, ja,  ka1)*aaa(6) + &
              Brfield(ia, ja1, ka1)*aaa(7) + Brfield(ia1, ja1, ka1)*aaa(8)

    out(2) =  Bzfield(ia, ja,  ka )*aaa(1) + Bzfield(ia1, ja,  ka )*aaa(2) + &
              Bzfield(ia, ja1, ka )*aaa(3) + Bzfield(ia1, ja1, ka )*aaa(4) + &
              Bzfield(ia, ja,  ka1)*aaa(5) + Bzfield(ia1, ja,  ka1)*aaa(6) + &
              Bzfield(ia, ja1, ka1)*aaa(7) + Bzfield(ia1, ja1, ka1)*aaa(8)

    out(3) =  Bphifield(ia, ja,  ka )*aaa(1) + Bphifield(ia1, ja,  ka )*aaa(2) + &
              Bphifield(ia, ja1, ka )*aaa(3) + Bphifield(ia1, ja1, ka )*aaa(4) + &
              Bphifield(ia, ja,  ka1)*aaa(5) + Bphifield(ia1, ja,  ka1)*aaa(6) + &
              Bphifield(ia, ja1, ka1)*aaa(7) + Bphifield(ia1, ja1, ka1)*aaa(8)
  end subroutine getField_3D

  subroutine getField_2D(rq, zq, phiq, out)
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
    real(8), intent(out)::out(3) !< Br, Bz, Bphi

    ! Interpolation coefficients.
    real(8)::aaa(4)
    integer::idx(4), ia, ia1, ja, ja1

    ! First, we have to check if the data point is out of range.
    if((rq .gt. grr%x1) .or.  (zq .gt. gzz%x1) .or. (zq .lt. gzz%x0) .or. (rq .lt. grr%x0))then
      out(1:3) = 0.0d0
      return
    end if

    call interpolate2D_coefficients(grr, gzz, rq, zq, aaa, idx)

    ia  = idx(1)
    ia1 = idx(2)
    ja  = idx(3)
    ja1 = idx(4)
    out(1) =  Brfield(ia, 1, ja)*aaa(1)   + Brfield(ia1, 1, ja)*aaa(2) &
            + Brfield(ia, 1, ja1)*aaa(3)  + Brfield(ia1, 1, ja1)*aaa(4)

    out(2) =  Bzfield(ia, 1, ja)*aaa(1)   + Bzfield(ia1, 1, ja)*aaa(2) &
            + Bzfield(ia, 1, ja1)*aaa(3)  + Bzfield(ia1, 1, ja1)*aaa(4)

    out(3) =  Bphifield(ia, 1, ja)*aaa(1)   + Bphifield(ia1, 1, ja)*aaa(2) &
            + Bphifield(ia, 1, ja1)*aaa(3)  + Bphifield(ia1, 1, ja1)*aaa(4)
  end subroutine getField_2D

  subroutine getField_0D(rq, zq, phiq, out)
    ! -----------------------------------------------------------------------
    ! GET THE FIELD of only 0D
    !> \brief Obtains Br, Bz, Bphi, interpolated in a given point,
    !! (R, z).
    ! Written by Jose Rueda
    ! -----------------------------------------------------------------------
    implicit none

    real(8), intent(in)::rq      !< Query point in R major to interpolate.
    real(8), intent(in)::zq      !< Query point in Z to interpolate.
    real(8), intent(in)::phiq    !< Query point in toroidal direciton.
    real(8), intent(out)::out(3) !< Br, Bz, Bphi

    out(1) =  Brfield(1, 1, 1)

    out(2) =  Bzfield(1, 1, 1)

    out(3) =  Bphifield(1, 1, 1)
  end subroutine getField_0D

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
    real(8), parameter :: eps = 1.0d-10

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

  subroutine checkCollision(particle)
    type(marker), intent(inout):: particle
    integer::iloop
    if (particle%kindOfCollision .eq. 1) then
      ! We already collided with the foil, check just the scintillator, which is
      ! supposed to be the last element of the geometry array
      call triangleRay(geometry(nGeomElements)%triangleNum, geometry(nGeomElements)%triangles, &
                       particle%position(:, istep), particle%position(:, istep + 1), &
                       particle%collision, particle%collision_point)
      if (particle%collision) then
        particle%kindOfCollision = geometry(nGeomElements)%kind
      endif
    else   ! We need to loop over all the plates
      plates: do iloop=1, nGeomElements
        call triangleRay(geometry(iloop)%triangleNum, geometry(iloop)%triangles, &
                         particle%position(:, istep), particle%position(:, istep + 1), &
                         particle%collision, particle%collision_point)
        if (part%collision) then
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
    real (8), dimension(3):: B, B1
    real (8), dimension(3)::field_data ! This contains the interpolated data.
    ! Note, 't' is a really bad notation for an auxiliary vector, but in most books
    ! is used, so I will maintain the name

    ! 1. We push the particle to n -> n+1/2
    r_plus_half = r0 + 0.5d0*v0*dt
    ! 2. We now have the r_{n+1/2}, let's evaluate the fields on that position.
    ! The fields, however, are stored in cylindrical coordinates, so the position
    ! must be firstly translated into that coordinate system.
    call cart2pol(r_plus_half, r_polar)
    call getField(r_polar(1), r_polar(2), r_polar(3), field_data)
    ! We have the fields in polar coordinates, but the Boris method requires the
    ! cartesian ones:
    call pol2cart_cov((/field_data(1),field_data(3), field_data(2) /), B, r_polar)
    ! 3. First half of the acceleration:
    vminus = v0

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
    v1 = vplus

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
    print*, 'reading: ', filename
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

  subroutine readGeometry(geomID, n, verb)
    ! Read the Geometry files. Note: No more than 9 files can be read. I could
    ! easily avoid this, but, are you going to define more than 9  geometric
    ! files??
    implicit none
    ! Dummy variables
    character (len=*), intent(in) :: geomID
    integer, intent(in) :: n


    logical, intent(in):: verb !< Write to console some info.

    ! Local variables
    integer :: kk,ll
    real (8), dimension(3):: vector1,vector2
    integer:: ierr, pinholeKind
    real(8), dimension(3) :: u1, u2, u3, rPin
    real(8) :: d1, d2
    character (len=1000) :: dummy_string, err_str, geometry_dir
    character (len=1) :: number


    ! Read namelist and configure the plates
    NAMELIST /ExtraGeometryParams/ u1, u2, u3, rPin, d1, d2, ps, rotation, pinholeKind

    ! Read the files
    geometry_dir = trim(SINPA_dir)//'Geometry/'

    ! read the plates
    allocate(geometry(n))
    do i = 1,n
      write(number, '(I1)') i
      dummy_string = trim(geometry_dir)//geomID//'/Element'//number//'.txt'
      call parseGeometry(trim(dummy_string),verb, geometry(i))
    enddo
    ! Read the pinhole system
    open (unit=60, file=trim(geometry_dir)//geomID//'/ExtraGeometryParams.txt',form='formatted',iostat=ierr)
    read(60, NML=ExtraGeometryParams, iostat=ierr)
    close(60)
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
      print*, 'Pinhole d1', pinhole%d1
      print*, 'Pinhole d2', pinhole%d2
    endif
  end subroutine readGeometry

!------------------------------------------------------------------------------
! SECTION 7: FIDASIM compatibility
!------------------------------------------------------------------------------
  subroutine readFIDASIM4Markers(filename)
    ! Dummy variables:
    character (*) :: filename !< Name of the file with the markers
    ! Local variables:
    real(4), dimension(:), allocatable:: dummy1D
    real(4), dimension(:, :), allocatable:: dummy2D
    print*, 'Reading markers info from: ', trim(filename)
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
    F4Markers%ipos = transpose(dble(dummy2D))
    read(60) dummy2D
    F4Markers%fpos = transpose(dble(dummy2D))
    read(60) dummy2D
    F4Markers%v = transpose(dble(dummy2D))
    read(60) dummy1D
    F4Markers%wght = dble(dummy1D)
    read(60) dummy1D
    F4Markers%kind = int(dummy1D)
    print*, F4Markers%counter, 'markers to be follwed'
  end subroutine readFIDASIM4Markers

 !------------------------------------------------------------------------------
 ! SECTION 8: Initialization of markers
 !------------------------------------------------------------------------------
 subroutine initMarker(index, vmod, dt_initial, xi)
   !----
   integer, intent(in):: index
   real(8), intent(in):: vmod
   real(8), intent(in):: dt_initial
   real(8), intent(in):: xi
   !---
   real(8):: distance
   ! Clean the marker position and trajectory
   part%position(:, :) = 0.0d0
   part%velocity(:, :) = 0.0d0
   part%collision    = .False.
   part%kindOfCollision = 9  ! Just initialise it to a dummy value
   part%dt = dt_initial
   part%qm    = Zin *qe / M /amu_to_kg
   part%xi = xi

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
      part%position(:, 1) = rPin + &
          pinhole%d1 * (-1+2*ran(1))*pinhole%u1 + &
          pinhole%d2 * (-1+2*ran(2))*pinhole%u2
    endif
    ! Initialise the random angle
    part%beta = minAngle + dAngle*ran(3)
    ! Init marker velocity
    if (FILDSIMmode) then
      print*, 'no implemented'
    else
      part%velocity(:, 1) = ((cos(xi)*pinhole%u3 + &
                              sin(xi)*pinhole%u1) * cos(part%beta) &
                           + sin(part%beta)*pinhole%u2) * vmod

      part%weight         = abs(sin(part%beta))
    endif
 end subroutine initMarker
end module sinpa_module
