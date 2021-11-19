!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ___  ____  _  _  ____   __       ___  _____  ____  ____
! / __)(_  _)( \( )(  _ \ /__\     / __)(  _  )(  _ \( ___)
! \__ \ _)(_  )  (  )___//(__)\   ( (__  )(_)(  )(_) ))__)
! (___/(____)(_)\_)(__) (__)(__)   \___)(_____)(____/(____)
!------------------------------------------------------------------------------
! -----------------------------------------------------------------------------
!  SINPA CODE
! -----------------------------------------------------------------------------
! TITLE         : SINPA
! PROJECT       : Synthetic INPA
! MODULE        : Main Core
! AFFILIATION   : University of Sevilla
!> \author Jose Rueda - Universidad de Sevilla
!> \date 01/10/2021
!> \version 0.0
!> \see https://gitlab.mpcdf.mpg.de/ruejo/sinpa
!
! DESCRIPTION:
!> \brief Main SINPA core
! -----------------------------------------------------------------------------
program sinpa
  !****************************************************************************
  ! --- Note on units to open files
  !       - 60: Temporal files: Inputs/output. Open, read(write) and close it
  !       - 61: File with the strike points for the mapping (signal)
  !       - 62: Colimator data
  !       - 63: Orbit file
  !****************************************************************************
  use sinpa_module
  implicit none
  !----------------------------------------------------------------------------
  !============================================================================
  ! ROUND 1: READING OF NAMELIST AND INITIALIZATION OF THE PROGRAM
  !============================================================================
  !----------------------------------------------------------------------------
  ! When the program it's called, the configuration file is passed as and input
  ! First of all we recovered the name of the configuration file
  call getarg(1, input_filename)

  ! Open and read the configuration namelist.
  open(unit=60, file=input_filename, form='formatted', iostat=ierr)
  read(60, NML=config, iostat=ierr)
  close(60)

  if (ierr /= 0) THEN
    print*,'Error when reading the input filename: ',input_filename
    print*,'Error in NAMELIST config: ',ierr
    stop
  end if
  ! Open and read the input_namelist. Even if all namelist are in the same file
  ! we open and close the namelist file each time because the python routine
  ! which generate the namelist can write the different namelist in different
  ! order in the file.
  ! Allocate the gyroradius and pitch arrays to be readed
  allocate (rl(nGyroradius))
  allocate (XI(nxi))
  ! Read the input namelist
  open (unit=60, file=input_filename, form='formatted', iostat=ierr)
  read(60, NML=inputParams, iostat = ierr)
  close(60)

  ! Read the specifict namelist depending on the diagnostic
  if (FILDSIMmode.eqv..False.) then
    open (unit=60, file=input_filename, form='formatted', iostat=ierr)
    read(60, NML=nbi_namelist, iostat = ierr)
    nbi%u = u
    nbi%p0 = p0
    close(60)
    ! Markers interaction (energy loss and scattering) parameters
    open (unit=60, file=input_filename, form='formatted', iostat=ierr)
    read(60, NML=markerinteractions, iostat=ierr)
    close(60)
  end if

  if (verbose) then
    print*,'------------------------------------------------------------------'
    print*,'Input parameters read from: ',trim(input_filename)
    print*, 'runid:', runID
    print*,'------------------------------------------------------------------'
  endif
  ! --- Set the geometry
  call readGeometry(trim(geomID), nGeomElements, verbose)
  if (verbose) then
    print*,'------------------------------------------------------------------'
    print*,'Setting up the geometry of the plates and scintillator'
    print*,'Geometry setup finished'
    print*,'------------------------------------------------------------------'
  endif

  ! --- Load the grid and magnetic field
  call parseField(trim(SINPA_dir)//'/runs/'//trim(runID)//'/inputs/field.bin',&
                  trim(SINPA_dir)//'/runs/'//trim(runID)//'/inputs/Efield.bin',&
                  verbose)
  call Bsystem()
  if (verbose) then
    print*, '-----------------------------------------------------------------'
    print*, 'Magnetic field at the pinhole: ', Bpinhole
    print*, 'Mod: ', BpinholeMod
    print*, '-----------------------------------------------------------------'
  endif

  ! --- Translate from FILDSIM pitch criteria to decent one
  if (FILDSIMmode) then
    ! Keep the crazy FILDSIM criteria
    XI = dble(IpBt) * cos(XI*pi/180.0d0)
  endif
  ! --- Caclualte the beta (gyrophases) for all X values
  call calculate_betas()
  !-----------------------------------------------------------------------------
  !=============================================================================
  ! ROUND 2: MAPPING
  !=============================================================================
  !-----------------------------------------------------------------------------
  mapeado: if (mapping) then
    if (verbose) then
      print*, '-----------------------------------------------------------------'
      print*,'Performing mapping'
    endif
    ! --- Allocate the necesary matrices
    if (FILDSIMmode) then
      allocate(Strike(12,nMap))            ! Strike points in the scint
      allocate(StrikeMap(9,nGyroradius * nxi))  ! Strike map
    else
      allocate(Strike(19,nMap))            ! Strike points in the scint
      allocate(StrikeMap(13,nGyroradius * nxi))  ! Strike map
    endif
    allocate(CollimatorStrikes(4,nMap))       ! Strike position on the coll
    ! Get the size (to later write it in the file)
    dummy_shape = shape(Strike)
    ! --- Open the files to save the data
    ! -- Strike points on the scintillator
    open(unit=61, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/StrikePoints.bin', access = 'stream', action='write', status='replace')
    ! Save the header of the file
    write(61) versionID1, versionID2, runID, nGyroradius, rL, nxi, XI, transfer(FILDSIMmode, 1), dummy_shape(1)
    ! -- Strike points on the collimator
    open(unit=62, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/CollimatorStrikePoints.bin', access = 'stream', action='write', status='replace')
    write(62) versionID1, versionID2, runID, nGyroradius, rL, nxi, XI, 4
    ! -- File to save the orbits
    if (saveOrbits) then
      open(unit=63, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
           '/results/Orbits.bin', access = 'stream', action='write', status='replace')
      write(63) versionID1, versionID2, runID
      if (saveOrbitLongMode) then
        write(63) 69
      else
        write(63) 6
      endif
    endif
    ! --- Allocate the particle:
    call omega(M, Zout, BpinholeMod, OmegaPart)  ! Gyrofrequency
    dt = 2 * pi / OmegaPart / nGyro
    part%n_t = int(maxT/dt)
    if (verbose) then
      print*, part%n_t, 'steps will be performed for each particle'
    endif
    allocate(part%position(3,part%n_t))
    allocate(part%velocity(3,part%n_t))
    Lenergies: do irl = 1, nGyroradius   ! Loop over energies (gyroradii)
      ! Get a dummy velocity modulus
      call rl2v0(rL(irl), M, Zout, BpinholeMod, v0)

      ! Set the dt to avoid doing unnnecesary steps for each markers
      if (Zin .gt. 0.1) then  ! if we launched ions
        dtNeutral = dt
      else
        dtNeutral = 0.04 / v0 ! if we launched neutrals
      endif
      LXI: do iXI = 1, nxi  ! Loop over the alpha angles
        ! -- Initialise all the counters
        cScintillator = 0
        cCollimator = 0
        cFoil = 0
        ! -- Clean the matrices with the complete information
        CollimatorStrikes(:,:) = 0.0d0
        Strike(:,:) = 0.0d0

        Lmarkers: do  imc = 1, nMap
          call initMarker(v0, dtneutral, XI(iXI), min_beta(iXI), delta_beta(iXI), rL(irl))

          tracking: do istep = 1, part%n_t-1
            call pushParticle(part%qm, part%position(:, istep), &
                              part%velocity(:, istep), &
                              part%position(:, istep + 1),&
                              part%velocity(:, istep + 1), part%dt)
            call checkCollision(part)
            ! If it has collided, do stuff
            if (part%collision) then
              if (part%kindOfCollision .eq. 0) then
                ! Collided with the collimator, dead marker
                cCollimator = cCollimator + 1
                CollimatorStrikes(1:3, cCollimator) = part%collision_point
                CollimatorStrikes(4, cCollimator) = part%weight


                if (saveOrbits) then
                  call random_number(rand_number)
                  if (rand_number .lt. saveRatio) then
                    ! Save the collision point in the trajectory for future use
                    part%position(:, istep+1) = part%collision_point
                    cOrb = cOrb + 1
                    call saveOrbit(63, part, istep+1, saveOrbitLongMode)
                  endif
                endif
                cycle Lmarkers
              elseif (part%kindOfCollision .eq. 2) then ! Scintillator collision
                cScintillator = cScintillator + 1 ! Update counter
                ! Save the common FILD and INPA variables:
                ! Store the other information of the marker
                Strike(1:3, cScintillator ) = part%collision_point ! f point
                Strike(4, cScintillator) = part%weight ! weight
                Strike(5, cScintillator) = part%beta ! Beta angle
                Strike(6:8, cScintillator) = part%position(:, 1) ! i pos
                Strike(9:11, cScintillator ) = &
                  MATMUL(rotation, part%collision_point - ps) ! f pos scint
                Strike(12, cScintillator) = &
                    180/pi*acos(sum(ScintNormal*part%velocity(:, istep))&
                    /norm2(part%velocity(:, istep))) ! incident angle
                ! Save INPA extra varialbes
                if (FILDSIMmode.eqv..False.) then
                  ! ToDo: Is the energy loss in the foil is included, change this
                  ! ratio vv0/v0 to the real velocity
                  call minimumDistanceLines(part%position(:, 1), &
                                            part%velocity(:, 1)/v0, nbi%p0,&
                                            nbi%u, dMin, posMin)

                  Strike(13:15,cScintillator) = posMin   ! Initial point (NBI)
                  Strike(16:18,cScintillator) = part%velocity(:, 1)  ! Initial velocity
                  Strike(19, cScintillator) = dMin ! Closer distance to NBI
                endif

                if (saveOrbits) then
                  call random_number(rand_number)
                  if (rand_number .lt. saveRatio) then
                    ! Save the collision point in the trajectory for future use
                    part%position(:, istep+1) = part%collision_point
                    cOrb = cOrb + 1
                    call saveOrbit(63, part, istep+1, saveOrbitLongMode)

                  endif
                endif
                cycle Lmarkers
              else ! we collided with the foil
                part%collision = .False.
                part%dt = dt   ! Good dt to evolve an ion
                part%qm = Zout * qe / M /amu_to_kg ! New charge after the foil
                part%position(:, istep + 1) = part%collision_point
                r0 = part%collision_point
                cFoil = cFoil + 1 ! Update the counter
              endif
            endif
          enddo tracking
          ! if we achieve this point, the particle has not collide, save the
          ! orbit for future analysis
          if (saveOrbits) then
            call random_number(rand_number)
            if (rand_number .lt. saveRatio) then
              ! Save the collision point in the trajectory for future use
              cOrb = cOrb + 1
              call saveOrbit(63, part, part%n_t, saveOrbitLongMode)
            endif
          endif
        enddo Lmarkers
        ! ! Rotate the data
        ! Strike(16:18, 1:cScintillator ) = &
        !   MATMUL(rotation, transpose(Strike(1:3, 1:cScintillator)) -ps)
        ! Write the markers in the file
        write(61) cScintillator, transpose(Strike(:, 1:cScintillator))
        write(62) cCollimator, transpose(CollimatorStrikes(:, 1:cCollimator))
        ! Print some information
        if (verbose) then
          print*, '-----'
          if (FILDSIMmode) then
            print*, 'Gyroradius:', rL(irl), ' Alpha:', 180.0*acos(XI(iXI)/IpBt)/pi
          else
            print*, 'Gyroradius:', rL(irl), ' Alpha:', XI(iXI)
          endif
          print*, 'Hitting Collimator', cCollimator
          print*, 'Hitting Scintillator', cScintillator
          print*, 'Wrong Neutrals', cWrongNeutral
          print*, 'Wrong Ions', cWrongIons
          print*, 'Foil', cFoil
        endif
        ! Save the strike map (average everything)
        ! ToDo: correct collimator factor for the interval or launched gyrophase
        StrikeMap(1, iXI + (irl - 1) * nxi) = rL(irl)  ! Gyroradius
        StrikeMap(2, iXI + (irl - 1) * nxi) = XI(iXI)  ! Pitch or NBI angle
        StrikeMap(3, iXI + (irl - 1) * nxi) = &   ! Strike point at Scint
          sum(Strike(9, :)) / cScintillator
        StrikeMap(4, iXI + (irl - 1) * nxi) = &
          sum(Strike(10, :)) / cScintillator
        StrikeMap(5, iXI + (irl - 1) * nxi) = &
          sum(Strike(11, :)) / cScintillator
        StrikeMap(6, iXI + (irl - 1) * nxi) = &
          sum(Strike(5, :)) / cScintillator   ! Average initial gyrophase (beta)
        StrikeMap(7, iXI + (irl - 1) * nxi) = cScintillator ! N strikes
        StrikeMap(8, iXI + (irl - 1) * nxi) = & ! Collimator factor
          100.0d0 * dble(cScintillator) / dble(nMap)
        StrikeMap(9, iXI + (irl - 1) * nxi) = &
          sum(Strike(12, :)) / cScintillator   ! Average incident angle
        if (FILDSIMmode.eqv..False.) then
          StrikeMap(10, iXI + (irl - 1) * nxi) = &  ! Birth position (NBI)
            sum(Strike(13, :)) / cScintillator
          StrikeMap(11, iXI + (irl - 1) * nxi) = &
            sum(Strike(14, :)) / cScintillator
          StrikeMap(12, iXI + (irl - 1) * nxi) = &
            sum(Strike(15, :)) / cScintillator
          StrikeMap(13, iXI + (irl - 1) * nxi) = &  ! Distance to NBI line
            sum(Strike(19, :)) / cScintillator
        endif
        ! De-allocate the variables
      enddo LXI
    enddo Lenergies
    ! Write the strike map
    print*,'Saving strike map in: '
    print*,trim(SINPA_dir)//'/runs/'//trim(runID)//'/results/StrikeMap.txt'
    open(unit=60, file=trim(SINPA_dir)//'runs/'//trim(runID)//&
         '/results/StrikeMap.txt', action='write', form='formatted')
    ! Write the header
    write(60,'(a,2x,a)') 'Simulation NAMELIST: ', input_filename
    write(60,'(a)') 'Dummy temporal line'
    if (FILDSIMmode) then
      write(60,'(a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a)') 'Gyroradius (cm)', &
        'Pitch-Angle (degree)', 'X (cm)', 'Y (cm)', 'Z (cm)',&
        'Average_initial_gyrophase','N_strike_points', 'Collimator_Factor (percent)', 'Average_incidence_angle'
      ! Revert the pitch to the criteria useed in old fildsim
      StrikeMap(2, :) = 180.0d0/pi*acos(StrikeMap(2, :)/dble(IpBt))
      write(60,'(9F14.6)') StrikeMap
    else
      write(60,'(a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a, 2x,a 2x,a)') 'Gyroradius (cm)', &
              'Alpha [rad]', 'X [cm]', 'Y [cm]', 'Z [cm]', &
              'X0 [cm]','Y0[cm]', 'Z0[cm]','d0[cm]','Collimator_Factor (%)', 'nMarkers'
      ! Write the data
      write(60,'(13F14.6)') StrikeMap
    endif
    ! Close all the openend files
    close(60)
    close(61)
    close(62)
    if (saveOrbits) then
      write(63) m
      write(63) cOrb  ! Write how many orbits we wrote in the file
      close(63)
    endif

    !De-allocate
    deallocate(StrikeMap)
    deallocate(part%position)
    deallocate(part%velocity)
    deallocate(Strike)
    deallocate(CollimatorStrikes)
  endif mapeado
  !-----------------------------------------------------------------------------
  !=============================================================================
  ! ROUND 3: FIDASIM signal
  !=============================================================================
  !-----------------------------------------------------------------------------
  signal_part: if (signal) then
    ! Fist load the FIDASIM simulation
    call readFIDASIM4Markers(trim(FIDASIMfolder)//"/npa.bin")
    ! --- Open the file to save the data
    ! -- Strike points on the scintillator
    open(unit=61, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/SignalStrikePoints.bin', access = 'stream', action='write')
    ! Save the header of the file
    write(61) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0, 15
    ! -- Strike points on the collimator
    open(unit=62, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/SignalCollimatorStrikePoints.bin', access = 'stream', action='write')
    write(62) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0,4
    ! -- File to save the orbits
    if (saveOrbits) then
      open(unit=63, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
           '/results/SignalOrbits.bin', access = 'stream', action='write')
      write(63) versionID1, versionID2, runID
    endif
    ! --- Allocate the particle:
    call omega(M, Zout, BpinholeMod, OmegaPart)  ! Gyrofrequency
    dt = 2 * pi / OmegaPart / nGyro
    part%n_t = int(maxT/dt)
    dtNeutral = 4.0E-7
    if (verbose) then
      print*, part%n_t, 'steps will be performed for each particle'
    endif
    allocate(part%position(3,part%n_t))
    allocate(part%velocity(3,part%n_t))
    ! --- Allocate the necesary matrix
    allocate(Strike(15,F4Markers%counter))            ! Strike points in the scint
    allocate(CollimatorStrikes(4,F4Markers%counter))       ! Strike position on the coll
    ! -- Initialise all the counters
    cScintillator = 0
    cCollimator = 0
    cWrongIons = 0
    cWrongNeutral = 0
    cFoil = 0
    markers: do i=1, F4Markers%counter
      ! Initialise the markers
      ! - Clean the marker data
      part%position(:, :) = 0.0d0
      part%velocity(:, :) = 0.0d0
      part%collision     = .False.
      part%kindOfCollision = 9  ! Just initialise it to a dummy value
      ! - FIDASIM data
      part%weight         = F4Markers%wght(i)
      part%position(:,1) = F4Markers%fpos(:, i)
      part%velocity(:,1) = F4Markers%v(:, i)
      ! - dt and charge
      part%dt    = dtNeutral
      part%qm    = Zin *qe / M /amu_to_kg
      signal_tracking: do istep = 1, part%n_t-1
        call pushParticle(part%qm, part%position(:, istep), &
                          part%velocity(:, istep), &
                          part%position(:, istep + 1),&
                          part%velocity(:, istep + 1), part%dt)
        call checkCollision(part)
        ! If it has collided, do stuff
        if (part%collision) then
          if (part%kindOfCollision .eq. 0) then
            ! Collided with the collimator, dead marker
            cCollimator = cCollimator + 1
            CollimatorStrikes(1:3, cCollimator) = part%collision_point
            CollimatorStrikes(4, cCollimator) = part%weight
            cycle markers
            if (saveOrbits) then
              call random_number(rand_number)
              if (rand_number .lt. saveRatio) then
                cOrb = cOrb + 1
                write(63) istep
                write(63) transpose(part%position(:, 1:istep))
              endif
            endif
          elseif (part%kindOfCollision .eq. 2) then ! Scintillator collision
            cScintillator = cScintillator + 1 ! Update counter
            Strike(1, cScintillator ) = dble(i) ! FIDASIM marker id
            Strike(2:4, cScintillator ) = part%collision_point ! f point
            Strike(5:7, cScintillator ) = &
              MATMUL(rotation, part%collision_point - ps) !f point in scintillator region
            Strike(8:10, cScintillator) = F4Markers%ipos(:, i) ! CX position
            Strike(11:13,cScintillator) = part%velocity(:, 1)   ! Initial velocity
            Strike(14,cScintillator) = part%weight ! weight
            Strike(15, cScintillator) = F4Markers%kind(i) ! Kind of signal
            if (saveOrbits) then
              call random_number(rand_number)
              if (rand_number .lt. saveRatio) then
                cOrb = cOrb + 1
                write(63) istep
                write(63) transpose(part%position(:, 1:istep))
              endif
            endif
            cycle markers
          else ! we collided with the foil
            part%collision = .False.
            part%dt = dt   ! Good dt to evolve an ion
            part%qm = Zout * qe / M /amu_to_kg ! New charge after the foil
            part%position(:, istep + 1) = part%collision_point
            r0 = part%collision_point
            cFoil = cFoil + 1 ! Update the counter
          endif
        endif
      enddo signal_tracking
    end do markers
    write(61) cScintillator, transpose(Strike(:, 1:cScintillator))
    write(62) cCollimator, transpose(CollimatorStrikes(:, 1:cCollimator))
    if (verbose) then
      print*, '-----'
      print*, 'Hitting Collimator', cCollimator
      print*, 'Hitting Scintillator', cScintillator
      print*, 'Wrong Neutrals', cWrongNeutral
      print*, 'Wrong Ions', cWrongIons
      print*, 'hitting foil', cFoil
    endif
    close(61)
    close(62)
    deallocate(part%position)
    deallocate(part%velocity)
    deallocate(Strike)
    deallocate(CollimatorStrikes)
  endif signal_part
end program sinpa
