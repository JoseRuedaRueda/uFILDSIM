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
!> \date 21/04/2022
!> \version 1.4
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
  !       - 64: Wrong markers file
  !       - 65: Strike position of self-shadowed markers
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
  open(unit=60, file=trim(input_filename), form='formatted', iostat=ierr)
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
  allocate (XI_input(nxi))
  ! Read the input namelist
  open (unit=60, file=input_filename, form='formatted', iostat=ierr)
  read(60, NML=inputParams, iostat = ierr)
  close(60)
  ! Sort the arrays
  call quicksort_nr(XI)
  call quicksort_nr(rl)
  ! Read the specifict namelist depending on the diagnostic
  if (FILDSIMmode.eqv..False.) then
    open (unit=60, file=input_filename, form='formatted', iostat=ierr)
    read(60, NML=nbi_namelist, iostat = ierr)
    close(60)
    nbi%u = u
    if (abs(norm2(u) - 1.0d0) .gt. 0.01) then
      stop 'NBI director vector is not well defined (not unitary)'
    endif
    nbi%p0 = p0
    ! Markers interaction (energy loss, tramisssion and scattering) parameters
    if (FoilElossModel .eq. 1) then
      allocate(FoilElossParameters(2))
      if (verbose) then
        print*,'Empirical energy loss will be applied'
      endif
    elseif (FoilElossModel .eq. 2) then
      allocate(FoilElossParameters(3))
      if (verbose) then
        print*,'SRIM energy loss will be applied'
      endif
    elseif (FoilElossModel .eq. 0) then
      if (verbose) then
        print*,'No model will be applied for carbon foil energy loss'
      endif
    else
      stop 'Interaction model not understood'
    endif
    if (ScintillatorYieldModel .eq. 1) then
      allocate(ScintillatorYieldParameters(1))
      if (verbose) then
        print*,'Just proportionality factor will be applied'
      endif
    elseif (ScintillatorYieldModel .eq. 2) then
      allocate(ScintillatorYieldParameters(2))
      if (verbose) then
        print*,'Birk model will be assumed'
      endif
    elseif (ScintillatorYieldModel .eq. 0) then
      if (verbose) then
        print*,'No model will be applied for scintillator efficiency'
      endif
    else
      stop 'Yield model not understood'
    endif
    if (FoilYieldModel .eq. 1) then
      allocate(FoilYieldParameters(1))
      if (verbose) then
        print*,'Just proportionality factor will be applied for foil transmission'
      endif
    elseif (FoilYieldModel .eq. 2) then
      allocate(FoilYieldParameters(4))
      if (verbose) then
        print*,'Empirical ionization yield applied'
      endif
    elseif (FoilYieldModel .eq. 0) then
      if (verbose) then
        print*,'No model will be applied for foil yield'
      endif
    else
      stop 'Yield model not understood'
    endif
    if ((FoilElossModel .ne. 0) .or. (FoilYieldModel .ne. 0) .or. (ScintillatorYieldModel .ne. 0))then
      open (unit=60, file=input_filename, form='formatted', iostat=ierr)
      read(60, NML=markerinteractions, iostat=ierr)
      close(60)
    endif

  end if

  if (verbose) then
    print*,'------------------------------------------------------------------'
    print*,'Input parameters read from: ',trim(input_filename)
    print*, 'runid:', runID
    print*,'------------------------------------------------------------------'
  endif
  ! --- Set the geometry
  call readGeometry(trim(GeomFolder), nGeomElements, verbose)
  if (verbose) then
    print*,'------------------------------------------------------------------'
    print*,'Setting up the geometry of the plates and scintillator'
    print*,'Geometry setup finished'
    print*,'------------------------------------------------------------------'
  endif

  ! --- Load the grid and magnetic field
  call parseField(trim(runFolder)//'/inputs/field.bin',&
                  trim(runFolder)//'/inputs/Efield.bin',&
                  verbose)
  call Bsystem()
  if (verbose) then
    print*, '-----------------------------------------------------------------'
    print*, 'Magnetic field at the pinhole: ', Bpinhole
    print*, 'Mod: ', BpinholeMod
    print*, '-----------------------------------------------------------------'
  endif

  ! --- Translate from FILDSIM pitch criteria to decent one
  XI_input = XI
  if (FILDSIMmode) then
    ! Keep the crazy FILDSIM criteria
    XI = dble(IpBt) * cos(XI_input*pi/180.0d0)
  endif
  ! --- Caclualte the beta (gyrophases) for all X values
  if (minAngle .le. 0.0d0) then
    minAngle = minAngle + 2.0d0*pi
  endif
  call calculate_betas(restrict_mode)
  ! --- Stablish the sign of the dt
  if (backtrace) then
    time_sign = -1.0d0
  else
    time_sign = 1.0d0
  endif
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
      allocate(StrikeMap(9,nGyroradius * nxi))  ! Strike map
      dummy_shape = 12    ! Get the size (to later write it in the file)
    else
      allocate(StrikeMap(13,nGyroradius * nxi))  ! Strike map
      dummy_shape = 19    ! Get the size (to later write it in the file)
    endif

    ! --- Open the files to save the data
    ! -- Strike points on the scintillator
    if (save_scintillator_strike_points) then
      open(unit=61, file=trim(runFolder)//&
           '/results/'//trim(runID)//'.spmap', access = 'stream', action='write', status='replace')
      ! Save the header of the file
      write(61) versionID1, versionID2, runID, nGyroradius, rL, nxi, XI_input, transfer(FILDSIMmode, 1), dummy_shape
    endif
    ! -- Strike points on the collimator
    if (save_collimator_strike_points) then
      open(unit=62, file=trim(runFolder)//&
           '/results/'//trim(runID)//'.spcmap', access = 'stream', action='write', status='replace')
      write(62) versionID1, versionID2, runID, nGyroradius, rL, nxi, XI_input, transfer(FILDSIMmode, 1), 4
    endif
    if (save_self_shadowing_collimator_strike_points) then
      open(unit=65, file=trim(runFolder)//&
           '/results/'//trim(runID)//'.spcself', access = 'stream', action='write', status='replace')
      write(65) versionID1, versionID2, runID, nGyroradius, rL, nxi, XI_input, transfer(FILDSIMmode, 1), 4
    endif
    ! -- File to save the orbits
    if (saveOrbits) then
      open(unit=63, file=trim(runFolder)//&
           '/results/'//trim(runID)//'.orb', access = 'stream', action='write', status='replace')
      write(63) versionID1, versionID2, runID
      if (saveOrbitLongMode) then
        write(63) 69
      else
        write(63) 6
      endif
    endif
    ! --  File to save the wrong markers
    if (save_wrong_markers_position) then
      open(unit=64, file=trim(runFolder)//&
           '/results/'//trim(runID)//'.wmmap', access = 'stream', action='write', status='replace')
      write(64) versionID1, versionID2, runID, nGyroradius, rL, nxi, XI_input, transfer(FILDSIMmode, 1), 4
    endif
    ! --- Allocate the particle:
    call cpu_time(t_initial_orbits)
    call omega(M, Zout, BpinholeMod, OmegaPart)  ! Gyrofrequency
    dt = 2 * pi / OmegaPart / nGyro * time_sign
    part%n_t = abs(int(maxT/dt))
    if (verbose) then
      print*, part%n_t, 'steps will be performed for each particle'
    endif
    allocate(part%position(3,part%n_t))
    allocate(part%velocity(3,part%n_t))
    if (self_shadowing) then
      backPart%n_t = nGyro * 3
      part%dt = dt
      allocate(backPart%position(3,backPart%n_t))
      allocate(backPart%velocity(3,backPart%n_t))
    endif
    Lenergies: do irl = 1, nGyroradius   ! Loop over energies (gyroradii)
      ! Get a dummy velocity modulus
      call rl2v0(rL(irl), M, Zout, BpinholeMod, v0)

      ! Set the dt to avoid doing unnnecesary steps for each markers
      if (Zin .gt. 0.1) then  ! if we launched ions
        dt1 = dt
      else
        dt1 = 0.04 / v0 * time_sign  ! if we launched neutrals
      endif
      LXI: do iXI = 1, nxi  ! Loop over the alpha angles
        ! -- Initialise all the counters
        cScintillator = 0
        cCollimator = 0
        backCollimator = 0
        cFoil = 0
        cWrong = 0
        ! -- Clean the matrices with the complete information

        nToLaunch = int(dble(nMap) * marker_factor(iXI) * (1.0d0 + n1/(rL(irl)-r1)**2))
        if (verbose) then
          print*,'---'
          print*, 'Following: ', nToLaunch, ' markers'
        endif
        if (FILDSIMmode) then
          allocate(Strike(12,nToLaunch+1))            ! Strike points in the scint
        else
          allocate(Strike(19,nToLaunch+1))            ! Strike points in the scint
        endif
        allocate(CollimatorStrikes(4,nToLaunch+1))       ! Strike position on the coll
        allocate(WrongMarkers(4,nToLaunch+1))       ! Strike position on the coll
        CollimatorStrikes(:,:) = 0.0d0
        Strike(:,:) = 0.0d0
        WrongMarkers(:,:) = 0.0d0
        if (self_shadowing) then
          allocate(BackCollimatorStrikes(4,nToLaunch+1))       ! Strike position on the coll
          BackCollimatorStrikes(:,:) = 0
        endif
        Lmarkers: do  imc = 1, nToLaunch
          call initMarker(v0, dt1, XI(iXI), min_beta(iXI), delta_beta(iXI), rL(irl))

          tracking: do istep = 1, part%n_t-1
            call pushParticle(part%qm, part%position(:, istep), &
                              part%velocity(:, istep), &
                              part%position(:, istep + 1),&
                              part%velocity(:, istep + 1), part%dt)
            call checkCollision(part, istep)
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
                incidentProjection = sum(ScintNormal*part%velocity(:, istep))&
                  /norm2(part%velocity(:, istep))
                if (incidentProjection .gt. 0) then
                  ! neglect marker
                  cycle Lmarkers
                endif
                ! If we reach this point, we have collided with the
                ! scintillator, so if there is the self_shadowing, now launch
                ! the particle backwards
                if (self_shadowing) then
                  ! Clean the backtracing particle
                  backPart%position(:, :) = 0.0d0
                  backPart%velocity(:, :) = 0.0d0
                  ! Set the initial position and velocity
                  backPart%position(:, 1) = part%position(:, 1)
                  backPart%velocity(:, 1) = part%velocity(:, 1)
                  ! Initialise the other varaibles
                  backPart%collision    = .False.
                  backPart%kindOfCollision = 9  ! Just initialise it to a dummy value
                  backPart%dt = -part%dt
                  backPart%qm = Zin *qe / M /amu_to_kg
                  ! Trace if and check collisions with the collimator
                  backtracking: do iiistep = 1, backPart%n_t-1
                    call pushParticle(backPart%qm, backPart%position(:, iiistep), &
                                      backPart%velocity(:, iiistep), &
                                      backPart%position(:, iiistep + 1),&
                                      backPart%velocity(:, iiistep + 1), backPart%dt)
                    call checkCollision(backPart, iiistep)
                    if (backPart%collision) then
                      if (backPart%kindOfCollision .eq. 0) then
                        ! Collided with the collimator, dead marker
                        backCollimator = backCollimator + 1
                        backCollimatorStrikes(1:3, backCollimator) = backPart%collision_point
                        backCollimatorStrikes(4, backCollimator) = backPart%weight
                        cycle Lmarkers

                      endif
                    end if
                  end do backtracking
                endif


                call yieldScintillator(part, istep)
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
                    180/pi*acos(incidentProjection) ! incident angle
                ! Save INPA extra varialbes
                if (FILDSIMmode.eqv..False.) then
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
              elseif (part%kindOfCollision .eq. 1) then ! we collided with the foil
                part%collision = .False.
                part%dt = dt   ! Good dt to evolve an ion
                part%qm = Zout * qe / M /amu_to_kg ! New charge after the foil
                part%position(:, istep + 1) = part%collision_point
                call foilInteraction(part,foilNormal, istep )
                r0 = part%collision_point
                cFoil = cFoil + 1 ! Update the counter
              else
                print*, 'Warning we collide with element of kind', part%kindOfCollision
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
          ! Save the final point of the marker
          cWrong = cWrong + 1
          WrongMarkers(1:3, cWrong) = part%position(:, istep)
          WrongMarkers(4, cWrong) = part%weight
        enddo Lmarkers
        ! ! Rotate the data
        ! Strike(16:18, 1:cScintillator ) = &
        !   MATMUL(rotation, transpose(Strike(1:3, 1:cScintillator)) -ps)
        ! Write the markers in the file
        if (save_scintillator_strike_points) then
          write(61) cScintillator, transpose(Strike(:, 1:cScintillator))
        endif
        if (save_collimator_strike_points) then
          write(62) cCollimator, transpose(CollimatorStrikes(:, 1:cCollimator))
        endif
        if ((save_self_shadowing_collimator_strike_points).and.(self_shadowing)) then
          write(65) backCollimator, transpose(backCollimatorStrikes(:, 1:backCollimator))
        endif
        if (save_wrong_markers_position) then
          write(64) cWrong, transpose(WrongMarkers(:, 1:cWrong))
        endif
        ! Print some information
        if (verbose) then
          if (FILDSIMmode) then
            ! The +0.001 is just a trick to avoid rounding issues
            print*, 'Gyroradius:', rL(irl), ' Pitch:', int(180.0*acos(XI(iXI)/IpBt)/pi+0.001)
          else
            print*, 'Gyroradius:', rL(irl), ' Alpha:', XI(iXI)
            print*, 'Hitting Foil', cFoil
          endif
          print*, 'Hitting Collimator', cCollimator
          print*, 'Hitting Scintillator', cScintillator
          print*, 'Not colliding', nToLaunch - cCollimator - cScintillator
          if (self_shadowing) then
            print*, 'Back colliding', backCollimator
          endif
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
          100.0d0 * dble(cScintillator) / dble(nToLaunch) * delta_beta(iXI) / 2.0d0 / pi
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
        deallocate(Strike)
        deallocate(CollimatorStrikes)
        deallocate(WrongMarkers)
        deallocate(backCollimatorStrikes)
      enddo LXI
    enddo Lenergies
    call cpu_time(t_final_orbits)
    print*,'It took (in seconds): ', t_final_orbits - t_initial_orbits
    ! Write the strike map
    print*,'Saving strike map in: '
    print*,trim(runFolder)//'/results/'//trim(runID)//'.map'
    open(unit=60, file=trim(runFolder)//'/results/'&
      //trim(runID)//'.map', action='write', form='formatted')
    ! Write the header
    write(60,'(a,2x,a)') 'Simulation NAMELIST: ', trim(input_filename)
    write(60,'(i2, i2, a)') versionID1, versionID2, ' # SINPA version'
    if (FILDSIMmode) then
      write(60,'(a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a)') 'Gyroradius (cm)', &
        'Pitch-Angle (degree)', 'X (m)', 'Y (m)', 'Z (m)',&
        'Average_initial_gyrophase','N_strike_points', 'Collimator_Factor (percent)', 'Average_incidence_angle'
      ! Revert the pitch to the criteria useed in old fildsim
      StrikeMap(2, :) = 180.0d0/pi*acos(StrikeMap(2, :)/dble(IpBt))
      write(60,'(9F14.6)') StrikeMap
    else
      write(60,'(a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a, 2x,a 2x,a)') 'Gyroradius [cm]', &
              'Alpha [rad]', 'X [m]', 'Y [m]', 'Z [m]', &
              'X0 [m]','Y0[m]', 'Z0[m]','d0[m]','Collimator_Factor (%)', 'nMarkers'
      ! Write the data
      write(60,'(13F14.6)') StrikeMap
    endif
    ! Close all the openend files
    close(60)
    if (save_scintillator_strike_points) then
      close(61)
    endif
    if (save_collimator_strike_points) then
      close(62)
    endif
    if (saveOrbits) then
      write(63) m
      write(63) cOrb  ! Write how many orbits we wrote in the file
      close(63)
    endif
    if (save_wrong_markers_position) then
      close(64)
    endif
    if (save_self_shadowing_collimator_strike_points) then
      close(65)
    endif

    !De-allocate
    deallocate(StrikeMap)
    deallocate(part%position)
    deallocate(part%velocity)

  endif mapeado
  !-----------------------------------------------------------------------------
  !=============================================================================
  ! ROUND 3: FIDASIM signal
  !=============================================================================
  !-----------------------------------------------------------------------------
  signal_part: if (signal) then
    if (verbose) then
      print*, '-----------------------------------------------------------------'
      print*,'Calculating signal'
    endif
    ! Fist load the markers simulation
    if (FILDSIMmode) then
      print*, "sorry not implemented"
    else
      call readFIDASIM4Markers(trim(FIDASIMfolder)//"/npa.bin", verbose)
    endif
    ! --- Open the file to save the data
    ! -- Strike points on the scintillator
    open(unit=61, file=trim(runFolder)//&
         '/results/'//trim(runID)//'.spsignal', access = 'stream', action='write')
    ! Save the header of the file
    write(61) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0, transfer(FILDSIMmode, 1), 18
    ! -- Strike points on the collimator
    open(unit=62, file=trim(runFolder)//&
         '/results/'//trim(runID)//'.spcsignal', access = 'stream', action='write')
    write(62) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0, transfer(FILDSIMmode, 1), 4
    ! -- File to save the orbits
    if (saveOrbits) then
      open(unit=63, file=trim(runFolder)//&
           '/results/'//trim(runID)//'.orbsignal', access = 'stream', action='write')
      write(63) versionID1, versionID2, runID
      if (saveOrbitLongMode) then
        write(63) 69
      else
        write(63) 6
      endif
    endif
    ! --  File to save the wrong markers
    if (save_wrong_markers_position) then
      open(unit=64, file=trim(runFolder)//&
           '/results/'//trim(runID)//'.wmsig', access = 'stream', action='write', status='replace')
      write(64) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0, transfer(FILDSIMmode, 1), 4
    endif
    ! --- Allocate the particle:
    call cpu_time(t_initial_orbits)
    call omega(M, Zout, BpinholeMod, OmegaPart)  ! Gyrofrequency
    dt = 2 * pi / OmegaPart / nGyro
    part%n_t = int(abs(maxT/dt))
    dt1 = 4.0E-7
    if (verbose) then
      print*, part%n_t, 'steps will be performed for each particle'
    endif
    allocate(part%position(3,part%n_t))
    allocate(part%velocity(3,part%n_t))
    ! --- Allocate the necesary matrix
    allocate(Strike(18,F4Markers%counter*nResampling))            ! Strike points in the scint
    allocate(CollimatorStrikes(4,F4Markers%counter*nResampling))       ! Strike position on the coll
    allocate(WrongMarkers(4,F4Markers%counter*nResampling))       ! Strike position on the coll
    CollimatorStrikes(:,:) = 0.0d0
    Strike(:,:) = 0.0d0
    WrongMarkers(:,:) = 0.0d0
    ! -- Initialise all the counters
    cScintillator = 0
    cCollimator = 0
    cWrongIons = 0
    cWrongNeutral = 0
    cWrong = 0
    cFoil = 0
    ! --- Get the scale
    ! See the normalization factor for the weight
    normalization_resample = max(dble(nResampling), 1.0d0)
    markers: do i=1, F4Markers%counter
      ! - Common data for all resampled markers
      part%energy0 = 0.5*sum(F4Markers%v(:, i)**2)*M/qe*amu_to_kg/1000.0
      resample: do j=1, nResampling
        ! Initialise the markers
        ! - Clean the marker data
        part%position(:, :) = 0.0d0
        part%velocity(:, :) = 0.0d0
        part%collision     = .False.
        part%kindOfCollision = 9  ! Just initialise it to a dummy value
        part%weight = F4Markers%wght(i) / normalization_resample
        part%cosalpha_foil = 0.0d0

        ! - FIDASIM data
        part%velocity(:, 1) = F4Markers%v(:, i)
        ! Resample the position (only valid for circular pinholes)
        pos1 = pinhole%r + pinhole%d1*(20.0d0*pinhole%u1 + 20.0d0*pinhole%u2)
        do while (sqrt(sum((pos1-pinhole%r)**2)) .gt.  pinhole%d1)
          call random_number(ran)
          pos1 = pinhole%r+ pinhole%d1*((-1+2*ran(2))*pinhole%u1+(-1+2*ran(3))*pinhole%u2)
        end do
        part%position(:, 1) = pos1

        ! - dt and charge
        part%dt    = dt1
        part%qm    = Zin * qe / M /amu_to_kg
        signal_tracking: do istep = 1, part%n_t-1
          call pushParticle(part%qm, part%position(:, istep), &
                            part%velocity(:, istep), &
                            part%position(:, istep + 1),&
                            part%velocity(:, istep + 1), part%dt)
          call checkCollision(part, istep)
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
              cycle resample

            elseif (part%kindOfCollision .eq. 2) then ! Scintillator collision
              call yieldScintillator(part, istep)
              cScintillator = cScintillator + 1 ! Update counter
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
              if (saveOrbits) then
                call random_number(rand_number)
                if (rand_number .lt. saveRatio) then
                  ! Save the collision point in the trajectory for future use
                  part%position(:, istep+1) = part%collision_point
                  cOrb = cOrb + 1
                  call saveOrbit(63, part, istep+1, saveOrbitLongMode)
                endif
              endif
              cycle resample
            elseif (part%kindOfCollision .eq. 1) then ! we collided with the foil
              part%collision = .False.
              part%dt = dt   ! Good dt to evolve an ion
              part%qm = Zout * qe / M /amu_to_kg ! New charge after the foil
              part%position(:, istep + 1) = part%collision_point
              call foilInteraction(part,foilNormal, istep )
              cFoil = cFoil + 1 ! Update the counter
            else
              print*, 'Warning we collide with element of kind', part%kindOfCollision
            endif
          endif
        enddo signal_tracking
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
        ! Save the final point of the marker
        cWrong = cWrong + 1
        WrongMarkers(1:3, cWrong) = part%position(:, istep)
        WrongMarkers(4, cWrong) = part%weight
      enddo resample
    end do markers
    write(61) cScintillator, transpose(Strike(:, 1:cScintillator))
    if (save_collimator_strike_points) then
      write(62) cCollimator, transpose(CollimatorStrikes(:, 1:cCollimator))
    endif
    if (save_wrong_markers_position) then
      write(64) cWrong, transpose(WrongMarkers(:, 1:cWrong))
    endif
    if (verbose) then
      print*, '-----'
      print*, 'Hitting Collimator', cCollimator
      print*, 'Hitting Scintillator', cScintillator
      print*, 'Not colliding Neutrals', cWrongNeutral
      print*, 'Not colliding Ions', cWrongIons
      print*, 'hitting foil', cFoil
    endif
    ! Write time and close the file
    write(61) F4Markers%time
    write(61) F4Markers%shot_number
    close(61)
    if (save_collimator_strike_points) then
      close(62)
    endif
    if (saveOrbits) then
      write(63) m
      write(63) cOrb  ! Write how many orbits we wrote in the file
      close(63)
    endif
    if (save_wrong_markers_position) then
      close(64)
    endif
    deallocate(part%position)
    deallocate(part%velocity)
    deallocate(Strike)
    deallocate(CollimatorStrikes)
  endif signal_part
end program sinpa
