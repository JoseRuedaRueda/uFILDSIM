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
  ! --- Note on units to open files
  !       - 60: Temporal files: Inputs/output. Open, read(write) and close it
  !       - 61: File with the strike points for the mapping (signal)
  !       - 62: Colimator data
  !       - 63: Orbit file

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
  allocate (rl(nGyroradius))
  allocate (alphas(nxi))
  open (unit=60, file=input_filename, form='formatted', iostat=ierr)
  read(60, NML=inputParams, iostat = ierr)
  close(60)

  ! Read the specifict namelist depending on the diagnostic
  if (FILDSIMmode.eqv..False.) then
    open (unit=60, file=input_filename, form='formatted', iostat=ierr)
    read(60, NML=nbi_namelist, iostat = ierr)
    nbi%u = u
    nbi%p0 = p0
    if (verbose) then
      print*, 'NBI namelist read.'
    endif
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
  call parseField(trim(SINPA_dir)//'/runs/'//trim(runID)//'/inputs/field.bin')
  ! Calculate the pinhole position in cylindrical coordinates
  rPinCyl(1) = sqrt(pinhole%r(1)**2 + pinhole%r(2)**2)
  rPinCyl(2) = atan2(pinhole%r(2), pinhole%r(1))
  rPinCyl(3) = pinhole%r(3)
  ! Get the magnetic field in the pinhole
  call getField(rPinCyl(1), rPinCyl(3), rPinCyl(2), Bpinhole)
  BpinholeMod = sqrt(Bpinhole(1)**2 + Bpinhole(2)**2 + Bpinhole(3)**2)
  if (verbose) then
    print*, '-----------------------------------------------------------------'
    print*, 'Magnetic field at the pinhole: ', Bpinhole
    print*, 'Mod: ', BpinholeMod
    print*, '-----------------------------------------------------------------'
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
    ! --- Open the files to save the data
    ! -- Strike points on the scintillator
    open(unit=61, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/StrikePoints.bin', access = 'stream', action='write')
    ! Save the header of the file
    write(61) versionID1, versionID2, runID, nGyroradius, rL, nxi, alphas, 18
    ! -- Strike points on the collimator
    open(unit=62, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/CollimatorStrikePoints.bin', access = 'stream', action='write')
    write(62) versionID1, versionID2, runID, nGyroradius, rL, nxi, alphas, 4
    ! -- File to save the orbits
    if (saveOrbits) then
      open(unit=63, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
           '/results/Orbits.bin', access = 'stream', action='write')
      write(63) versionID1, versionID2, runID
    endif

    ! --- Allocate the necesary matrices
    allocate(StrikeMap(11,nGyroradius * nxi))  ! Strike map
    allocate(CollimatorStrikes(4,nMap))       ! Strike position on the coll
    allocate(Strike(18,nMap))            ! Strike points in the scint

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
        dtNeutral = 4./v0 ! if we launched neutrals
      endif
      Lalphas: do ialpha = 1, nxi  ! Loop over the alpha angles
        ! -- Initialise all the counters
        cScintillator = 0
        cCollimator = 0
        cWrongIons = 0
        cWrongNeutral = 0
        cFoil = 0
        ! -- Clean the matrices with the complete information
        CollimatorStrikes(:,:) = 0.0d0
        Strike(:,:) = 0.0d0
        !
        Lmarkers: do  imc = 1, nMap
          call initMarker(imc, v0, dtneutral, alphas(ialpha))

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
                cycle Lmarkers
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
                ! See from where this markers was coming
                ! ToDo: Is the energy loss in the foil is included, change this
                ! ratio vv0/v0 to the real velocity
                call minimumDistanceLines(part%position(:, 1), part%velocity(:, 1)/v0, nbi%p0, nbi%u, dMin, posMin)
                ! Store the information of the marker
                Strike(1:3, cScintillator ) = part%collision_point ! f point
                Strike(4, cScintillator) = part%weight ! weight
                Strike(5:7,cScintillator) = posMin   ! Initial point (NBI)
                Strike(8:10,cScintillator) = part%velocity(:, 1)  ! Initial velocity
                Strike(11, cScintillator) = dMin ! Closer distance to NBI
                Strike(12, cScintillator) = part%beta ! Beta angle
                Strike(13:15, cScintillator) = part%position(:, 1)
                Strike(16:18, cScintillator ) = &
                  MATMUL(rotation, part%collision_point - ps)
                if (saveOrbits) then
                  call random_number(rand_number)
                  if (rand_number .lt. saveRatio) then
                    cOrb = cOrb + 1
                    write(63) istep
                    write(63) transpose(part%position(:, 1:istep))
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
          print*, 'Gyroradius:', rL(irl), ' Alpha:', alphas(ialpha)
          print*, 'Hitting Collimator', cCollimator
          print*, 'Hitting Scintillator', cScintillator
          print*, 'Wrong Neutrals', cWrongNeutral
          print*, 'Wrong Ions', cWrongIons
          print*, 'Foil', cFoil
        endif
        ! Average the position of the markers
        StrikeMap(1, ialpha + (irl - 1) * nxi) = rL(irl)
        StrikeMap(2, ialpha + (irl - 1) * nxi) = alphas(ialpha)
        StrikeMap(3, ialpha + (irl - 1) * nxi) = &   ! Strike point at Scint
          sum(Strike(:, 16)) / cScintillator
        StrikeMap(4, ialpha + (irl - 1) * nxi) = &
          sum(Strike(:, 17)) / cScintillator
        StrikeMap(5, ialpha + (irl - 1) * nxi) = &
          sum(Strike(:, 18)) / cScintillator
        StrikeMap(6, ialpha + (irl - 1) * nxi) = &  ! Birth position (NBI)
          sum(Strike(:, 5)) / cScintillator
        StrikeMap(7, ialpha + (irl - 1) * nxi) = &
          sum(Strike(:, 6)) / cScintillator
        StrikeMap(8, ialpha + (irl - 1) * nxi) = &
          sum(Strike(:, 7)) / cScintillator
        StrikeMap(9, ialpha + (irl - 1) * nxi) = &  ! Distance to NBI line
          sum(Strike(:, 11)) / cScintillator
        StrikeMap(10, ialpha + (irl - 1) * nxi) = & ! Collimator factor
          dble(cScintillator) / dble(nMap)            ! and striking ions
        StrikeMap(11, ialpha + (irl - 1) * nxi) = cScintillator
        ! De-allocate the variables
      enddo Lalphas
    enddo Lenergies
    ! Write the strike map
    print*,'Saving strike map in: '
    print*,trim(SINPA_dir)//'/runs/'//trim(runID)//'/results/StrikeMap.txt'
    open(unit=60, file=trim(SINPA_dir)//'runs/'//trim(runID)//&
         '/results/StrikeMap.txt', action='write', form='formatted')
    ! Write the header
    write(60,'(a,2x,a)') 'Simulation NAMELIST: ', input_filename
    write(60,'(a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a, 2x,a 2x,a)') 'Gyroradius (cm)', &
            'Alpha [rad]', 'X [cm]', 'Y [cm]', 'Z [cm]', &
            'X0 [cm]','Y0[cm]', 'Z0[cm]','d0[cm]','Collimator_Factor (%)', 'nMarkers'
    ! Write the data
    write(60,'(11f14.5)') StrikeMap
    ! Close all the openend files
    close(60)
    close(61)
    close(62)
    if (saveOrbits) then
      write(63) cOrb  ! Write how many orbits we wrote in the file
    endif
    close(63)
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
  ! signal_part: if (signal) then
  !   ! Fist load the FIDASIM simulation
  !   call readFIDASIM4Markers(trim(FIDASIMfolder)//"/npa.bin")
  !   ! --- Open the file to save the data
  !   ! -- Strike points on the scintillator
  !   open(unit=61, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
  !        '/results/SignalStrikePoints.bin', access = 'stream', action='write')
  !   ! Save the header of the file
  !   write(61) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0, 15
  !   ! -- Strike points on the collimator
  !   open(unit=62, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
  !        '/results/SignalCollimatorStrikePoints.bin', access = 'stream', action='write')
  !   write(62) versionID1, versionID2, runID, 1, 0.0d0, 1, 0.0d0,4
  !   ! -- File to save the orbits
  !   if (saveOrbits) then
  !     open(unit=63, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
  !          '/results/SignalOrbits.bin', access = 'stream', action='write')
  !     write(63) versionID1, versionID2, runID
  !   endif
  !   ! --- Allocate the particle:
  !   call omega(M, Zout, BpinholeMod, OmegaPart)  ! Gyrofrequency
  !   dt = 2 * pi / OmegaPart / nGyro
  !   part%n_t = int(maxT/dt)
  !   dtNeutral = 4.0E-7
  !   if (verbose) then
  !     print*, part%n_t, 'steps will be performed for each particle'
  !   endif
  !   allocate(part%position(3,part%n_t))
  !   allocate(part%velocity(3,part%n_t))
  !   ! --- Allocate the necesary matrix
  !   allocate(Strike(15,F4Markers%counter))            ! Strike points in the scint
  !   allocate(CollimatorStrikes(4,F4Markers%counter))       ! Strike position on the coll
  !   ! -- Initialise all the counters
  !   cScintillator = 0
  !   cCollimator = 0
  !   cWrongIons = 0
  !   cWrongNeutral = 0
  !   cEnter = 0
  !   cFoil = 0
  !   markers: do i=1, F4Markers%counter
  !     ! Initialise the markers
  !     part%weight         = F4Markers%wght(i)
  !     part%position(:, :) = 0.0d0
  !     part%velocity(:, :) = 0.0d0
  !     part%position(:,1) = F4Markers%fpos(:, i)
  !     part%velocity(:,1) = F4Markers%v(:, i)
  !     part%collision1     = .False.
  !     part%collision2     = .False.
  !     part%dt    = dtNeutral
  !     part%qm    = Zin *qe / M /amu_to_kg
  !     ! -- Collimator and carbon foil evolution
  !     ! The neutral dt is quite large, so do a maximum of 10 steps, if
  !     ! there is not collision for that, there is a problem
  !     signal_neutral: do istep = 1, 10
  !       ! Push the particle
  !       call pushParticle(part%qm, part%position(:, istep), &
  !                         part%velocity(:, istep), &
  !                         part%position(:, istep + 1),&
  !                         part%velocity(:, istep + 1), part%dt)
  !       ! Check collision with the collimator
  !       call triangleRay(collimator%triangleNum, collimator%triangles, &
  !                        part%position(:, istep), part%position(:, istep + 1), &
  !                        part%collision1, part%collision_point1)
  !       ! if collide with the collimator, there is not much to be done:
  !       if (part%collision1) then
  !           cCollimator = cCollimator + 1
  !           CollimatorStrikes(1:3, cCollimator) = part%collision_point1
  !           CollimatorStrikes(4, cCollimator) = part%weight
  !           if (saveOrbits) then
  !             call random_number(rand_number)
  !             if (rand_number .lt. saveRatio) then
  !               cOrb = cOrb + 1
  !               write(63) istep
  !               write(63) transpose(part%position(:, 1:istep))
  !             endif
  !           endif
  !           cycle markers
  !       endif
  !       ! Check collision with the foil
  !       call triangleRay(foil%triangleNum, foil%triangles, &
  !                        part%position(:, istep), part%position(:, istep + 1), &
  !                        part%collision1, part%collision_point1)
  !       ! if the particle collided, this loops need to end
  !       if (part%collision1) then
  !         ! Update particle variables
  !         part%dt = dt   ! Good dt to evolve an ion
  !         part%qm = Zout *qe / M /amu_to_kg ! New charge after the foil
  !         part%position(:, istep + 1) = part%collision_point1
  !         cFoil = cFoil + 1 ! Update the counter
  !         iistep = istep ! Save the step to eneter in the next loop
  !         if (saveOrbits) then
  !           call random_number(rand_number)
  !           if (rand_number .lt. saveRatio) then
  !             cOrb = cOrb + 1
  !             write(63) istep
  !             write(63) transpose(part%position(:, 1:istep))
  !           endif
  !         endif
  !         exit signal_neutral
  !       endif
  !     enddo signal_neutral
  !     ! If there is no collision, we got a problem
  !     if (.not. part%collision1) then
  !         ! ToDo Save fails
  !         cWrongNeutral = cWrongNeutral + 1
  !         cycle markers
  !     end if
  !     ! Now start the ion tracking part
  !     signal_charged: do istep = iistep + 1 , part%n_t - 1
  !       ! Push the particle
  !       call pushParticle(part%qm, part%position(:, istep), &
  !                         part%velocity(:, istep), &
  !                         part%position(:, istep + 1),&
  !                         part%velocity(:, istep + 1), part%dt)
  !       ! Check collision with the scintillator
  !       call triangleRay(scintillator%triangleNum, scintillator%triangles, &
  !                        part%position(:, istep), part%position(:, istep + 1), &
  !                        part%collision2, part%collision_point2)
  !       if (part%collision2) then
  !         cScintillator = cScintillator + 1 ! Update counter
  !         Strike(1, cScintillator ) = dble(i) ! FIDASIM marker id
  !         Strike(2:4, cScintillator ) = part%collision_point2 ! f point
  !         Strike(5:7, cScintillator ) = &
  !           MATMUL(rotation, part%collision_point2 - ps) !f point in scintillator region
  !         Strike(8:10, cScintillator) = F4Markers%ipos(:, i) ! CX position
  !         Strike(11:13,cScintillator) = part%velocity(:, 1)   ! Initial velocity
  !         Strike(14,cScintillator) = part%weight ! weight
  !         Strike(15, cScintillator) = F4Markers%kind(i) ! Kind of signal
  !
  !         if (saveOrbits) then
  !           call random_number(rand_number)
  !           if (rand_number .lt. saveRatio) then
  !             cOrb = cOrb + 1
  !             write(63) istep
  !             write(63) transpose(part%position(:, 1:istep))
  !           endif
  !         endif
  !         cycle markers
  !       endif
  !     enddo signal_charged
  !   end do markers
  !   write(61) cScintillator, transpose(Strike(:, 1:cScintillator))
  !   write(62) cCollimator, transpose(CollimatorStrikes(:, 1:cCollimator))
  !   if (verbose) then
  !     print*, '-----'
  !     print*, 'Hitting Collimator', cCollimator
  !     print*, 'Hitting Scintillator', cScintillator
  !     print*, 'Wrong Neutrals', cWrongNeutral
  !     print*, 'Wrong Ions', cWrongIons
  !   endif
  !   close(61)
  !   close(62)
  !   deallocate(part%position)
  !   deallocate(part%velocity)
  !   deallocate(Strike)
  !   deallocate(CollimatorStrikes)
  ! endif signal_part
end program sinpa
