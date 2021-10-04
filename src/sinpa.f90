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
  !       - 61: File with the strike points for the mapping
  !       - 62: Colimator data
  !       - 63: Orbits reaching scintillator

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

  ! Open and read the configuration namelist
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
  allocate (alphas(nAlpha))
  open (unit=60, file=input_filename, form='formatted', iostat=ierr)
  read(60, NML=inputParams, iostat = ierr)
  close(60)

  open (unit=60, file=input_filename, form='formatted', iostat=ierr)
  read(60, NML=nbi_namelist, iostat = ierr)
  print*, nbi%u, nbi%p0
  close(60)


  open (unit=60, file=input_filename, form='formatted', iostat=ierr)
  read(60, NML=particlesFoil, iostat=ierr)
  close(60)

  if (verbose) then
    print*,'------------------------------------------------------------------'
    print*,'Input parameters read from: ',input_filename
    print*, 'runid:', runID
    print*,'------------------------------------------------------------------'
  endif
  ! --- Set the geometry
  call readGeometry(trim(geomID), scintillator, collimator, foil, .TRUE., &
                    u1, u2, u3, rHead, pinRadius)
  if (verbose) then
    print*,'------------------------------------------------------------------'
    print*,'Setting up the geometry of the plates and scintillator'
    print*,'Geometry setup finished'
    print*, ps, rotation
    print*,'------------------------------------------------------------------'
  endif

  ! --- Load the grid and magnetic field
  call parseField(trim(SINPA_dir)//'/runs/'//trim(runID)//'/inputs/field.bin')
  ! Calculate the pinhole position in cylindrical coordinates
  rHeadCyl(1) = sqrt(rHead(1)**2 + rHead(2)**2)
  rHeadCyl(2) = atan2(rHead(2), rHead(1))
  rHeadCyl(3) = rHead(3)
  ! Get the magnetic field in the pinhole
  call getField(rHeadCyl(1), rHeadCyl(3), rHeadCyl(2), Bpinhole)
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
    write(61) versionID1, versionID2, runID, nGyroradius, rL, nAlpha, alphas, 18
    ! -- Strike points on the collimator
    open(unit=62, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/CollimatorStrikePoints.bin', access = 'stream', action='write')
    write(62) versionID1, versionID2, runID, nGyroradius, rL, nAlpha, alphas, 4
    ! -- File to save the orbits
    if (saveOrbits) then
      open(unit=63, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
           '/results/obits.bin', access = 'stream', action='write')
      write(63) versionID1, versionID2, runID
    endif

    ! --- Allocate the necesary matrices
    allocate(StrikeMap(11,nGyroradius * nAlpha))  ! Strike map
    allocate(CollimatorStrikes(4,nMap))       ! Strike position on the coll
    allocate(MappingData(18,nMap))            ! Strike points in the scint

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
      Lalphas: do ialpha = 1, nAlpha  ! Loop over the alpha angles
        ! -- Initialise all the counters
        cScintillator = 0
        cCollimator = 0
        cWrongIons = 0
        cWrongNeutral = 0
        cEnter = 0
        cFoil = 0
        ! -- Clean the matrices with the complete information
        CollimatorStrikes(:,:) = 0.0
        MappingData(:,:) = 0.0

        Lmarkers: do  imc = 1, nMap
          ! -- Get the initial position of the marker
          ! We launch a marker in the pinhole area. For that we use a basic
          ! rejection method to distribute them uniformly along the pinhole
          call random_number(ran)
          r0 = rHead + pinRadius * ((-1+2*ran(2))*u1+(-1+2*ran(3))*u2)
          if (sqrt(sum((r0 - rHead)**2)).gt.pinRadius) then
            cycle Lmarkers
          end if
          ! If we have reached this point, we have entered, so update the counter
          cEnter = cEnter + 1
          ! -- Get the initial velocity of the marker
          ! Define the velocity vector of the particles
          beta = beta0*(-1 + 2*ran(1))
          vv0 = ((cos(alphas(ialpha))*u3 + sin(alphas(ialpha))*u1) * cos(beta) &
                 + sin(beta)*u2) * v0
          ! Initiallise the marker
          part%weight         = sin(beta)
          part%position(:, :) = 0.0
          part%velocity(:, :) = 0.0
          part%position(:,1) = r0
          part%velocity(:,1) = vv0
          part%collision1     = .False.
          part%collision2     = .False.
          part%dt    = dtNeutral
          part%qm    = Zin *qe / M /amu_to_kg
          part%alpha = alphas(ialpha)
          part%beta  = beta
          ! -- Collimator and carbon foil evolution
          ! The neutral dt is quite large, so do a maximum of 10 steps, if
          ! there is not collision for that, there is a problem
          neutral: do istep = 1, 10
            ! Push the particle
            call pushParticle(part%qm, part%position(:, istep), &
                              part%velocity(:, istep), &
                              part%position(:, istep + 1),&
                              part%velocity(:, istep + 1), part%dt)
            ! Check collision with the collimator
            call triangleRay(collimator%triangleNum, collimator%triangles, &
                             part%position(:, istep), part%position(:, istep + 1), &
                             part%collision1, part%collision_point1)
            ! if collide with the collimator, there is not much to be done:
            if (part%collision1) then
                cCollimator = cCollimator + 1
                CollimatorStrikes(1:3, cCollimator) = part%collision_point1
                CollimatorStrikes(4, cCollimator) = part%weight
                cycle Lmarkers
            endif
            ! Check collision with the foil
            call triangleRay(foil%triangleNum, foil%triangles, &
                             part%position(:, istep), part%position(:, istep + 1), &
                             part%collision1, part%collision_point1)
            ! if the particle collided, this loops need to end
            if (part%collision1) then
              ! Update particle variables
              part%dt = dt   ! Good dt to evolve an ion
              part%qm = Zout *qe / M /amu_to_kg ! New charge after the foil
              part%position(:, istep + 1) = part%collision_point1
              cFoil = cFoil + 1 ! Update the counter
              iistep = istep ! Save the step to eneter in the next loop
              exit neutral
            endif
          enddo neutral
          ! If there is no collision, we got a problem
          if (.not. part%collision1) then
              ! ToDo Save fails
              cWrongNeutral = cWrongNeutral + 1
              cycle Lmarkers
          end if
          ! Now start the ion tracking part
          charged: do istep = iistep + 1 , part%n_t - 1
            ! Push the particle
            call pushParticle(part%qm, part%position(:, istep), &
                              part%velocity(:, istep), &
                              part%position(:, istep + 1),&
                              part%velocity(:, istep + 1), part%dt)
            ! Check collision with the scintillator
            call triangleRay(scintillator%triangleNum, scintillator%triangles, &
                             part%position(:, istep), part%position(:, istep + 1), &
                             part%collision2, part%collision_point2)
            if (part%collision2) then
              cScintillator = cScintillator + 1 ! Update counter
              ! See from where this markers was coming
              ! ToDo: Is the energy loss in the foil is included, change this
              ! ratio vv0/v0 to the real velocity
              call minimumDistanceLines(r0, vv0/v0, nbi%p0, nbi%u, dMin, posMin)
              ! Store the information of the marker
              MappingData(1:3, cScintillator ) = part%collision_point2 ! f point
              MappingData(4:6,cScintillator) = posMin   ! Initial point (NBI)
              MappingData(7:9,cScintillator) = vv0  ! Initial velocity
              MappingData(10, cScintillator) = dMin ! Closer distance to NBI
              MappingData(11, cScintillator) = beta ! Beta angle
              MappingData(12:14, cScintillator) = part%collision_point1 ! Ionization point
              MappingData(15, cScintillator) = part%weight ! weight
              MappingData(16:18, cScintillator ) = &
                MATMUL(rotation, part%collision_point2 - ps)

              ! Save the orbits if we need it: Note, this way is cheating a bit
              ! because only orbits which collide with the scintillator are
              ! saved, so the orbits saved can't be used to check why some
              ! orbit has not collide... need to think if I should change this
              ! If you are reading this and you are not me and have an idea,
              ! please send it via email to jrrueda@us.es or ruejo@ipp.mpg.de
              if (saveOrbits) then
                call random_number(rand_number)
                if (rand_number .lt. saveRatio) then
                  cOrb = cOrb + 1
                  write(63) istep
                  write(63) transpose(part%position(:, 1:istep))
                endif
              endif
              cycle Lmarkers
            endif
          enddo charged
          ! If we reached this point and there is no collision, problem!
          if (.not. part%collision2) then
            cWrongIons = cWrongIons + 1
          endif
        enddo Lmarkers
        ! ! Rotate the data
        ! MappingData(16:18, 1:cScintillator ) = &
        !   MATMUL(rotation, transpose(MappingData(1:3, 1:cScintillator)) -ps)
        ! Write the markers in the file
        write(61) cScintillator, transpose(MappingData(:, 1:cScintillator))
        write(62) cCollimator, transpose(CollimatorStrikes(:, 1:cCollimator))
        ! Print some information
        if (verbose) then
          print*, '-----'
          print*, 'Gyroradius:', rL(irl), ' Alpha:', alphas(ialpha)
          print*, 'Hitting Collimator', cCollimator
          print*, 'Hitting Scintillator', cScintillator
          print*, 'Wrong Neutrals', cWrongNeutral
          print*, 'Wrong Ions', cWrongIons
        endif
        ! Average the position of the markers
        StrikeMap(1, ialpha + (irl - 1) * nAlpha) = rL(irl)
        StrikeMap(2, ialpha + (irl - 1) * nAlpha) = alphas(ialpha)
        StrikeMap(3, ialpha + (irl - 1) * nAlpha) = &   ! Strike point at Scint
          sum(MappingData(:, 1)) / cScintillator
        StrikeMap(4, ialpha + (irl - 1) * nAlpha) = &
          sum(MappingData(:, 2)) / cScintillator
        StrikeMap(5, ialpha + (irl - 1) * nAlpha) = &
          sum(MappingData(:, 3)) / cScintillator
        StrikeMap(6, ialpha + (irl - 1) * nAlpha) = &  ! Birth position (NBI)
          sum(MappingData(:, 4)) / cScintillator
        StrikeMap(7, ialpha + (irl - 1) * nAlpha) = &
          sum(MappingData(:, 5)) / cScintillator
        StrikeMap(8, ialpha + (irl - 1) * nAlpha) = &
          sum(MappingData(:, 6)) / cScintillator
        StrikeMap(9, ialpha + (irl - 1) * nAlpha) = &  ! Distance to NBI line
          sum(MappingData(:, 10)) / cScintillator
        StrikeMap(10, ialpha + (irl - 1) * nAlpha) = & ! Collimator factor
          dble(cScintillator) /dble(cEnter)            ! and striking ions
        StrikeMap(11, ialpha + (irl - 1) * nAlpha) = cScintillator
        ! De-allocate the variables
      enddo Lalphas
    enddo Lenergies
    ! Write the strike map
    print*,'Saving strike map in: '
    print*,trim(SINPA_dir)//'/runs/'//trim(runID)//'/results/StrikeMap.txt'
    open(unit=60, file=trim(SINPA_dir)//'/runs/'//trim(runID)//&
         '/results/StrikeMap.txt', action='write', form='formatted')
    ! Write the header
    write(60,'(a,2x,a)') 'Simulation NAMELIST: ', input_filename
    write(60,'(a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a,2x,a, 2x,a 2x,a)') 'Gyroradius (cm)', &
            'Alpha [rad]', 'X [cm]', 'Y [cm]', 'Z [cm]', &
            'X0 [cm]','Y0[cm]', 'Z0[cm]','d0[cm]','Collimator_Factor (%)', 'nMarkers'
    ! Write the data
    write(60,'(11f14.7)') StrikeMap
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
    deallocate(MappingData)
    deallocate(CollimatorStrikes)
  endif mapeado
end program sinpa
