MODULE genericParams
    !====================================================================
    !Define variables that are jointly used by the global optimization
    !as well as in the objective function.
    !Note, parameters that are specific objective function .
    !should not be here. For them use a separate module.
    !====================================================================
    USE nrtype
    USE stateControl
    IMPLICIT NONE

    INTEGER(I4B) :: p_default_alg = 0 !default algorithm.
          !  0-bobyqa (DFNLS), 1-amoeba (simplex), 2-dfpmin (BFGS Quasi-Newton method)
    INTEGER(I4B) :: p_nx, p_nmom
    INTEGER(I4B) :: p_iprint=2, p_maxeval, p_ninterppt ! BOBYQ variables
    INTEGER(I4B):: p_qr_ndraw, p_maxpoints,p_legitimate
    INTEGER(I4B) :: p_searchType, p_lotteryPoints
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_wspace
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_init
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_range, p_bound
    REAL(DP) :: p_fvalmax
    INTEGER(I4B), PARAMETER  :: p_nsmplx=100
        ! initialize the seed for random numbers
    INTEGER(I4B), PARAMETER :: p_SEED=314159265
    REAL(DP), parameter :: p_tolf_amoeba= 1.0d-4    ! Amoeba variables
    INTEGER(I4B), parameter :: p_itmax_amoeba=1000  ! Amoeba variables
    REAL(DP), PARAMETER :: gtol_DFPMIN = 1.0d-06			      ! DFPMIN PARAMETER
    INTEGER,  PARAMETER :: itmax_DFPMIN = 1000      ! DFPMIN PARAMETER
    REAL(DP), parameter :: neg_miss=-10000000.0_dp
    REAL(DP), parameter :: pos_miss= 10000000.0_dp

CONTAINS
    SUBROUTINE initialize(option,seqNo,config)
        !Initialize the instance. This includes parsing the configuration file,
        !setting parameter values, etc.
! 0ptions:
!       -1 = exit state - Stop all instances.
!        0 = cold start - The first invocation of the program, that will set up all the parallel helper files.
!                    Should always be the parameter when this is the first attempt at solving the problem.
!        1 = warm start - after one process is already running, all helper programs should be invoked using
!                    a warm start.
!        2 = update number of Sobol points - update the sobol point parameters over which to search
!                    as well as the local optimizations from these sobol points, but assume everything
!                    else in the config file has not been changed.
!        3 = update number of local minimizations - Increase the number of local minimizations but keep
!                    everything else same in the config file. This option uses the results from previously
!                    found local minimums.
!        4 = Diagnostic - Just evaluate the objective functions once for given initial guess in the config file with
!                    diagnostic option.
!        5 = Local Minimization - Run local minimization once for given initial guess in the config file.


        INTEGER(I4B), INTENT(in) :: option
        INTEGER(I4B), INTENT(out) :: seqNo
        CHARACTER(LEN=25), INTENT(IN), OPTIONAL :: config

        LOGICAL :: THERE_EXISTS,isUpdate
        INTEGER(I4B) :: i, n, fileDesc, openStat
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Old_sobolFnVal, Old_SearchResults
        INTEGER(I4B) :: numsolved,old_p_qr_ndraw,legitSobol,old_p_maxpoints
        CHARACTER(LEN=1000) :: errorString

        call random_seed(size = n)
        allocate(seed(n))
        seed(1)=123456
        call random_seed(put = seed)
        deallocate(seed)

        IF(OPTION==4 .or. OPTION==5) THEN
              !not updating, just parse the config file and run the objective value
              !for one set of parameters (option 4).
              !or run the given local minimization around the given point (option 5).
              isUpdate= .FALSE.
              call parseConfig(isUpdate, config)
              close(41)
              open(UNIT=41, FILE='seqNo.dat', STATUS='replace'); write(41,*) 0 ; close(41)
              IF(OPTION==5) THEN
                INQUIRE(FILE='searchResults.dat', EXIST=THERE_EXISTS)
                IF(THERE_EXISTS .eqv. .FALSE.) THEN
                  open(UNIT=41, FILE='searchResults.dat', STATUS='replace'); close(41)
                ENDIF

                INQUIRE(FILE='searchStart.dat', EXIST=THERE_EXISTS)
                IF(THERE_EXISTS .eqv. .FALSE.) THEN
                  open(UNIT=41, FILE='searchStart.dat', STATUS='replace'); close(41)
                ENDIF
              ENDIF
              RETURN
        END IF


        !if config is present,
        IF (PRESENT(config)) THEN

            IF (OPTION == 0) THEN
              ! if not updating, then reset all files
                  !delete any of the files managing the state
              isUpdate= .FALSE.
              close(41)
              open(UNIT=41, FILE='seqNo.dat', STATUS='replace'); write(41,*) 0 ; close(41)
              open(UNIT=41, FILE='lastSobol.dat', STATUS='replace');  write(41,*) 0 ; close(41)
              open(UNIT=41, FILE='legitSobol.dat', STATUS='replace');  write(41,*) 0 ; close(41)
              open(UNIT=41, FILE='lastParam.dat', STATUS='replace');  write(41,*) 0 ; close(41)
              open(UNIT=41, FILE='state.dat', STATUS='replace');  write(41,*) 0; close(41)
              open(UNIT=41, FILE='logFile.txt', STATUS='replace'); close(41)
              open(UNIT=41, FILE='internalConfig.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='sobol.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='sobolFnVal.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='missingSobol.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='x_starts.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='searchResults.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='searchStart.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='FinalResults.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='FinalStart.dat', STATUS='replace'); close(41)

            ELSEIF(OPTION==2) THEN
              !if updating Sobol points, then don't reset state, lastSobol, previously
              !calculated function values, or previously used parameters
              isUpdate= .TRUE.
              close(41)
              open(UNIT=41, FILE='missingSobol.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='seqNo.dat', STATUS='replace'); write(41,*) 0 ; close(41)
              open(UNIT=41, FILE='lastParam.dat', STATUS='replace');  write(41,*) 0 ; close(41)
              open(UNIT=41, FILE='x_starts.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='searchResults.dat', STATUS='replace');  close(41)
              open(UNIT=41, FILE='searchStart.dat', STATUS='replace');  close(41)
              open(UNIT=41, FILE='FinalResults.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='FinalStart.dat', STATUS='replace'); close(41)

              ! if we are updating, then store old values first
              call parseConfig(.FALSE., 'internalConfig.dat')
              old_p_qr_ndraw = p_qr_ndraw

            ELSEIF(OPTION==3) THEN
              !if updating number of local points, only reset  the following files
              isUpdate= .TRUE.
              close(41)
              open(UNIT=41, FILE='seqNo.dat', STATUS='replace'); write(41,*) 0 ; close(41)
              open(UNIT=41, FILE='x_starts.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='FinalResults.dat', STATUS='replace'); close(41)
              open(UNIT=41, FILE='FinalStart.dat', STATUS='replace'); close(41)

              ! if we are updating, then store old values first
              call parseConfig(.FALSE., 'internalConfig.dat')
              old_p_qr_ndraw = p_qr_ndraw
              old_p_maxpoints = p_maxpoints
            END IF

            call parseConfig(isUpdate, config)

            !print to internal config file so we have values for warm starts.
            call myopen(UNIT=fileDesc, FILE='internalConfig.dat', STATUS='replace', IOSTAT=openStat, ACTION='write')
            write(fileDesc,*) p_nx
            write(fileDesc,*) p_nmom
            write(fileDesc,*) p_maxeval
            write(fileDesc,*) p_qr_ndraw
            write(fileDesc,*) p_legitimate
            write(fileDesc,*) p_fvalmax
            write(fileDesc,*) p_maxpoints-1
            write(fileDesc,*) p_searchType
            write(fileDesc,*) p_lotteryPoints
            DO i=1,p_nx
                write(fileDesc,*) p_range(i,1),p_range(i,2)
            END DO
            DO i=1,p_nx
                write(fileDesc,*) p_init(i)
            END DO
            DO i=1,p_nx
                write(fileDesc,*) p_bound(i,1),p_bound(i,2)
            END DO
            call myclose(fileDesc)

            IF (OPTION==2) THEN  ! Update sobol points
                ! Then decide which sobol point to start with or
                ! to move on to local minimizations.
                INQUIRE(FILE='sobolFnVal.dat', EXIST=THERE_EXISTS)
                IF(THERE_EXISTS) THEN
                  allocate(Old_sobolFnVal(2*old_p_qr_ndraw,p_nx+2))
                  call myread2(old_sobolFnVal,'sobolFnVal.dat',numsolved)
                  numsolved=maxval(INT(old_sobolFnVal(1:numsolved,1)))
                  deallocate(Old_sobolFnVal)
                  legitSobol=getNumber('legitSobol.dat')
                  print*, "new # sobol= ",p_qr_ndraw
                  print*, "old # sobol= ",old_p_qr_ndraw
                  print*, "solved # sobol= ",numsolved
                  print*, "legitimate Sobol #= ",legitSobol
                  IF (numsolved >= p_qr_ndraw .or. legitSobol >= p_legitimate) THEN
                    write(errorString,*) "Update sobol: numsolved > p_qr_ndraw or legitSobol > p_legitimate ", &
                    NEW_LINE('A'), old_p_qr_ndraw, p_qr_ndraw, legitSobol , p_legitimate, NEW_LINE('A'), &
                    " Change p_qr_ndraw to ", max(numsolved+1,p_qr_ndraw)," or p_legitimate to ", max(legitSobol+1,p_legitimate)
                    call exitState(errorString);
                  ENDIF
                  CALL setState(numsolved, 'lastSobol.dat')
                ELSE
                  CALL setState(0, 'lastSobol.dat')
                ENDIF
            ENDIF

            IF (OPTION==3) THEN  ! Update number of local searches

                ! Decide  which sobol point to start the local minimization with.

                if(old_p_qr_ndraw .NE. p_qr_ndraw) THEN
                  print*, "new # sobol= ",p_qr_ndraw
                  print*, "old # sobol= ",old_p_qr_ndraw
                  write(errorString,*) "Update LocalMinimizations: old_p_qr_ndraw .NE. p_qr_ndraws", &
                  NEW_LINE('A'), "Run with option=2, update sobol points first "
                  call exitState(errorString);
                ENDIF

                INQUIRE(FILE='sobolFnVal.dat', EXIST=THERE_EXISTS)
                IF(THERE_EXISTS) THEN
                  allocate(Old_sobolFnVal(2*old_p_qr_ndraw,p_nx+2))
                  call myread2(old_sobolFnVal,'sobolFnVal.dat',numsolved)
                  numsolved=maxval(INT(old_sobolFnVal(1:numsolved,1)))
                  deallocate(Old_sobolFnVal)
                  print*, "solved # sobol= ",numsolved
                  print*, "new # local mimizations = ",p_maxpoints
                  IF (numsolved < p_maxpoints) THEN
                    write(errorString,*) "Update LocalMinimizations: numsolved < p_maxpoints", &
                    NEW_LINE('A'), "Run with option=2, update sobol points first!"
                    call exitState(errorString);
                  ENDIF
                ELSE
                  write(errorString,*) "Update LocalMinimizations:  sobolFnVal.dat cannot be found!"
                  call exitState(errorString);
                ENDIF

                INQUIRE(FILE='searchResults.dat', EXIST=THERE_EXISTS)
                IF(THERE_EXISTS) THEN
                  ALLOCATE(Old_SearchResults(2*old_p_maxpoints,p_nx+4))
                  call myread2(Old_SearchResults,'searchResults.dat',numsolved)
                  numsolved=maxval(INT(Old_SearchResults(1:numsolved,2)))
                  deallocate(Old_SearchResults)
                  print*, "new # local mimizations = ",p_maxpoints
                  print*, "old # local mimizations = ",old_p_maxpoints
                  print*, "solved # local mimizations= ",numsolved
                  IF (numsolved > p_maxpoints) THEN
                    write(errorString,*) "Update LocalMinimizations: numsolved > p_maxpoints", &
                    NEW_LINE('A'), old_p_maxpoints, p_maxpoints, &
                    " Change p_maxpoints to >", numsolved+1
                    call exitState(errorString);
                  ENDIF
                  CALL setState(numsolved, 'lastParam.dat')
                ELSE
                  CALL setState(0, 'lastParam.dat')
                ENDIF
            ENDIF

            if(p_qr_ndraw==0 .or. OPTION==3) THEN
              CALL setState(6, 'state.dat')
            ELSE
              CALL setState(1, 'state.dat')
            endif
        ELSE
            !Otherwise, read from internal config file
            call parseConfig(.FALSE., 'internalConfig.dat')
        END IF
        seqNo = getNextNumber('seqNo.dat')
    END SUBROUTINE initialize

    SUBROUTINE parseConfig(isUpdate, configFile)
        !parse the configuration file
        LOGICAL, INTENT(IN) :: isUpdate
        CHARACTER(LEN=*), INTENT(IN) :: configFile
        INTEGER(I4B) :: openStat, i
        INTEGER(I4B) :: nx, nmom, maxeval, ninterppt,legitimate
        REAL(DP) :: fvalmax
        INTEGER(I4B) :: qr_ndraw,maxpoints
        INTEGER(I4B) :: searchType, lotteryPoints, fileDesc
        CHARACTER(LEN=100) :: line
        LOGICAL :: parsedLine

        call myopen(UNIT=fileDesc, FILE=configFile, STATUS='old', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            print *,"<genericParams.parseConfig()> Error: unable to open file:",configFile
            stop 0
        END IF

        !Parse the basic parameters in the first line
        parsedLine = .FALSE.
        DO WHILE (parsedLine .eqv. .FALSE.)
            read(fileDesc,'(A100)',END=10) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.

            read(line,*) nx
            read(fileDesc,*,END=10) nmom
            read(fileDesc,*,END=10) maxeval
            read(fileDesc,*,END=10) qr_ndraw
            read(fileDesc,*,END=10) legitimate
            read(fileDesc,*,END=10) fvalMax
            read(fileDesc,*,END=10) maxpoints
            read(fileDesc,*,END=10) searchType
            read(fileDesc,*,END=10) lotteryPoints
        END DO

        IF (maxeval < 0) THEN
            maxeval = 40*(nx+1)
        END IF

        ninterppt = 2*nx+1

        IF (qr_ndraw < 0) THEN
            qr_ndraw = 500
        END IF

        IF (legitimate < 0) THEN ! minimum # of legitimate draws is not binding.
            legitimate = qr_ndraw + 10000
        END IF

        IF (fvalMax < 0.0_DP) THEN ! All objective values are legitimate, no matter how big.
            fvalMax = 100000000.0_DP
        END IF

        IF (qr_ndraw == 0) THEN
            maxpoints = 0

        ELSEIF (maxpoints < 0 ) THEN
            maxpoints=min(200,qr_ndraw)

        ELSE IF(maxpoints > qr_ndraw) THEN
            print *, "<genericParams:parseConfig> Error in config file. maxpoints > # sobol points"
            stop 0
        END IF

        IF (lotteryPoints < 0) THEN
            lotteryPoints = maxpoints
        END IF

        IF (isUpdate) THEN
            !if updating, don't allow changes in nx or nmom
            IF (nx /= p_nx) THEN
                print *, "<genericParams:parseConfig> Cannot update the number of parameters.",nx,p_nx
                stop 0
            END IF

            IF (nmom /= p_nmom) THEN
                print *, "<genericParams:parseConfig> Cannot update the number of moments."
                stop 0
            END IF
        ELSE
            IF (allocated(p_range) .eqv. .FALSE.) THEN
                allocate(p_range(nx,2))
            END IF
            IF (allocated(p_init) .eqv. .FALSE.) THEN
                allocate(p_init(nx))
            END IF
            IF (allocated(p_bound) .eqv. .FALSE.) THEN
                allocate(p_bound(nx,2))
            END IF
        END IF

        p_nx=nx
        p_nmom=nmom
        p_maxeval=maxeval
        p_qr_ndraw=qr_ndraw
        p_legitimate = legitimate
        p_fvalmax = fvalMax
        p_maxpoints = maxpoints+1 !we add one to account for the initial guess
        p_ninterppt = ninterppt
        p_searchType = searchType
        p_lotteryPoints = lotteryPoints

        ! now parse the range for the sobol points
        parsedLine = .FALSE.
        DO WHILE ( (parsedLine .eqv. .FALSE.))
            read(fileDesc,'(A100)',END=11) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.

            ! We are not allowed to change the parameter ranges except for warm start
            IF(.not. isUpdate) read(line,*) p_range(1,1), p_range(1,2)

            DO i=2,p_nx
                IF(.not. isUpdate) THEN
                read(fileDesc,*,END=11) p_range(i,1), p_range(i,2)
                ELSE
                  read(fileDesc,'(A100)',END=11) line
                ENDIF
            END DO
        END DO

        ! now parse the initial guess for the parameters
        ! initial guess can be updated during cold starts phase.
        parsedLine = .FALSE.
        DO WHILE ( (parsedLine .eqv. .FALSE.))
            read(fileDesc,'(A100)',END=12) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.
            ! We will always change the initial parameter values
            read(line,*) p_init(1)

            DO i=2,p_nx
                read(fileDesc,*,END=12) p_init(i)
            END DO
        END DO

        ! now parse the bounds for the parameters
        parsedLine = .FALSE.
        DO WHILE ( (parsedLine .eqv. .FALSE.))
            read(fileDesc,'(A100)',END=13) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.
            IF(.not. isUpdate)  read(line,*) p_bound(1,1), p_bound(1,2)

            DO i=2,p_nx
              IF(.not. isUpdate) THEN
                read(fileDesc,*,END=13) p_bound(i,1), p_bound(i,2)
              ELSE
                read(fileDesc,'(A100)',END=11) line
              ENDIF
            END DO
        END DO

        call myclose(fileDesc)

        IF(allocated(p_wspace)) deallocate(p_wspace)
        allocate(p_wspace((p_ninterppt+5)*(p_ninterppt+p_nx)+3*p_nx*(p_nx+5)/2))

        return

10      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read first line of parameters"
11      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read range of values for all parameters"
12      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read initial guesses for all parameters"
13      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read bounds for all parameters"
!100     format(i10)

    END SUBROUTINE parseConfig
END MODULE genericParams
