PROGRAM GlobalSearch

    ! The main program for the TikTak global optimization algorithm. This algorithm
    !	evolved out of Fatih Guvenen's joint projects with Tony Smith, Serdar Ozkan,
    ! Fatih Karahan, Tatjana Kleineberg, and Antoine Arnaud. This version of the
    ! code was written by Arun Kandanchatha and Serdar Ozkan.
    !
    ! The code is written as a pseudo state machine, and
    ! multiple instances can be run at once. The main driver (which is
    ! the initially the cold start) will set states for all other
    ! instances of the program.
    !
    ! The states are as follows:
    !     -1: exit state - informs all processes to terminate
    !      1: initial state. Make sures that all variables, parameters, etc are Initialized.
    !           Prev States: None
    !           Next States: 2 (main driver), 3 (others)
    !      2: This state generates the sobol points.
    !         Only the "main driver" program can be in this state.
    !           Prev States: 1
    !           Next States: 3
    !      3: This state begins solving the objective function at the sobol points.
    !           Prev States: 1, 2
    !           Next States: 4 or 5
    !      4: This state finds if there are any missing sobol points.
    !         If not then this state sorts the sobol points by minimum objective function value.
    !         Only the "main driver" program can be in this state.
    !           Prev States: 3
    !           Next States: 5 or 6
    !      5: This state sorts the sobol points by minimum objective function value
    !         Only the "main driver" program can be in this state.
    !           Prev States: 4
    !           Next States: 6
    !      6: This state tries to minimize the objective function starting at
    !         the smallest calculated sobol value, and iterating through the
    !         number of points specified in the config file.
    !           Prev States: 4 or 5
    !           Next States: 7
    !      7: All instances: Check if any local minimization points have been missed.
    !         If so, run local minimization around them. After this, the program
    !         continues with taking the minimum objective value found, and restarting
    !         the minimization at this point one last time.
    !           Prev State: 5
    !           Next States: Exit

    USE nrtype
    USE genericParams
    USE stateControl
    USE utilities
    USE minimize
    USE OBJECTIVE, only : objFun, dfovec, obj_initialize, diagnostic

    IMPLICIT NONE

    ! command line options
    INTEGER(I4B) :: option      ! running the program with option=-1 (exit all instances ), 0 (cold start),
                                ! =1 (warm start), =2 (update sobol points), =3 (update local minimizations)
                                ! =4 (running diagnostics for given objective value parameters),
                                ! =5 (running local minimization around the given initial guess).

    LOGICAL :: isWarm             ! if this is a warm start or not
    LOGICAL :: updateSobolPoints  ! if this instance invocation is to update the sobol points
    LOGICAL :: updateLocalSearch  ! if this instance invocation is to update the number of local searches
    LOGICAL :: runDiagnostics     ! if this instance invocation is to compute the objective value
                                  ! once for given initial guess.
    LOGICAL :: runLocalMin        ! if this instance invocation is to run the given local minimization
                                      ! algorithm once around the given point.
    INTEGER(I4B) :: alg           ! the default algorithm this instance will use for
                                  ! minimizing the objective function

    ! temporary variables
    LOGICAL :: args_missing     ! if program arguments are missing
    LOGICAL :: complete         ! if the number of points to evaluate (either sobol or
                                ! minimization) has been completed
    CHARACTER(LEN=1000) :: errorString  ! variable used for writing error messages
    INTEGER(I4B) :: temp, temp2, i, driver, newDriver, currentState
    CHARACTER(LEN=25) :: arg1, config
    CHARACTER(LEN=1) :: arg2
    REAL(DP), DIMENSION(:), ALLOCATABLE  :: qrdraw
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sobol_trial
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: x_starts
    REAL(DP) :: itratio, rhobeg, rhoend,fval
    ! variables for file manipulation
    INTEGER(I4B) :: fileDesc, iostat, ierr,complete_status

    ! variables that will be used globally in the program
    INTEGER(I4B) :: LeadTerm ! 1 if if initial terminal when isWarm=0 or updateSobolPoints=1
    INTEGER(I4B) :: seqNo ! Unique number of the instance. We set it only once in the beginning.

    LeadTerm=0
    alg = p_default_alg

    call parseCommandLine(args_missing)
    if (args_missing) STOP

    IF (runDiagnostics .eqv. .TRUE.) THEN
      ! We are evaluating the objective value once for given parameter values in the config file.
      ! We run the simulation with diagnostic=1 so that it can produce moments other than the targets.
      diagnostic = 1
      print*,"Running diagnostics for the initial guess:"
      call initialize(option,seqNo,config)
      call obj_initialize
      write(*,'(A9,/,200(f12.6,/))') "p_init = ",(p_init(i), i=1,p_nx)
      fval=objFun(p_init)
      write(*,'(A18,200(f12.6,/))') "Objective value = ",fval
      stop
    ENDIF

    IF (runLocalMin .eqv. .TRUE.) THEN
      ! We are evaluating the objective value once for given parameter values in the config file.
      ! We run the simulation with diagnostic=1 so that it can produce moments other than the targets.
      print*,"Running local minimization around the initial guess:"
      call initialize(option,seqNo,config)
      call obj_initialize
      write(*,'(A9,/,200(f12.6,/))') "p_init = ",(p_init(i), i=1,p_nx)

      !save the initial point
      call myopen(UNIT=fileDesc, FILE='searchStart.dat', STATUS='unknown', IOSTAT=iostat, ACTION='write',position='append')
      write(fileDesc,'(2i10, 200f40.20)') -1, -1, p_init
      call myclose(fileDesc)


      SELECT CASE (alg)
      !minimize according to the algorithm
          CASE (0)
              !search
              itratio=0.6
              IF (itratio>0.5) THEN
                  rhobeg  = (minval(p_bound(:,2)-p_bound(:,1))/2.5_DP)/(4*itratio)
                  rhoend  = 1.0D-3/(4*itratio)
              ELSE
                  rhobeg = minval(p_bound(:,2)-p_bound(:,1))/2.50_DP
                  rhoend  = 1.0D-3
              END IF
              call bobyqa_h(p_nx,p_ninterppt,p_init,p_bound(:,1),p_bound(:,2),rhobeg,rhoend,p_iprint,p_maxeval,p_wspace,p_nmom)
          CASE (1)
              call runAmoeba(p_init,p_tolf_amoeba,p_itmax_amoeba)
          CASE (2)
                CALL EST_dfpmin(p_init,fval,objFun,1,itmax_DFPMIN,gtol_DFPMIN)
          CASE DEFAULT
              write(errorString,*) "<main:search> Error, unknown search method: ",alg
              call writeToLog(errorString)
              stop 0
      END SELECT
      fval = objFun(p_init)
      write(*,'(A9,/,200(f12.6,/))') "Solution = ",(p_init(i), i=1,p_nx)
      write(*,'(A18,200(f12.6,/))') "Objective value = ",fval

      !save the results
      call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=iostat, ACTION='write',position='append')
      write(fileDesc,'(3i10, 200f40.20)') -1, -1, -1, fval, p_init
      call myclose(fileDesc)

      stop
    ENDIF

    IF ((isWarm .eqv. .FALSE.) .or. (updateSobolPoints .eqv. .TRUE.)) THEN
        ! If this instance is true cold start or the one that updates sobol points,
        ! then it is the leader/initial
        LeadTerm=1
        CALL setState(0, 'state.dat')
        !If cold start or update, then initialize with config file (mandatory parameter)
        call initialize(option,seqNo,config)
        ! We set the proper state in initialize depending on whether warm/cold start or update
    ELSE
        ! we are doing a warm start. Let's make sure that all the cold state initialization
        ! is complete before trying to get a new sequence number.
        call waitState(1)
        call initialize(option,seqNo)
    END IF
    call obj_initialize

    ! now, let's go through the state machine
    currentState = 0
    DO WHILE(currentState .NE. p_exitState)

        ! Get the current state.
        currentState = getState()
        print*,"currentState",currentState

        ! perform the action as per the current state. The details of each state
        ! are given at the top of this file.
        SELECT CASE (currentState)
            CASE (p_exitState)
                cycle
            CASE (1)
              IF (LeadTerm==1) THEN
                  CALL setState(2, 'state.dat')
                  cycle
              ELSE
                  ! not driving. Let's wait until points are generated so we can start
                  ! evaluating the functions at those points.
                  call waitState(3)
                  cycle
              END IF
            CASE (2)
              IF (LeadTerm==1) THEN
                call setupSobol()
                CALL setState(3, 'state.dat')
                cycle
              ELSE
                !Otherwise, let's wait until we can start solving.
                call waitState(3)
                cycle
              END IF
            CASE (3)
              !Now everyone should be solving the function at the
              !sobol points.
              call solveAtSobolPoints(seqno,complete)
              IF (complete) THEN
                  ! If we have completed solving the objective function
                  ! at each sobol point, then check  if everything is
                  ! done.
                  IF(LeadTerm==1) THEN
                      call setState(4,'state.dat')
                      cycle
                  ELSE
                      call waitState(5)
                      cycle
                  END IF
              ELSE
                write(errorString,*) "EXIT STATE: COMPLETE=FALSE from solveAtSobolPoints."
                call exitState(trim(errorString))
              END IF
            CASE (4)
              ! Lead terminal finds the missing sobol points and write their identifiers to.
              ! missingSobol.dat. The other terminals wait until this is finished.
              ! If there are no missing sobol points lead terminal prepares for local minimization.
              IF(LeadTerm==1) THEN
                  call setupMissingSobol(complete_status)
                  IF (complete_status==1) THEN
                      ! Lead terminal finished finding missing sobol points.
                      ! There are missing sobol points so we run solveMissingSobol
                      call setState(5, 'state.dat')
                      cycle
                  ELSEIF(complete_status==2) THEN
                    ! Lead terminal did not find missing sobol points.
                    ! We skip to chooseSobol
                    write(errorString, *) " There is no missing sobol points or  &
                                there are enough legitimate sobols."
                    call writeToLog(errorString); print*, trim(errorString)
                    call setState(6, 'state.dat')
                    cycle
                  ELSE
                    write(errorString,*) "EXIT STATE: COMPLETE=FALSE from setupMissingSobol."
                    call exitState(trim(errorString))
                  END IF
              ELSE
                  call waitState(5)
                  cycle
              END IF
            CASE (5)
              ! We compute missing sobol points.
              call solveMissingSobol(seqno,complete)
              IF (complete) THEN
                  ! we have finished computing all sobol points (including the missing).
                  ! Lead terminal sort the sobol points - the algorithm will use the best n points
                  ! as specified in the configuration file.
                  ! Other instances wait until lead terminal is done with picking those points.
                  IF(LeadTerm==1) THEN
                    write(errorString, *) "<chooseSobol> Now we finished all sobol points (including the missing)."
                    call writeToLog(trim(errorString)); print*,trim(errorString)
                    call setState(6, 'state.dat')
                    cycle
                  ELSE
                    call waitState(6)
                    cycle
                  END IF
              ELSE
                write(errorString,*) "EXIT STATE: COMPLETE=FALSE from solveMissingSobol."
                call exitState(trim(errorString))
              END IF
            CASE (6)
              IF(LeadTerm==1) THEN
                write(errorString, *) "<chooseSobol> Now we choose the sobol points"
                call writeToLog(trim(errorString)); print*,trim(errorString)
                call chooseSobol(seqno)
                call setState(7, 'state.dat')
                cycle
              ELSE
                call waitState(7)
                cycle
              END IF
            CASE (7)
              call LocalMinimizations(seqno,alg,complete)
              IF (complete) THEN
                  ! If we have minimized at every point, look for any
                  ! missing points
                  IF(LeadTerm==1) THEN
                    call setState(8, 'state.dat')
                    cycle
                  ELSE
                    call waitState(8)
                    cycle
                  END IF
              ELSE
                write(errorString,*) "EXIT STATE: COMPLETE=FALSE from LocalMinimizations."
                call exitState(trim(errorString))
              END IF
            CASE (8)
                ! all instances search for the missing points.
                call SLEEP(10)
                call findMissingSearch(COMPLETE)
                IF (COMPLETE) THEN
                  ! There are no missing local minimizations
                  ! perform the last optimization
                  call lastSearch(seqNo)
                  EXIT
                ELSE
                  write(errorString,*) "EXIT STATE: COMPLETE=FALSE from findMissingSearch."
                  call exitState(trim(errorString))
                END IF
            CASE DEFAULT
              write(errorString, *) "<main> : Error, unknown state: ", currentState
              call writeToLog(trim(errorString)); print*,trim(errorString)
              stop 0
        END SELECT
        currentState = getState()
    END DO
    write(errorString, *) seqNo," has no more tasks."
    call writeToLog(trim(errorString)); print*, trim(errorString)
contains
    SUBROUTINE solveAtSobolPoints(seqNo, complete)
        !This routine solves the given function at the sobol point
        !indicated in lastSobol.dat If we have solved all points,
        !then complete is returned as true.
        implicit none
        INTEGER(I4B), INTENT(IN) :: seqNo
        LOGICAL, INTENT(OUT) :: complete
        INTEGER(I4B) :: whichPoint, openStat,numrows,legitSobol
        REAL(DP) :: fval

        IF(.not. allocated(sobol_trial)) THEN
          allocate(sobol_trial(p_qr_ndraw,p_nx))
          call myread2(sobol_trial,'sobol.dat',numrows)
          IF (numrows .NE. p_qr_ndraw) THEN
              ! There is an error
              write(errorString,*) "<genericSearch.solveAtSobolPoints()> Error: ", numrows,"  is less than",p_qr_ndraw
              print*, trim(errorString); call writeToLog(trim(errorString))
          END IF
        ENDIF

        LeadTerm = 0 ! while sobol points are executed no instance is leader.

        whichPoint = getNextNumber('lastSobol.dat');  print*, "sobol point", whichPoint
        legitSobol = getNumber('legitSobol.dat');  print*, "legitSobol=", legitSobol

        do while (whichPoint <= p_qr_ndraw .and. legitSobol<p_legitimate)
          IF (whichPoint < 0) THEN
              ! If the number from lastSobol.dat is negative an exit instance is run!
              complete = .FALSE.
              RETURN
          ELSE
              write(errorString, 7001) "Sequence# ", seqNo," is solving sobol point ",whichPoint
              call writeToLog(errorString)
              fval=objFun(sobol_trial(whichPoint,:))
              if(fval<p_fvalmax) then
                legitSobol = getNextNumber('legitSobol.dat')
              ELSE
                legitSobol = getNumber('legitSobol.dat')
                IF(legitSobol==p_legitimate) legitSobol = p_legitimate + 1
              ENDIF
              if(mod(legitSobol,10)==1) THEN
                write(errorString, 7002) "Sequence# ", seqNo," has found ", &
                legitSobol , "legitimate sobol points"
                call writeToLog(errorString)
              ENDIF
              call myopen(unit=fileDesc, file='sobolFnVal.dat', STATUS='old', &
              IOSTAT=openStat, ACTION='write', position='append')
              write(fileDesc,7000) whichPoint, fval, sobol_trial(whichPoint,:)
              call myclose(fileDesc)
          END IF

          IF (whichPoint == p_qr_ndraw .or. legitSobol==p_legitimate) THEN
            ! We make the instance that executes the last sobol point to be leader.
            LeadTerm=1
            write(errorString, *) seqNo," has become the leader in solveAtSobolPoints"
            call writeToLog(trim(errorString)); print*,trim(errorString)
          ENDIF

          whichPoint = getNextNumber('lastSobol.dat');
          if(mod(whichPoint,10)==1) print*, "sobol point",whichPoint
        END DO
        complete = .TRUE.

7000   format(i8, 200f40.20)
7001   format(A11,i5,A25,i6)
7002   format(A11,i5,A25,i6,A25)

    END SUBROUTINE solveAtSobolPoints

    SUBROUTINE LocalMinimizations(seqNo, algor, complete)
        !This routine searches for a minimum at the next
        !point, as given by lastParam.dat, using the algorithm
        !specified in algor. If we have gone through all the
        !parameter guesses, complete is set to TRUE. Otherwise,
        !it is set to .FALSE.
        implicit none
        INTEGER(I4B), INTENT(IN) :: seqNo, algor
        LOGICAL, INTENT(OUT) :: complete
        INTEGER(I4B) :: i,k, whichPoint, lotPoint
        REAL(DP), DIMENSION(p_nx) :: evalParam

        IF(LeadTerm==0) THEN
          ! We read the starting points once for the efficiency.
          ! We already have these starting points for the initial terminal.
          allocate(x_starts(p_maxpoints,p_nx+1))
          call myread2(x_starts,'x_starts.dat')
        ENDIF
        LeadTerm = 0 ! while sobol points are executed no instance is leader.

        !get which point we want to evaluate
        whichPoint = getNextNumber('lastParam.dat')
        print*,"whichPoint in LocalMinimizations=", whichPoint
        IF (whichPoint == 1) THEN
            i=0
            call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=ioStat,&
                ACTION='write',position='append')
            write(fileDesc,74521) seqNo, i,i, x_starts(1,:)
            call myclose(fileDesc)

            call myopen(UNIT=fileDesc, FILE='searchStart.dat', STATUS='unknown', IOSTAT=ioStat, &
                ACTION='write',position='append')
            write(fileDesc,7452) seqNo, i, x_starts(1,2:)
            call myclose(fileDesc)
        END IF

        do while (whichPoint <= p_maxpoints)
          IF (whichPoint == p_maxpoints) THEN
            ! We make the instance that executes the last sobol point to be leader.
            LeadTerm=1
            write(errorString, *) seqNo," has become the leader in search"
            call writeToLog(errorString); print*,trim(errorString)
          ENDIF

          IF (whichPoint < 0) THEN
              ! If the number from lastSobol.dat is negative an exit instance is run!
              complete = .FALSE.
              RETURN
          ELSE
            ! get parameter value based on previous solutions. Note there are multiple search types,
            !   0: the best point
            !   1: a lottery using n points (where n is given in the config file)
            evalParam = getModifiedParam(whichPoint, x_starts(whichPoint,2:), p_searchType, lotPoint)
            write(errorString, 5452) "Sequence # ", seqNo," is searching using sobol point ", &
                whichPoint," and best point ",lotPoint
            call writeToLog(errorString);  print*, trim(errorString)

            !We now have the point at which to solve. So, finish the search
            call completeSearch(seqNo, algor, whichPoint, evalParam)
          END IF
          whichPoint = getNextNumber('lastParam.dat')
          print*,"whichPoint in LocalMinimizations=", whichPoint
        END DO
        complete = .TRUE.

        deallocate(x_starts)
        RETURN
5452    format(A11, i3, A32,i4,A16,i4)
74521    format(3i10, 200f40.20)
7452    format(2i10, 200f40.20)

    END SUBROUTINE LocalMinimizations

    SUBROUTINE completeSearch (seqNo, algor, whichPoint, evalParam)
        !This subroutine completes the search given a specific point at which to evaluate.
        INTEGER(I4B), INTENT(IN) :: seqNo, algor, whichPoint
        INTEGER(I4B) :: openStat,nwrite=1
        REAL(DP), DIMENSION(p_nx), INTENT(INOUT) :: evalParam
        REAL(DP) :: fn_val
        REAL(DP), DIMENSION(p_nx+4) :: currentBest

        !save the initial point
        call myopen(UNIT=fileDesc, FILE='searchStart.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,7453) seqNo, whichPoint, evalParam
        call myclose(fileDesc)

        SELECT CASE (algor)
        !minimize according to the algorithm
            CASE (0)
                !search
                itratio=REAL(whichPoint)/REAL(p_maxpoints)
                IF (itratio>0.5) THEN
                    rhobeg  = (minval(p_bound(:,2)-p_bound(:,1))/2.5_DP)/(4*itratio)
                    rhoend  = 1.0D-3/(4*itratio)
                ELSE
                    rhobeg = minval(p_bound(:,2)-p_bound(:,1))/2.50_DP
                    rhoend  = 1.0D-3
                END IF
                call bobyqa_h(p_nx,p_ninterppt,evalParam,p_bound(:,1),p_bound(:,2),rhobeg,rhoend,p_iprint,p_maxeval,p_wspace,p_nmom)
            CASE (1)
                call runAmoeba(evalParam,p_tolf_amoeba,p_itmax_amoeba)
            CASE (2)
                  CALL EST_dfpmin(evalParam,fn_val,objFun,nwrite,itmax_DFPMIN,gtol_DFPMIN)
            CASE DEFAULT
                write(errorString,*) "<main:search> Error, unknown search method: ",algor
                call writeToLog(errorString)
                stop 0
        END SELECT
        fn_val = objFun(evalParam)

        !output result
        currentBest = getBestLine(openStat)
        i = INT(currentBest(3))
        IF (fn_val < currentBest(4)) THEN
            i = i+1
        END IF

        !save the results
        call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,74531) seqNo, whichPoint, i, fn_val, evalParam
        call myclose(fileDesc)

74531    format(3i10, 200f40.20)
7453    format(2i10, 200f40.20)

    END SUBROUTINE completeSearch

    SUBROUTINE runAmoeba(amoeba_pt,tolf,itmax)
        !This subroutine executes the amoeba search. It takes the point passed in
        !as a parameter, and generates a simplex centered on it. The points of the
        !simplex are proportional to the number of sobol points generated and the
        !dimension of the parameter space. For example, if we have 2 parameters and
        !100 sobol points, then the simplex distance is approximately 1/10 of the
        !range of each parameter.
        use simplex, only : simplex_coordinates2
        IMPLICIT NONE
        REAL(DP), intent(in) :: tolf
        INTEGER(I4B), intent(in) :: itmax
        INTEGER(I4B) :: ia
        REAL(DP), DIMENSION(p_nx), intent(inout) :: amoeba_pt
        REAL(DP), DIMENSION(p_nx,p_nx+1) ::  x
        REAL(DP), DIMENSION(p_nx+1,p_nx) :: xT
        REAL(DP), DIMENSION(p_nx+1) :: fnVals
        REAL(DP), DIMENSION(p_nx) :: temp

        call simplex_coordinates2 ( p_nx, x )

        do ia=1,p_nx+1
            temp = amoeba_pt+x(1:p_nx,ia)*(p_range(:,2)-p_range(:,1))/(DBLE(p_maxpoints)**(1.0_dp/p_nx))
            temp=max(min(temp,p_bound(:,2)),p_bound(:,1))
            fnVals(ia)=objFun(temp)
            xT(ia,:)=temp
        end do

        ia=itmax
        call amoeba(xT,fnVals,tolf,objFun,ia)
        ia=minloc(fnVals,1)
        amoeba_pt = xT(ia,:)
    END SUBROUTINE runAmoeba

    SUBROUTINE lastSearch(seqNo)
        !This routine solves the given function for the last time,
        !taking as the initial point the best point we have found
        !so far. It always uses the same algorithm, bobyqa.
        implicit none
        INTEGER(I4B), INTENT(IN) :: seqNo
        INTEGER(I4B) :: i, openStat,whichPoint,nwrite=1
        REAL(DP) :: evalParam(p_nx),temp(p_nx+4)
        REAL(DP) :: rhobeg, rhoend, fn_val

        temp=getBestLine(whichPoint)
        evalParam=temp(5:)
        fn_val=temp(4)

        i=-1
        call myopen(UNIT=fileDesc, FILE='FinalStart.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,7451) seqNo, whichPoint, evalParam
        call myclose(fileDesc)

        write(errorString, 7452) " "
        call writeToLog(errorString); print*,trim(errorString)
        write(errorString, 7452) "The best final point is ", temp
        call writeToLog(errorString); print*,trim(errorString)
        ! write(errorString, 7453) "Instance# ", seqno, " is executing the final search at with max times: ",p_maxeval
        write(errorString, 7453) "Instance# ", seqno, " is executing the final search at with max times: ",itmax_DFPMIN
        call writeToLog(errorString);  print*,trim(errorString)

        CALL EST_dfpmin(evalParam,fn_val,objFun,nwrite,itmax_DFPMIN,gtol_DFPMIN)

        ! rhobeg  = minval(p_bound(:,2)-p_bound(:,1))/10.0_DP
        ! rhoend  = 1.0D-3/4.0_dp
        ! call bobyqa_h(p_nx,p_ninterppt,evalParam,p_bound(:,1),p_bound(:,2),rhobeg,rhoend,p_iprint,p_maxeval,p_wspace,p_nmom)

        ! fn_val = objFun(evalParam)

        call myopen(UNIT=fileDesc, FILE='FinalResults.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,7451) seqNo, whichPoint, fn_val, evalParam
        call myclose(fileDesc)

        return

7451    format(2i10, 200f40.20)
7452    format(A25, 200f12.6)
7453    format(A10,i4,A54,i5)

    END SUBROUTINE lastSearch

    FUNCTION getModifiedParam(whichPoint, evalPoint, whichMethod, lotPoint) RESULT (y)
        !This routine returns the next parameter at which to evaluate. It combines
        !the next sobol point with a previously evaluated point specified by whichMethod.
        !Method 1: simply get the best point, and scale to that
        !Method 2: each point has a probability of being selected, based
        !          on rank. Pick one and then scale to that
        INTEGER(I4B), INTENT(IN) :: whichPoint
        REAL(DP), DIMENSION(p_nx), INTENT(IN) :: evalPoint
        INTEGER(I4B), INTENT(IN) :: whichMethod
        INTEGER(I4B), INTENT(OUT) :: lotPoint
        REAL(DP), DIMENSION(p_nx) :: y
        REAL(DP), DIMENSION(p_nx) :: basePoint
        REAL(DP) :: w11
        REAL(DP), DIMENSION(2*p_maxpoints,p_nx+4) :: sortedPoints
        INTEGER(I4B) :: count

        !sort the points and find out how many we have
        call getSortedResults(sortedPoints,count)
        IF(count<=1) THEN ! if we haven't finished any local minization yet just use evalPoint
          y= evalPoint
          lotPoint=0
          return
        ENDIF

        SELECT CASE (whichMethod)
            CASE (0)
                basePoint = getBestPoint(lotPoint)
            CASE (1)
                basePoint = getLotteryPoint(lotPoint)
            CASE DEFAULT
                write(errorString, *) "<genericSearch:getModifiedParam> : Error, unknown selectionMethod: ", whichMethod
                call writeToLog(errorString); print*,trim(errorString)
                stop 0
        END SELECT

        w11=min(max(0.10, (sqrt(real(whichPoint)/real(p_maxpoints)))),0.95)
        y= w11*basePoint+(1.0_DP-w11)*evalPoint
    END FUNCTION getModifiedParam

    FUNCTION getLotteryPoint(lotPointSelect) RESULT (y)
        !This routine returns a parameter where the probability of that parameter is
        !proportional to the inverse of its rank amongst all sorted points.
        INTEGER(I4B), INTENT(OUT) :: lotPointSelect
        REAL(DP), DIMENSION(p_nx) :: y
        REAL(DP), DIMENSION(2*p_maxpoints,p_nx+4) :: sortedPoints
        INTEGER(I4B) :: numPoints, sumVal, i,nonmiss
        REAL(DP), DIMENSION(p_maxpoints) :: probOfPoint
        REAL(DP) :: sum2

        !sort the points and find out how many we have
        call getSortedResults(sortedPoints,numPoints)
        IF(numPoints .eq. 0) THEN
            y = p_init
            lotPointSelect=0
            return
        END IF
        !but, only use the number of points we are allowed to as specified
        !in the configuration file
        numPoints = min(numPoints, p_lotteryPoints)

        sumVal = (numPoints * (numPoints+1)) / 2

        !calculate the probability of each point
        probOfPoint = 0.0_dp
        DO i=1,numPoints
            probOfPoint(i) = DBLE(sumVal)/DBLE(i)
        END DO
        sum2 = sum(probOfPoint)
        probOfPoint = probOfPoint/sum2

        DO i=2,numPoints
            probOfPoint(i) = probOfPoint(i-1)+probOfPoint(i)
        END DO

        !now, let's get a random number and find the appropriate point
        call random_number(sum2)
        DO i=1,numPoints
            IF (sum2 < probOfPoint(i)) THEN
                exit
            END IF
        END DO
        y=sortedPoints(i,4:)
        lotPointSelect=INT(sortedPoints(i,2))
    END FUNCTION getLotteryPoint

    FUNCTION getBestPoint(lotPoint) RESULT (y)
        !This routine returns the best set of parameters we have solved so far, as
        !stored in searchSaved.dat
        INTEGER(I4B), INTENT(OUT) :: lotPoint
        REAL(DP), DIMENSION(p_nx) :: y
        REAL(DP), DIMENSION(p_nx+4) :: y2

        y2=getBestLine(lotPoint)
        y=y2(5:)
        print*,"y", y

    END FUNCTION getBestPoint

    FUNCTION getBestLine(lotPoint) RESULT(y)
        !This function returns the best line that we have calculateod
        !so far (in terms of the minimum function value)
        INTEGER(I4B), INTENT(OUT) :: lotPoint
        REAL(DP), DIMENSION(p_nx+4) :: y
        INTEGER(I4B) :: count
        REAL(DP), DIMENSION(2*p_maxpoints,p_nx+4) :: sortedPoints

        !sort the points and find out how many we have
        call getSortedResults(sortedPoints,count)
        IF(count==0) THEN
          y(1:3) = 0.0_DP
          y(4) = pos_miss
          y(5:) = p_init
          lotPoint = 0
          return
        ELSE
          y=sortedPoints(1,:)
          lotPoint=INT(sortedPoints(1,2))
          return
        END IF
    END FUNCTION getBestLine

    SUBROUTINE getSortedResults(y,nonmiss)
        !This function sorts all the results so far, and returns them in an array.
        !We use this to determine the top points for the lottery selection.
        REAL(DP), DIMENSION(2*p_maxpoints,p_nx+4),INTENT(OUT) :: y
        INTEGER, INTENT(OUT) :: nonmiss
        INTEGER(I4B), DIMENSION(2*p_maxpoints) :: fval_init_index
        INTEGER(I4B) :: i
        REAL(DP), DIMENSION(2*p_maxpoints,p_nx+4) :: SearchRes

        y = 10.0_dp*pos_miss

        call myread2(SearchRes,'searchResults.dat', nonmiss)
        IF(nonmiss>0) THEN
          !Sort the function values. fval_init_index(1) gives the smallest
          !element in fval, fval_init_index(2) gives the next smallest, etc.
          CALL indexx(nonmiss,SearchRes(1:nonmiss,4),fval_init_index(1:nonmiss))
          y(1:nonmiss,:) = SearchRes(fval_init_index(1:nonmiss),:)
        ENDIF
        return
    END SUBROUTINE getSortedResults

    SUBROUTINE setupSobol()
      !If this process is able to set the state, then it must be the "driver"
      !So generate sobol points, etc.
      INTEGER(I4B) :: i

      allocate(qrdraw(p_nx), sobol_trial(p_qr_ndraw,p_nx))
      call insobl(p_nx, p_qr_ndraw)
      DO i = 1,p_qr_ndraw
          CALL I4_SOBOL(p_nx,qrdraw)
          sobol_trial(i,:) = p_range(:,1)+qrdraw*(p_range(:,2)-p_range(:,1))
      END DO
      call mywrite2(sobol_trial,'sobol.dat')
      deallocate(qrdraw)
    END SUBROUTINE setupSobol

    SUBROUTINE chooseSobol(seqno)
        !This routine returns the best set of parameters we have solved so far, as
        !stored in searchSaved.dat
        !This is only called by the leader terminal/instance
        INTEGER(I4B), INTENT(IN) :: seqno
        INTEGER(I4B) :: i,j,numrows,legitSobol
        LOGICAL :: allFinished
        INTEGER  :: num_pos_sob ! (Hopefully) Maximum number of rows in 'sobolFnVal.dat'
        INTEGER(I4B),ALLOCATABLE, DIMENSION(:) :: fval_init_index
        REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: sobolFnVal

        allocate(x_starts(p_maxpoints,p_nx+1))

        x_starts(1,1) = objFun(p_init)
        x_starts(1,2:p_nx+1) = p_init

        IF(p_qr_ndraw>0) THEN
          num_pos_sob=2*p_qr_ndraw ! (Hopefully) Maximum number of rows in 'sobolFnVal.dat'

          ALLOCATE(sobolFnVal(num_pos_sob,p_nx+2),fval_init_index(num_pos_sob))
          call myread2(sobolFnVal,'sobolFnVal.dat',numrows)

          IF(numrows>=num_pos_sob) THEN
            ! There are more sobol evaluations in 'sobolFnVal.dat' than assumed number num_pos_sob
            ! Increase the number num_pos_sob!
            write(errorString, *) "<main:chooseSobol> Error. More evaluations in 'sobolFnVal.dat", &
             numrows ,"than assumed num_pos_sob ",num_pos_sob, "Increase num_pos_sob and run update (option 2)!"

            call exitState(errorString)
          ENDIF

          legitSobol=getNumber('legitSobol.dat')
          if(numrows<p_qr_ndraw .and. legitSobol<p_legitimate) THEN
            ! Didn't read enough values, which shouldn't happen
            write(errorString, *) "<main:chooseSobol> Error. Could read ", numrows, &
            " sobol evaluations. Less than ",p_qr_ndraw, " and legitSobol<p_legitimate", &
            legitSobol, p_legitimate
            call exitState(errorString)
          endif

          !Sort the function values. fval_init_index(1) gives the smallest
          !element in fval, fval_init_index(2) gives the next smallest, etc.
          CALL indexx(numrows,sobolFnVal(1:numrows,2),fval_init_index(1:numrows))

          x_starts(2,:)=sobolFnVal(fval_init_index(1),2:)
          j=3
          DO i=2,numrows
          IF(INT(sobolFnVal(fval_init_index(i),1)) .NE. INT(sobolFnVal(fval_init_index(i-1),1))) THEN
             x_starts(j,:)=sobolFnVal(fval_init_index(i),2:)
             j=j+1
          ENDIF
          IF(j>p_maxpoints) EXIT
          END DO

          DEALLOCATE(sobolFnVal,fval_init_index)
        ENDIF

        call mywrite2(x_starts,'x_starts.dat')
        return
    END SUBROUTINE chooseSobol

    SUBROUTINE setupMissingSobol(complete_status)
        !Called by the main instance. Waits until function values have been derived
        !for all sobol points. If some points are missing (often due to a warm start
        !instance having been killed), find those points and prepare for .
        !This is only run by the leader terminal/instance.
        INTEGER, INTENT(OUT) :: complete_status
        INTEGER(I4B) :: openStat, i,missing,legitSobol,numrows,nummiss,numsobol
        LOGICAL, DIMENSION(p_qr_ndraw) :: solvedPoints
        CHARACTER(LEN=1000) :: errorString
        REAL(DP) :: fval,missingSobol(p_qr_ndraw,1)


        open(UNIT=41, FILE='legitSobolMiss.dat', STATUS='replace');  write(41,*) 0 ; close(41)
        complete_status=0

        !we aren't actually guaranteed that all points are complete. Just that we have tried everything. So
        !find out which ones are missing, and have all processes go back and solve them
        solvedPoints = .FALSE.
        legitSobol = 0
        numsobol = 0

        call myopen(UNIT=fileDesc, FILE='sobolFnVal.dat', STATUS='unknown', IOSTAT=openStat, ACTION='read')
        DO
            read(fileDesc,270, END=10) i, fval
            IF (i > p_qr_ndraw) THEN
                !This shouldn't happen. let's note the error and stop.
                write(errorString, *) seqNo, " found point: ",i,"greater than max: ",p_qr_ndraw
                call exitState(errorString)
            END IF
            if(fval<p_fvalmax .and. solvedPoints(i) == .FALSE.) legitSobol = legitSobol +1
            if(solvedPoints(i) == .FALSE.) numsobol = numsobol +1
            solvedPoints(i) = .TRUE.
        END DO
10      call myclose(fileDesc)

        missing=count(solvedPoints==.FALSE.)

        IF(missing>0 .and. legitSobol<p_legitimate) THEN
          !Add the sobol points that need to be solved.
          call myopen(unit=fileDesc, file='missingSobol.dat', STATUS='replace', IOSTAT=openStat, ACTION='write')
          DO i=1,p_qr_ndraw
            IF(solvedPoints(i)) THEN
                cycle
            END IF
            write(fileDesc,271) i
          END DO
          call myclose(fileDesc)
          call setState(legitSobol,'legitSobolMiss.dat')
          if(p_legitimate < p_qr_ndraw) THEN
            call myread2(missingSobol,'missingSobol.dat',numrows)
            nummiss=min(numrows,INT(4.0_DP*numsobol/legitSobol)*(p_legitimate-legitSobol))
            call mywrite2(missingSobol(1:nummiss,:),'missingSobol.dat')

            write(errorString,*) "<setupMissingSobol> There are ", nummiss, " missing sobols to be solved."
            print*, trim(errorString); call writeToLog(trim(errorString))
          endif
          complete_status=1
        ELSE
          call myopen(unit=fileDesc, file='missingSobol.dat', STATUS='replace', IOSTAT=openStat, ACTION='write')
          call myclose(fileDesc)
          complete_status=2
        ENDIF
7000    format(i8, 200f40.20)
270     format(i10,f40.20)
271     format(i10)

    END SUBROUTINE setupMissingSobol

    SUBROUTINE solveMissingSobol(seqno,allFinished)
        !Called by the main instance. Waits until function values have been derived
        !for all sobol points. If some points are missing (often due to a warm start
        !instance having been killed), solve for it.
        !This is only run by the leader terminal/instance.
        LOGICAL, INTENT(OUT) :: allFinished
        INTEGER(I4B), INTENT(IN) :: seqno
        INTEGER(I4B) :: openStat, i,missing,numrows,legitSobol
        CHARACTER(LEN=1000) :: errorString
        REAL(DP) :: fval
        REAL(DP) :: missingSobol(p_qr_ndraw,1)

        LeadTerm=0
        allFinished = .TRUE.
        !we aren't actually guaranteed that all points are complete. Just that we have tried everything. So
        !find out which ones are missing, and have all processes go back and solve them

        IF (.not. allocated(sobol_trial)) THEN
          allocate(sobol_trial(p_qr_ndraw,p_nx))
          call myread2(sobol_trial,'sobol.dat',numrows)
          IF (numrows .NE. p_qr_ndraw) THEN
              ! There is an error
              write(errorString,*) "<solveMissingSobol()> Error: ", numrows,"  is less than",p_qr_ndraw
              print*, trim(errorString); call writeToLog(trim(errorString))
              STOP
          END IF
        ENDIF

        legitSobol=getNumber('legitSobolMiss.dat')

        DO
          call myread2(missingSobol,'missingSobol.dat',missing)
          print*, "missing", missing
          if(missing>0) THEN
            if(missing>1) THEN
              call mywrite2(missingSobol(2:missing,:),'missingSobol.dat')
            ELSE
              call myopen(unit=fileDesc, file='missingSobol.dat', STATUS='replace', IOSTAT=openStat, ACTION='write')
              call myclose(fileDesc)
            ENDIF
            write(errorString, *) seqNo," solving missing sobol point ",INT(missingSobol(1,1))
            call writeToLog(errorString); print*, trim(errorString);
            fval=objFun(sobol_trial(INT(missingSobol(1,1)),:))

            IF(fval<p_fvalmax) THEN
              legitSobol = getNextNumber('legitSobolMiss.dat')
            ELSE
              legitSobol = getNumber('legitSobolMiss.dat')
              IF(legitSobol==p_legitimate) legitSobol = p_legitimate + 1
            ENDIF

            call myopen(unit=fileDesc, file='sobolFnVal.dat', STATUS='old', IOSTAT=openStat, ACTION='write', position='append')
            write(fileDesc,7000) INT(missingSobol(1,1)), fval, sobol_trial(INT(missingSobol(1,1)),:)
            call myclose(fileDesc)
            if(missing==1 .or. legitSobol==p_legitimate) THEN
              LeadTerm=1
              allFinished = .TRUE.
              EXIT
            ENDIF
          ENDIF
          IF(missing==0 .or. legitSobol > p_legitimate) THEN
            LeadTerm=0
            allFinished = .TRUE.
            EXIT
          ENDIF
        END DO
        deallocate(sobol_trial)

7000    format(i8, 200f40.20)
270     format(i10)

    END SUBROUTINE solveMissingSobol

    SUBROUTINE findMissingSearch(COMPLETE)
        !Called by the main instance. Waits until function values have been derived
        !for all sobol points. If some points are missing (often due to a warm start
        !instance having been killed), solve for it.
        LOGICAL, INTENT(OUT) :: COMPLETE
        INTEGER(I4B) :: openStat, lastP, i,nummiss,lotPoint,seqn
        LOGICAL :: solvedPoints(p_maxpoints)
        REAL(DP) :: evalParam(p_nx)

        allocate(x_starts(p_maxpoints,p_nx+1))
        call myread2(x_starts,'x_starts.dat')

        solvedPoints = .FALSE.
        do while(any(solvedPoints == .FALSE.))

          call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=openStat, ACTION='read')
          DO
            READ (fileDesc,7460, END=10) seqn, lastP
            IF (lastP > p_maxpoints) THEN
                !This shouldn't happen. let's note the error and stop.
                write(errorString, *) seqNo, " found point: ",lastP,"greater than max: ",p_maxpoints
                call myclose(fileDesc)
                call exitState(errorString)
            END IF
            IF (lastP>0) solvedPoints(lastP) = .TRUE.
          END DO
  10      call myclose(fileDesc)

          nummiss=count(solvedPoints == .FALSE.)

          IF(nummiss==0) EXIT ! No missing local minimizations.

          !There are missing local minimizations.
          write(errorString, *) "Missing local search: ", nummiss, " missing local searches."
          call writeToLog(errorString);  print*,errorString

          DO i=1,p_maxpoints
              IF(solvedPoints(i))THEN
                  cycle
              END IF
              evalParam = getModifiedParam(p_maxpoints, x_starts(i,2:), p_searchType, lotPoint)
              write(errorString, *) seqNo," searching using sobol point ",i," and best point ",lotPoint
              call writeToLog(errorString);  print*,errorString

              !We now have the point at which to solve. So, finish the search
              call completeSearch(seqNo, alg, i, evalParam)
              solvedPoints(i) = .TRUE.
              exit
          END DO
          solvedPoints = .FALSE.
        ENDDO

        COMPLETE = .TRUE.

7460    format(3i10, 200f40.20)
    END SUBROUTINE findMissingSearch

    SUBROUTINE parseCommandLine(args_missing)
        !As the name suggests, parse the command line
        LOGICAL, INTENT(OUT) :: args_missing

        args_missing=.FALSE.
        temp=COMMAND_ARGUMENT_COUNT()

        IF (temp == 0) THEN
            print *,"Invoke with:"
            print *,"               ./genericSearch <-1|0|1|2|3|4> <configfile> <(a)moeba | (d)fpmin | (b)obyqa>"
            print *,"where"
            print *,"              -1 = stop all instances"
            print *,"               0 = cold start"
            print *,"               1 = warm start"
            print *,"               2 = update number of Sobol points"
            print *,"               3 = update number of local minimizationss"
            print *,"               4 = run diagnostics for given parameters"
            print *,"               5 = run local minimization for given parameters"
            print *,"           configfile is mandatory for all but warm start"
            print *,"           a,d,b = local minimization algorithm. Optional, default is (b) bobyqa"
            args_missing = .TRUE.
            return
        END IF

        IF (temp .GE. 1) THEN
            call GET_COMMAND_ARGUMENT(1, arg1, ierr)
            read (arg1,*) temp2
            isWarm = .FALSE.
            updateSobolPoints = .FALSE.
            updateLocalSearch = .FALSE.
            runDiagnostics = .FALSE.
            runLocalMin = .FALSE.
            option = temp2
            IF (temp2 == -1) THEN
                write(errorString, *) "ending all starts."
                call exitState(errorString)
            ELSE IF (temp2 == 1) THEN
                isWarm = .TRUE.
            ELSE IF (temp2 == 2) THEN
                updateSobolPoints = .TRUE.
              ELSE IF (temp2 == 3) THEN
                  updateLocalSearch = .TRUE.
            ELSE IF (temp2 == 4) THEN
                runDiagnostics = .TRUE.
            ELSE IF (temp2 == 5) THEN
                runLocalMin = .TRUE.
            END IF
        END IF

        IF ( (temp == 1) .AND. (isWarm .EQV. .FALSE.)) THEN
            print *, "<main> Error. Must specify a config file if not running warm start!"
            stop 0
        END IF

        IF (temp .GE. 2) THEN
            IF (isWarm) THEN
                CALL GET_COMMAND_ARGUMENT(2, arg2, ierr)
                SELECT CASE (arg2)
                    CASE ('a')
                        alg = 1   ! AMOEBA (Simplex algorithm)
                    CASE ('d')
                        alg = 2   ! DFPMIN (BFGS Quasi-Newton method)
                    CASE ('b')
                        alg = 0   ! BOBYQA (DFNLS algorithm)
                    CASE default
                        alg = p_default_alg
                END SELECT
            ELSE
                call GET_COMMAND_ARGUMENT(2, config, ierr)
                alg = p_default_alg
            END IF
        END IF

        IF (temp .GE. 3) THEN
            IF (isWarm) THEN
                print *,"<main>: Error. Cannot specify a config file for warm starts."
                STOP 0
            ELSE
                call GET_COMMAND_ARGUMENT(3, arg2, ierr)
                SELECT CASE (arg2)
                    CASE ('a')
                        alg = 1   ! AMOEBA (Simplex algorithm)
                    CASE ('d')
                        alg = 2   ! DFPMIN (BFGS Quasi-Newton method)
                    CASE ('b')
                        alg = 0   ! BOBYQA (DFNLS algorithm)
                    CASE default
                        alg = p_default_alg
                END SELECT
            END IF
        END IF

        IF (temp .GE. 4) THEN
            print *, "<main>: Error. Unknown parameters."
            STOP 0
        END IF
    END SUBROUTINE parseCommandLine

END PROGRAM GlobalSearch
