!=================================================================
! 1. STATECONTROL: This module does all the file io stuff that manages
!                  control of states.
!=================================================================
MODULE stateControl
    USE nrtype
    IMPLICIT NONE

    INTEGER(I4B), parameter :: p_exitState=-1
CONTAINS

    SUBROUTINE myopen(unit, file, status, iostat, action, SHARE,position)
      !open a file. Locks files if they are opened for any purpose other
      !than reading. Locked files cannot be read/written until the
      !lock is removed.

      INTEGER(I4B), INTENT(OUT) :: unit
      CHARACTER(LEN=*), INTENT(IN) :: file, status, action,SHARE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: position
      INTEGER(I4B), INTENT(OUT) :: iostat
      CHARACTER(LEN=200) :: msg
      LOGICAL :: done,exists
      INTEGER(I4B), SAVE :: fileUnit = 200
      INTEGER(I4B) :: waitMax

      fileUnit = 200+mod(fileUnit+2,100)
      unit = fileUnit

      waitMax = 1000
      iostat=1
      DO WHILE(iostat > 0)
        IF (present(position)) THEN
            open(UNIT=unit,FILE=file,STATUS=status,IOSTAT=iostat, &
            ACTION=action,SHARE = SHARE,POSITION=position, err=111)
        ELSE
            open(UNIT=unit,FILE=file,STATUS=status,IOSTAT=iostat, &
            ACTION=action,SHARE = SHARE,err=111)
        END IF
        EXIT
111     waitMax = waitMax - 1
        IF(waitMax < 0) THEN
            msg = trim(file)//'. Cannot be accesed!'
            msg = 'Waited too long for; now EXIT '//trim(msg)
            call exitState(msg);
        END IF
        if(mod(waitMax,5)==0) call SLEEP(1)
      ENDDO

    END SUBROUTINE myopen

    SUBROUTINE myclose(unit)
        !closes a file. Removes the lock on any locked file.
        INTEGER(I4B), INTENT(IN) :: unit
        INTEGER(I4B) :: closeStat,waitMax
        CHARACTER(LEN=100) :: filename
        LOGICAL :: exists
        CHARACTER(LEN=200) :: msg

        waitMax = 30
        closeStat=1
        DO WHILE(closeStat > 0)
          close(UNIT=unit,IOSTAT=closeStat,err=121)
          EXIT
121       waitMax = waitMax - 1
          IF(waitMax < 0) THEN
              msg = ' Something is wrong with  myclose! '
              msg = ' Waited too long for; now EXIT. '//trim(msg)
              call writeToLog(msg); print*, msg
              INQUIRE(UNIT=unit, EXIST=exists, NAME=filename)
              print*, "closeStat= ", closeStat
              print*, "Unit#= ",  unit
              print*, "Exist? ", exists
              print*, "filename= ", filename
              stop
          END IF
          call SLEEP(1)
        ENDDO

        RETURN

    END SUBROUTINE myclose

    FUNCTION getNextNumber(file) RESULT (y)
      ! This function returns a the next number specified in the file.
      CHARACTER(LEN=*), INTENT(IN) :: file
      INTEGER(I4B) :: openStat, y,fileDesc,ios,wait

      wait = 0; ios=1
      DO WHILE(ios.NE.0)
         call myopen(UNIT=fileDesc, FILE=file, STATUS='unknown', IOSTAT=openstat, &
          ACTION='read',SHARE='DENYRW')
         READ(fileDesc,*,IOSTAT=ios) y
         call myclose(fileDesc)
         wait = wait + 1
         IF(ios.NE.0 .AND. mod(wait,10)==0) PRINT*,'getNextNumber, read: ',file,' ios= ', ios
      ENDDO
      y = y+1
      wait = 0; ios=1;
      DO WHILE(ios.NE.0)
        call myopen(UNIT=fileDesc, FILE=file, STATUS='replace', IOSTAT = openStat, &
         ACTION='write',SHARE='DENYRW')
        write(fileDesc, *,IOSTAT=ios) y
        call myclose(fileDesc)
        wait = wait + 1
        IF(ios.NE.0 .AND. mod(wait,10)==0) PRINT*,'getNextNumber, write: ',file,' ios= ', ios
      ENDDO
      return
    END FUNCTION getNextNumber

    FUNCTION getNumber(file) RESULT (y)
      ! This function returns the number specified in the file.
      CHARACTER(LEN=*), INTENT(IN) :: file
      INTEGER(I4B) :: openStat, y,fileDesc,ios,wait

      wait = 0; ios=1;
      DO WHILE(ios.NE.0)
         call myopen(UNIT=fileDesc, FILE=file, STATUS='unknown', IOSTAT=openstat, &
         ACTION='read',SHARE='DENYRW')
         READ(fileDesc,*,IOSTAT=ios) y
         call myclose(fileDesc)
         wait = wait + 1
         IF(ios.NE.0 .AND. mod(wait,10)==0) PRINT*,'getNumber, read: ',file,' ios= ', ios
      ENDDO
      return
    END FUNCTION getNumber

    SUBROUTINE setState(s,file)
      !Sets the the value s to the specified file.
      INTEGER(I4B), INTENT(IN) :: s
      CHARACTER(LEN=*), INTENT(IN) :: file
      INTEGER(I4B) :: openStat,fileDesc,ios,wait
      CHARACTER(LEN=1000) :: errorString

      wait = 0; ios=1;
      DO WHILE(ios.NE.0)
        call myopen(UNIT=fileDesc, FILE=file, STATUS='replace', IOSTAT=openStat, &
         ACTION='write',SHARE='DENYRW')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            write(errorString,*) "Error: unable to open the file ", file, " Cannot set state to ", s
            call exitState(errorString)
        ENDIF
        write(fileDesc, *,IOSTAT=ios) s
        call myclose(fileDesc)
        wait = wait + 1
        IF(ios.NE.0 .AND. mod(wait,10)==0) PRINT*,'setState: ',file,' s= ',s,' ios= ', ios
      ENDDO
    END SUBROUTINE setState

    FUNCTION getState() RESULT (y)
        !returns the current state value
        INTEGER(I4B) :: y, openStat,fileDesc,ios,wait

        wait = 0; ios=1;
        DO WHILE(ios.NE.0)
          call myopen(UNIT=fileDesc, FILE='state.dat', STATUS='old', IOSTAT=openStat, &
           ACTION='read',SHARE='DENYRW')
          IF (openStat /= 0) THEN
              ! file doesn't exist. This is an error
              print *,"<stateControl.getState()> Error: unable to open state file."
              stop 0
          END IF
          read(fileDesc,*,IOSTAT=ios) y
          call myclose(fileDesc)
          wait = wait + 1
          IF(ios.NE.0 .AND. mod(wait,10)==0) PRINT*,'getState:, ios= ', ios
        ENDDO
    END FUNCTION getState

    SUBROUTINE waitState(s)
        !instance waits until the current state is the one requested.
        INTEGER(I4B), INTENT(IN) :: s
        INTEGER(I4B) :: currentstate
        CHARACTER(LEN=1000) :: errorString

        currentstate = getState()
        DO WHILE ( (currentstate < s))
          IF (currentstate == p_exitState) then
            write(errorString,*) "ExitState: called for exit "
            call exitState(errorString)
          endif
          call SLEEP(30)
          currentstate = getState()
        END DO
    END SUBROUTINE waitState

    SUBROUTINE writeToLog(str)
        !writes a string to the log file
        CHARACTER(LEN=*), INTENT(IN) :: str
        INTEGER(I4B) ::  iostat,fileDesc,ios,wait

        wait = 0; ios=1;
        DO WHILE(ios.NE.0)
          call myopen(unit=fileDesc, file='logFile.txt',status='unknown',IOSTAT=iostat, &
           ACTION='write',SHARE='DENYRW',position="append")
          write(fileDesc,*,IOSTAT=ios) trim(str)
          call myclose(fileDesc)
          wait = wait + 1
          IF(ios.NE.0 .AND. mod(wait,10)==0)  PRINT*,'writeToLog: ', trim(str)
        ENDDO
    END SUBROUTINE writeToLog

    SUBROUTINE exitState(str)
        !instance waits until the current state is the one requested.
        CHARACTER(LEN=*), INTENT(IN) :: str

        call setState(p_exitState,'state.dat')
        call setState(-10000000,'lastSobol.dat')
        call setState(-10000000,'lastParam.dat')
        call setState(-10000000,'seqNo.dat')

        call writeToLog(str); print*,trim(str)

        stop 0

    END SUBROUTINE exitState

    SUBROUTINE myread2(myarray2, filename2, numrows)
    IMPLICIT NONE

    !dummy arguments
    REAL(DP), DIMENSION(:,:), INTENT(OUT):: myarray2
    CHARACTER(LEN =*), INTENT(IN) :: filename2
    INTEGER, INTENT(OUT), OPTIONAL ::  numrows
    !local
    INTEGER :: count,openStat,fileDesc,wait
    INTEGER, DIMENSION(2) :: fileDIM
    CHARACTER(LEN =1000) :: errorString

    fileDIM = SHAPE(myarray2)
    call myopen(UNIT=fileDesc, FILE=filename2, STATUS='old', IOSTAT=openStat, &
    ACTION='read',SHARE='DENYRW')
    IF (openStat /= 0) THEN
        ! file doesn't exist. This is an error
        call myclose(fileDesc)
        print*,"<myread2> Error: unable to open file: ", filename2
        write(errorString,*)"<myread2> Error: unable to open file: ", filename2
        call writeToLog(errorString)
        stop 0
    END IF

    count=1; wait = 0;
    DO WHILE (count <= fileDIM(1))
        ! read until we run out of lines. Possible that file is really large (by mistake), so
        ! stop if we've reached the maximum number of iterations
        read(fileDesc,*,IOSTAT=openStat,end=45) myarray2(count,:)
        IF(openStat /= 0) THEN
          wait = wait + 1
          IF(mod(wait,10)==0) PRINT*,'myread2, read(fileDesc,*,IOSTAT=openStat,end=45)', openStat
          openStat=0
       ELSE
        count=count+1
       ENDIF
    END DO
45  call myclose(fileDesc)
    If(PRESENT(numrows)) numrows=count-1
    RETURN
  END SUBROUTINE  myread2

  subroutine mywrite2(myarray2, filename2)
    IMPLICIT NONE

    !dummy arguments
    REAL(DP), DIMENSION(:,:), INTENT(IN):: myarray2
    CHARACTER(LEN =*), INTENT(IN) :: filename2

    !local
    INTEGER, DIMENSION(2) :: fileDIM
    INTEGER:: i,j,openStat,fileDesc

    fileDIM= SHAPE(myarray2)
    if (fileDIM(2)>1000) then
      write (*,*) "ERROR THE SECOND DIMENSION OF FILE >1000, FORMAT STATEMENT CURRENTLY USES CANNOT HANDLE THIS SIZE"
    end if

    call myopen(UNIT=fileDesc, FILE=filename2, STATUS='replace', IOSTAT=openStat, &
     ACTION='write',SHARE='DENYRW')
    DO i=1,fileDIM(1)
      WRITE(fileDesc,100) (myarray2(i,j), j=1,fileDIM(2))
    END DO
    call myclose(fileDesc)

100 FORMAT(1x,1000(3x,F30.15))

  end subroutine  mywrite2

END MODULE stateControl
