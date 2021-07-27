! OBJECTIVE.f90 contains the simulation and objective
!               (model-data distance) function
! This file must contain:
!   function objFun : used by amoeba (Nelder Mead)
!   subroutine dfovec: used by dfnls
!   other subroutines are specific to our income process function
!
! Variables from genericParams (values read from config.txt file) can be used here:
!   p_nx is the dimension of the problem
!   p_nmom is the number of moments
!   p_bound is the range of the search space
! These three parameters come from config.txt
! --------------------------------------------------------------------
! --------------------------------------------------------------------

MODULE OBJECTIVE
  USE genericParams, ONLY : p_nx,p_nmom,p_bound,p_fvalmax
  USE nrtype

  IMPLICIT NONE
  PRIVATE
  INTEGER :: diagnostic=0 ! default is 0 during the estimation,
                          ! only switch to 1 if diagnostic moments are estimated
  REAL(DP) :: penalty0, penaltyc, max_penalty
  REAL(DP), PARAMETER :: myeps=EPSILON(1.0_DP)
  PUBLIC objFun, dfovec, obj_initialize, GetData, getPenalty,diagnostic

CONTAINS

  ! call to GetData, possibly other subroutines, to load data
  SUBROUTINE obj_initialize
    ! this routine is called in the main program. it is used to load data in memory or other operations before minimization.
    IMPLICIT NONE
    CALL GetData
    penalty0=p_fvalmax; penaltyc=10*p_fvalmax
    max_penalty=1000.0_DP*penalty0
  END SUBROUTINE obj_initialize

  SUBROUTINE GetData
    ! this routine process all the data (reading data files etc...)
    USE nrtype
    USE genericParams

    ! this routine can be used to load data
    !....
    !....

  END SUBROUTINE GetData


  FUNCTION getPenalty(theta)
    ! this function adds a very high penalty when x is out of bounds
    USE genericParams
    IMPLICIT NONE

    REAL(DP),DIMENSION(p_nx),INTENT(IN) :: theta
    REAL(DP),DIMENSION(p_nx)            :: pencons,temp
    REAL(DP)                            :: penalty, getPenalty
    INTEGER :: ii,ind(1)

    getPenalty = 0.0_DP
    pencons=penaltyc

    temp=pencons*(MAX(0.0_DP,p_bound(:,1)-theta)/MAX(0.1_DP,dabs(p_bound(:,1))))**2 + &
             pencons*(MAX(0.0_DP,theta-p_bound(:,2))/MAX(0.1_DP,dabs(p_bound(:,2))))**2
    penalty=sum(temp)

    IF(penalty>myeps) getPenalty = min(max_penalty,penalty0 + penalty)
    IF(diagnostic==1 .and. penalty>myeps) THEN
      ind = MINLOC(temp,temp>myeps)
      WRITE(*,'(A25,f9.2,A12,I5,A8,f9.2)') "Bound penalty:",penalty," on index:",ind(1)," value:",theta(ind(1))
    ENDIF

  END FUNCTION getPenalty

  SUBROUTINE dfovec(np, nm, x, Fnout)
    ! This routine computes a np-dimensional vector
    ! (e.g. a vector of scaled distances between data moment and simulated moment.)
    ! As an example here, we compute the np-dimensional griewank function. Call its value val
    ! we then return the nm-dimensional vector Fnout = (val, val, val,..., val)

    USE nrtype
    USE genericParams
    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: np, nm
    REAL(DP), DIMENSION(np), INTENT(IN)  :: x
    REAL(DP), DIMENSION(nm),INTENT(OUT) :: Fnout
    REAL(DP) :: v_err

    ! USER: MODIFY THESE VARIABLES AS NEEDED
      REAL(DP) :: val1, val2
      INTEGER :: ii
    ! END OF USER MODIFICATION FOR VARIABLES

    v_err=getPenalty(x)
    IF(v_err>myeps) THEN
      Fnout = v_err/DSQRT(DBLE(nm))
    ELSE

    ! USER: MODIFY THIS PART AS NEEDED
      ! Fnout computation: using griewank as example
      val1 = (x(1)**2.0)/200.0
      do ii = 2, np
        val1 = val1 + (x(ii)**2.0)/200.0
      enddo
      val2 = cos(x(1))
      do ii = 2,np
          val2=val2*cos(x(ii)/sqrt(ii+0.0D0))
      enddo
      v_err = val1 - val2 + 1.0D0
      v_err = v_err + 1.0D0 !for relative criterion
      call sleep(3)
      ! using the same value for all the "moments"
      FORALL(ii=1:nm) Fnout(ii) = v_err
    ! END OF USER MODIFICATION
    ENDIF
    if(diagnostic==1) then
      ! Calculate detailed moments or other statistics that are not needed in the estimation
      print*,"v_err=",v_err
    endif

  END SUBROUTINE dfovec

  REAL(DP) FUNCTION objFun(theta)
    ! this function is the objective function. it calls dfovec which returns a vector Fout
    ! each component of the vector Fout is the distance between data moment and simulated moment, multiplied by a weight.
    ! this function computes the norm of Fout and adds a penalty if the points computed is out of the search space.
    ! -------------------------------------------------------------------------

    IMPLICIT NONE

    REAL(DP),DIMENSION(:),INTENT(IN) :: theta
    REAL(DP)                            :: objFun0(1,1)
    REAL(DP),DIMENSION(p_nmom)          :: Fout

    CALL dfovec(p_nx,p_nmom,theta,Fout)

    objFun0=DOT_PRODUCT(Fout,Fout)
    objFun = objFun0(1,1)+getPenalty(theta)

  END FUNCTION objFun

END MODULE objective
