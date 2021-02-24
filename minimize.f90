module minimize
    !This module contains all the minimization routines: amoeba and bobyq
    use nrtype
  REAL(DP) :: eps=epsilon(1.0_DP)

contains
    SUBROUTINE amoeba(p,y,ftol,func,iter,ITMAX_USER,global_fval)
        USE utilities, ONLY : assert_eq,imaxloc,iminloc,nrerror,swap
        IMPLICIT NONE
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p ! vertices. If we have n vertices, then we must be
                                                         ! in n-1 dimensional space (we need one extra vertex
                                                         ! than dimensions. For each row, the n-1 vector
                                                         ! specifies the vertex
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: y   ! "func" evaluated at the n vertices provided in "p"
        REAL(DP), INTENT(IN) :: ftol                 ! tolerance for amoeba
        INTERFACE
          REAL(DP)  FUNCTION func(x)
                use genericParams
                IMPLICIT NONE
                REAL(DP), DIMENSION(:), INTENT(IN) :: x
!                REAL(DP) :: func
            END FUNCTION func
        END INTERFACE
        INTEGER(I4B), INTENT(OUT) :: iter             ! the number of iterations required
        REAL(DP), INTENT(IN), OPTIONAL :: global_fval ! the best global objective value so far
        INTEGER(I4B), INTENT(IN),  OPTIONAL :: ITMAX_USER
        INTEGER(I4B) :: ITMAX=5000
        REAL(dp), PARAMETER :: TINY=1.0e-10
        INTEGER(I4B) :: ihi,ndim
        REAL(dp), DIMENSION(size(p,2)) :: psum
        if(PRESENT(ITMAX_USER))  THEN
          ITMAX=ITMAX_USER
        else
          ITMAX=5000
        ENDIF
        call amoeba_private

        CONTAINS
        !BL
        SUBROUTINE amoeba_private
            IMPLICIT NONE
            INTEGER(I4B) :: i,ilo,inhi,j
            REAL(dp) :: rtol,ysave,ytry,ytmp
            REAL(DP) :: ylo_pre,improv
            INTEGER :: iter_pre

            ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
            iter=0
            psum(:)=sum(p(:,:),dim=1)
            ihi=imaxloc(y(:))
            ylo_pre=y(ihi)
            iter_pre=0
            do
              ilo=iminloc(y(:))
              ihi=imaxloc(y(:))
              ytmp=y(ihi)
              y(ihi)=y(ilo)
              inhi=imaxloc(y(:))
              y(ihi)=ytmp
              rtol=2.0_DP*dabs(y(ihi)-y(ilo))/(dabs(y(ihi))+dabs(y(ilo))+TINY)
              if(mod(iter,100)==0) THEN
                improv=ylo_pre-y(ilo)
                ylo_pre=y(ilo)
                WRITE (*,*) 'in ameoba,  iteration =', iter
                100       FORMAT(i2,2i10,5i12)
                WRITE(*,*)     ' ya(ilo)     ya(ihi)      rtol     improvement'
                WRITE(6,101) y(ilo),y(ihi),rtol,improv
                101       FORMAT(5f12.6)

                IF(PRESENT(global_fval) .AND. iter>=400 .AND. improv>0.0_DP) THEN
                   IF((y(ilo)-improv*DBLE(ITMAX-iter)/DBLE(iter-iter_pre)) > 1.02_DP*global_fval) THEN
                      PRINT*,'NOT ENOUGH IMPROVEMENT'
                      RETURN
                   ENDIF
                ENDIF
                iter_pre=iter
              ENDIF
              ilo=iminloc(y(:))
              ihi=imaxloc(y(:))
              ytmp=y(ihi)
              y(ihi)=y(ilo)
              inhi=imaxloc(y(:))
              y(ihi)=ytmp
              rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
              if (rtol < ftol) then
                  call swap(y(1),y(ilo))
                  call swap(p(1,:),p(ilo,:))
                  RETURN
              end if
              !if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
              if (iter >= ITMAX) then
                  write(*,*) 'reached max number of evaluations in amoeba; returning best result so far'
                  call swap(y(1),y(ilo))
                  call swap(p(1,:),p(ilo,:))
                  RETURN
              end if
              ytry=amotry(-1.0_dp)
              iter=iter+1
              if (ytry <= y(ilo)) then
                  ytry=amotry(2.0_dp)
                  iter=iter+1
              else if (ytry >= y(inhi)) then
                  ysave=y(ihi)
                  ytry=amotry(0.5_dp)
                  iter=iter+1
                  if (ytry >= ysave) then
                      p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
                      do i=1,ndim+1
                          if (i /= ilo) y(i)=func(p(i,:))
                      end do
                      iter=iter+ndim
                      psum(:)=sum(p(:,:),dim=1)
                  end if
              end if
            end do
        END SUBROUTINE amoeba_private
        !BL
        FUNCTION amotry(fac)
            IMPLICIT NONE
            REAL(dp), INTENT(IN) :: fac
            REAL(dp) :: amotry
            REAL(dp) :: fac1,fac2,ytry
            REAL(dp), DIMENSION(size(p,2)) :: ptry
            fac1=(1.0_dp-fac)/ndim
            fac2=fac1-fac
            ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
            ytry=func(ptry)
            if (ytry < y(ihi)) then
                y(ihi)=ytry
                psum(:)=psum(:)-p(ihi,:)+ptry(:)
                p(ihi,:)=ptry(:)
            end if
            amotry=ytry
        END FUNCTION amotry
    END SUBROUTINE amoeba

    SUBROUTINE amebsa(p,y,pb,yb,ftol,func,iter,temptr)
        USE utilities, ONLY : assert_eq,imaxloc,iminloc,swap
        IMPLICIT NONE
        INTEGER(I4B), INTENT(INOUT) :: iter
        REAL(DP), INTENT(INOUT) :: yb
        REAL(DP), INTENT(IN) :: ftol,temptr
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: y,pb
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p
        INTERFACE
            FUNCTION func(x)
                USE nrtype
                IMPLICIT NONE
                REAL(DP), DIMENSION(:), INTENT(IN) :: x
                REAL(DP) :: func
            END FUNCTION func
        END INTERFACE
        INTEGER(I4B), PARAMETER :: NMAX=200
        INTEGER(I4B) :: ihi,ndim
        REAL(DP) :: yhi
        REAL(DP), DIMENSION(size(p,2)) :: psum
        call amebsa_private
    CONTAINS
        !BL
        SUBROUTINE amebsa_private
            INTEGER(I4B) :: i,ilo,inhi
            REAL(DP) :: rtol,ylo,ynhi,ysave,ytry
            REAL(DP), DIMENSION(size(y)) :: yt,harvest
            ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,size(pb),'amebsa')
            psum(:)=sum(p(:,:),dim=1)
            do
                call random_number(harvest)
                yt(:)=y(:)-temptr*log(harvest)
                ilo=iminloc(yt(:))
                ylo=yt(ilo)
                ihi=imaxloc(yt(:))
                yhi=yt(ihi)
                yt(ihi)=ylo
                inhi=imaxloc(yt(:))
                ynhi=yt(inhi)
                rtol=2.0_dp*abs(yhi-ylo)/(abs(yhi)+abs(ylo))
                if (rtol < ftol .or. iter < 0) then
                    call swap(y(1),y(ilo))
                    call swap(p(1,:),p(ilo,:))
                    RETURN
                end if
                ytry=amotsa(-1.0_dp)
                iter=iter-1
                if (ytry <= ylo) then
                    ytry=amotsa(2.0_dp)
                    iter=iter-1
                else if (ytry >= ynhi) then
                    ysave=yhi
                    ytry=amotsa(0.5_dp)
                    iter=iter-1
                    if (ytry >= ysave) then
                        p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
                        do i=1,ndim+1
                            if (i /= ilo) y(i)=func(p(i,:))
                        end do
                        iter=iter-ndim
                        psum(:)=sum(p(:,:),dim=1)
                    end if
                end if
            end do
        END SUBROUTINE amebsa_private
        !BL
        FUNCTION amotsa(fac)
            IMPLICIT NONE
            REAL(DP), INTENT(IN) :: fac
            REAL(DP) :: amotsa
            REAL(DP) :: fac1,fac2,yflu,ytry,harv
            REAL(DP), DIMENSION(size(p,2)) :: ptry
            fac1=(1.0_dp-fac)/ndim
            fac2=fac1-fac
            ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
            ytry=func(ptry)
            if (ytry <= yb) then
                pb(:)=ptry(:)
                yb=ytry
            end if
            call random_number(harv)
            yflu=ytry+temptr*log(harv)
            if (yflu < yhi) then
                y(ihi)=ytry
                yhi=yflu
                psum(:)=psum(:)-p(ihi,:)+ptry(:)
                p(ihi,:)=ptry(:)
            end if
            amotsa=yflu
        END FUNCTION amotsa
    END SUBROUTINE amebsa

!  Important Notice:
!   These algorithms are modifications and based on the Software BOBYQA, authored by M. J. D. Powell,
!   to minimize sum of squares with bound constraints by taking advantage of the problem structure.
!   i.e. Min  F(x) := Sum_{i=1}^{mv}  v_err_i(x)^2, s.t. xl <= x <= xu, x \in R^n,
!   where v_err(x) : R^n \to R^{mv} is a vector function.
!   This subroutine seeks the least value of sum of the squres of the components of v_err(x)
!   by combing trust region method and Levenberg-Marquardt method
!
!   References:
!
!   1.  M. J. D. Powell, The NEWUOA software for unconstrained optimization without derivatives,
!       DAMTP 2004/ NA 05
!   2.  M. J. D. Powell, The BOBYQA algorithm for bound constrained optimization without derivatives,
!       DAMTP 2009/ NA 06
!   3.  H. Zhang, A. R. CONN, AND K. SCHEINBERG, A DERIVATIVE-FREE ALGORITHM FOR THE LEAST-SQUARE
!       MINIMIZATION, technical report, 2009
!
!      -----------------------------------------------------------------
!      | This program is free software; you can redistribute it and/or  |
!      |modify it under the terms of the GNU General Public License as  |
!      |published by the Free Software Foundation; either version 2 of  |
!      |the License, or (at your option) any later version.             |
!      |This program is distributed in the hope that it will be useful, |
!      |but WITHOUT ANY WARRANTY; without even the implied warranty of  |
!      |MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   |
!      |GNU General Public License for more details.                    |
!      |                                                                |
!      |You should have received a copy of the GNU General Publi!       |
!      |License along with this program; if not, write to the Free      |
!      |Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, |
!      |MA  02110-1301  USA                                             |
!      -----------------------------------------------------------------|

      SUBROUTINE BOBYQA_H (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT, MAXFUN,W,mv)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),W(*)

!     N must be set to the number of variables and must be at least two.
!     mv must be set to the lengh of the vector function v_err(x):  R^n \to R^{mv}.
!     The maximum number variables in this codes: nmax =  100
!     The maximum lengh of the vector function v_err(x): mmax = 400
!     If n > 100 or m > 400, the parameter nmax and mmax need to be creased in
!     subroutine  BOBYQB_H, PRELIM_H, RESCUE_H and TRSBOX_H

!     NPT is the number of interpolation conditions. Its value must be in the
!     interval [N+2,2N+1].  Recommended: NPT = 2*N+1
!
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!     will be changed to the values that give the least calculated F= Sum_{i=1}^{mv} v_err_i(x)^2..
!
!     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
!       bounds, respectively, on X(I). The construction of quadratic models
!       requires XL(I) to be strictly less than XU(I) for each I. Further,
!       the contribution to a model from changes to the I-th variable is
!       damaged severely by rounding errors if XU(I)-XL(I) is too small.

!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while RHOEND should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG. Default: RHOBEG = 1.0, RHOEND = 10^{-8}

!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.

!     MAXFUN must be set to an upper bound on the number of calls of subroutine
!     dfovec(n, mv, x, v_err) which provides the values of the vector function v_err(x).
!     Here: n, mv, x \in R^n are input, v_err \in R^{mv} are output.
!     Default:  MAXFUN= 400(n+1), i.e 400 (simplex) gradients for reasonable accuracy.
!               MAXFUN= infinity, to let the algorithm explore the lowest function value
!                       as much as it could.

!     The array W will be used for working space. Its length must be at least
!       (NPT+5)*(NPT+N)+3*N*(N+5)/2.

!     SUBROUTINE dfovec(n, mv, x, v_err) must be provided by the user.
!     It must provide the values of the vector function v_err(x) : R^n to R^{mv}
!     at the variables X(1),X(2),...,X(N), which are generated automatically in
!     a way that satisfies the bounds given in XL and XU.

!     Return if the value of NPT is unacceptable.

      NP=N+1
      IF (NPT < N+2 .OR. NPT > 2*N+1) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in','[N+2, 2N+1]')
          GO TO 40
      END IF

!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.

      NDIM=NPT+N
      IXB=1
      IXP=IXB+N
      IFV=IXP+N*NPT
      IXO=IFV+NPT
      IGO=IXO+N
      IHQ=IGO+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ISL=IZMAT+NPT*(NPT-NP)
      ISU=ISL+N
      IXN=ISU+N
      IXA=IXN+N
      ID=IXA+N
      IVL=ID+N
      IW=IVL+NDIM

!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.

      ZERO=0.0_DP
      DO 30 J=1,N
      TEMP=XU(J)-XL(J)
      IF (TEMP < RHOBEG+RHOBEG) THEN
          PRINT 20
   20     FORMAT (/4X,'Return from BOBYQA_H because one of the',&
     &     ' differences XU(I)-XL(I)'/6X,' is less than 2*RHOBEG.')
          GO TO 40
      END IF
      JSL=ISL+J-1
      JSU=JSL+N
      W(JSL)=XL(J)-X(J)
      W(JSU)=XU(J)-X(J)
      if (w(jsl) >= -rhobeg) then
          if (w(jsl) >= zero) then
              x(j)=xl(j)
              w(jsl)=zero
              w(jsu)=temp
          else
              x(j)=xl(j)+rhobeg
              w(jsl)=-rhobeg
              w(jsu)=dmax1(xu(j)-x(j),rhobeg)
          end if
      else if (w(jsu) <= rhobeg) then
          if (w(jsu) <= zero) then
              x(j)=xu(j)
              w(jsl)=-temp
              w(jsu)=zero
          else
              x(j)=xu(j)-rhobeg
              w(jsl)=dmin1(xl(j)-x(j),-rhobeg)
              w(jsu)=rhobeg
          end if
      end if

   30 CONTINUE

!     Make the call of BOBYQB_H.

      CALL BOBYQB_H (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),&
       W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),&
      NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW),mv)
   40 RETURN
      END subroutine


!
!   Important Notice:
!   This NEWUOB_H are modifications and based on the subroutine BOBYQB in the software
!   BOBYQA, authored by M. J. D. Powell.
!
      SUBROUTINE BOBYQB_H (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,&
           MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,&
           SL,SU,XNEW,XALT,D,VLAG,W,mv)
      use OBJECTIVE, only: dfovec, getPenalty
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),&
       XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),&
       SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)


      integer nmax, mmax, nptmax
      parameter (nmax = 100, mmax=3000, nptmax=2*nmax+1)
      dimension GQV(mmax,nmax),HQV(mmax,(nmax+1)*nmax/2),&
     &          PQV(mmax,nptmax), PQVold(mmax),&
     &          FVAL_V(mmax,nptmax),v_sum(nptmax),&
     &          WV(mmax,nmax), v_sumpq(mmax), &
     &          v_err(mmax), v_beg(mmax), v_temp(mmax),&
     &          v_diff(mmax), v_vquad(mmax),&
     &          HD1(nmax)
      logical model_update, opt_update
      integer v_itest(mmax)
  !   external dfovec


      if (n>nmax) then
        print *, "in bobyqb_h.f increase the dimension  &
     &            nmax to be at least", n
        stop
      endif
      if (mv>mmax) then
        print *, "in bobyqb_h.f increase the dimension &
     &            mmax to be at least", mv
        stop
      endif

!     Set some constants.

      model_update = .true.
      opt_update = .true.
      HALF=0.5_DP
      ONE=1.0_DP
      TEN=10.0_DP
      TENTH=0.1_DP
      TWO=2.0_DP
      ZERO=0.0_DP
      NP=N+1
      NPTM=NPT-NP
      NH=(N*NP)/2

!     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, with the corresponding values of
!     of NF and KOPT, which are the number of calls of CALFUN so far and the
!     index of the interpolation point at the trust region centre. Then the
!     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
!     less than NPT. GOPT will be updated if KOPT is different from KBASE.

      CALL PRELIM_H (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT,&
       FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,&
       mv,HQV,GQV,PQV,FVAL_V,v_err,v_beg,FBEG)
      XOPTSQ=ZERO
      DO 10 I=1,N
      XOPT(I)=XPT(KOPT,I)
   10 XOPTSQ=XOPTSQ+XOPT(I)**2
      FSAVE=FVAL(1)
      IF (NF < NPT) THEN
          IF (IPRINT > 0) PRINT 390
          GOTO 720
      END IF
      KBASE=1

!     Complete the settings that are required for the iterative procedure.

      RHO=RHOBEG
      DELTA=RHO
      NRESC=NF
      NTRITS=0
      DIFFA=ZERO
      DIFFB=ZERO
      do 15 m1=1,mv
   15 v_itest(m1)=0
      NFSAV=NF

!     Update GOPT if necessary before the first iteration and after each
!     call of RESCUE that makes a call of CALFUN.

   20 IF (KOPT .NE. KBASE) THEN
          model_update = .true.
          IH=0
          DO 30 J=1,N
          DO 30 I=1,J
          IH=IH+1
          IF (I < J) then
             do 24 m1=1,mv
               GQV(m1,J)=GQV(m1,J)+HQV(m1,IH)*XOPT(I)
   24        continue
          endif
          do 26 m1=1,mv
             GQV(m1,I)=GQV(m1,I)+HQV(m1,IH)*XOPT(J)
   26     continue
   30     continue
          IF (NF > NPT) THEN
              DO 50 K=1,NPT
              TEMP=ZERO
              DO 40 J=1,N
   40         TEMP=TEMP+XPT(K,J)*XOPT(J)
              do 45 m1=1,mv
                 t=PQV(m1,k)*TEMP
                 do 45 i=1,n
   45            GQV(m1,i)=GQV(m1,i)+t*XPT(k,i)
   50         continue
          END IF
      END IF

!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.

   60 CALL TRSBOX_H (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,&
      W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN,&
       mv,HQV,GQV,PQV,FVAL_V,KOPT,HD1,v_temp,vquad1,&
       model_update,opt_update)
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM < HALF*RHO) THEN
          NTRITS=-1
          DISTSQ=(TEN*RHO)**2
          IF (NF <= NFSAV+2) GOTO 650

!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.

          ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
          FRHOSQ=0.125_DP*RHO*RHO
          IF (CRVMIN > ZERO .AND. ERRBIG > FRHOSQ*CRVMIN)    GOTO 650
          BDTOL=ERRBIG/RHO
          DO 80 J=1,N
          BDTEST=BDTOL
          IF (dabs(XNEW(J)-SL(J)) <= eps) BDTEST=W(J)
          IF (dabs(XNEW(J)-SU(J)) <= eps) BDTEST=-W(J)
          IF (BDTEST < BDTOL) THEN
              CURV=HQ((J+J*J)/2)
              BDTEST=BDTEST+HALF*CURV*RHO
              IF (BDTEST < BDTOL) GOTO 650
          END IF
   80     CONTINUE
          GOTO 680
      END IF
      NTRITS=NTRITS+1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     zero. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.

!   90 IF (DSQ <= 1.0D-3*XOPTSQ) THEN
   90 IF (DSQ <= 1.0D-1*XOPTSQ) THEN
          model_update = .true.
          FRACSQ=0.25_DP*XOPTSQ
          do 95 m1=1,mv
   95     v_sumpq(m1)=ZERO
          DO 110 K=1,NPT
          do 98 m1=1,mv
   98     v_sumpq(m1)= v_sumpq(m1)+PQV(m1,K)
          SUM=-HALF*XOPTSQ
          DO 100 I=1,N
  100     SUM=SUM+XPT(K,I)*XOPT(I)
          W(NPT+K)=SUM
          TEMP=FRACSQ-HALF*SUM
          DO 110 I=1,N
          W(I)=BMAT(K,I)
          VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
          IP=NPT+I
          DO 110 J=1,I
  110     BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)

!     Then the revisions of BMAT that depend on ZMAT are calculated.

          DO 150 JJ=1,NPTM
          SUMZ=ZERO
          SUMW=ZERO
          DO 120 K=1,NPT
          SUMZ=SUMZ+ZMAT(K,JJ)
          VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
  120     SUMW=SUMW+VLAG(K)
          DO 140 J=1,N
          SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
          DO 130 K=1,NPT
  130     SUM=SUM+VLAG(K)*XPT(K,J)
          W(J)=SUM
          DO 140 K=1,NPT
  140     BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
          DO 150 I=1,N
          IP=I+NPT
          TEMP=W(I)
          DO 150 J=1,I
  150     BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)

!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.

          IH=0
          DO 170 J=1,N
          do 152 m1=1,mv
  152     WV(m1,J)=-HALF*v_sumpq(m1)*XOPT(J)
          DO 160 K=1,NPT
          do 154 m1=1,mv
  154     WV(m1,J)=WV(m1,J)+PQV(m1,K)*XPT(K,J)
  160     XPT(K,J)=XPT(K,J)-XOPT(J)
          DO 170 I=1,J
          IH=IH+1
          do 165 m1=1,mv
            HQV(m1,IH)=HQV(m1,IH)+WV(m1,I)*XOPT(J) &
                      +XOPT(I)*WV(m1,J)
  165     continue
  170     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 180 I=1,N
          XBASE(I)=XBASE(I)+XOPT(I)
          XNEW(I)=XNEW(I)-XOPT(I)
          SL(I)=SL(I)-XOPT(I)
          SU(I)=SU(I)-XOPT(I)
  180     XOPT(I)=ZERO
          XOPTSQ=ZERO
      END IF
      IF (NTRITS .EQ. 0) GOTO 210
      GOTO 230

!     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
!     more expensive than the previous shift, because new matrices BMAT and
!     ZMAT are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. Therefore RESCUE is called
!     only if rounding errors have reduced by at least a factor of two the
!     denominator of the formula for updating the H matrix. It provides a
!     useful safeguard, but is not invoked in most applications of BOBYQA.

  190 NFSAV=NF
      KBASE=KOPT
      CALL RESCUE_H (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL,&
       XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT,&
       VLAG,W,W(N+NP),W(NDIM+NP),&
       mv,HQV,GQV,PQV,FVAL_V,WV,v_err,v_sumpq,v_vquad,v_diff,v_temp)
      model_update = .true.

!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.

      XOPTSQ=ZERO
      IF (KOPT .NE. KBASE) THEN
          DO 200 I=1,N
          XOPT(I)=XPT(KOPT,I)
  200     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      IF (NF < 0) THEN
          NF=MAXFUN
          IF (IPRINT > 0) PRINT 390
          GOTO 720
      END IF
      NRESC=NF
      IF (NFSAV < NF) THEN
          NFSAV=NF
          GOTO 20
      END IF
      IF (NTRITS > 0) GOTO 60

!     Pick two alternative vectors of variables, relative to XBASE, that
!     are suitable as new positions of the KNEW-th interpolation point.
!     Firstly, XNEW is set to the point on a line through XOPT and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
!     and SU bounds. Secondly, XALT is set to the best feasible point on
!     a constrained version of the Cauchy step of the KNEW-th Lagrange
!     function, the corresponding value of the square of this function
!     being returned in CAUCHY. The choice between these alternatives is
!     going to be made when the denominator is calculated.


  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,&
      KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
      DO 220 I=1,N
  220 D(I)=XNEW(I)-XOPT(I)

!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.

  230 DO 250 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 240 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  240 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
      VLAG(K)=SUM
  250 W(NPT+K)=SUMA
      BETA=ZERO
      DO 270 JJ=1,NPTM
      SUM=ZERO
      DO 260 K=1,NPT
  260 SUM=SUM+ZMAT(K,JJ)*W(K)
      BETA=BETA-SUM*SUM
      DO 270 K=1,NPT
  270 VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
      DSQ=ZERO
      BSUM=ZERO
      DX=ZERO
      DO 300 J=1,N
      DSQ=DSQ+D(J)**2
      SUM=ZERO
      DO 280 K=1,NPT
  280 SUM=SUM+W(K)*BMAT(K,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 290 I=1,N
  290 SUM=SUM+BMAT(JP,I)*D(I)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  300 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
!     If NTRITS is zero, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.

      IF (NTRITS .EQ. 0) THEN
          DENOM=VLAG(KNEW)**2+ALPHA*BETA
          IF (DENOM < CAUCHY .AND. CAUCHY > ZERO) THEN
              DO 310 I=1,N
              XNEW(I)=XALT(I)
  310         D(I)=XNEW(I)-XOPT(I)
              CAUCHY=ZERO
              GO TO 230
          END IF
          IF (DENOM <= HALF*VLAG(KNEW)**2) THEN
              IF (NF > NRESC) GOTO 190
              IF (IPRINT > 0) PRINT 320
  320         FORMAT (/5X,'Return from BOBYQA_H because of much',&
               ' cancellation in a denominator.')
              GOTO 720
          END IF

!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.

      ELSE
          DELSQ=DELTA*DELTA
          SCADEN=ZERO
          BIGLSQ=ZERO
          KNEW=0
          DO 350 K=1,NPT
          IF (K .EQ. KOPT) GOTO 350
          HDIAG=ZERO
          DO 330 JJ=1,NPTM
  330     HDIAG=HDIAG+ZMAT(K,JJ)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          DISTSQ=ZERO
          DO 340 J=1,N
  340     DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
          TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
          IF (TEMP*DEN > SCADEN) THEN
              SCADEN=TEMP*DEN
              KNEW=K
              DENOM=DEN
          END IF
          BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
  350     CONTINUE
          IF (SCADEN <= HALF*BIGLSQ) THEN
              IF (NF > NRESC) GOTO 190
              IF (IPRINT > 0) PRINT 320
              GOTO 720
          END IF
      END IF

!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.

!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.

  360 DO 380 I=1,N
      X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
      IF (dabs(XNEW(I)-SL(I)) <= eps) X(I)=XL(I)
      IF (dabs(XNEW(I)-SU(I)) <= eps) X(I)=XU(I)
  380 CONTINUE
      IF (NF >= MAXFUN) THEN
          IF (IPRINT > 0) PRINT 390
  390     FORMAT (/4X,'Return from BOBYQA_H because CALFUN has been',&
          ' called MAXFUN times.')
          GOTO 720
      END IF
      NF=NF+1
      IF (IPRINT .EQ. 3) THEN
         PRINT 398, (X(I),I=1,N)
  398      FORMAT (/4X,'Before the call to CALFUN',&
            '    The corresponding X is:'/(2X,5D15.6))

      END IF

!     dfovec(n, mv, x, v_err) provides the values of the vector function v_err(x): R^n \to R^{mv}.
!     Here: n, mv, x \in R^n are input, v_err \in R^{mv} are output.
      call dfovec(n, mv, x, v_err)

!     f_value(mv,v_err,F) provides the value of the sum of the squres of the components of v_err(x)
!     i.e. F = sum_{i=1}^{mv} v_err_i (x)^2
      call f_value(mv,v_err,F)

	  !LS: Add penalty to result
	  F = F + getPenalty(x)

      IF (IPRINT .EQ. 3) THEN
          PRINT 400, NF,F,(X(I),I=1,N)
  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,&
            '    The corresponding X is:'/(2X,5D15.6))
      END IF

      if(F<=dmax1(1.d-12,1.d-20*FBEG)) then
         print *, ""
         print *, "   Return: F<=dmax1(1.d-12,1.d-20*FBEG)"
         go to 720
      endif

      IF (NTRITS .EQ. -1) THEN
          FSAVE=F
          GOTO 720
      END IF

!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.

!      VQUAD=ZERO
      do 405 m1=1,mv
  405 v_vquad(m1)=ZERO
      IH=0
      DO 410 J=1,N
      do 408 m1=1,mv
  408 v_vquad(m1)=v_vquad(m1)+D(J)*GQV(m1,J)
      DO 410 I=1,J
      IH=IH+1
      TEMP=D(I)*D(J)
      IF (I .EQ. J) TEMP=HALF*TEMP
      do 410 m1=1,mv
  410 v_vquad(m1)=v_vquad(m1)+HQV(m1,IH)*TEMP
      DO 420 K=1,NPT
      do 420 m1=1,mv
  420 v_vquad(m1)=v_vquad(m1)+HALF*PQV(m1,K)*W(NPT+K)**2
      do 425 m1=1,mv
  425 v_diff(m1)=v_err(m1)-FVAL_V(m1,KOPT)-v_vquad(m1)
      if (NTRITS <= 0) then
         DO 426 I=1,N
  426    HD1(I) = zero
         IH=0
         DO 427 J=1,N
         DO 427 I=1,J
         IH=IH+1
         IF (I < J) HD1(J)=HD1(J)+HQ(IH)*D(I)
  427    HD1(I)=HD1(I)+HQ(IH)*D(J)
         vquad1 = zero
         do 428 i=1,n
  428    vquad1 = vquad1 + D(i)*(GOPT(i)+HALF*HD1(i))
      endif
      FOPT=FVAL(KOPT)
      DIFF=F-FOPT-VQUAD1
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM > RHO) NFSAV=NF

!     Pick the next value of DELTA after a trust region step.

      IF (NTRITS > 0) THEN
          IF (VQUAD1 >= ZERO) THEN
              IF (IPRINT > 0) PRINT 430
  430         FORMAT (/4X,'Return from BOBYQA_H because a trust',&
               ' region step has failed to reduce Q.')
              GOTO 720
          END IF
          RATIO=(F-FOPT)/VQUAD1
          IF (RATIO <= TENTH) THEN
              DELTA=DMIN1(HALF*DELTA,DNORM)
          ELSE IF (RATIO <= 0.7_DP ) THEN
              DELTA=DMAX1(HALF*DELTA,DNORM)
          ELSE
!              DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
              DELTA=DMIN1(DMAX1(2._DP*DELTA,4._DP*DNORM),1.d10)
          END IF
          IF (DELTA <= 1.5_DP*RHO) DELTA=RHO

!     Recalculate KNEW and DENOM if the new F is less than FOPT.

          IF (F < FOPT) THEN
              KSAV=KNEW
              DENSAV=DENOM
              DELSQ=DELTA*DELTA
              SCADEN=ZERO
              BIGLSQ=ZERO
              KNEW=0
              DO 460 K=1,NPT
              HDIAG=ZERO
              DO 440 JJ=1,NPTM
  440         HDIAG=HDIAG+ZMAT(K,JJ)**2
              DEN=BETA*HDIAG+VLAG(K)**2
              DISTSQ=ZERO
              DO 450 J=1,N
  450         DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
              TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
              IF (TEMP*DEN > SCADEN) THEN
                  SCADEN=TEMP*DEN
                  KNEW=K
                  DENOM=DEN
              END IF
  460         BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
              IF (SCADEN <= HALF*BIGLSQ) THEN
                  KNEW=KSAV
                  DENOM=DENSAV
              END IF
          END IF
      END IF

!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.

      CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
      model_update = .true.
      IH=0
      do 465 m1=1,mv
        PQVold(m1) = PQV(m1,KNEW)
        PQV(m1,KNEW)=ZERO
  465 continue
      DO 470 I=1,N
      do 468 m1=1,mv
        v_temp(m1)=PQVold(m1)*XPT(KNEW,I)
  468 continue
      DO 470 J=1,I
      IH=IH+1
      do 470 m1=1,mv
        HQV(m1,IH)=HQV(m1,IH)+v_temp(m1)*XPT(KNEW,J)
  470 continue
      DO 480 JJ=1,NPTM
      do 475 m1=1,mv
         v_temp(m1)=v_diff(m1)*ZMAT(KNEW,JJ)
  475 continue
      DO 480 K=1,NPT
      do 480 m1=1,mv
         PQV(m1,K)=PQV(m1,K)+v_temp(m1)*ZMAT(K,JJ)
  480 continue

!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.

      FVAL(KNEW)=F
      do 485 m1=1,mv
         FVAL_V(m1,KNEW)=v_err(m1)
  485 continue
      DO 490 I=1,N
      XPT(KNEW,I)=XNEW(I)
  490 W(I)=BMAT(KNEW,I)
      DO 520 K=1,NPT
      SUMA=ZERO
      DO 500 JJ=1,NPTM
  500 SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
      SUMB=ZERO
      DO 510 J=1,N
  510 SUMB=SUMB+XPT(K,J)*XOPT(J)
      TEMP=SUMA*SUMB
      DO 520 I=1,N
  520 W(I)=W(I)+TEMP*XPT(K,I)
      DO 530 I=1,N
      do 530 m1=1,mv
         GQV(m1,I)=GQV(m1,I)+v_diff(m1)*W(I)
  530 continue

!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.

      IF (F < FOPT) THEN
          opt_update = .true.
          model_update = .true.
          KOPT=KNEW
          XOPTSQ=ZERO
          IH=0
          DO 540 J=1,N
          XOPT(J)=XNEW(J)
          XOPTSQ=XOPTSQ+XOPT(J)**2
          DO 540 I=1,J
          IH=IH+1
          IF (I < J) then
             do 535 m1=1,mv
  535        GQV(m1,J)=GQV(m1,J)+HQV(m1,IH)*D(I)
          endif
          do 540 m1=1,mv
             GQV(m1,I)=GQV(m1,I)+HQV(m1,IH)*D(J)
  540     continue
          DO 560 K=1,NPT
          TEMP=ZERO
          DO 550 J=1,N
  550     TEMP=TEMP+XPT(K,J)*D(J)
          do 555 m1=1,mv
  555     v_temp(m1)=PQV(m1,K)*TEMP
          DO 560 I=1,N
          do 560 m1=1,mv
             GQV(m1,I)=GQV(m1,I)+v_temp(m1)*XPT(K,I)
  560     continue
      END IF

!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.

      IF (NTRITS > 0) THEN
        do 568 k=1,NPT
          t=ZERO
          do 565 j=1,N
  565     t=t+XPT(k,j)*XOPT(j)
  568   v_sum(k) = t

        do 645 m1=1,mv
          DO 570 K=1,NPT
          VLAG(K)=FVAL_V(m1,K)-FVAL_V(m1,KOPT)
  570     W(K)=ZERO
          DO 590 J=1,NPTM
          SUM=ZERO
          DO 580 K=1,NPT
  580     SUM=SUM+ZMAT(K,J)*VLAG(K)
          DO 590 K=1,NPT
  590     W(K)=W(K)+SUM*ZMAT(K,J)
          DO 610 K=1,NPT
          W(K+NPT)=W(K)
  610     W(K)=v_sum(K)*W(K)
          GQSQ=ZERO
          GISQ=ZERO
          DO 630 I=1,N
          SUM=ZERO
          DO 620 K=1,NPT
  620     SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
          IF (dabs(XOPT(I)-SL(I)) <= eps) THEN
              GQSQ=GQSQ+DMIN1(ZERO,GQV(m1,I))**2
              GISQ=GISQ+DMIN1(ZERO,SUM)**2
          ELSE IF (dabs(XOPT(I)-SU(I)) <= eps) THEN
              GQSQ=GQSQ+DMAX1(ZERO,GQV(m1,I))**2
              GISQ=GISQ+DMAX1(ZERO,SUM)**2
          ELSE
              GQSQ=GQSQ+GQV(m1,I)**2
              GISQ=GISQ+SUM*SUM
          END IF
  630     VLAG(NPT+I)=SUM

!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.

          v_itest(m1)=v_itest(m1)+1
          if (GQSQ < TEN*GISQ) v_itest(m1)=0
          if (v_itest(m1) >= 3) then
              model_update = .true.
              do 640 i=1,MAX0(NPT,NH)
              if (i < n) GQV(m1,i)=VLAG(NPT+i)
              IF (i <= NPT) PQV(m1,i)=W(NPT+i)
              IF (i <= NH) HQV(m1,i)=ZERO
              v_itest(m1)=0
  640         CONTINUE
          endif
  645   continue
      END IF

!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.

      IF (NTRITS .EQ. 0) GOTO 60
      IF (F <= FOPT+TENTH*VQUAD1) GOTO 60

!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.

      DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
  650 KNEW=0
      DO 670 K=1,NPT
      SUM=ZERO
      DO 660 J=1,N
  660 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM > DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  670 CONTINUE

!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.

      IF (KNEW > 0) THEN
          DIST=DSQRT(DISTSQ)
          IF (NTRITS .EQ. -1) THEN
              DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
              IF (DELTA <= 1.5_DP*RHO) DELTA=RHO
          END IF
          NTRITS=0
          ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
          DSQ=ADELT*ADELT
          GOTO 90
      END IF
      IF (NTRITS .EQ. -1) GOTO 680
      IF (RATIO > ZERO) GOTO 60
      IF (DMAX1(DELTA,DNORM) > RHO) GOTO 60

!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.

  680 IF (RHO > RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO <= 16.0_DP) THEN
              RHO=RHOEND
          ELSE IF (RATIO <= 250.0_DP) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT >= 2) THEN
              IF (IPRINT >= 3) PRINT 690
  690         FORMAT (5X)
              PRINT 700, RHO,NF
  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',&
               ' function values =',I6)
              PRINT 710, FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
  710         FORMAT (4X,'Least value of F =',1PD23.15,9X,&
               'The corresponding X is:'/(2X,5D15.6))
          END IF
          NTRITS=0
          NFSAV=NF
          GOTO 60
      END IF
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
      IF (NTRITS .EQ. -1) GOTO 360
  720 IF (FVAL(KOPT) <= FSAVE) THEN
          DO 730 I=1,N
          X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
          IF (dabs(XOPT(I)-SL(I)) <= eps) X(I)=XL(I)
          IF (dabs(XOPT(I)-SU(I)) <= eps) X(I)=XU(I)
  730     CONTINUE
          F=FVAL(KOPT)
          do 735 m1=1,mv
  735     v_err(m1)=FVAL_V(m1,KOPT)
      END IF
      IF (IPRINT >= 1) THEN
          PRINT 740, NF
  740     FORMAT (/4X,'At the return from BOBYQA_H',5X,&
           'Number of function values =',I6)
          PRINT 710, F,(X(I),I=1,N)
      END IF
      RETURN
      END subroutine

      subroutine f_value(mv,v_err,F)
      integer mv
      double precision v_err(*), F
      integer m1

      F=0._DP
      do 5 m1=1,mv
    5 F = F + v_err(m1)**2

      return
      end subroutine


      SUBROUTINE ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,&
       KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL,W)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),&
       SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)

!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XNEW will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XALT also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all components of GLAG (see below) are zero.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to zero if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.

!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.

      HALF=0.5_DP
      ONE=1.0_DP
      ZERO=0.0_DP
      CONST=ONE+DSQRT(2.0_DP)
      DO 10 K=1,NPT
   10 HCOL(K)=ZERO
      DO 20 J=1,NPT-N-1
      TEMP=ZMAT(KNEW,J)
      DO 20 K=1,NPT
   20 HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
      ALPHA=HCOL(KNEW)
      HA=HALF*ALPHA

!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.

      DO 30 I=1,N
   30 GLAG(I)=BMAT(KNEW,I)
      DO 50 K=1,NPT
      TEMP=ZERO
      DO 40 J=1,N
   40 TEMP=TEMP+XPT(K,J)*XOPT(J)
      TEMP=HCOL(K)*TEMP
      DO 50 I=1,N
   50 GLAG(I)=GLAG(I)+TEMP*XPT(K,I)

!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.

      PRESAV=ZERO
      DO 80 K=1,NPT
      IF (K .EQ. KOPT) GOTO 80
      DDERIV=ZERO
      DISTSQ=ZERO
      DO 60 I=1,N
      TEMP=XPT(K,I)-XOPT(I)
      DDERIV=DDERIV+GLAG(I)*TEMP
   60 DISTSQ=DISTSQ+TEMP*TEMP
      SUBD=ADELT/DSQRT(DISTSQ)
      SLBD=-SUBD
      ILBD=0
      IUBD=0
      SUMIN=DMIN1(ONE,SUBD)

!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.

      DO 70 I=1,N
      TEMP=XPT(K,I)-XOPT(I)
      IF (TEMP > ZERO) THEN
          IF (SLBD*TEMP < SL(I)-XOPT(I)) THEN
              SLBD=(SL(I)-XOPT(I))/TEMP
              ILBD=-I
          END IF
          IF (SUBD*TEMP > SU(I)-XOPT(I)) THEN
              SUBD=DMAX1(SUMIN,(SU(I)-XOPT(I))/TEMP)
              IUBD=I
          END IF
      ELSE IF (TEMP < ZERO) THEN
          IF (SLBD*TEMP > SU(I)-XOPT(I)) THEN
              SLBD=(SU(I)-XOPT(I))/TEMP
              ILBD=I
          END IF
          IF (SUBD*TEMP < SL(I)-XOPT(I)) THEN
              SUBD=DMAX1(SUMIN,(SL(I)-XOPT(I))/TEMP)
              IUBD=-I
          END IF
      END IF
   70 CONTINUE

!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.

      IF (K .EQ. KNEW) THEN
          DIFF=DDERIV-ONE
          STEP=SLBD
          VLAG=SLBD*(DDERIV-SLBD*DIFF)
          ISBD=ILBD
          TEMP=SUBD*(DDERIV-SUBD*DIFF)
          IF (DABS(TEMP) > DABS(VLAG)) THEN
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          END IF
          TEMPD=HALF*DDERIV
          TEMPA=TEMPD-DIFF*SLBD
          TEMPB=TEMPD-DIFF*SUBD
          IF (TEMPA*TEMPB < ZERO) THEN
              TEMP=TEMPD*TEMPD/DIFF
              IF (DABS(TEMP) > DABS(VLAG)) THEN
                  STEP=TEMPD/DIFF
                  VLAG=TEMP
                  ISBD=0
              END IF
          END IF

!     Search along each of the other lines through XOPT and another point.

      ELSE
          STEP=SLBD
          VLAG=SLBD*(ONE-SLBD)
          ISBD=ILBD
          TEMP=SUBD*(ONE-SUBD)
          IF (DABS(TEMP) > DABS(VLAG)) THEN
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          END IF
          IF (SUBD > HALF) THEN
              IF (DABS(VLAG) < 0.25_DP) THEN
                  STEP=HALF
                  VLAG=0.25_DP
                  ISBD=0
              END IF
          END IF
          VLAG=VLAG*DDERIV
      END IF

!     Calculate PREDSQ for the current line search and maintain PRESAV.

      TEMP=STEP*(ONE-STEP)*DISTSQ
      PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
      IF (PREDSQ > PRESAV) THEN
          PRESAV=PREDSQ
          KSAV=K
          STPSAV=STEP
          IBDSAV=ISBD
      END IF
   80 CONTINUE

!     Construct XNEW in a way that satisfies the bound constraints exactly.

      DO 90 I=1,N
      TEMP=XOPT(I)+STPSAV*(XPT(KSAV,I)-XOPT(I))
   90 XNEW(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
      IF (IBDSAV < 0) XNEW(-IBDSAV)=SL(-IBDSAV)
      IF (IBDSAV > 0) XNEW(IBDSAV)=SU(IBDSAV)

!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.

      BIGSTP=ADELT+ADELT
      IFLAG=0
  100 WFIXSQ=ZERO
      GGFREE=ZERO
      DO 110 I=1,N
      W(I)=ZERO
      TEMPA=DMIN1(XOPT(I)-SL(I),GLAG(I))
      TEMPB=DMAX1(XOPT(I)-SU(I),GLAG(I))
      IF (TEMPA > ZERO .OR. TEMPB < ZERO) THEN
          W(I)=BIGSTP
          GGFREE=GGFREE+GLAG(I)**2
      END IF
  110 CONTINUE
      IF (dabs(GGFREE) <= eps) THEN
          CAUCHY=ZERO
          GOTO 200
      END IF

!     Investigate whether more components of W can be fixed.

  120 TEMP=ADELT*ADELT-WFIXSQ
      IF (TEMP > ZERO) THEN
          WSQSAV=WFIXSQ
          STEP=DSQRT(TEMP/GGFREE)
          GGFREE=ZERO
          DO 130 I=1,N
          IF (dabs(W(I)-BIGSTP) <= eps) THEN
              TEMP=XOPT(I)-STEP*GLAG(I)
              IF (TEMP <= SL(I)) THEN
                  W(I)=SL(I)-XOPT(I)
                  WFIXSQ=WFIXSQ+W(I)**2
              ELSE IF (TEMP >= SU(I)) THEN
                  W(I)=SU(I)-XOPT(I)
                  WFIXSQ=WFIXSQ+W(I)**2
              ELSE
                  GGFREE=GGFREE+GLAG(I)**2
              END IF
          END IF
  130     CONTINUE
          IF (WFIXSQ > WSQSAV .AND. GGFREE > ZERO) GOTO 120
      END IF

!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.

      GW=ZERO
      DO 140 I=1,N
      IF (dabs(W(I)-BIGSTP) <= eps) THEN
          W(I)=-STEP*GLAG(I)
          XALT(I)=DMAX1(SL(I),DMIN1(SU(I),XOPT(I)+W(I)))
      ELSE IF (dabs(W(I)) <= eps) THEN
          XALT(I)=XOPT(I)
      ELSE IF (GLAG(I) > ZERO) THEN
          XALT(I)=SL(I)
      ELSE
          XALT(I)=SU(I)
      END IF
  140 GW=GW+GLAG(I)*W(I)

!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.

      CURV=ZERO
      DO 160 K=1,NPT
      TEMP=ZERO
      DO 150 J=1,N
  150 TEMP=TEMP+XPT(K,J)*W(J)
  160 CURV=CURV+HCOL(K)*TEMP*TEMP
      IF (IFLAG .EQ. 1) CURV=-CURV
      IF (CURV > -GW .AND. CURV < -CONST*GW) THEN
          SCALE=-GW/CURV
          DO 170 I=1,N
          TEMP=XOPT(I)+SCALE*W(I)
  170     XALT(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
          CAUCHY=(HALF*GW*SCALE)**2
      ELSE
          CAUCHY=(GW+HALF*CURV)**2
      END IF

!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.

      IF (IFLAG .EQ. 0) THEN
          DO 180 I=1,N
          GLAG(I)=-GLAG(I)
  180     W(N+I)=XALT(I)
          CSAVE=CAUCHY
          IFLAG=1
          GOTO 100
      END IF
      IF (CSAVE > CAUCHY) THEN
          DO 190 I=1,N
  190     XALT(I)=W(N+I)
          CAUCHY=CSAVE
      END IF
  200 RETURN
      END subroutine




!   Important Notice:
!   This PRELIM_H are modifications and based on the subroutine PRELIM in the software
!   BOBYQA, authored by M. J. D. Powell.

      SUBROUTINE PRELIM_H (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,&
       XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,&
       mv,HQV,GQV,PQV,FVAL_V,v_err,v_beg,FBEG)
      use OBJECTIVE, only: dfovec, getPenalty

      IMPLICIT double precision (A-H,O-Z)
      parameter (nmax = 100, mmax=3000)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),&
       HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
      dimension HQV(mmax,*),GQV(mmax,*),PQV(mmax,*),FVAL_V(mmax,*)
      dimension v_err(*),v_beg(*)

!external dfovec

!     Set some constants.

      HALF=0.5_DP
      ONE=1.0_DP
      TWO=2.0_DP
      ZERO=0.0_DP
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      NP=N+1

!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.

      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 IH=1,(N*NP)/2
      do 25 m1=1, mv
   25 HQV(m1,IH)=zero
   30 HQ(IH)=ZERO
      DO 40 K=1,NPT
      do 35 m1=1, mv
   35 PQV(m1,K)=zero
      PQ(K)=ZERO
      DO 40 J=1,NPT-NP
   40 ZMAT(K,J)=ZERO

!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).

      NF=0
   50 NFM=NF
      NFX=NF-N
      NF=NF+1
      IF (NFM <= 2*N) THEN
          IF (NFM >= 1 .AND. NFM <= N) THEN
              STEPA=RHOBEG
              IF (dabs(SU(NFM)) <= eps) STEPA=-STEPA
              XPT(NF,NFM)=STEPA
          ELSE IF (NFM > N) THEN
              STEPA=XPT(NF-N,NFX)
              STEPB=-RHOBEG
              IF (dabs(SL(NFX)) <= eps) STEPB=DMIN1(TWO*RHOBEG,SU(NFX))
              IF (dabs(SU(NFX)) <= eps) STEPB=DMAX1(-TWO*RHOBEG,SL(NFX))
              XPT(NF,NFX)=STEPB
          END IF
      END IF
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
      DO 60 J=1,N
      X(J)=DMIN1(DMAX1(XL(J),XBASE(J)+XPT(NF,J)),XU(J))
      IF (dabs(XPT(NF,J)-SL(J)) <= eps) X(J)=XL(J)
      IF (dabs(XPT(NF,J)-SU(J)) <= eps) X(J)=XU(J)
   60 CONTINUE
!
!     dfovec(n, mv, x, v_err) provides the values of the vector function v_err(x): R^n \to R^{mv}.
!     Here: n, mv, x \in R^n are input, v_err \in R^{mv} are output.
      call dfovec(n, mv, x, v_err)
!
!     f_value(mv,v_err,F) provides the value of the sum of the squres of the components of v_err(x)
!     i.e. F = sum_{i=1}^{mv} v_err_i (x)^2
      call f_value(mv,v_err,F)

	  !LS: Add penalty to result
	  F = F + getPenalty(x)

      IF (IPRINT .EQ. 3) THEN
          PRINT 70, NF,F,(X(I),I=1,N)
   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,&
            '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(NF)=F
      do 72 m1=1,mv
   72 FVAL_V(m1,NF)=v_err(m1)
      IF (NF .EQ. 1) THEN
          FBEG=F
          KOPT=1
          do 75 m1=1,mv
   75     v_beg(m1) = v_err(m1)
      ELSE IF (F < FVAL(KOPT)) THEN
          KOPT=NF
      END IF
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1.
!
      IF (NF <= 2*N+1) THEN
          IF (NF >= 2 .AND. NF <= N+1) THEN
              do 78 m1=1,mv
   78         GQV(m1,NFM)=(v_err(m1)-v_beg(m1))/STEPA
              IF (NPT < NF+N) THEN
                  BMAT(1,NFM)=-ONE/STEPA
                  BMAT(NF,NFM)=ONE/STEPA
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NF >= N+2) THEN
              IH=(NFX*(NFX+1))/2
              DIFF=STEPB-STEPA
              do 79 m1=1, mv
                TEMP=(v_err(m1)-v_beg(m1))/STEPB
                HQV(m1,IH)=TWO*(TEMP-GQV(m1,NFX))/DIFF
                GQV(m1,NFX)=(GQV(m1,NFX)*STEPB-TEMP*STEPA)/DIFF
   79         continue
              IF (STEPA*STEPB < ZERO) THEN
                  IF (F < FVAL(NF-N)) THEN
                      FVAL(NF)=FVAL(NF-N)
                      FVAL(NF-N)=F
                      do 80 m1=1, mv
                         FVAL_V(m1,NF)=FVAL_V(m1,NF-N)
                         FVAL_V(m1,NF-N)=v_err(m1)
   80                 continue
                      IF (KOPT .EQ. NF) KOPT=NF-N
                      XPT(NF-N,NFX)=STEPB
                      XPT(NF,NFX)=STEPA
                  END IF
              END IF
              BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
              BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
              BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
              ZMAT(1,NFX)=DSQRT(TWO)/(STEPA*STEPB)
              ZMAT(NF,NFX)=DSQRT(HALF)/RHOSQ
              ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
          END IF
      END IF
      IF (NF < NPT .AND. NF < MAXFUN) GOTO 50
      RETURN
      END subroutine



!   Important Notice:
!   This RESCUE_H are modifications and based on the subroutine RESCUE in the software
!   BOBYQA, authored by M. J. D. Powell.

      SUBROUTINE RESCUE_H (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,&
       FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,&
       KOPT,VLAG,PTSAUX,PTSID,W,&
       mv,HQV,GQV,PQV,FVAL_V,WV,v_err,v_sumpq,v_vquad,v_diff,v_temp)
      use OBJECTIVE, only: dfovec, getPenalty

      IMPLICIT double precision (A-H,O-Z)
      parameter (nmax = 100, mmax=3000, nptmax=2*nmax+1)
      DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),&
       GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),&
       VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
      dimension HQV(mmax,*),GQV(mmax,*),PQV(mmax,*),FVAL_V(mmax,*),&
       WV(mmax,*),v_err(*),v_sumpq(*),v_vquad(*),v_diff(*),v_temp(*)
      dimension v_fbase(mmax)

!   external dfovec

!     Set some constants.

      HALF=0.5_DP
      ONE=1.0_DP
      ZERO=0.0_DP
      NP=N+1
      SFRAC=HALF/DFLOAT(NP)
      NPTM=NPT-NP

!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to zero. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.

!      SUMPQ=ZERO
      do 5 m1=1,mv
    5 v_sumpq(m1)=ZERO
      WINC=ZERO
      DO 20 K=1,NPT
      DISTSQ=ZERO
      DO 10 J=1,N
      XPT(K,J)=XPT(K,J)-XOPT(J)
   10 DISTSQ=DISTSQ+XPT(K,J)**2
      do 15 m1=1,mv
   15 v_sumpq(m1)= v_sumpq(m1)+PQV(m1,k)
      W(NDIM+K)=DISTSQ
      WINC=DMAX1(WINC,DISTSQ)
      DO 20 J=1,NPTM
   20 ZMAT(K,J)=ZERO

!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.

      IH=0
      DO 40 J=1,N
      do 25 m1=1,mv
   25 WV(m1,J)=HALF*v_sumpq(m1)*XOPT(J)
      DO 30 K=1,NPT
      do 30 m1=1,mv
   30 WV(m1,J)=WV(m1,J)+PQV(m1,K)*XPT(K,J)
      DO 40 I=1,J
      IH=IH+1
      do 40 m1=1,mv
         HQV(m1,IH)=HQV(m1,IH)+WV(m1,I)*XOPT(J) &
     &              +WV(m1,J)*XOPT(I)
   40 continue

!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.

      DO 50 J=1,N
      XBASE(J)=XBASE(J)+XOPT(J)
      SL(J)=SL(J)-XOPT(J)
      SU(J)=SU(J)-XOPT(J)
      XOPT(J)=ZERO
      PTSAUX(1,J)=DMIN1(DELTA,SU(J))
      PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
      IF (PTSAUX(1,J)+PTSAUX(2,J) < ZERO) THEN
          TEMP=PTSAUX(1,J)
          PTSAUX(1,J)=PTSAUX(2,J)
          PTSAUX(2,J)=TEMP
      END IF
      IF (DABS(PTSAUX(2,J)) < HALF*DABS(PTSAUX(1,J))) THEN
          PTSAUX(2,J)=HALF*PTSAUX(1,J)
      END IF
      DO 50 I=1,NDIM
   50 BMAT(I,J)=ZERO
      FBASE=FVAL(KOPT)
      do 55 m1=1,mv
   55 v_fbase(m1)=FVAL_V(m1,KOPT)

!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.

      PTSID(1)=SFRAC
      DO 60 J=1,N
      JP=J+1
      JPN=JP+N
      PTSID(JP)=DFLOAT(J)+SFRAC
      IF (JPN <= NPT) THEN
          PTSID(JPN)=DFLOAT(J)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
          BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
          BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
          BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
          ZMAT(1,J)=DSQRT(2.0_DP)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
          ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
          ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
      ELSE
          BMAT(1,J)=-ONE/PTSAUX(1,J)
          BMAT(JP,J)=ONE/PTSAUX(1,J)
          BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
      END IF
   60 CONTINUE

!     Set any remaining identifiers with their nonzero elements of ZMAT.

      IF (NPT >= N+NP) THEN
          DO 70 K=2*NP,NPT
          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
          IP=K-NP-IW*N
          IQ=IP+IW
          IF (IQ > N) IQ=IQ-N
          PTSID(K)=DFLOAT(IP)+DFLOAT(IQ)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
          ZMAT(1,K-NP)=TEMP
          ZMAT(IP+1,K-NP)=-TEMP
          ZMAT(IQ+1,K-NP)=-TEMP
   70     ZMAT(K,K-NP)=TEMP
      END IF
      NREM=NPT
      KOLD=1
      KNEW=KOPT

!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).

   80 DO 90 J=1,N
      TEMP=BMAT(KOLD,J)
      BMAT(KOLD,J)=BMAT(KNEW,J)
   90 BMAT(KNEW,J)=TEMP
      DO 100 J=1,NPTM
      TEMP=ZMAT(KOLD,J)
      ZMAT(KOLD,J)=ZMAT(KNEW,J)
  100 ZMAT(KNEW,J)=TEMP
      PTSID(KOLD)=PTSID(KNEW)
      PTSID(KNEW)=ZERO
      W(NDIM+KNEW)=ZERO
      NREM=NREM-1
      IF (KNEW .NE. KOPT) THEN
          TEMP=VLAG(KOLD)
          VLAG(KOLD)=VLAG(KNEW)
          VLAG(KNEW)=TEMP

!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.

          CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
          IF (NREM .EQ. 0) GOTO 350
          DO 110 K=1,NPT
  110     W(NDIM+K)=DABS(W(NDIM+K))
      END IF

!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.

  120 DSQMIN=ZERO
      DO 130 K=1,NPT
      IF (W(NDIM+K) > ZERO) THEN
          IF (dabs(DSQMIN) <= eps .OR. W(NDIM+K) < DSQMIN) THEN
              KNEW=K
              DSQMIN=W(NDIM+K)
          END IF
      END IF
  130 CONTINUE
      IF (dabs(DSQMIN) <= eps) GOTO 260

!     Form the W-vector of the chosen original interpolation point.

      DO 140 J=1,N
  140 W(NPT+J)=XPT(KNEW,J)
      DO 160 K=1,NPT
      SUM=ZERO
      IF (K .EQ. KOPT) THEN
          CONTINUE
      ELSE IF (dabs(PTSID(K)) <= eps) THEN
          DO 150 J=1,N
  150     SUM=SUM+W(NPT+J)*XPT(K,J)
      ELSE
          IP=PTSID(K)
          IF (IP > 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IF (IQ > 0) THEN
              IW=1
              IF (IP .EQ. 0) IW=2
              SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
          END IF
      END IF
  160 W(K)=HALF*SUM*SUM

!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.

      DO 180 K=1,NPT
      SUM=ZERO
      DO 170 J=1,N
  170 SUM=SUM+BMAT(K,J)*W(NPT+J)
  180 VLAG(K)=SUM
      BETA=ZERO
      DO 200 J=1,NPTM
      SUM=ZERO
      DO 190 K=1,NPT
  190 SUM=SUM+ZMAT(K,J)*W(K)
      BETA=BETA-SUM*SUM
      DO 200 K=1,NPT
  200 VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
      BSUM=ZERO
      DISTSQ=ZERO
      DO 230 J=1,N
      SUM=ZERO
      DO 210 K=1,NPT
  210 SUM=SUM+BMAT(K,J)*W(K)
      JP=J+NPT
      BSUM=BSUM+SUM*W(JP)
      DO 220 IP=NPT+1,NDIM
  220 SUM=SUM+BMAT(IP,J)*W(IP)
      BSUM=BSUM+SUM*W(JP)
      VLAG(JP)=SUM
  230 DISTSQ=DISTSQ+XPT(KNEW,J)**2
      BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE

!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.

      DENOM=ZERO
      VLMXSQ=ZERO
      DO 250 K=1,NPT
      IF (dabs(PTSID(K)) > eps) THEN
          HDIAG=ZERO
          DO 240 J=1,NPTM
  240     HDIAG=HDIAG+ZMAT(K,J)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          IF (DEN > DENOM) THEN
              KOLD=K
              DENOM=DEN
          END IF
      END IF
  250 VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
      IF (DENOM <= 1.0D-2*VLMXSQ) THEN
          W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
          GOTO 120
      END IF
      GOTO 80

!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.

  260 DO 340 KPT=1,NPT
      IF (dabs(PTSID(KPT)) <= eps) GOTO 340
      IF (NF >= MAXFUN) THEN
          NF=-1
          GOTO 350
      END IF
      IH=0
      DO 270 J=1,N
      W(J)=XPT(KPT,J)
      XPT(KPT,J)=ZERO
      do 265 m1=1,mv
  265 v_temp(m1)=PQV(m1,KPT)*W(J)
      DO 270 I=1,J
      IH=IH+1
      do 270 m1=1,mv
  270 HQV(m1,IH)=HQV(m1,IH)+v_temp(m1)*W(I)
      do 275 m1=1,mv
  275 PQV(m1,KPT)=ZERO
      IP=PTSID(KPT)
      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
      IF (IP > 0) THEN
          XP=PTSAUX(1,IP)
          XPT(KPT,IP)=XP
      END IF
      IF (IQ > 0) THEN
          XQ=PTSAUX(1,IQ)
          IF (IP .EQ. 0) XQ=PTSAUX(2,IQ)
          XPT(KPT,IQ)=XQ
      END IF

!     Set VQUAD1 to the value of the current model at the new point.

      do 276 m1=1,mv
  276 v_vquad(m1)=v_fbase(m1)
      IF (IP > 0) THEN
          IHP=(IP+IP*IP)/2
          do 277 m1=1,mv
  277     v_vquad(m1)=v_vquad(m1)+XP*(GQV(m1,IP)+HALF*XP*HQV(m1,IHP))
      END IF
      IF (IQ > 0) THEN
          IHQ=(IQ+IQ*IQ)/2
          do 278 m1=1,mv
  278     v_vquad(m1)=v_vquad(m1)+XQ*(GQV(m1,IQ)+HALF*XQ*HQV(m1,IHQ))
          IF (IP > 0) THEN
              IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
              do 279 m1=1,mv
  279         v_vquad(m1)=v_vquad(m1)+XP*XQ*HQV(m1,IW)
          END IF
      END IF
      DO 280 K=1,NPT
      TEMP=ZERO
      IF (IP > 0) TEMP=TEMP+XP*XPT(K,IP)
      IF (IQ > 0) TEMP=TEMP+XQ*XPT(K,IQ)
      do 280 m1=1,mv
  280 v_vquad(m1)=v_vquad(m1)+HALF*PQV(m1,K)*TEMP*TEMP

!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.

      DO 290 I=1,N
      W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
      IF (dabs(XPT(KPT,I)-SL(I)) <= eps) W(I)=XL(I)
      IF (dabs(XPT(KPT,I)-SU(I)) <= eps) W(I)=XU(I)
  290 CONTINUE
      NF=NF+1
!
!     dfovec(n, mv, x, v_err) provides the values of the vector function v_err(x): R^n \to R^{mv}.
!     Here: n, mv, x \in R^n are input, v_err \in R^{mv} are output.
      call dfovec(n, mv, W, v_err)
!
!     f_value(mv,v_err,F) provides the value of the sum of the squres of the components of v_err(x)
!     i.e. F = sum_{i=1}^{mv} v_err_i (x)^2
      call f_value(mv,v_err,F)

	  !LS: Add penalty to result
	  F = F + getPenalty(W)

      IF (IPRINT .EQ. 3) THEN
          PRINT 300, NF,F,(W(I),I=1,N)
  300     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,&
         '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(KPT)=F
      do 305 m1=1,mv
  305 FVAL_V(m1,KPT)=v_err(m1)
      IF (F < FVAL(KOPT)) KOPT=KPT
      do 308 m1=1,mv
  308 v_diff(m1)=v_err(m1)-v_vquad(m1)
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
      DO 310 I=1,N
      do 310 m1=1,mv
  310 GQV(m1,I)=GQV(m1,I)+v_diff(m1)*BMAT(KPT,I)
      DO 330 K=1,NPT
      SUM=ZERO
      DO 320 J=1,NPTM
  320 SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
      do 322 m1=1,mv
  322 v_temp(m1)=v_diff(m1)*SUM
      IF (dabs(PTSID(K)) <= eps) THEN
          do 324 m1=1,mv
  324     PQV(m1,K)=PQV(m1,K)*v_temp(m1)
      ELSE
          IP=PTSID(K)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IHQ=(IQ*IQ+IQ)/2
          IF (IP .EQ. 0) THEN
              do 326 m1=1,mv
  326         HQV(m1,IHQ)=HQV(m1,IHQ)+v_temp(m1)*PTSAUX(2,IQ)**2
          ELSE
              IHP=(IP*IP+IP)/2
              HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
              IF (IQ > 0) THEN
                IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                do 328 m1=1,mv
                HQV(m1,IHQ)=HQV(m1,IHQ)+v_temp(m1)*PTSAUX(1,IQ)**2
  328           HQV(m1,IW)=HQV(m1,IW) &
                         +v_temp(m1)*PTSAUX(1,IP)*PTSAUX(1,IQ)
              END IF
          END IF
      END IF
  330 CONTINUE
      PTSID(KPT)=ZERO
  340 CONTINUE
  350 RETURN
      END subroutine




!   Important Notice:
!   This TRSBOX_H are modifications and based on the subroutine TRSBOX in the software
!   BOBYQA, authored by M. J. D. Powell.

      SUBROUTINE TRSBOX_H (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,&
       XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN,&
       mv,HQV,GQV,PQV,FVAL_V,KOPT,HD1,v_temp,vquad,&
       model_update,opt_update)
      IMPLICIT double precision (A-H,O-Z)
      parameter (nmax = 100, mmax=3000)
      DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),&
       XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
      dimension HQV(mmax,*),GQV(mmax,*),PQV(mmax,*),FVAL_V(mmax,*)
      dimension HD1(*),v_temp(*)
      logical model_update, opt_update
      integer cases

      if (n>nmax) then
        print *, "in trsbox_h.f increase the dimension &
     &            nmax to be at least", n
        stop
      endif
      if (mv>mmax) then
        print *, "in trsbox_h.f increase the dimension &
     &            mmax to be at least", mv
        stop
      endif

      if ((.not.model_update).and.(.not.opt_update)) go to 8

        model_update = .false.
        opt_update = .false.
        voptmax = 0._DP
        do 5 m1=1,mv
        t=FVAL_V(m1,KOPT)
        voptmax = dmax1(voptmax, dabs(t))
    5   v_temp(m1)=t

! Use the gradient at xopt to formulate \sum_i (2*f_i \nabla f_i) = 2 J^t m(x_opt)

        gnorm2 = 0._DP
        do i=1,n
          GOPT(i) = 0._DP
          do m1=1,mv
            GOPT(i) = GOPT(i) + v_temp(m1)*GQV(m1,i)
          enddo
          GOPT(i) = 2._DP*GOPT(i)
          gnorm2 = gnorm2 + GOPT(i)**2
       enddo

! Calculate the explicite Hessian.

       f_base = 0._DP
       call f_value(mv,v_temp,f_base)
       if (gnorm2>=1._DP) then
! xopt is far awary from a stationary point
          cases = 1
       elseif (f_base<=dsqrt(gnorm2)) then
! xopt is close to a stationary point and zero residue
          cases = 2
       else
! xopt is close to a stationary point and nonzero residue
          cases = 3
       endif



       IH=0
       do j=1,n
         do i=1,j
           IH=IH+1
           if (cases.eq.1) then
             t1 = 0._DP
             do m1=1,mv
               t1 = t1+GQV(m1,i)*GQV(m1,j)
             enddo
             HQ(IH) = 2._DP*t1
           elseif (cases.eq.2) then
             t1 = 0._DP
             do m1=1,mv
               t1 = t1+GQV(m1,i)*GQV(m1,j)
             enddo
             HQ(IH) = 2._DP*t1
             if (i.eq.j) HQ(IH) = HQ(IH) + 1.d-3*voptmax
           else
             t1 = 0._DP
             do m1=1,mv
               t2 = 0._DP
               do k=1,NPT
                 t2 = t2 + XPT(k,i)*PQV(m1,k)*XPT(k,j)
               enddo
               t2 = t2 + HQV(m1,IH)
               t1 = t1+(GQV(m1,i)*GQV(m1,j)+v_temp(m1)*t2)
             enddo
             HQ(IH) = 2._DP*t1
           endif
         enddo
       enddo

!     A version of the truncated conjugate gradient is applied. If a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. If
!     the trust region boundary is reached, then further changes may be made
!     to D, each one being in the two dimensional space that is spanned
!     by the current D and the gradient of Q at XOPT+D, staying on the trust
!     region boundary. Termination occurs when the reduction in Q seems to
!     be close to the greatest reduction that can be achieved.

!     Set some constants.

    8 HALF=0.5_DP
      ONE=1.0_DP
      ONEMIN=-1.0_DP
      ZERO=0.0_DP

!     The sign of GOPT(I) gives the sign of the change to the I-th variable
!     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
!     or not to fix the I-th variable at one of its bounds initially, with
!     NACT being set to the number of fixed variables. D and GNEW are also
!     set for the first iteration. DELSQ is the upper bound on the sum of
!     squares of the free variables. QRED is the reduction in Q so far.

      ITERC=0
      NACT=0
      DO 10 I=1,N
      XBDI(I)=ZERO
      IF (XOPT(I) <= SL(I)) THEN
          IF (GOPT(I) >= ZERO) XBDI(I)=ONEMIN
      ELSE IF (XOPT(I) >= SU(I)) THEN
          IF (GOPT(I) <= ZERO) XBDI(I)=ONE
      END IF
      IF (dabs(XBDI(I)) > eps) NACT=NACT+1
      D(I)=ZERO
   10 GNEW(I)=GOPT(I)
      DELSQ=DELTA*DELTA
      QRED=ZERO
      CRVMIN=ONEMIN

!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.

   20 BETA=ZERO
   30 STEPSQ=ZERO
      DO 40 I=1,N
      IF (dabs(XBDI(I)) > eps) THEN
          S(I)=ZERO
      ELSE IF (dabs(BETA) <= eps) THEN
          S(I)=-GNEW(I)
      ELSE
          S(I)=BETA*S(I)-GNEW(I)
      END IF
   40 STEPSQ=STEPSQ+S(I)**2
      IF (dabs(STEPSQ) <= eps) GOTO 190
      IF (dabs(BETA) <= eps) THEN
          GREDSQ=STEPSQ
          ITERMAX=ITERC+N-NACT
      END IF
!      save original gredsq0
      if (ITERC.eq.0) gredsq0=GREDSQ
      if (GREDSQ<=dmin1(1.0D-6*gredsq0,1.d-18)) goto 190
      if (GREDSQ<= 1.0D-14*gredsq0) goto 190
      IF (GREDSQ*DELSQ <= dmin1(1.0D-6*QRED*QRED,1.d-18)) GOTO 190
      IF (GREDSQ*DELSQ <= 1.0D-14*QRED*QRED) GOTO 190

!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.

      GOTO 210
   50 RESID=DELSQ
      DS=ZERO
      SHS=ZERO
      DO 60 I=1,N
      IF (dabs(XBDI(I)) <= eps) THEN
          RESID=RESID-D(I)**2
          DS=DS+S(I)*D(I)
          SHS=SHS+S(I)*HS(I)
      END IF
   60 CONTINUE
      IF (RESID <= ZERO) GOTO 90
      TEMP=DSQRT(STEPSQ*RESID+DS*DS)
      IF (DS < ZERO) THEN
          BLEN=(TEMP-DS)/STEPSQ
      ELSE
          BLEN=RESID/(TEMP+DS)
      END IF
      STPLEN=BLEN
      IF (SHS > ZERO) THEN
          STPLEN=DMIN1(BLEN,GREDSQ/SHS)
      END IF
      if (STPLEN<=1.e-30) goto 190

!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.

      IACT=0
      DO 70 I=1,N
      IF (dabs(S(I)) > eps) THEN
          XSUM=XOPT(I)+D(I)
          IF (S(I) > ZERO) THEN
              TEMP=(SU(I)-XSUM)/S(I)
          ELSE
              TEMP=(SL(I)-XSUM)/S(I)
          END IF
          IF (TEMP < STPLEN) THEN
              STPLEN=TEMP
              IACT=I
          END IF
      END IF
   70 CONTINUE

!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.

      SDEC=ZERO
      IF (STPLEN > ZERO) THEN
          ITERC=ITERC+1
          TEMP=SHS/STEPSQ
          IF (IACT .EQ. 0 .AND. TEMP > ZERO) THEN
              CRVMIN=DMIN1(CRVMIN,TEMP)
              IF (dabs(CRVMIN-ONEMIN) <= eps) CRVMIN=TEMP
          END IF
          GGSAV=GREDSQ
          GREDSQ=ZERO
          DO 80 I=1,N
          GNEW(I)=GNEW(I)+STPLEN*HS(I)
          IF (dabs(XBDI(I)) <= eps) GREDSQ=GREDSQ+GNEW(I)**2
   80     D(I)=D(I)+STPLEN*S(I)
          SDEC=DMAX1(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
          QRED=QRED+SDEC
      END IF

!     Restart the conjugate gradient method if it has hit a new bound.

      IF (IACT > 0) THEN
          NACT=NACT+1
          XBDI(IACT)=ONE
          IF (S(IACT) < ZERO) XBDI(IACT)=ONEMIN
          DELSQ=DELSQ-D(IACT)**2
          IF (DELSQ <= ZERO) GOTO 90
          GOTO 20
      END IF

!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.

      IF (STPLEN < BLEN) THEN
          IF (ITERC .EQ. ITERMAX) GOTO 190
          If (SDEC <= 1.D-6*QRED) GOTO 190
          BETA=GREDSQ/GGSAV
          GOTO 30
      END IF
   90 CRVMIN=ZERO

!     Prepare for the alternative iteration by calculating some scalars
!     and by multiplying the reduced D by the second derivative matrix of
!     Q, where S holds the reduced D in the call of GGMULT.

  100 IF (NACT >= N-1) GOTO 190
      DREDSQ=ZERO
      DREDG=ZERO
      GREDSQ=ZERO
      DO 110 I=1,N
      IF (dabs(XBDI(I)) <= eps) THEN
          DREDSQ=DREDSQ+D(I)**2
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
          S(I)=D(I)
      ELSE
          S(I)=ZERO
      END IF
  110 CONTINUE
      ITCSAV=ITERC
      GOTO 210

!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.

  120 ITERC=ITERC+1
      TEMP=GREDSQ*DREDSQ-DREDG*DREDG
      IF (TEMP <= 1.0D-4*QRED*QRED) GOTO 190
      TEMP=DSQRT(TEMP)
      DO 130 I=1,N
      IF (dabs(XBDI(I)) <= eps) THEN
          S(I)=(DREDG*D(I)-DREDSQ*GNEW(I))/TEMP
      ELSE
          S(I)=ZERO
      END IF
  130 CONTINUE
      SREDG=-TEMP

!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.

      ANGBD=ONE
      IACT=0
      DO 140 I=1,N
      IF (dabs(XBDI(I)) <= eps) THEN
          TEMPA=XOPT(I)+D(I)-SL(I)
          TEMPB=SU(I)-XOPT(I)-D(I)
          IF (TEMPA <= ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONEMIN
              GOTO 100
          ELSE IF (TEMPB <= ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONE
              GOTO 100
          END IF
          RATIO=ONE
          SSQ=D(I)**2+S(I)**2
          TEMP=SSQ-(XOPT(I)-SL(I))**2
          IF (TEMP > ZERO) THEN
              TEMP=DSQRT(TEMP)-S(I)
              IF (ANGBD*TEMP > TEMPA) THEN
                  ANGBD=TEMPA/TEMP
                  IACT=I
                  XSAV=ONEMIN
              END IF
          END IF
          TEMP=SSQ-(SU(I)-XOPT(I))**2
          IF (TEMP > ZERO) THEN
              TEMP=DSQRT(TEMP)+S(I)
              IF (ANGBD*TEMP > TEMPB) THEN
                  ANGBD=TEMPB/TEMP
                  IACT=I
                  XSAV=ONE
              END IF
          END IF
      END IF
  140 CONTINUE
!
!     Calculate HHD and some curvatures for the alternative iteration.

      GOTO 210
  150 SHS=ZERO
      DHS=ZERO
      DHD=ZERO
      DO 160 I=1,N
      IF (dabs(XBDI(I)) <= eps) THEN
          SHS=SHS+S(I)*HS(I)
          DHS=DHS+D(I)*HS(I)
          DHD=DHD+D(I)*HRED(I)
      END IF
  160 CONTINUE
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
!     the alternative iteration.
!
      REDMAX=ZERO
      ISAV=0
      REDSAV=ZERO
      IU=17.0_DP*ANGBD+3.1_DP
      DO 170 I=1,IU
      ANGT=ANGBD*DFLOAT(I)/DFLOAT(IU)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (REDNEW > REDMAX) THEN
          REDMAX=REDNEW
          ISAV=I
          RDPREV=REDSAV
      ELSE IF (I .EQ. ISAV+1) THEN
          RDNEXT=REDNEW
      END IF
  170 REDSAV=REDNEW
!
!     Return if the redu!tion is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and cal!ulate SDE!.
!
      IF (ISAV .EQ. 0) GOTO 190
      IF (ISAV < IU) THEN
          TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
          ANGT=ANGBD*(DFLOAT(ISAV)+HALF*TEMP)/DFLOAT(IU)
      END IF
      CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (SDEC <= ZERO) GOTO 190
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
      DREDG=ZERO
      GREDSQ=ZERO
      DO 180 I=1,N
      GNEW(I)=GNEW(I)+(CTH-ONE)*HRED(I)+STH*HS(I)
      IF (dabs(XBDI(I)) <= eps) THEN
          D(I)=CTH*D(I)+STH*S(I)
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
      END IF
  180 HRED(I)=CTH*HRED(I)+STH*HS(I)
      QRED=QRED+SDEC
      IF (IACT > 0 .AND. ISAV .EQ. IU) THEN
          NACT=NACT+1
          XBDI(IACT)=XSAV
          GOTO 100
      END IF
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
      IF (SDEC > 0.01_DP*QRED) GOTO 120
  190 DSQ=ZERO
      DO 200 I=1,N
      XNEW(I)=DMAX1(DMIN1(XOPT(I)+D(I),SU(I)),SL(I))
      IF (dabs(XBDI(I)-ONEMIN) <= eps) XNEW(I)=SL(I)
      IF (dabs(XBDI(I)-ONE) <= eps) XNEW(I)=SU(I)
      D(I)=XNEW(I)-XOPT(I)
  200 DSQ=DSQ+D(I)**2

      do 202 i=1,n
  202 HD1(i) = 0._DP
      IH=0
      DO 204 J=1,N
      DO 204 I=1,J
      IH=IH+1
      IF (I < J) HD1(J)=HD1(J)+HQ(IH)*D(I)
  204 HD1(I)=HD1(I)+HQ(IH)*D(J)
      vquad= 0._DP
      do 163 i=1,n
  163 vquad = vquad + D(i)*(GOPT(i)+0.5_DP*HD1(i))
      if (vquad>zero) then
         print *," Warning: the TR subproblem was not well solved!"
         t = 0._DP
         do i=1,n
           t = t + D(i)**2
         enddo
         print *, " vquad=", vquad, " Stepsize=",dsqrt(t)
!         if (dsqrt(t)>=0.5_DP*DELTA) stop
      endif
      RETURN
!
!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
  210 IH=0
      DO 220 J=1,N
      HS(J)=ZERO
      DO 220 I=1,J
      IH=IH+1
      IF (I < J) HS(J)=HS(J)+HQ(IH)*S(I)
  220 HS(I)=HS(I)+HQ(IH)*S(J)

      IF (dabs(CRVMIN) > eps) GOTO 50
      IF (ITERC > ITCSAV) GOTO 150
      DO 260 I=1,N
  260 HRED(I)=HS(I)
      GOTO 120
      END subroutine



      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,&
       KNEW,W)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)

!     The arrays BMAT and ZMAT are updated, as required by the new position
!     of the interpolation point that has the index KNEW. The vector VLAG has
!     N+NPT components, set on entry to the first NPT and last N components
!     of the product Hw in equation (4.11) of the Powell (2006) paper on
!     NEWUOA. Further, BETA is set on entry to the value of the parameter
!     with that name, and DENOM is set to the denominator of the updating
!     formula. Elements of ZMAT may be treated as zero if their moduli are
!     at most ZTEST. The first NDIM elements of W are used for working space.

!     Set some constants.

      ONE=1.0_DP
      ZERO=0.0_DP
      NPTM=NPT-N-1
      ZTEST=ZERO
      DO 10 K=1,NPT
      DO 10 J=1,NPTM
   10 ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
      ZTEST=1.0D-20*ZTEST

!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.

      JL=1
      DO 30 J=2,NPTM
      IF (DABS(ZMAT(KNEW,J)) > ZTEST) THEN
          TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,1)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 20 I=1,NPT
          TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
   20     ZMAT(I,1)=TEMP
      END IF
      ZMAT(KNEW,J)=ZERO
   30 CONTINUE

!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.

      DO 40 I=1,NPT
      W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
   40 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      VLAG(KNEW)=VLAG(KNEW)-ONE

!     Complete the updating of ZMAT.

      TEMP=DSQRT(DENOM)
      TEMPB=ZMAT(KNEW,1)/TEMP
      TEMPA=TAU/TEMP
      DO 50 I=1,NPT
   50 ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)

!     Finally, update the matrix BMAT.

      DO 60 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 60 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I > NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   60 CONTINUE
      RETURN
      END subroutine


      !=======================================================================
      ! The DFPMIN procedure minimizes a user-written function Func of two or more independent variables using the
      ! Broyden-Fletcher-Goldfarb-Shanno variant of the Davidon-Fletcher-Powell method, using its gradient.
      ! DFPMIN is based on the routine dfpmin described in section 10.7 of Numerical Recipes in C:
      ! The Art of Scientific Computing (Second Edition), published by Cambridge University Press.

      SUBROUTINE EST_dfpmin(p,fret,func, nwrite,itmax,gtol)
        USE UTILITIES, ONLY : imaxloc,iminloc, outerprod

        !   use nr, only : lnsrch
        IMPLICIT NONE
        REAL(dp), DIMENSION(:), INTENT(inout) :: p
        REAL(dp), INTENT(OUT) :: fret
        INTEGER, INTENT(IN)::nwrite,  itmax
        REAL(dp), INTENT(IN) :: gtol
        INTERFACE
           FUNCTION func(x)
             USE UTILITIES, ONLY: DP
             IMPLICIT NONE
             REAL(DP), DIMENSION(:), INTENT(IN) :: x
             REAL(DP) :: func
           END FUNCTION func
        END INTERFACE

        REAL(dp), PARAMETER :: stpmx=100.0_dp,eps=EPSILON(p),tolx=4.0_dp*eps, xinc=1.0D-05
        INTEGER(i4b) :: its,i,k,nn,iter
        LOGICAL(lgt) :: check
        REAL(dp) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi, se1avg, se2avg,coravg,se1obs,se2obs,corobs
        REAL(dp), DIMENSION(SIZE(p)) :: dg,g,hdg,pnew,xi
        REAL(dp), DIMENSION(SIZE(p),SIZE(p)) :: hessin


        !==================================start here ======================
        WRITE(*,*) 'dfpmin is running..'
      498 FORMAT(1i6,f20.10)
      332 FORMAT(1i6,20f20.10)
        nn=SIZE(p)
        fp=func(p)
        g=dograd(p,xinc,func)
        !  write(*,*) 'after dograd is:',g

        IF (nwrite == 1) THEN
           its = 0
           WRITE(6,498) its,fp
           WRITE(6,332) its,(p(i),i=1,nn)
           WRITE(6,332) its,(g(i),i=1,nn)
        ENDIF

        !  pause
        hessin=0.0_DP
        DO i=1,nn
           hessin(i,i)=1.0_dp
        END DO

        xi=-g
        stpmax=stpmx*MAX(dsqrt(DOT_PRODUCT(p,p)),DBLE(SIZE(p)))
        DO its=1,itmax
           iter=its
           ! write (*,*) 'just before lnsrch'
           ! write(*,*) size(p), size(g),nn
           CALL lnsrch(func,p,   fp,  g,xi,pnew,fret,stpmax,check,nn)
           fp=fret
           xi=pnew-p
           p=pnew
           IF (MAXVAL(dabs(xi)/MAX(dabs(p),1.0_dp)) < tolx) RETURN
           dg=g
           g=dograd(p,xinc,func)
           den=MAX(fret,1.0_dp)
           IF (MAXVAL(dabs(g)*MAX(dabs(p),1.0_dp)/den) < gtol) RETURN
           dg=g-dg
           hdg=MATMUL(hessin,dg)
           fac=DOT_PRODUCT(dg,xi)
           fae=DOT_PRODUCT(dg,hdg)
           sumdg=DOT_PRODUCT(dg,dg)
           sumxi=DOT_PRODUCT(xi,xi)
           IF (fac > dsqrt(eps*sumdg*sumxi)) THEN
              fac=1.0_dp/fac
              fad=1.0_dp/fae
              dg=fac*xi-fad*hdg
              hessin=hessin+fac*outerprod(xi,xi)-&
                   fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
           END IF

           IF (nwrite == 1 .AND. MOD(its,5)==0) THEN
              WRITE(6,498) its,fp
              WRITE(6,332) its,(p(i),i=1,nn)
              WRITE(6,332) its,(g(i),i=1,nn)

           ENDIF

           xi=-MATMUL(hessin,g)
        END DO

      CONTAINS

        FUNCTION dograd(x,xinc,func)
          REAL(dp), DIMENSION(:), INTENT(IN):: x
          REAL(dp), DIMENSION(SIZE(x)):: dograd
          REAL(dp), INTENT(IN):: xinc
          INTERFACE
             FUNCTION func(x)
               USE UTILITIES, ONLY: DP
               IMPLICIT NONE
               REAL(DP), DIMENSION(:), INTENT(IN) :: x
               REAL(DP) :: func
             END FUNCTION func
          END INTERFACE

          INTEGER, PARAMETER:: nmax=100
          INTEGER:: n,i,j
          REAL(dp), DIMENSION(SIZE(x)):: dh,  yy, grad
          REAL(dp), DIMENSION(SIZE(x),SIZE(x)):: ee
          REAL(dp):: tempf2, tempf1, tolera

          n=SIZE(x)
          !	    write (*,*) 'inside dograd', n,x
          tolera=10.0_dp*xinc

          IF (n .GT. nmax) PAUSE 'nmax exceeded in dograd'
          DO i = 1,n
             IF (dabs(x(i)) >= tolera) THEN
                dh(i) = x(i)*xinc
             ELSE
                dh(i) = xinc
             ENDIF
          END DO

          ee(1:n,1:n)= 0.0_dp
          DO i = 1,n
             ee(i,i) = dh(i)
          END DO


          DO i = 1,n
             DO j = 1,n
                yy(j) = x(j) - ee(i,j)
             END DO
             tempf1 = func(yy)
             !            write (*,*) 'inside dograd, yy1', yy, tempf1
             DO j = 1,n
                yy(j) = x(j) + ee(i,j)
             END DO
             tempf2 = func(yy)
             grad(i) = (tempf2-tempf1)/(2.0d+00*dh(i))
             !   	     write (*,*) 'inside dograd, n, x', n,x, yy
             !		     write (*,*) 'inside dograd, yy2', yy, tempf2
          END DO

          dograd=grad
          RETURN
        END FUNCTION dograd

        !========================================================================
        SUBROUTINE lnsrch(func,xold,fold,g,pdir,x,f,stpmax,check,nest)
          USE UTILITIES, ONLY : assert_eq

          IMPLICIT NONE
          INTERFACE
             FUNCTION func(xx)
               USE UTILITIES, ONLY: DP
               IMPLICIT NONE
               REAL(DP), DIMENSION(:), INTENT(IN) :: xx
               REAL(DP) :: func
             END FUNCTION func
          END INTERFACE
          INTEGER(i4B), INTENT(IN):: nest
          REAL(dp), DIMENSION(:), INTENT(IN) :: xold,g
          REAL(dp), DIMENSION(:), INTENT(inout) :: pdir
          REAL(dp), INTENT(IN) :: fold,stpmax
          REAL(dp), DIMENSION(:), INTENT(OUT) :: x
          REAL(dp), INTENT(OUT) :: f
          LOGICAL(lgt), INTENT(OUT) :: check


          REAL(dp), PARAMETER :: alf=1.0e-4_dp,tolx=EPSILON(x)
          INTEGER(i4b) :: ndum
          REAL(dp) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
          !    write(*,*)'inside linesearch', nest
          !    write(*,*) size(g), size(pdir),  size(xold)
          !    write(*,*)'poin1'
          !
          !  	write(*,*) xold
          !    write(*,*)'point2'
          !	write(*,*) pdir
          !    write(*,*)'point3'
          !    write(*,*) g
          !  	332 format(8f15.8)

          !   ndum=assert_eq((/size(g),size(pdir),size(x),size(xold)/),'lnsrch')
          ndum=nest
          check=.FALSE.
          pabs=dsqrt(DOT_PRODUCT(pdir(:),pdir(:)))
          IF (pabs > stpmax) pdir(:)=pdir(:)*stpmax/pabs
          slope=DOT_PRODUCT(g,pdir)
          IF (slope >= 0.0) THEN
             WRITE (*,*) 'roundoff problem in lnsrch.. quitting'
             RETURN
          ENDIF

          alamin=tolx/MAXVAL(dabs(pdir(:))/MAX(dabs(xold(:)),1.0_dp))
          alam=1.0_DP
          !    write(*,*) 'alamin is',alamin
          !    write(*,*) 'alam is', alam

          DO
             x(:)=xold(:)+alam*pdir(:)
             !        write(*,*) 'x is'
             !		write(*,332)  x
             !		write(*,*) 'pdir is'
             !		write(*,332)  pdir
             !		write(*,*) 'xold is'
             !		write(*,332)  xold

             f=func(x)
             IF (alam < alamin) THEN
                x(:)=xold(:)
                check=.TRUE.
                RETURN
             ELSE IF (f <= fold+alf*alam*slope) THEN
                RETURN
             ELSE
                IF (alam == 1.0) THEN
                   tmplam=-slope/(2.0_dp*(f-fold-slope))
                ELSE
                   rhs1=f-fold-alam*slope
                   rhs2=f2-fold-alam2*slope
                   a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
                   b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                        (alam-alam2)
                   IF (a == 0.0) THEN
                      tmplam=-slope/(2.0_dp*b)
                   ELSE
                      disc=b*b-3.0_dp*a*slope
                      IF (disc < 0.0) THEN
                         tmplam=0.5_dp*alam
                      ELSE IF (b <= 0.0) THEN
                         tmplam=(-b+dsqrt(disc))/(3.0_dp*a)
                      ELSE
                         tmplam=-slope/(b+dsqrt(disc))
                      END IF
                   END IF
                   IF (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
                END IF
             END IF
             alam2=alam
             f2=f
             alam=MAX(tmplam,0.1_dp*alam)
          END DO
        END SUBROUTINE lnsrch

      END SUBROUTINE EST_dfpmin

END MODULE minimize
