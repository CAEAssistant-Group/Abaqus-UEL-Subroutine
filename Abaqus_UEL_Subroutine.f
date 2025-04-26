C      ######################################################################
C      #################      CAE Assistant Company          ################
C      ##############         CAEassistant-com              #############
C      ###########   Copy right by CAE Assistant Company    ###############
C      ######################################################################

C      ######################################################################
C      ######################################################################
C      CAE Assisitant Services: 
C      Toturial Packages,Consultancy,Articles,Q&A,Video Gallery,Online Course
C      ######################################################################
C      Need help with your project? 
C      You can get initial free consultation from (our website caeassistant com)
C      ###################################################################### 	

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

C      DOUBLE PERCISION   AK,F,WEIGHT,SHAP,GP
      DIMENSION WEIGHT(2,1),GP(2,1),AK(8,8),F(8,1),UU(8,1),UN(9,9)
      DIMENSION SHAP(2,8),DE(3,3),XX(4,2),B(3,8),AJACOBIAN(2,2),VN(9,9) 
      DIMENSION AINVJ(2,2),AA(8,3),CC(8,8),AB(12,1),AC(1,8),BB(1,1) 
      DIMENSION DSHAP0(2,4),DSHAP1(2,4),AMASS(8,8),SS1(3,1),STRESS(3,1) 
c******************************************************
      PARAMETER (NDIM=2, NDOF=3, NDI=2, NSHR=1, NNODEMAX=4,
     1 NTENS=6, NINPT=4, NSVINT=3, NELEMENT=4)
     
C
c      DATA WGHT /ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE/      
c******************************************************      
      REAL*8 UVAR
      INTEGER JELEM,K1,KINTK
      COMMON/KUSER/UVAR(NELEMENT,NSVINT,NINPT)
c     NINPT=NUMBER OF INTEGRATION POINTS PER ELEMENT
C     NSVINT=NUMBER OF STATE VARIBLE PER INTEGRATION POINT
C     NELEMENT=NUMBER OF ELEMENTS      
c******************************************************      
      

C
      AK=0.D0  
      NG=2.D0
      WEIGHT(1,1)=1.D0
      WEIGHT(2,1)=1.D0
      GP(1,1)=-0.57735D0
      GP(2,1)=0.57735D0


      RO=PROPS(1)
      EE=PROPS(2)
      ANU=PROPS(3)
      
      DO K1=1,NDOFEL
        RHS(K1,1)=0.
          DO K2=1,NDOFEL
             AMATRX(K2,K1)=0.
          END DO
      END DO
        
      DE(1,1)=EE/(1.0-ANU**2)
      DE(1,2)=EE*ANU/(1.0-ANU**2)
      DE(2,1)=EE*ANU/(1.0-ANU**2)
      DE(2,2)=EE/(1.0-ANU**2)
      DE(3,3)=EE/(2*(1.0+ANU))
        
      k = 0
      DO R = 1, NDOFEL
        k = k + 1
        UU(R,1) = U(k)
      ENDDO
        
      DO II=1,NG       
        X=GP(II,1)
            DO JJ=1,NG
C======================================================
                IF(II==1.0 .AND. JJ==1.0)THEN
                   KINTK=1.0
                ELSEIF(II==2.0 .AND. JJ==1.0)THEN
                   KINTK=2.0
                ELSEIF(II==1.0 .AND. JJ==2.0)THEN
                   KINTK=3.0
                ELSEIF(II==2.0 .AND. JJ==2.0)THEN
                   KINTK=4.0
                END IF
C=====================================================
                Y=GP(JJ,1)
                SHAP=0.D0
                SHAP(1,1)=(1.-X)*(1.-Y)/4
                SHAP(2,2)=(1.-X)*(1.-Y)/4
                SHAP(1,3)=(1.+X)*(1.-Y)/4
                SHAP(2,4)=(1.+X)*(1.-Y)/4
                SHAP(1,5)=(1.+X)*(1.+Y)/4
                SHAP(2,6)=(1.+X)*(1.+Y)/4
                SHAP(1,7)=(1.-X)*(1.+Y)/4
                SHAP(2,8)=(1.-X)*(1.+Y)/4
      
      
                DSHAP0(1,1)=(-1)*(1.-Y)/4
                DSHAP0(1,2)=(1)*(1.-Y)/4
                DSHAP0(1,3)=(1)*(1.+Y)/4
                DSHAP0(1,4)=(-1)*(1.+Y)/4
                DSHAP0(2,1)=(1.-X)*(-1)/4
                DSHAP0(2,2)=(1.+X)*(-1)/4
                DSHAP0(2,3)=(1.+X)*(1)/4
                DSHAP0(2,4)=(1.-X)*(1)/4
        
					!Hidden 100 Line(s)
                    END DO    
        
                STRESS=0.D0
                !!!! stress calculation
                SS1=MATMUL(B,UU)
                STRESS=MATMUL(DE,SS1)
        
C                    DO K1=1,3
C                      SVARS(K1)=0.0
C                  END DO        
c***********************************
c                    PRINT*,'KINTK',KINTK
c                    PRINT*,'STRESS', STRESS
c                    PRINT*,'JELEM', JELEM               
                    DO K1=1,NSVINT
                       SVARS(NSVINT*(KINTK-1)+K1)=STRESS(K1,1)
                       UVAR(JELEM,K1,KINTK) = SVARS(NSVINT*(KINTK-1)+K1)
c                       PRINT*,'K1',K1
c                      PRINT*,'STRESS(K1,1)',STRESS(K1,1)
c                       PRINT*,'SVARS(K1)',SVARS(NSVINT*(KINTK-1)+K1)
c                      PRINT*,'UVAR(K1)',UVAR(JELEM,K1,KINTK)
                    END DO 
c***********************************                
        
               AA=MATMUL(TRANSPOSE(B),DE)
               CC=MATMUL(AA,B)
                !Hidden 20 Line(s)
        END DO
        END DO
        
        
        IF (LFLAGS(3).EQ.4) THEN
            AMATRX=AMASS
        END IF
        
        IF (LFLAGS(3).EQ.1) THEN
C       Normal implicit time incrementation procedure.
            !Hidden Line(s)
C       *STATIC
                AMATRX=AK

        
			!Hidden 10 Line(s)
                    END DO
            END IF
        END IF
c********************************************************        

c********************************************************              
      RETURN
      END
      
c*********************************************************      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV)
      PARAMETER(NINT=4, NSTV=3, NELEMENT=4, ELEMOFFSET=4)
		!Hidden 15 Lines
      END DO
      RETURN
      END
c*********************************************************      
      
  
      

      
      
