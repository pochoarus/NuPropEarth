*! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
*! TAUSIC: TAU-lepton SImulation Code (for sea water)
*! supported by V.A. Kudryavtsev, University of Sheffield
*!
*!
*! Last modification by Vitaly Kudryavtsev, November 2004.
*!
*! Currently two random number generators RANLUX and RANMAR are used in TAUSIC.
*! Do not forget to initialize random number generators 
*! using RLUXGO and RMARIN (see CERNLIB documentation)
*! Cross-sections are calculated up to 10^20 eV muon energy.
*!
*! For any problems please contact:
*! 
*! v.kudryavtsev@sheffield.ac.uk
*! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
*!
*! TAUSIC is a 3D muon propagation code:
*! it takes into account energy losses and deflections
*! due to breemstrahlung, pair production, inelastic scattering
*! and ionisation.
*! It simulates also the angular and lateral deflections due to
*! multiple scattering
*! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

      subroutine initialize_tausic_sw(ic)
      implicit real*8 (a-h,o-z)

      REAL*8 n_a,n_z,n_rho,n_lambda

      integer*4 j0
      real*4 v2(100),cs1(100,111,4),cs2(100,111,4)
      real*4 emulo(121,6),emc(121)
      real*4 ema(61),va(30),anga(50),dang(50,30,61)

      COMMON/sig_tau_sw/v2,cs1,cs2,j0
      COMMON/coenlo_tau_sw/emulo,emc
      COMMON/ang_tau_sw/ema,va,anga,dang
      COMMON/rock2_tau_sw/n_a,n_z,n_rho,n_lambda

* Rock description
      n_z=6.8442     !Z, 11.00 for standard rock, 6.6 for pure water
      n_a=12.393     !A, 22.00 for standard rock, 11.89 for pure water
      n_rho=1.03     !density, 2.65 for standard rock
      n_lambda=36.08 !rad. length, 26.48 for standard rock
      rho=n_rho

*       reading of the muon energy losses
*      PRINT *,'tausic: loading energy losses table...'
      if(ic.eq.1) then
      OPEN(UNIT=4,file='tausic-eloss-sw-bs-2004.dat',
     *form='formatted',status='old')
      end if
      if(ic.eq.2) then
      OPEN(UNIT=4,file='tausic-eloss-sw-f2-2004.dat',
     *form='formatted',status='old')
      end if
      READ(4,99)emulo
   99 FORMAT(24(5E12.5/),E12.5/)
      CLOSE(4)

*       reading of the arrays with the normalized 
*	(and integrated) double differential
*       cross-section of muon inelastic scattering
*       log10(E_mu(GeV))= 0.0 - 6.0 with a step of 0.1 (ema(61))
*       log10(v)= -0.1 - -3.0 with a step of 0.1 (va(30))
*       log10(theta(rad))= 0.4 - -4.5 with a step of 0.1 (anga(50))
*       normalized double differential cross-section dang(50,30,61)
*       was calculated using the program sign_ang_arr.for with
*       the subroutine sign_ang.for
*      PRINT *,'tausic: loading inel. scatt. cross section...'
      if(ic.eq.1) then
      OPEN(UNIT=3,file='tausic-double-diff-seawater.dat',
     *  FORM='formatted',STATUS='old')
      end if
      if(ic.eq.2) then
      OPEN(UNIT=3,file='tausic-double-diff-seawater-f2.dat',
     *  FORM='formatted',STATUS='old')
      end if
      READ(3,122)ema,va,anga,dang
  122 FORMAT(5e15.7)
      CLOSE(3)

C       J0 - THE MINIMAL VALUE OF RELATIVE ENERGY LOSSES (E/EMU);
C               STARTING FROM THIS VALUE ALL KINDS OF TAU ENERGY
C               LOSSES EXCEPT FOR IONISATION ARE TREATED AS THE
C               THE FLUCTUATED ONES, FOR EXAMPLE: J0=201 MEANS
C               THAT FROM V=E/EMU=10**(-10) UP TO V=1 THE LOSSES
C               ARE TREATED AS THE FLUCTUATED ONES.
      J0=61

      DO I=1,100
        V2(I)=0.0-0.05*(I-1)
      END DO

      DO IL=1,121
        EMC(IL)=-1.0+0.1*(IL-1)
      END DO

C       READING THE INTEGRAL CROSS-SECTION
*      PRINT *,'tausic: loading integral cross sections...'
      if(ic.eq.1) then
      OPEN(UNIT=1,file='tausic-cross-sections-sw-bs-2004.dat',
     *  FORM='formatted',STATUS='old')
      end if
      if(ic.eq.2) then
      OPEN(UNIT=1,file='tausic-cross-sections-sw-f2-2004.dat',
     *  FORM='formatted',STATUS='old')
      end if
      READ(1,3)CS2
    3 FORMAT(4(111(20(5E14.6/)/)/)/)
      CLOSE(UNIT=1)


C       FILLING THE ARRAY OF INTEGRAL CROSS-SECTIONS FOR THE SIMULATION

      do i=1,4
        do k=1,111
          do j=1,j0
             if(cs2(j0,k,i).gt.0.) then
      cs1(j,k,i)=cs2(j,k,i)/cs2(j0,k,i)
             end if
          end do
        end do
      end do

      return
      end

***********************************************************************

      subroutine tau_transport_sw
     *     (x0,y0,z0,cx0,cy0,cz0,Emuin0,depth0,tmu0,rho,
     *     ms_flag,d_flag,xf,yf,zf,cxf,cyf,czf,Emuout,depthf,tmuf,
     *     idec,itflag,tauti,tautf)
      implicit real*8 (a-h,o-z)
      REAL*8 n_a,n_z,n_rho,n_lambda,mmu
      real*4 yfl
      INTEGER*4 ms_flag,d_flag
      COMMON/rock2_tau_sw/n_a,n_z,n_rho,n_lambda
      parameter (pi=3.141592654d0)
      parameter (mmu=1.777d0)
      emuin=emuin0
      tautf = 0.d0
      n_rho = rho
      idec = 0
      ttau = 0.d0

                                !input and output parameters:
        emu_f=0.d0                !final tau energy

        depth=depth0*n_rho      !depth along z-axis in g/cm^2, input
        path_max=1000000.*100.  !maximal path of tau, input
        emu=emuin               !initial tau energy, input

        theta=0.d0                !angle with respect to z-axis, input and output
        phi=0.d0                  !angle with respect to x-axis, input and output
        dr=0.d0                   !tau path in the sea water, output
        x=0.d0                    !deviation along the axes x,y, perpendicular
        y=0.d0                    !to z-axis, x,y-axes should be defined,
                                !input and output
        z=0.d0                    !z-coordinate
        t=0.d0                    !Initial tau path in the sea water, input
                                !emu_f - final energy of tau

        emu_0=emu
        theta0=theta
        phi0=phi

        if(itflag.eq.0) then
           taudec0 = 2.91d-4      ! nsec - the mean life time of tau
           call ranlux(yfl,1)
           taudec = -taudec0 * alog(yfl) ! the life time of tau
        end if
        if(itflag.gt.0) then
           taudec = tauti
        end if
        depth_tau = depth
c        if(depth_dec.lt.depth) then
c           depth_tau = depth_dec
c           idec = 1
c        end if
         
        CALL tausic_sw(emu,depth,emu_f,x,y,z,t,theta,phi,dr,
     *       path_max,ms_flag,d_flag,idec,taudec,ttau)

        gamma=(Emu_0+Emu_f)/2./mmu
        beta=dsqrt(gamma*gamma-1.)/gamma
        vmu=beta*29.9792458
        tmu=t/vmu/n_rho
        tmuf=tmu0+tmu

        if(idec.eq.0) then
           tautf = taudec - ttau
        end if

        if(idec.eq.1) then
           trat = taudec / tmu * gamma
           if(trat.gt.1.) trat = 1.
           emu_f = emu_0 + trat * (emu_f-emu_0)
           tmuf = tmu0 + taudec * gamma
           x = x * trat
           y = y * trat
           z = z * trat
           theta = theta * trat
           dr = dr * trat
        end if
        depthf = z / rho

        
        thetac=dacos(cz0)
        IF(cx0.NE.0.) THEN
         phic=datan(cy0/cx0)
        IF(cx0.LT.0.) phic=pi+phic
        END IF
        IF(cx0.EQ.0.) THEN
        phic=pi/2.
        IF(cy0.LT.0.) phic=pi*3./2.
        END IF
        IF(phic.GT.pi*2.) phic=-pi*2.+phic
        IF(phic.LT.0.) phic=phic+pi*2.

        call coord_transform_tau_sw(x1,y1,z1,x,y,z,thetac,phic)
        xf=x1/n_rho+x0
        yf=y1/n_rho+y0
        zf=z1/n_rho+z0
        theta0=thetac
        phi0=phic
        call angle_transform_tau_sw(theta0,phi0,theta,phi)
        czf=dcos(theta0)
        cxf=dsin(theta0)*dcos(phi0)
        cyf=dsin(theta0)*dsin(phi0)
        emuout=emu_f                   !if emu_f is less than 0.106 the muon
                                       ! was stopped
        if(emuout.le.mmu+0.001) emuin0=0.d0

        return
        end

*! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
*! MUSIC: MUon SImulation Code
*! by P. Antonioli, C. Ghetti, E.V. Korolkova, 
*! V.A. Kudryavtsev, G. Sartorelli
*!
*! For any problems please contact:
*! 
*! antonioli@bo.infn.it
*! v.kudryavtsev@sheffield.ac.uk
*! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

*! MUSIC is a 3D muon propagation code:
*! it takes into account energy losses and deflections
*! due breemstrahlung, pair production, inelastic scattering
*! and ionisation.
*! It simulates also the angular and lateral deflections due to
*! multiple scattering

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE TAUSIC_SW(Emu,zf,Emu_f,x,y,z,t0,theta,phi,dr,t10,
     *     ms_flag,d_flag,idec,taudec,ttau)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C       V2(100) - ARRAY OF LOG OF RELATIVE ENERGY TRANSFERS (E/EMU) FOR WHICH
C               THE MUON INTERACTION INTEGRAL CROSS-SECTIONS WERE CALCULATED:
C               FROM 0.0 TO -5.0 WITH THE STEP OF -0.05. THE INTEGRATION STARTS
C               FROM FROM V2=ALOG10(E/EMU)=0.0 AND ENDS AT V2(I).
C       CS1(100,111,4) - ARRAY OF NORMALIZED INTEGRAL CROSS-SECTIONS FOR V2(50),
C               LOG OF MUON ENERGY (111) (FROM 0.0 (1 GeV) TO 11.0 (100 EeV)
C               WITH THE STEP OF 0.1) AND 4 PROCESSES TAKING INTO ACCOUNT
C               (BREEMSTRAHLUNG, PAIR PRODUCTION, INELASTIC SCATTERING
C               AND IONISATION)

	implicit real*8 (a-h,o-z)
      REAL*8 mmu,n_a,n_z,n_rho,n_lambda
      real*4 yfl
      INTEGER*4 err
      INTEGER*4 ms_flag,d_flag
      PARAMETER (mmu=1.777d0)
c      PARAMETER (ms_flag=1)   !to switch off mult.scatt. put ms_flag=0
c      PARAMETER (d_flag=1)    !to switch off scatt. due to other processes
c                              ! put d_flag=0

      real*4 v2(100),cs1(100,111,4),cs2(100,111,4)
      integer*4 j0
      real*4 emulo(121,6),emc(121)
      COMMON/sig_tau_sw/v2,cs1,cs2,j0
      COMMON/coenlo_tau_sw/emulo,emc
      COMMON/rock2_tau_sw/n_a,n_z,n_rho,n_lambda

*       ema(61) - array of log of muon energies (GeV), 0.0, 0.1,...6.0.
*       va(30) - array of log of relative energy transfers, 0.0,-0.1,...-3.0.
*       anga(50) - array of log of scattering angles (rad), 0.4,0.3,...-4.5.
*       dang(50,30,61) - integral probability for muon scattering at an angle
*       larger than anga(i) for fixed values of ema(k) and va(j).

      real*4 ema(61),va(30),anga(50),dang(50,30,61)
      COMMON/ang_tau_sw/ema,va,anga,dang

      em=emu
      t1=t0
      err=0

    9 em1=EM
      J=dlog10(em1)*10.+2-10
      IF(J.GE.111) J=111
      J1=J-1

C       IF THE MUON ENERGY BECOME LESS THAN 1 GEV THERE IS NO
C       SIMULATION OF FLUCTUATED ENERGY LOSSES

      IF(em1.LE.10.) GO TO 7

C       SIMULATION OF THE INTERACTION POINT

      call ranlux(yfl,1)

      FP=1./(CS2(J0,J,1)+CS2(J0,J,2)+CS2(J0,J,3)+CS2(J0,J,4))
      FP1=1./(CS2(J0,J1,1)+CS2(J0,J1,2)+CS2(J0,J1,3)+CS2(J0,J1,4))
      FP=(FP-FP1)/(EMC(J+10)-EMC(J1+10))*(DLOG10(EM1)-EMC(J1+10))+FP1
      T=-FP*ALOG(YFL)
      T3=T1+T
      z1=T*dcos(theta)

   15 IF(z+z1.LE.zf) GO TO 10

*       If the next interaction point is deeper than the level of observation
*       we will calculate the continuous energy losses and multiple scattering
*       up to the level of observation

      IF(t3.GT.t10) THEN
        em1=mmu
        dr=t10
        go to 16
      END IF
      t2=(zf-z)/dabs(dcos(theta))
      dz=t2*dcos(theta)

C       CALCULATION OF THE CONTINUOUS ENERGY LOSSES:
C       DUE TO IONISATION AND OTHER PROCESSES FOR V=E/EMU<10**(-3)

      CALL vdem_tau_sw(em1,B0,BS)
      emf=em1-(em1*BS+B0)*T2
      IF(emf.LE.mmu) emf=mmu
      emi=10.**((dlog10(emf)+dlog10(em1))/2.)
      CALL vdem_tau_sw(emi,b0,bs)
      CALL vem_tau_sw(em1,b0,bs,t2,em3)
      IF(em3.LE.mmu) THEN
        dr=t1-t0+em1/b0
        em1=mmu
        go to 16
      END IF
      theta2=0.
      phi2=0.
      deltax=0.
      deltay=0.
      deltar=0.
      dx1=0.
      dy1=0.
      dz1=dz

*       simulation of the multiple scattering
      IF(ms_flag.EQ.1) THEN
        CALL multiple_tau_sw
     *(t2,em1,em3,b0,bs,theta2,phi2,deltax,deltay,deltar,err)
cc      IF (err.ne.0) WRITE (*,3057) err,em1,em3
 3057   FORMAT(1x,'%MUSIC-W, multiple error status ',i1,2(1x,f10.5),
     &    ' (standard call)')
      ENDIF

*       transformation of the coordinates
      CALL coord_transform_tau_sw
     *     (dx1,dy1,dz1,deltax,deltay,dz,theta,phi)
*       transformation of the angles
      IF(theta2.GT.0.) THEN
        CALL angle_transform_tau_sw(theta,phi,theta2,phi2)
      END IF

      gamma = (em1 + em3) / 2. / mmu
      beta=dsqrt(gamma*gamma-1.)/gamma
      vmu=beta*29.9792458
      ttau = ttau + t2 / vmu / n_rho / gamma
      em1=EM3
      emu_f = em3
      IF(em3.LE.mmu) go to 16
      NM=1
      T=T-T2+deltar
      z=z+dz1
      x=x+dx1
      y=y+dy1
      t1=t1+t2+deltar
      t0 = t1
      dr = t1
*      print *,taudec,ttau
      if(taudec.le.ttau) then
         idec = 1
         go to 16
      end if

      GO TO 11

   10 CONTINUE

*       Calculation of the continuous energy losses of tau up to the next
*       interaction point
      CALL vdem_tau_sw(em1,B0,BS)
      emf=em1-(em1*BS+B0)*T
      IF(emf.LE.mmu) emf=mmu
      emi=10.**((dlog10(emf)+dlog10(em1))/2.)
      CALL vdem_tau_sw(emi,B0,BS)
      CALL vem_tau_sw(em1,B0,BS,T,EM2)
      IF(em2.LE.mmu) THEN
        dr=t1-t0+em1/b0
        em1=mmu
        go to 16
      END IF
      theta2=0.
      phi2=0.
      deltax=0.
      deltay=0.
      deltar=0.
      dx1=0.
      dy1=0.
      dz1=z1

      IF(ms_flag.EQ.1) THEN
        CALL multiple_tau_sw           !multiple scattering between 2 interactions
     *(t,em1,em2,b0,bs,theta2,phi2,deltax,deltay,deltar,err)
 3058   FORMAT(1x,'%MUSIC-W, multiple error status ',i1,2(1x,f10.5),
     &    ' (between two int. call)')
      END IF

*       transformation of the muon deflection to the base reference frame
      CALL coord_transform_tau_sw
     *     (dx1,dy1,dz1,deltax,deltay,z1,theta,phi)

      gamma = (em1 + em2) / 2. / mmu
      beta=dsqrt(gamma*gamma-1.)/gamma
      vmu=beta*29.9792458
      ttau = ttau + t / vmu / n_rho / gamma
      em1=em2
      emu_f = em2
      IF(em1.LE.mmu) go to 16
      T1=T3+deltar
      z=z+dz1
      x=x+dx1
      y=y+dy1

      IF(theta2.GT.0.) THEN
*       transformation of the deflection angles to the base reference frame
        CALL angle_transform_tau_sw(theta,phi,theta2,phi2)
      END IF

*      print *,taudec,ttau
      t0 = t1
      dr = t1
      if(taudec.le.ttau) then
         idec = 1
         go to 16
      end if

      if(z.gt.zf) go to 11
      IF(em1.LE.10.) GO TO 7

      call ranlux(yfl,1)
      J=dlog10(em1)*10.+2-10
      J2=dlog10(em1)*10.+0.5+1-10
      IF(J.GE.111) J=111
      IF(J2.GE.111) J2=111
      J1=J-1

C       SIMULATION OF THE TYPE OF MUON INTERACTION
      css=cs2(j0,j,1)+cs2(j0,j,2)+cs2(j0,j,3)+cs2(j0,j,4)
      css1=cs2(j0,j1,1)+cs2(j0,j1,2)+cs2(j0,j1,3)+cs2(j0,j1,4)
      css=(css-css1)/(EMC(J+10)-EMC(J1+10))*(DLOG10(EM1)-EMC(J1+10))+
     *     css1
      al1=((cs2(j0,j,1)-cs2(j0,j1,1))/(EMC(J+10)-EMC(J1+10))*
     *     (DLOG10(EM1)-EMC(J1+10))+cs2(j0,j1,1))/css
      al2=((cs2(j0,j,2)-cs2(j0,j1,2))/(EMC(J+10)-EMC(J1+10))*
     *     (DLOG10(EM1)-EMC(J1+10))+cs2(j0,j1,2))/css
      al3=((cs2(j0,j,3)-cs2(j0,j1,3))/(EMC(J+10)-EMC(J1+10))*
     *     (DLOG10(EM1)-EMC(J1+10))+cs2(j0,j1,3))/css
      IF(yfl.GT.AL1) GO TO 20
      ip=1
      GO TO 23
   20 IF(YFL.GT.AL1+AL2) GO TO 21
      IP=2
      GO TO 23
   21 IF(YFL.GT.AL1+AL2+AL3) GO TO 22
      IP=3
      GO TO 23
   22 IP=4
   23 CONTINUE

C       SIMULATION OF THE FRACTION OF MUON ENERGY TRANSFERED
C       TO THE SECONDARY PARTICLE
      call ranlux(yfl,1)
      CALL vcs_tau_sw(j2,yfl,ip,v1)

      em=em1
c       simulation of the deflection angle of muon
      theta1=0.
      phi1=0.
      IF(ip.LE.3.AND.ip.GE.1.AND.d_flag.EQ.1) THEN
        CALL mu_scatt_tau_sw(ip,em,v1,theta1,phi1)
      END IF

      Demf=V1*em1
      EM=EM-Demf
      em1=em
c       transformation of the deflection angles to the base reference frame
      IF(theta1.GT.0.) THEN
        CALL angle_transform_tau_sw(theta,phi,theta1,phi1)
      END IF

      IF(em1.LE.mmu) go to 16
      GO TO 9
    7 CONTINUE

C       IF THE TAU ENERGY IS LESS THAN 10 GEV ONLY THE IONISATION
C       LOSSES ARE CALCULATED
      CALL vdem_tau_sw(em1,b0,bs)
      z4=zf-z
      t4=z4/dabs(dcos(theta))
      emf=em1-(em1*bs+b0)*t4
      IF(emf.LE.mmu) emf=mmu
      emi=10.**((dlog10(emf)+dlog10(em1))/2.)
      CALL vdem_tau_sw(emi,b0,bs)
      CALL vem_tau_sw(em1,b0,bs,t4,em3)
      IF(em3.LE.mmu) THEN
         dr=t1-t0+em1/b0
         em1=mmu
         go to 16
      END IF
      theta2=0.
      phi2=0.
      deltax=0.
      deltay=0.
      deltar=0.
      dx1=0.
      dy1=0.
      dz1=z4

*       simulation of the multiple scattering up to the level of observation
      IF(ms_flag.EQ.1) THEN
        CALL multiple_tau_sw
     *(t4,em1,em3,b0,bs,theta2,phi2,deltax,deltay,deltar,err)
 3059   FORMAT(1x,'%MUSIC-W, multiple error status ',i1,2(1x,f10.5),
     &    ' (after int. call)')
      END IF

*       transformation of the coordinates
      CALL coord_transform_tau_sw
     *     (dx1,dy1,dz1,deltax,deltay,z4,theta,phi)

      gamma = (em1 + em3) / 2. / mmu
      beta=dsqrt(gamma*gamma-1.)/gamma
      vmu=beta*29.9792458
      ttau = ttau + t4 / vmu / n_rho / gamma
      em1=em3
      emu_f = em3
      IF(em1.LE.mmu) GO TO 16
      nm=2
      z=z+dz1
      x=x+dx1
      y=y+dy1
      t1=t1+t4+deltar

*       transformation of the angles
      IF(theta2.GT.0.) THEN
        CALL angle_transform_tau_sw(theta,phi,theta2,phi2)
      END IF

      t0 = t1
      dr = t1
      if(taudec.le.ttau) then
         idec = 1
         go to 16
      end if

   11 CONTINUE

      dr=t1-t0
      emu_f=em1

   16 CONTINUE

      IF(em1.LE.mmu) emu_f=mmu

      RETURN
      END


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE vcs_tau_sw(j,yfl,ip,v1)
C       CALCULATION OF THE FRACTION OF MUON ENERGY TRANSFERRED
C       TO THE SECONDARY PARTICLE
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       implicit real*8 (a-h,o-z)
      real*4 v2(100),cs1(100,111,4),cs2(100,111,4)
      integer*4 j0
      COMMON/sig_tau_sw/v2,cs1,cs2,j0
      real*4 yfl
      i=1
    1 i=i+1
      if(yfl.gt.cs1(i,j,ip)) go to 1
      v3=(v2(i)-v2(i-1))/(cs1(i,j,ip)-cs1(i-1,j,ip))*
     *  (yfl-cs1(i-1,j,ip))+v2(i-1)
      v1=10.**v3
      RETURN
      END


*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE vdem_tau_sw(EM,B0,BS)
C       CALCULATION OF THE IONISATION LOSSES PER G/CM**2 AND
C       THE RELATIVE ENERGY LOSSES DUE TO OTHER PROCESSES FOR
C       V=E/EMU<10**(-3)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit real*8 (a-h,o-z)
      real*4 f(121,6),em0(121)
      COMMON/coenlo_tau_sw/f,em0

C       EM0 - ARRAY OF TAU ENERGIES FOR WHICH THE RELATIVE ENERGY LOSSES
C               FOR V<10**(-3) HAVE BEEN ALREADY CALCULATED

C       F - ARRAY OF RELATIVE ENERGY LOSSES FOR V<10**(-3)
C       NOTE: IN THIS VARIANT OF PROGRAM THE ARRAY FOR RELATIVE
C       ENERGY LOSSES FOR V<10**(-3) IS USED

C       IF THE TAU ENERGY IS LESS THAN 1 GEV ONLY THE IONISATION
C       LOSSES ARE CALCULATED

C       CALCULATION OF THE RELATIVE ENERGY LOSSES DUE TO BREEMSTRAHLUNG,
C       PAIR PRODUCTION AND INELASTIC SCATTERING FOR V<10**(-3)

      e1=dlog10(em)
      i=1
    4 i=i+1
      if(i.gt.121) go to 6
      if(e1.gt.em0(i)) go to 4
    6 if(i.gt.121) i=121
      f1=(f(i,5)-f(i-1,5))/(em0(i)-em0(i-1))*(e1-em0(i-1))+f(i-1,5)
      bs=f1
      if(bs.le.0.) bs=0.
      if(e1.le.0.) bs=0.
      f2=(f(i,6)-f(i-1,6))/(em0(i)-em0(i-1))*(e1-em0(i-1))+f(i-1,6)
      b0=f2
      IF(e1.LT.0.4) b0=f(4,6)
      if(b0.le.0.) b0=0.

      RETURN
      END

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE vem_tau_sw(EM0,B0,BS,T,EM)
C       CALCULATION OF THE CONTINUOUS ENERGY LOSSES
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit real*8 (a-h,o-z)
      IF(bs*t.LE.1.e-6) GO TO 1
      em=em0*dexp(-bs*t)-b0/bs*(1.-dexp(-bs*t))
      return
    1 em=em0-b0*t-bs*t*em0
      RETURN
      END

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE defl_tau_sw(em,ip,v,theta,phi)
*       Simulation of the muon scattering angles due to
*       inelastic scattering
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      implicit real*8 (a-h,o-z)
      real*4 ema(61),va(30),anga(50),dang(50,30,61)
      COMMON/ang_tau_sw/ema,va,anga,dang
      REAL*8 ang(2,2)
      real*4 yfl
	parameter (pi=3.141592654d0)	
      IF(ip.EQ.3) THEN
        call ranlux(yfl,1)
        IF(em.LE.1.) RETURN
        eml=dlog10(em)
        k1=10.*eml+1.
        if(k1.gt.60) go to 3
        IF(k1.LT.1) k1=1
        IF(k1.GT.60) k1=60
        k2=k1+1
        IF(k2.LT.2) k2=2
        IF(k2.GT.61) k2=61
        vl=dlog10(v)
        m1=-vl*10.
        if(m1.gt.29) go to 3
        IF(m1.LT.1) m1=1
        IF(m1.GT.29) m1=29
        m2=m1+1
        IF(m2.LT.2) m2=2
        IF(m2.GT.30) m2=30
        DO kk=1,2
          IF(kk.EQ.1) k=k1
          IF(kk.EQ.2) k=k2
          DO mm=1,2
            IF(mm.EQ.1) m=m1
            IF(mm.EQ.2) m=m2
            j=1
    1       j=j+1
            IF(j.GT.49) THEN
              ang(mm,kk)=-5.
              go to 2
            END IF
            IF(yfl.GT.dang(j,m,k)) go to 1
            ang(mm,kk)=(anga(j)-anga(j-1))/(dang(j,m,k)-dang(j-1,m,k))*
     *        (yfl-dang(j-1,m,k))+anga(j-1)
    2       CONTINUE
          END DO
        END DO
        ang1=(ang(2,1)-ang(1,1))/(va(m2)-va(m1))*(vl-va(m1))+ang(1,1)
        ang2=(ang(2,2)-ang(1,2))/(va(m2)-va(m1))*(vl-va(m1))+ang(1,2)
        angt=(ang2-ang1)/(ema(k2)-ema(k1))*(eml-ema(k1))+ang1
        IF(angt.GT.-4.5) theta=10.**angt
        IF(angt.LE.-4.5) theta=0.
        go to 4
 3      theta=0.
 4      continue
	call ranlux(yfl,1)
        phi=2.*pi*yfl
        if(theta.eq.0.) phi=0.
      END IF
      RETURN
      END

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE mu_scatt_tau_sw(jp,E,v,theta,phi)
      IMPLICIT NONE
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Input:
      INTEGER*4 jp
!               * jp = 1 : bremsstrahlung
!               * jp = 2 : pair-production
!               * jp = 3 : nuclear
      REAL*8 E    !initial energy of the muon
      REAL*8 v        !fractional energy loss of muon
* Output:
      REAL*8 theta,phi
      REAL*4 yfl
* theta respect to parent muon direction,always positive!!!
* phi orthogonal to plane defined by parent muon direction
* theta phi in deg
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      REAL*8 mmu,pi,n_z,n_a,n_rho,n_lambda,mel,conv
      PARAMETER (pi=3.141592654d0)
      PARAMETER (conv=57.29577951d0)
      PARAMETER (mmu=1.777d0)
      PARAMETER (mel=0.511d-3)

      real*4 ee,vv

      REAL*8 eoutmu,p2q
      COMMON/rock2_tau_sw/n_a,n_z,n_rho,n_lambda

* Bremsstrahlung declarations
      REAL*8 ak1,ak2,ak3,ak4,ak5,an,the1,the2,the3,the_mean,at


* Pair prod. declarations
      REAL*8 a_min

* Bremsstrahlung: following Van Ginneken

      IF (jp.EQ.1) THEN
         ak1=0.245*E**(-1./2.)          ! as for proton
         ak3=1.28*E**(-0.92)            ! as for proton
         ak4=0.86*E**(-0.85)            ! as for proton
         an=0.81*E**(1./3.)/(E**(1./3.)+2.0)
         ak2=0.059/E*n_z**(-0.25)
        at=min(ak1*dsqrt(v),ak2)
        the1=max(at,ak3*v)
        the2=ak4*v**(1.+an)*(1.-v)**(-an)
        ak5=ak4*0.5**(1.+an)*0.5**(-an)/0.5**(-0.5)
        the3=ak5/dsqrt(1.-v)
        IF(v.LE.0.5) the_mean=the1
        IF(v.GT.0.5.AND.the2.LT.0.2) the_mean=the2
        IF(v.GT.0.5.AND.the2.GE.0.2) the_mean=the3
        the_mean=the_mean*the_mean
   10   call ranlux(yfl,1)
        theta=-the_mean*alog(yfl)
        theta=dsqrt(theta)
        IF(theta.GT.pi) go to 10
        call ranlux(yfl,1)
        phi=(2*pi)*(yfl)                   !uniform between 0 and 2*pi
c        theta=theta*(conv)              !from rad to deg
c        phi=(phi)*(conv)

      ELSE IF (jp.EQ.2) THEN
        ee=E
        vv=v
         a_min=1.06e-3*sqrt(sqrt(vv))+0.090*vv/           ! as for proton
     *        (vv+0.35)
         the_mean=(2.3+dlog(E))/E*(1.-v)**(-0.1*dlog(E))* ! as for proton
     *        (v-2.*mel/E)**2/v/v*a_min
        the_mean=the_mean*the_mean
   11   call ranlux(yfl,1)
        theta=-the_mean*alog(yfl)
        theta=dsqrt(theta)
        IF(theta.GT.pi) go to 11
        call ranlux(yfl,1)
        phi=(2*pi)*(yfl)                   !uniform between 0 and 2*pi
c        theta=theta*(conv)              !from rad to deg
c        phi=(phi)*(conv)

* Inelastic scattering

      ELSE IF (jp.EQ.3) THEN
        CALL defl_tau_sw(E,jp,v,theta,phi)

      ELSE
        PRINT *,'Invalid jp argument on mu_scat',jp
        STOP
      ENDIF
      RETURN
      END

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE multiple_tau_sw
     *(x0,einmu,eoutmu,alpha,beta,theta,phi,deltax,deltay,deltar,err)
      IMPLICIT NONE
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Multiple scattering: following
* P.Lipari and T. Stanev, Phys. Rev. D, 44, 3543 (1991) and
* L. Highland, Nucl. Instr. and Methods, 129, 497 (1975)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Input:
      REAL*8 x0        !in km w.e.
      REAL*8 einmu     !initial energy of the tau
      REAL*8 eoutmu    !final energy of the tau
      REAL*8 alpha     !ionization energy loss of tau
      REAL*8 beta      !relative (continuous) energy losses of tau
                       !due to other processes
* Output:
      REAL*8 theta
      REAL*8 phi
      REAL*8 deltax,deltay
      REAL*8 deltar             !approximate increase of a tau path
      INTEGER*4 err             !error flag (err=1) indicates that
                                !einmu or eoutmu is less than mmu
                                !error flag (err=2) indicates that
                                !integral along path is LE.0

* theta respect to parent tau direction
* phi orthogonal to plane defined by parent tau direction
* theta phi in deg
* deltax,deltay - lateral displacement of the tau in g/cm^2
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Multiple scattering declarations
      REAL*8 mmu,pi,mel,conv
      PARAMETER (pi=3.141592654d0)
      PARAMETER (conv=57.29577951d0)
      PARAMETER (mmu=1.777d0)
      PARAMETER (mel=0.511d-3)


      REAL*8 ems,sigma,sigt,xx,yy,ro
      REAL*8 sigmat,sigmax,thetax,thetay,zz,sz,step,sz1
      REAL*8 p1,g1,bet1,p2,g2,bet2
      REAL*8 n_a,n_z,n_rho,n_lambda
      REAL*4 sigmatt,sigmaxx,roo,thetaxx,thetayy,xxx,yyy
      REAL*4 vv(2,2),cc(2,2),vcx(2)

C      PARAMETER (ems=0.0124)                 !in Gev (cfr Highland)
C      PARAMETER (ems=0.014)                  !in Gev (cfr Lipari)
      PARAMETER (ems=0.015)                  !in Gev (060697)

      COMMON/rock2_tau_sw/n_a,n_z,n_rho,n_lambda

      err=0
      IF(einmu.LE.mmu.OR.eoutmu.LE.mmu) THEN
        err=1
        RETURN
      END IF
      p1=dsqrt(einmu*einmu-mmu*mmu)
      g1=einmu/mmu
      bet1=dsqrt(g1*g1-1.)/g1
      p2=dsqrt(eoutmu*eoutmu-mmu*mmu)
      g2=eoutmu/mmu
      bet2=dsqrt(g2*g2-1.)/g2

      IF(beta.GT.0.) THEN
        sigt=1./alpha*(1./p2/bet2-1./p1/bet1)+beta/alpha/alpha*
     *    (dlog(dabs(1.+alpha/beta/p1/bet1))-
     *    dlog(dabs(1.+alpha/beta/p2/bet2)))
      END IF
      IF(beta.EQ.0.) THEN
        sigt=1./alpha*(1./p2/bet2-1./p1/bet1)
      END IF
      IF(sigt.LE.0.) THEN
        sigt=1./p1/bet1/p1/bet1
      END IF
      sigma=sigt*ems*ems/n_lambda
      sigmat=dsqrt(sigma)              !sigma(theta_x)=sigma(theta_y)
      IF(einmu.EQ.eoutmu) THEN                        !for test programs
        sigmat=dsqrt(ems*ems*x0/einmu/einmu/n_lambda)    !
      END IF						!
      sz=0.
      step=x0/100.
      zz=0.
    1 zz=zz+step/2.
      IF(beta.GT.0.) THEN
        sz1=p1*bet1*dexp(-beta*zz)-alpha/beta*(1.-dexp(-beta*zz))
      END IF
      IF(beta.EQ.0.) THEN
        sz1=p1*bet1-alpha*zz
      END IF
      IF(sz1.LE.mmu) THEN
        sz=0.
        go to 3
      END IF
      sz=sz+zz*zz*step/sz1/sz1
      zz=zz+step/2.
      IF(zz.LT.x0) go to 1
    3 CONTINUE
      IF(sz.LE.0.) THEN
        sz=x0*x0*x0/3./p1/bet1/p1/bet1
      END IF
      IF(sz.LE.0.) THEN
        err=2
        RETURN
      END IF
      sigmax=dsqrt(sz*ems*ems/n_lambda)        !sigma(x)=sigma(y)
      ro=sqrt(3.)/2.                  !correlation coefficient
    2 CONTINUE
      sigmatt=sigmat
      sigmaxx=sigmax
      roo=ro
c      CALL normco(thetaxx,xxx,0.,0.,sigmatt,sigmaxx,roo) !simulation of thetax,
c      CALL normco(thetayy,yyy,0.,0.,sigmatt,sigmaxx,roo) !thetay,xx,yy
	vv(1,1)=sigmatt*sigmatt
	vv(2,2)=sigmaxx*sigmaxx
	vv(1,2)=roo*sigmatt*sigmaxx
	vv(2,1)=roo*sigmatt*sigmaxx
	call corset(vv,cc,2)
	call corgen(cc,vcx,2)
	thetaxx=vcx(1)
	xxx=vcx(2)
	call corgen(cc,vcx,2)
	thetayy=vcx(1)
	yyy=vcx(2)
c      thetax=conv*thetaxx
c      thetay=conv*thetayy
      thetax=thetaxx
      thetay=thetayy
cc      IF(dabs(thetax).GE.pi/2..OR.dabs(thetay).GE.pi/2.) go to 2
      theta=datan(dsqrt(dtan(thetax)*dtan(thetax)+     !evaluation of theta
     *  dtan(thetay)*dtan(thetay)))                     !using thetax,thetay
      IF(dsin(thetax).NE.0.) THEN
        phi=datan(dtan(thetay)/dtan(thetax))            !evaluation of phi
        IF(dtan(thetax).LT.0.) phi=pi+phi             !using thetax,thetay
      END IF
      IF(dsin(thetax).EQ.0.) THEN
        phi=pi/2.
        IF(dtan(thetay).LT.0.) phi=pi*3./2.
      END IF
      IF(phi.GT.pi*2.) phi=-pi*2.+phi
      IF(phi.LT.0.) phi=phi+pi*2.
      deltax=xxx
      deltay=yyy
      deltar=sigmat*sigmat*x0

      RETURN
      END

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	SUBROUTINE coord_transform_tau_sw(x,y,z,x1,y1,z1,theta1,phi1)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IMPLICIT real*8 (a-h,o-z)
	parameter (pi=3.141592654d0)	
	phi2=phi1
	phi=phi2
	theta=theta1
	t11=dcos(theta)*dcos(phi)*dcos(phi)+dsin(phi)*dsin(phi)
	t21=dcos(theta)*dcos(phi)*dsin(phi)-dcos(phi)*dsin(phi)
	t31=-dsin(theta)*dcos(phi)
	t12=dcos(theta)*dcos(phi)*dsin(phi)-dcos(phi)*dsin(phi)
	t22=dcos(theta)*dsin(phi)*dsin(phi)+dcos(phi)*dcos(phi)
	t32=-dsin(theta)*dsin(phi)
	t13=dsin(theta)*dcos(phi)
	t23=dsin(theta)*dsin(phi)
	t33=dcos(theta)
	z=t31*x1+t32*y1+t33*z1
	x=t11*x1+t12*y1+t13*z1
	y=t21*x1+t22*y1+t23*z1
	RETURN
	END

*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	SUBROUTINE angle_transform_tau_sw(theta,phi,theta1,phi1)
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	IMPLICIT real*8 (a-h,o-z)
	parameter (pi=3.141592654d0)
	if(theta.eq.0.d0) then
	theta=theta1
	phi=phi1
	return
	end if
	t11=dcos(theta)*dcos(phi)*dcos(phi)+dsin(phi)*dsin(phi)
	t21=dcos(theta)*dcos(phi)*dsin(phi)-dcos(phi)*dsin(phi)
	t31=-dsin(theta)*dcos(phi)
	t12=dcos(theta)*dcos(phi)*dsin(phi)-dcos(phi)*dsin(phi)
	t22=dcos(theta)*dsin(phi)*dsin(phi)+dcos(phi)*dcos(phi)
	t32=-dsin(theta)*dsin(phi)
	t13=dsin(theta)*dcos(phi)
	t23=dsin(theta)*dsin(phi)
	t33=dcos(theta)
	x2=dsin(theta1)*dcos(phi1)
	y2=dsin(theta1)*dsin(phi1)
	z2=dcos(theta1)
	z21=t31*x2+t32*y2+t33*z2
	x21=t11*x2+t12*y2+t13*z2
	y21=t21*x2+t22*y2+t23*z2
	if(z21.gt.1.d0) z21=1.d0
	if(z21.lt.-1.d0) z21=-1.d0
	theta=dacos(z21)
	IF(x21.NE.0.d0) THEN
	phi=datan(y21/x21)
	IF(x21.LT.0.d0) phi=pi+phi
	END IF
	IF(x21.EQ.0.d0) THEN
	phi=pi/2.
	IF(y21.LT.0.d0) phi=1.5*pi
	END IF
	IF(phi.GT.2.*pi) phi=-2.*pi+phi
	IF(phi.LT.0.d0) phi=phi+2.*pi
	RETURN
	END
