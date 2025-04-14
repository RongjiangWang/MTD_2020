      program mtdmain
      implicit none
c
      integer*4 nb,mb
      parameter(nb=17,mb=35)
      integer*4 i,j,ics,idcmp,iexpo,ierr
      real*8 mxx,myy,mzz,mxy,myz,mzx
      real*8 xmxx,xmyy,xmzz
      real*8 omxx,omyy,omzz,omxy,omyz,omzx
      real*8 smxx,smyy,smzz,smxy,smyz,smzx
      real*8 st1,di1,ra1,st2,di2,ra2
      real*8 m0,mw,miso,mten,mdc,poisson,mscale
      real*8 eig(3),plg(3),azm(3)
      character*1 beachball(nb,mb,2)
c
      real*8 eps
      data eps/1.0d-12/
c
      write(*,'(a)')'*************************************************'
     &            //'***********************'
      write(*,'(a)')'*       Moment Tensor Decomposition (MTD)'
      write(*,'(a)')'*       into isotropic, tensile/CLVD and DC parts'
      write(*,'(a)')'*       by'
      write(*,'(a)')'*       Rongjiang Wang'
      write(*,'(a)')'*       Helmholtz Centre Potsdam'
      write(*,'(a)')'*       GFZ German Research Centre for Geosciences'
      write(*,'(a)')'*       Potsdam, last updated on April 24, 2020'
      write(*,'(a)')'*       wang@gfz-potsdam.de'
      write(*,'(a)')'*************************************************'
     &            //'***********************'
      write(*,'(a)')' Inputs:'
      write(*,'(a)')' ======='
      write(*,'(a)')' Conventions of coordinate system'
      write(*,'(a)')'   0 = Spherical'
     &            //' (1/2/3 = Radial/co-latitude Theta/longitude Phi)'
      write(*,'(a)')'   1 = Aki'//char(39)//'s'
     &      //' (1/2/3 = North/East/Down)'
      write(*,'(a)')'   2 = Geodetic (1/2/3 = East/North/Up)'
      write(*,'(a,$)')' Choose the convention of coordinate'
     &              //' system (0, 1 or 2): '
      read(*,*)ics
      write(*,'(a,$)')' Type the scale of moment in [Nm]: '
      read(*,*)mscale
      write(*,'(a)')' Type the 6 moment tensor elements'
      if(ics.eq.0)then
        write(*,'(a,$)')' Mrr, Mtt, Mpp, Mrt, Mrp, Mtp = '
        read(*,*)mzz,mxx,myy,mzx,myz,mxy
        myz=-myz
        mxy=-mxy
      else if(ics.eq.1)then
        write(*,'(a,$)')' Mnn, Mee, Mdd, Mne, Mnd, Med = '
        read(*,*)mxx,myy,mzz,mxy,mzx,myz
      else if(ics.eq.2)then
        write(*,'(a,$)')' Mee, Mnn, Muu, Men, Meu, Mnu = '
        read(*,*)myy,mxx,mzz,mxy,myz,mzx
        myz=-myz
        mzx=-mzx
      else
        stop ' Wrong choice of convention of coordinate system!'
      endif
      mxx=mxx*mscale
      myy=myy*mscale
      mzz=mzz*mscale
      mxy=mxy*mscale
      myz=myz*mscale
      mzx=mzx*mscale
      write(*,'(a,$)')' Choose option tensile (1) or CLVD (2): '
      read(*,*)idcmp
      if(idcmp.eq.1)then
        write(*,'(a,$)')' Type Poisson ratio of medium'
     &                //' (> 0 and < 0.5): '
        read(*,*)poisson
        if(poisson.le.0.d0.or.poisson.ge.0.5d0)then
          stop ' Wrong Poisson ratio!'
        endif
      else
c
c       tensile => CLVD if lambda/mue is set to -2/3 (Poisson = -1)
c
        poisson=-1.d0
      endif
      
      call mtdcmp(mxx,myy,mzz,mxy,myz,mzx,poisson,
     &            miso,mten,mdc,eig,plg,azm,
     &            st1,di1,ra1,st2,di2,ra2,
     &            beachball,nb,mb,ierr)
c
      m0=dsqrt(0.5d0*(mxx**2+myy**2+mzz**2)+mxy**2+myz**2+mzx**2)
      mw=(dlog10(m0)-9.1d0)/1.5d0
c
      write(*,'(a)')'*************************************************'
     &            //'***********************'
      write(*,'(a)')' Outputs:'
      write(*,'(a)')' ========'
      write(*,'(a,E12.4,a,f5.1,a)')'        Total seismic moment = ',m0
     &                     ,' Nm (Mw = ',mw,')'
      write(*,'(a)')'        ----------------------------------'
     &            //'-----------------'
      if(idcmp.eq.1)then
        write(*,'(a,E12.4,a)')'        (Residual) isotropic = ',
     &                      miso,' Nm'

        write(*,'(a,E12.4,a)')'                     Tensile = ',
     &                      mten,' Nm'
      else
        miso=(mxx+myy+mzz)*dsqrt(1.5d0)/3.d0
        write(*,'(a,E12.4,a)')'                   Isotropic = ',
     &                      miso,' Nm'
        write(*,'(a,E12.4,a)')'                        CLVD = ',
     &                      mten,' Nm'
      endif
      write(*,'(a,E12.4,a)')'               Double-couple = ',mdc,' Nm'
      write(*,'(a)')'        ----------------------------------'
     &            //'-----------------'
      write(*,'(a)')' Principal axes'
      write(*,'(a,E12.4,2(a,f7.2))')'   T      Val[Nm] = ',eig(1),
     &              '  Plg[deg] = ',plg(1),'  Azm[deg] = ',azm(1)
      write(*,'(a,E12.4,2(a,f7.2))')'   N                ',eig(3),
     &              '             ',plg(3),'             ',azm(3)
      write(*,'(a,E12.4,2(a,f7.2))')'   P                ',eig(2),
     &              '             ',plg(2),'             ',azm(2)
      write(*,'(a)')'        ---------------------------------------'
      if(idcmp.eq.1)then
        write(*,'(a)')'        NP1 (fault with slip and tensile)'
        write(*,'(a,3f9.2)')'           strike, dip and rake [deg] = '
     &        //'     ',st1,di1,ra1
        write(*,'(a)')'        NP2 (      conjugate nodal plane)'
        write(*,'(a,3f9.2)')'           strike, dip and rake [deg] = '
     &        //'     ',st2,di2,ra2
      else
        write(*,'(a)')'        NP1 (fault with slip and CLVD)'
        write(*,'(a,3f9.2)')'           strike, dip and rake [deg] = '
     &        //'     ',st2,di2,ra2
        write(*,'(a)')'        NP2 (   conjugate nodal plane)'
        write(*,'(a,3f9.2)')'           strike, dip and rake [deg] = '
     &        //'     ',st1,di1,ra1
      endif
c
c     beachball plot
c
      write(*,'(a)')' '
      write(*,'(a)')'        Upper and lower hemispherical beachball'
      write(*,'(a)')'        ---------------------------------------'
      write(*,'(a)')' '
      do i=1,nb
        write(*,'(a,$)')' '
        do j=1,mb
          write(*,'(a1,$)')beachball(i,j,1)
        enddo
        write(*,'(a,$)')'  '
        do j=1,mb
          write(*,'(a1,$)')beachball(i,j,2)
        enddo
        write(*,'(a)')' '
      enddo
      write(*,'(a)')' '
c
c     test
c
      miso=miso/dsqrt(1.5d0)
c
      call mtop(st1,di1,poisson,mten,omxx,omyy,omzz,omxy,omyz,omzx,ierr)
      call mtdc(st1,di1,ra1,mdc,smxx,smyy,smzz,smxy,smyz,smzx,ierr)
c
      mxx=miso+omxx+smxx
      myy=miso+omyy+smyy
      mzz=miso+omzz+smzz
      mxy=omxy+smxy
      myz=omyz+smyz
      mzx=omzx+smzx
c
      mscale=dmax1(dabs(mxx),dabs(myy),dabs(mzz),
     &             dabs(mxy),dabs(myz),dabs(mzx))
      iexpo=max0(0,idint(dlog10(mscale)))
      mscale=10.d0**iexpo
c
      miso=miso/mscale
c
      omxx=omxx/mscale
      omyy=omyy/mscale
      omzz=omzz/mscale
      omxy=omxy/mscale
      omyz=omyz/mscale
      omzx=omzx/mscale
c
      smxx=smxx/mscale
      smyy=smyy/mscale
      smzz=smzz/mscale
      smxy=smxy/mscale
      smyz=smyz/mscale
      smzx=smzx/mscale
c
      mxx=mxx/mscale
      myy=myy/mscale
      mzz=mzz/mscale
      mxy=mxy/mscale
      myz=myz/mscale
      mzx=mzx/mscale
c
      write(*,'(a)')'=================================================='
     &            //'======================'
      write(*,'(a,$)')' Decomposed moment tensors (scale = 10**'
      if(iexpo.ge.10.and.iexpo.le.99)then
        write(*,'(i2,$)')iexpo
      else if(iexpo.ge.0.and.iexpo.le.9)then
        write(*,'(a1,i1,$)')'0',iexpo
      else
      endif
      write(*,'(a)')' Nm)'
      write(*,'(a)')'=================================================='
     &            //'======================'
c
      write(*,'(a,$)')'   Mechanism'
      if(ics.eq.0)then
        write(*,'(a)')'       Mrr       Mtt       Mpp'
     &              //'       Mrt       Mrp       Mtp'
        write(*,'(a,6f10.6)')'   Isotropic',miso,miso,miso,
     &                       0.d0,0.d0,0.d0
        if(idcmp.eq.1)then
          write(*,'(a,$)')'     Tensile'
        else
          write(*,'(a,$)')'        CLVD'
        endif
        write(*,'(6f10.6)')omzz,omxx,omyy,omzx,-omyz,-omxy
        
        write(*,'(a,6f10.6)')'          DC',
     &                       smzz,smxx,smyy,smzx,-smyz,-smxy
        write(*,'(a,6f10.6)')'       Total',mzz,mxx,myy,mzx,-myz,-mxy
      else if(ics.eq.1)then
        write(*,'(a)')'       Mnn       Mee       Mdd'
     &              //'       Mne       Mnd       Med'
        write(*,'(a,6f10.6)')'   Isotropic',miso,miso,miso,
     &                       0.d0,0.d0,0.d0
        if(idcmp.eq.1)then
          write(*,'(a,$)')'     Tensile'
        else
          write(*,'(a,$)')'        CLVD'
        endif
        write(*,'(6f10.6)')omxx,omyy,omzz,omxy,omzx,omyz
        
        write(*,'(a,6f10.6)')'          DC',
     &                       smxx,smyy,smzz,smxy,smzx,smyz
        write(*,'(a,6f10.6)')'       Total',mxx,myy,mzz,mxy,mzx,myz
      else if(ics.eq.2)then
        write(*,'(a)')'       Mee       Mnn       Muu'
     &              //'       Men       Meu       Mnu'
        write(*,'(a,6f10.6)')'   Isotropic',miso,miso,miso,
     &                       0.d0,0.d0,0.d0
        if(idcmp.eq.1)then
          write(*,'(a,$)')'     Tensile'
        else
          write(*,'(a,$)')'        CLVD'
        endif
        write(*,'(6f10.6)')omyy,omxx,omzz,omxy,-omyz,-omzx
        
        write(*,'(a,6f10.6)')'          DC',
     &                       smyy,smxx,smzz,smxy,-smyz,-smzx
        write(*,'(a,6f10.6)')'       Total',myy,mxx,mzz,mxy,-myz,-mzx
      endif
      write(*,'(a)')'*************************************************'
     &            //'***********************'
      stop
      end
      
      