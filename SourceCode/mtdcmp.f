      subroutine mtdcmp(mxx,myy,mzz,mxy,myz,mzx,poisson,
     &            miso,mten,mdc,eig,plg,azm,
     &            st1,di1,ra1,st2,di2,ra2,
     &            beachball,nb,mb,ierr)
      implicit none
c
c     by
c     Rongjiang Wang
c     Helmholtz Centre Potsdam
c     GFZ German Research Centre for Geosciences
c     Potsdam, April 20, 2020
c
c     This subroutine is used to decompose full moment tensor into the tensile,
c     double-couple and (residual) isotropic parts
c
c     Inputs:
c     mxx,myy,mzz,mxy,myz,mzx = full moment tensor
c     x = northward, y = eastward, z = downward
c
      real*8 mxx,myy,mzz,mxy,myz,mzx,poisson
c
c     Outputs:
c     miso,mten,mdc = (residual) isotropic, tensile and double-couple moment
c     eig,plg,azm = eigenvalue, plunge and azimuth of tension, nodal and compression axes
c     st1,di1,ra1[deg] = strike, dieigp and rake [deg] of the 1. double-couple
c                        focal mechanism with the opening (closing) component
c     st2,di2,ra2[deg] = strike, dip and rake [deg] of the 2. double-couple
c                        focal mechanism.
c     beach(nb,mb,2) = upper and lower hemispherical beachball
c
c     ierr = 0 number of errors (> 0 means eigabnormal return)
c
      integer*4 ierr,nb,mb
      real*8 miso,mten,mdc
      real*8 eig(3),plg(3),azm(3)
      real*8 st1,di1,ra1,st2,di2,ra2
      character*1 beachball(nb,mb,2)
c
c     local memories:
c
      integer*4 i,j,k,l,xt,yt,xp,yp,xn,yn,ibb
      real*8 pi,rad2deg
      real*8 b,c,d,eig1,eig2,eig3,alpha,beta,theta
      real*8 am,det1,det2,det3,tam,pam,nam
      real*8 x,y,dx,dy,azi,tkf
      real*8 st(2),di(2),ra(2),rst(3),rdi(3)
      real*8 mt(3,3),r(3,3),ns(3,2),ts(3,2),swap(3)
c
      real*8 eps
      data eps/1.0d-12/
c
      miso=0.d0
      mten=0.d0
      mdc=0.d0
c
      st1=0.d0
      di1=0.d0
      ra1=0.d0
c
      st2=0.d0
      di2=0.d0
      ra2=0.d0
c
      do j=1,3
        plg(j)=0.d0
        azm(j)=0.d0
        azm(j)=0.d0
      enddo
c
      if(poisson.eq.0.d0.or.poisson.eq.0.5d0)then
        ierr=1
        return
      endif
c
      ierr=0
      pi=4.d0*datan(1.d0)
      rad2deg=180.d0/pi
c
      if(mxy.eq.0.d0.and.myz.eq.0.d0.and.mzx.eq.0.d0)then
        eig(1)=mxx
        eig(2)=myy
        eig(3)=mzz
      else
        am=dsqrt(0.5d0*(mxx**2+myy**2+mzz**2)
     &                 +mxy**2+myz**2+mzx**2)
        b=-(mxx+myy+mzz)/am
        c=(mxx*myy+myy*mzz+mzz*mxx-mxy**2-myz**2-mzx**2)/am**2
        d=(mxx*myz**2+myy*mzx**2+mzz*mxy**2
     &    -2.d0*mxy*myz*mzx-mxx*myy*mzz)/am**3
        call rootm3(b,c,d,eig)
        do j=1,3
          eig(j)=eig(j)*am
        enddo
      endif
c
      eig1=dmax1(eig(1),eig(2),eig(3))
      eig2=dmin1(eig(1),eig(2),eig(3))
      eig3=eig(1)+eig(2)+eig(3)-eig1-eig2
c
      eig(1)=eig1
      eig(2)=eig2
      eig(3)=eig3
c
      if(dabs(eig1-eig2).le.eps*dsqrt(eig1**2+eig2**2+eig3**2))then
c
c       (nearly) isotropic stress state
c
        miso=eig(1)*dsqrt(1.5d0)
        if(eig(1).lt.0.d0)then
          eig(1)=0.d0
        else
          eig(2)=0.d0
        endif
        eig(3)=0.d0
c
        dx=pi/dble(nb)
        dy=pi/dble(mb)
        do i=1,nb
          x=0.5d0*pi-(dble(i)-0.5d0)*dx
          do j=1,mb
            y=(dble(j)-0.5d0)*dy-0.5d0*pi
            if(dsqrt(x*x+y*y).gt.0.5d0*pi)then
              beachball(i,j,1)=' '
              beachball(i,j,2)=' '
            else if(miso.ge.0.d0)then
              beachball(i,j,1)='#'
              beachball(i,j,2)='#'
            else
              beachball(i,j,1)='-'
              beachball(i,j,2)='-'
            endif
          enddo
        enddo
        return
      endif
c
c     alpha = lambda/mue = poisson/(0.5-poisson)
c
      alpha=poisson/(0.5d0-poisson)
c
      miso=(1.d0+alpha)*eig3-0.5d0*alpha*(eig1+eig2)
      mten=0.5d0*(eig1+eig2)-eig3
      mdc =dsqrt(dabs((eig1-eig3)*(eig3-eig2)))
c
c     determine eigenvectors of the max. and min. eigenvalueseig (principal stresses)
c
      do j=1,2
        if(dabs(eig(j)-eig(3)).le.
     &     eps*dsqrt(eig1**2+eig2**2+eig3**2))then
c
c         two of the three principal stresses are identical
c         -> more than a pair of optimal oriented fault planes
c         here any arbitrary pair is chosen
c
          det1=dmax1(dabs(mxx-eig(j)),dabs(mxy),dabs(mzx))
          det2=dmax1(dabs(mxy),dabs(myy-eig(j)),dabs(myz))
          det3=dmax1(dabs(mzx),dabs(myz),dabs(mzz-eig(j)))
          am=dmax1(det1,det2,det3)
          if(det1.ge.am)then
            if(dabs(mxx-eig(j)).ge.am)then
              r(1,j)=-(mxy+mzx)
              r(2,j)=mxx-eig(j)
              r(3,j)=mxx-eig(j)
            else if(dabs(mxy).ge.am)then
              r(1,j)=mxy
              r(2,j)=-(mxx-eig(j)+mzx)
              r(3,j)=mxy
            else
              r(1,j)=mzx
              r(2,j)=-(mxx-eig(j)+mxy)
              r(3,j)=mxy
            endif
          else if(det2.ge.am)then
            if(dabs(mxy).ge.am)then
              r(1,j)=-(myy-eig(j)+myz)
              r(2,j)=mxy
              r(3,j)=mxy
            else if(dabs(myy-eig(j)).ge.am)then
              r(1,j)=myy-eig(j)
              r(2,j)=-(mxy+myz)
              r(3,j)=myy-eig(j)
            else
              r(1,j)=myz
              r(2,j)=myz
              r(3,j)=-(myy-eig(j)+mxy)
            endif
          else
            if(dabs(mzx).ge.am)then
              r(1,j)=-(mzz-eig(j)+myz)
              r(2,j)=mzx
              r(3,j)=mzx
            else if(dabs(myz).ge.am)then
              r(1,j)=myz
              r(2,j)=-(mzz-eig(j)+mzx)
              r(3,j)=myz
            else
              r(1,j)=mzz-eig(j)
              r(2,j)=mzz-eig(j)
              r(3,j)=-(mzx+myz)
            endif
          endif
        else
          det1=myz*myz-(myy-eig(j))*(mzz-eig(j))
          det2=mzx*mzx-(mxx-eig(j))*(mzz-eig(j))
          det3=mxy*mxy-(mxx-eig(j))*(myy-eig(j))
          am=dmax1(dabs(det1),dabs(det2),dabs(det3))
          if(dabs(det1).ge.am)then
            r(1,j)=det1
            r(2,j)=(mzz-eig(j))*mxy-myz*mzx
            r(3,j)=(myy-eig(j))*mzx-myz*mxy
          else if(dabs(det2).ge.am)then
            r(1,j)=(mzz-eig(j))*mxy-mzx*myz
            r(2,j)=det2
            r(3,j)=(mxx-eig(j))*myz-mzx*mxy
          else
            r(1,j)=(myy-eig(j))*mzx-mxy*myz
            r(2,j)=(mxx-eig(j))*myz-mxy*mzx
            r(3,j)=det3
          endif
        endif
        am=dsqrt(r(1,j)**2+r(2,j)**2+r(3,j)**2)
        do i=1,3
          r(i,j)=r(i,j)/am
        enddo
      enddo
c
      beta=(eig1+eig2-2.d0*eig3)/(eig1-eig2)
c
      if(beta.le.0.d0)then
        theta=0.5d0*dacos(-beta)
      else
        theta=0.5d0*(pi-dacos(beta))
      endif
c
c     determine the two optimal fault-plane normals, ns,
c     their unit vectors in the rake direction, ts.
c
      do i=1,3
        ts(i,1)=r(i,1)*dcos(theta)+r(i,2)*dsin(theta)
        ns(i,1)=r(i,1)*dsin(theta)-r(i,2)*dcos(theta)
c
        ts(i,2)=r(i,1)*dsin(theta)-r(i,2)*dcos(theta)
        ns(i,2)=r(i,1)*dcos(theta)+r(i,2)*dsin(theta)
      enddo
c
      do j=1,2
        if(ns(3,j).gt.0.d0)then
          do i=1,3
            ns(i,j)=-ns(i,j)
            ts(i,j)=-ts(i,j)
          enddo
        endif
      enddo
c
c     determine strike, dip and rake
c
      do j=1,2
        st(j)=datan2(ns(2,j),ns(1,j))-0.5d0*pi
        di(j)=datan2(dsqrt(ns(1,j)**2+ns(2,j)**2),dabs(ns(3,j)))
c
c       rst = unit vector in the strike direction
c       rdi = unit vector in the down-dip direction
c
        rst(1)=dcos(st(j))
        rst(2)=dsin(st(j))
        rst(3)=0.d0
c
        rdi(1)=dcos(di(j))*dcos(st(j)+0.5d0*pi)
        rdi(2)=dcos(di(j))*dsin(st(j)+0.5d0*pi)
        rdi(3)=dsin(di(j))
c
        ra(j)=datan2(-(ts(1,j)*rdi(1)+ts(2,j)*rdi(2)+ts(3,j)*rdi(3)),
     &                (ts(1,j)*rst(1)+ts(2,j)*rst(2)+ts(3,j)*rst(3)))
      enddo
c
      st1=dmod(st(1)*rad2deg+360.d0,360.d0)
      di1=di(1)*rad2deg
      ra1=dmod(ra(1)*rad2deg+360.d0,360.d0)
      if(ra1.gt.180.d0)ra1=ra1-360.d0
c
      st2=dmod(st(2)*rad2deg+360.d0,360.d0)
      di2=di(2)*rad2deg
      ra2=dmod(ra(2)*rad2deg+360.d0,360.d0)
      if(ra2.gt.180.d0)ra2=ra2-360.d0
c
c     moment amplitude
c
      miso=dsqrt(1.5d0)*miso
      mten=dsqrt(2.d0+2.d0*alpha+1.5d0*alpha**2)*mten
c
      r(1,3)= r(2,1)*r(3,2)-r(3,1)*r(2,2)
      r(2,3)=-r(1,1)*r(3,2)+r(3,1)*r(1,2)
      r(3,3)= r(1,1)*r(2,2)-r(2,1)*r(1,2)
      am=dsqrt(r(1,3)**2+r(2,3)**2+r(3,3)**2)
      do i=1,3
        r(i,3)=r(i,3)/am
      enddo
c
      do j=1,3
        if(r(3,j).lt.0.d0)then
          do i=1,3
            r(i,j)=-r(i,j)
          enddo
        endif
        plg(j)=datan2(dabs(r(3,j)),dsqrt(r(1,j)**2+r(2,j)**2))*rad2deg
        azm(j)=datan2(r(2,j),r(1,j))
        azm(j)=dmod(azm(j)*rad2deg+360.d0,360.d0)
      enddo
c
c     beachball representation
c
      mt(1,1)=mxx
      mt(1,2)=mxy
      mt(1,3)=mzx
      mt(2,1)=mxy
      mt(2,2)=myy
      mt(2,3)=myz
      mt(3,1)=mzx
      mt(3,2)=myz
      mt(3,3)=mzz
c
      dx=pi/dble(nb)
      dy=pi/dble(mb)
      tam=5.d0
      nam=5.d0
      pam=5.d0
      xt=0
      yt=0
      xp=0
      yp=0
      xn=0
      yn=0
      do i=1,nb
        x=0.5d0*pi-(dble(i)-0.5d0)*dx
        do j=1,mb
          y=(dble(j)-0.5d0)*dy-0.5d0*pi
          if(dsqrt(x*x+y*y).gt.0.5d0*pi)then
            beachball(i,j,2)=' '
          else
            azi=datan2(y,x)
            tkf=dsqrt(x*x+y*y)
            rdi(1)=dsin(tkf)*dcos(azi)
            rdi(2)=dsin(tkf)*dsin(azi)
            rdi(3)=(-1.d0)**ibb*dcos(tkf)
c
            am=0.d0
            do k=1,3
              am=am+(r(k,1)-rdi(k))**2
            enddo
            if(am.le.tam)then
              xt=i
              yt=j
              tam=am
            endif
c
            am=0.d0
            do k=1,3
              am=am+(r(k,2)-rdi(k))**2
            enddo
            if(am.le.pam)then
              xp=i
              yp=j
              pam=am
            endif
c
            am=0.d0
            do k=1,3
              am=am+(r(k,3)-rdi(k))**2
            enddo
            if(am.le.nam)then
              xn=i
              yn=j
              nam=am
            endif
c
            am=0.d0
            do k=1,3
              swap(k)=0.d0
              do l=1,3
                swap(k)=swap(k)+mt(k,l)*rdi(l)
              enddo
              am=am+swap(k)*rdi(k)
            enddo
            if(am.gt.0.d0)then
              beachball(i,j,2)='#'
            else
              beachball(i,j,2)='-'
            endif
          endif
        enddo
      enddo
c
      if(xt*yt.gt.0)then
        beachball(xt,yt,2)='T'
        if(yt.gt.1 )beachball(xt,yt-1,2)='('
        if(yt.lt.mb)beachball(xt,yt+1,2)=')'
      endif
c
      if(xp*yp.gt.0)then
        beachball(xp,yp,2)='P'
        if(yp.gt.1 )beachball(xp,yp-1,2)='('
        if(yp.lt.mb)beachball(xp,yp+1,2)=')'
      endif
c
      if(xn*yn.gt.0)then
        beachball(xn,yn,2)='N'
        if(yn.gt.1 )beachball(xn,yn-1,2)='('
        if(yn.lt.mb)beachball(xn,yn+1,2)=')'
      endif
c
      k=0
      do i=nb,1,-1
        k=k+1
        l=0
        do j=mb,1,-1
          l=l+1
          beachball(k,l,1)=beachball(i,j,2)
          if(beachball(k,l,1).eq.'(')then
            beachball(k,l,1)=')'
          else if(beachball(k,l,1).eq.')')then
            beachball(k,l,1)='('
          endif
        enddo
      enddo
c
      return
      end