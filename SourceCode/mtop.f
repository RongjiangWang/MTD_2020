      subroutine mtop(strike,dip,poisson,m0,
     &                mxx,myy,mzz,mxy,myz,mzx,ierr)
      implicit none
c
c     Inputs:
c     fault strike and dip angles in degree, poisson ratio and seismic moment
c
      real*8 strike,dip,poisson,m0
c
c     Return:
c     unified opening moment tensor: mxx,myy,mzz,mxy,myz,mzx
c     (x = northward, y = eastward, z = downward)
c     ierr = number of errors (> 0 means abnormal return)
c
      integer*4 ierr
      real*8 mxx,myy,mzz,mxy,myz,mzx
c
c     local memory
c
      real*8 st,di
      real*8 norm,alpha,pi,deg2rad
      real*8 sss,sss2,ss2s,css,cs2s,css2
      real*8 ssd,ssd2,ss2d,csd,cs2d,csd2
c
      if(poisson.eq.0.5d0)then
        mxx=0.d0
        myy=0.d0
        mzz=0.d0
        mxy=0.d0
        myz=0.d0
        mzx=0.d0
c
        ierr=1
        return
      endif
c
      ierr=0
      pi=4.d0*datan(1.d0)
      deg2rad=pi/180.d0
      alpha=poisson/(0.5d0-poisson)
c
      st=strike*deg2rad
      di=dip*deg2rad
c
      sss=dsin(st)
      sss2=sss*sss
      ss2s=dsin(2.d0*st)
c
      css=dcos(st)
      css2=css*css
      cs2s=dcos(2.d0*st)
c
      ssd=dsin(di)
      ssd2=ssd*ssd
      ss2d=dsin(2.d0*di)
c
      csd=dcos(di)
      csd2=csd*csd
      cs2d=dcos(2.d0*di)
c
      norm=m0/dsqrt(2.d0+2.d0*alpha+1.5d0*alpha**2)
c
      mxx=(alpha+2.d0*sss2*ssd2)*norm
      myy=(alpha+2.d0*css2*ssd2)*norm
      mzz=(alpha+2.d0*csd2)*norm
      mxy=-ss2s*ssd2*norm
      myz=-css*ss2d*norm
      mzx=sss*ss2d*norm
c
      return
      end