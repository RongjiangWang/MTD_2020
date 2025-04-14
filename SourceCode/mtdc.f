      subroutine mtdc(strike,dip,rake,m0,mxx,myy,mzz,mxy,myz,mzx,ierr)
      implicit none
c
c     Inputs:
c     fault strike and dip, and slip rake angles in degree, and seismic moment
c
      real*8 strike,dip,rake,m0
c
c     Return:
c     unified double-couple moment tensor: mxx,myy,mzz,mxy,myz,mzx
c     (x = northward, y = eastward, z = downward)
c     ierr = number of errors (> 0 means abnormal return)
c
      integer*4 ierr
      real*8 mxx,myy,mzz,mxy,myz,mzx
c
c     local memory
c
      real*8 st,di,ra,pi,deg2rad
      real*8 sss,sss2,ss2s,css,cs2s,css2
      real*8 ssd,ssd2,ss2d,csd,cs2d,csd2
      real*8 ssr,ssr2,ss2r,csr,csr2,cs2r
c
      ierr=0
      pi=4.d0*datan(1.d0)
      deg2rad=pi/180.d0
c
      st=strike*deg2rad
      di=dip*deg2rad
      ra=rake*deg2rad
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
      ssr=dsin(ra)
      ssr2=ssr*ssr
      ss2r=dsin(2.d0*ra)
c
      csr=dcos(ra)
      csr2=csr*csr
      cs2r=dcos(2.d0*ra)
c
      mxx=(-ssd*csr*ss2s-ss2d*ssr*sss2)*m0
      myy=(ssd*csr*ss2s-ss2d*ssr*css2)*m0
      mzz=-(mxx+myy)
      mxy=(ssd*csr*cs2s+0.5d0*ss2d*ssr*ss2s)*m0
      myz=(-csd*csr*sss+cs2d*ssr*css)*m0
      mzx=(-csd*csr*css-cs2d*ssr*sss)*m0
c
      return
      end