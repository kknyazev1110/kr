      implicit real*8(a-b,d-h,o-z),integer*4(i-n),complex*16(c)
C propnharm.h
      parameter(nx=160)
      parameter(ny=16)
      parameter(nphi=17)
      parameter(nyx=nx*ny)
      parameter(nxyphi=nx*ny*nphi)
      dimension r(0:nx+1), y(1:ny), tetta(1:ny),
     &   x(0:nx+1), rpr(0:nx+1), cpropall(0:nphi,nyx,nyx), phi(1:nphi),
     & w(1:2),etta(1:2),tau(1:2),npuls(1:2),t_pau(1:2),
     & cpsi0(1:nyx),cpsi(1:nxyphi),
     & dip(1:3,0:82000), bz(1:3,0:82000),
     & cd(1:3), cbz(1:3),
     & cua(1:nyx,1:2),cub(1:nyx,1:2),
     & cnpsi(1:nphi),cnpsim(1:nyx,-(nphi-1)/2:(nphi-1)/2),
     & ddz(0:63),ddx(0:63),ddy(0:63),
     & daz(0:63),dax(0:63),day(0:63)
      allocatable :: cprop(:,:),a(:,:)
     

      
