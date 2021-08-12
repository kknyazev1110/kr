      implicit real*8(a-b,d-h,o-z),integer*4(i-n),complex*16(c)
      parameter(nx=160)
      parameter(ny=24)
      parameter(nphi=17)
      parameter(nxy=nx*ny)
      parameter(nxyphi=nx*ny*nphi)
      dimension x(0:nx+1), y(0:ny+1), xi(0:nx+1), eta(0:ny+1), 
     &          derxi(0:nx+1), dereta(0:ny+1),
     &          dx(1:nx+1,1:nx+1), dy(1:ny,1:ny),
     &          dermatx(1:nx+1,1:nx+1), dermaty(1:ny,1:ny),
     &          tii(1:nx,1:nx), tjj(1:ny,1:ny),
     &          tkk(1:nphi,1:nphi), h0(1:nxy,1:nxy), 
     &          val(1:nxy), cpsi0(1:nxy)
      allocatable :: work(:),ce(:,:),cpsim(:,:),cprop(:,:),cp(:,:),
     &   cvec(:,:),z(:,:)
