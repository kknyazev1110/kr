      implicit real*8(a-b,d-h,o-z),integer*4(i-n),complex*16(c)
C kinmat.h
      parameter(nx=160)
      parameter(ny=16)
      parameter(nphi=17)
      parameter(nxy=nx*ny)
      parameter(nxyphi=nx*ny*nphi)
C nx - MAXIMUM NUMBER OF GRID POINTS IN PSEUDORADIAL  GRID
C ny - MAXIMUM NUMBER OF GRID POINTS IN PSEUDOANGULAR GRID
      dimension r(0:nx+1), y(1:ny), tetta(1:ny),
     &   x(0:nx+1), rpr(0:nx+1),
     &   tii(1:nx,1:nx), tjj0(1:ny,1:ny),tjj1(1:ny,1:ny),
     &   cpsi0(1:nxy),val(1:nxy),h0(1:nxy,1:nxy)
      allocatable :: work(:),ce(:,:),cpsim(:,:),cprop(:,:),cp(:,:),
     &   cvec(:,:),a(:,:)


