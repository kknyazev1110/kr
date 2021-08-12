      program eigen
      
      include 'eigen.h'
      
      pi=4.d0*datan(1.d0) 
      n=10
      m=20

      open(4,file='fieldparams.dat',status='old')
      read(4,*)ndt
      read(4,*)ndw
      read(4,*)o_au
      read(4,*)f_eau
      read(4,*)s_l
      read(4,*)alpha
      read(4,*)Rm
      read(4,*)Rb
      read(4,*)R0
      read(4,*)a
      read(4,*)za

      lr10=nxy*4      
      open(unit=10,file='kinmat.mat',form='unformatted',
     &    access='direct',recl=lr10,status='replace')
      write(*,*)a
c  LEGENDRE-GAUSS-RADAU for xi
      call radau(nx,x(1))
      x(0)=-1.d0
      x(nx+1)=1.d0
      
      do i=0,nx+1
      call xtoxi(x(i),Rm,Rb,a,res)
      xi(i)=res
c      write(*,*)x(i), xi(i)
      end do
      
      do i2=1, nx+1
       do i1=1, nx+1
       if(i1.eq.i2) then
       dx(i1,i2) = -1.d0/(2.d0*(1.d0+x(i1)))
       else
       dx(i1,i2) = 1.d0/(x(i1)-x(i2))
       endif
       end do
      end do
      dx(nx+1,nx+1)=nx*(nx+2.d0)/4.d0
      write(*,*)'---------'


c  LEGENDRE-GAUSS for eta
      call gauss(ny,y(1))
      y(0)=-1.d0
      y(ny+1)=1.d0

      do j=0, ny+1
      call ytoeta(y(j),pi,res)
      eta(j)=res
c      write(*,*)y(j), eta(j)
      end do

      do j1=1, ny
       do j2=1, ny
       if (j1.eq.j2) then
       dy(j1,j2)=y(j1)/(1.d0-y(j1)*y(j1))
       else
       dy=1.d0/(y(j1)-y(j2))
       end if
       end do
      end do

      do i=0,nx+1
      call xider(x(i),Rm,Rb,a,res)
      derxi(i)=res
c      write(*,*)derxi(i)
      end do
      do j=0, ny+1
      call etader(y(j),pi,res)
      dereta(j)=res
      end do


c  KINETIC ENERGY MATRIX COMPONENTS 
c   T_ii'
      do i1=1, nx
       do i2=1, nx
       s=0.d0
        do k=1, nx
        s=s+(xi(k)**2-1.d0)/(derxi(k)*(1.d0+x(k)))*dx(k,i1)*dx(k,i2)
        end do
       tii(i1,i2)=s*dsqrt((1.d0+x(i1))/derxi(i1))*
     &     dsqrt((1.d0+x(i2))/derxi(i2))  
       end do
      end do
      write(*,*)tii(1,1),tii(1,8),tii(8,1)
      

c   T_jj'
      do j1=1, ny
       do j2=1, ny
       s=0.d0
        do k=1, ny
        s=s+(1.d0-eta(k)**2)/(dereta(k)*(1.d0-y(k)**2))*
     &    dy(k,j1)*dy(k,j2)
        end do
       tjj(j1,j2)=s*dsqrt((1.d0-y(j1)**2)/dereta(j1))*
     &     dsqrt((1.d0-y(j2)**2)/derxi(j2))
       end do
      end do
      write(*,*)tjj(1,1),tjj(1,8),tjj(8,1)

c   T_kk'
      do k1=1, nphi
       do k2=1, nphi
       s=0.d0
        do n=1, (nphi-1)/2
        s=s+n*n*dcos(2.d0*pi*n*(k1-k2)/nphi)
        end do
       tkk(k1,k2)=2.d0*s/nphi
       end do
      end do
      write(*,*)tkk(1,1),tkk(1,8),tkk(8,1)


c     COMPUTING OF H_0 MATRIX
      do k1=1, nphi
      do j1=1, nx
      do i1=1, ny
       i=((i1-1)*ny+j1-1)*nphi+k1
       do k2=1, nphi
       do j2=1, nx
       do i2=1, ny
        j=((i2-1)*ny+j2-1)*nphi+k2
        if ((j1.eq.j2).and.(k1.eq.k2)) then
         h0(i,j) = h0(i,j)+tii(i1,i2)/dsqrt((xi(i1)**2-eta(j1)**2)*
     &      (xi(i2)**2-eta(j2)**2))
        end if
        if ((i1.eq.i2).and.(k1.eq.k2)) then
         h0(i,j) = h0(i,j)+tjj(j1,j2)/dsqrt((xi(i1)**2-eta(j1)**2)*
     &      (xi(i2)**2-eta(j2)**2))
        end if
        if ((i1.eq.i2).and.(j1.eq.j2)) then
         h0(i,j) = h0(i,j)+tkk(k1,k2)/((xi(i1)**2-1.d0)*
     &      (1.d0-eta(j2)**2))
        end if                                         
        if ((i1.eq.i2).and.(j1.eq.j2).and.(k1.eq.k2)) then
         h0(i,j)=h0(i,j)-2.d0*za*xi(i1)/(a*(xi(i1)**2-eta(j1)**2))
        end if
       end do
       end do
       end do
      end do
      end do
      end do
      

c   COMPUTING EIGENVALUES AND VECTORS

      allocate(cprop(nxyphi,nxyphi))
      do i=1,nxyphi
       do j=1,nxyphi
       cprop(i,j)=i+j
       enddo
      end do
      lwmax=99939520
      lda=nxyphi
c      lwork=-1
      lwork=lwmax
      allocate(work(lwmax))
      write(*,*)'start dsyev'
      call dsyev('V', 'U', nxyphi, cprop, lda, val, work, lwork, info)
c      lwork=min(lwmax,int(work(1)))
      write(*,*)info
      write(*,*)lwork,work(1)
      write(*,*)val(1)
      deallocate(work)
      deallocate(cprop)






      end
c----------------------------------------------------------            
      subroutine derLezh(n,x,res)
      implicit real*8(a-h,o-z),integer*4(i-n)
      call Lezh(n-1,x,res)
      z1=res
      call Lezh(n,x,res)
      z2=res      
      res=dble(n)/(1.d0-x**2)*(z1-x*z2)
      return
      end
c---------------------------------------------------------
      subroutine Lezh(n,x,c)
      implicit real*8(a-h,o-z),integer*4(i-n) 
      if (n.eq.0) then 
       c=1.d0
      else 
      if (n.eq.1) then
       c=x
      else 
      a=1.d0
      b=x
      do i=2,n
       c=b*x+dble(i-1)*(b*x-a)/dble(i)
       a=b
       b=c
      end do
      end if
      end if
      return
      end 
c---------------------------------------------------------
      subroutine xtoxi(x,Rm,Rb,a,res)
      implicit real*8(a-h,o-z),integer*4(i-n) 
      frac = Rm/(Rb-a)
      frac = 4.d0*frac
      res = Rm*(1.d0+x)**2/a/(1.d0-x+frac)  + 1.d0 
c       write(*,*)x,  Rm*(1.d0+x)**2/a/(1.d0-x+frac) 
!      res = res +1.d0
      return
      end 
c---------------------------------------------------------
      subroutine ytoeta(y,pi,res)
      implicit real*8(a-h,o-z),integer*4(i-n) 
      arg=(1.d0+y)*pi/2.d0
      res = -dcos(arg)
      return
      end       
c---------------------------------------------------------
      subroutine xider(x,Rm,Rb,a,res)
      implicit real*8(a-h,o-z),integer*4(i-n) 
      frac = Rm/(Rb-a)
      frac = 4.d0*frac
      fracn = frac + 1.d0 - x
      res = (1.d0 + x)*(2.d0*fracn+ 1.d0)/(a*fracn**2)
c      write(*,*)fracn, a, fracn**2
      return
      end 
c---------------------------------------------------------
      subroutine etader(y,pi,res)
      implicit real*8(a-h,o-z),integer*4(i-n) 
      arg = (1.d0+y)*pi/2.d0
      res = pi*dsin(arg)/2.d0
      return
      end 
c--------------------------------------------------------
      subroutine tqli(n,np,z,d,e)
      implicit real(16) (a-h,o-z)
      dimension z(np,np),d(np),e(np)
      if (n.gt.1) then
        do i=2,n
          e(i-1)=e(i)
        end do
        e(n)=0.d0
        do 15 l=1,n
          iter=0
1         do 12 m=l,n-1
            dd=abs(d(m))+abs(d(m+1))
            if (abs(e(m))+dd.eq.dd) go to 2
12        end do
          m=n
2         if(m.ne.l)then
            if(iter.ge.50)write(*,*) 'too many iterations ',iter
            iter=iter+1
            g=(d(l+1)-d(l))/(2.0_16*e(l))
            r=sqrt(g**2+1.0_16)
            if(abs(r+g).lt.abs(r-g)) r=-r
            g=d(m)-d(l)+e(l)/(g+r)
            s=1.0_16
            c=1.0_16
            p=0.0_16
            do 14 i=m-1,l,-1
              f=s*e(i)
              b=c*e(i)
              if(abs(f).ge.abs(g))then
                c=g/f
                r=sqrt(c**2+1.0_16)
                e(i+1)=f*r
                s=1.0_16/r
               c=c*s
              else
                s=f/g
                r=sqrt(s**2+1.0_16)
                e(i+1)=g*r
                c=1.0_16/r
                s=s*c
              endif
              g=d(i+1)-p
              r=(d(i)-g)*s+2.0_16*c*b
              p=s*r
              d(i+1)=g+p
              g=c*r-b
              do 13 k=1,n
                f=z(k,i+1)
                z(k,i+1)=s*z(k,i)+c*f
                z(k,i)=c*z(k,i)-s*f
13            end do
14          end do
            d(l)=d(l)-p
           e(l)=g
            e(m)=0.0_16
            go to 1
          endif
15      end do
      endif
      return
      end            
c-------------------------------------------------------------      
      subroutine gauss(n,dd)
      implicit real(16) (a-h,o-z)
      real(8)  dd
      dimension d(n+1),e(n+1),z(n+1,n+1), dd(n)
      
      np=n+1
      
      do i=1, n
       di=float(i)
       d(i)=0.0_16
       e(i+1)=0.5_16*di/sqrt((di-0.5_16)*(di+0.5_16))
      end do
      
      call tqli(n,np,z,d,e)
      do i=1, n-1
       do j=i+1, n
       if (d(j).lt.d(i)) then
        d2=d(j)
        d(j)=d(i)
        d(i)=d2
       endif
       enddo
      dd(i)=d(i)
      end do
      dd(n)=d(n)
      return
      end
c--------------------------------------------------------
      subroutine radau(n,dd)
      implicit real(16) (a-h,o-z)
      real(8)  dd
      dimension d(n+1),e(n+1),z(n+1,n+1), dd(n)
      
      np=n+1      
      
      do i=1, n
       di=float(i)
       d(i)=-0.25_16/(di-0.5_16)/(di+0.5_16)
       e(i+1)=0.5_16*sqrt(di*(di+1.0_16))/(di+0.5_16)
      end do
      
      call tqli(n,np,z,d,e)
      do i=1, n-1
       do j=i+1, n
       if (d(j).lt.d(i)) then
        d2=d(j)
        d(j)=d(i)
        d(i)=d2
       endif
       enddo
      dd(i)=d(i)
      end do
      dd(n)=d(n)
      return
      end 
      

