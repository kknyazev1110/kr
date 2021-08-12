      program disret
      implicit real*8(a-b,d-h,o-z),integer*4(i-n),complex*16(c)
      parameter(nx=72)
      parameter(ny=32)
      dimension x(0:nx+1), y(1:ny), xi(0:nx+1), eta(1:ny), y1(0:ny+1),
     &          gry(0:ny)
      allocatable :: z1(:,:)

      pi=4.d0*datan(1.d0)
      m=30
       
c     LEGENDRE-GAUSS for eta
      allocate(z1(1:ny,0:m))
      do i=1,ny
       z1(ny-i+1,0)=dcos(pi*(4*i+1)/(4*ny+6))
      end do
      do i=1,ny
       do k=1,m 
c         call Lezh(ny-1,z(i,k-1),res)
c         ylezh=res
         call Lezh(ny,z1(i,k-1),res) 
         wlezh=res
         call derLezh(ny,z1(i,k-1),res)
         prlezh=res
         z1(i,k)=z1(i,k-1)-wlezh/prlezh
       end do  ! do k=1,m 
      end do  ! do i=1,n

      do i=1,ny
       write(*,*) z1(i,m)
      end do      
      write(*,*)'-----------'
       do j=1,ny
       y1(j)=z1(j,m) 
       call ytoeta(y(j),pi,res)	
       eta(j)=res!y(j)
       write (*,*) y1(j), eta(j)
      end do 
      deallocate(z1)     
      write(*,*)'-----------'
      
      call gauss(ny,y1(1))      
      y1(0)=-1.d0
      y1(ny+1)=1.d0
      
      do j=0, ny+1
       write(*,*)y1(j)
      end do
      
      
      
      
      
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
      subroutine ytoeta(y,pi,res)
      implicit real*8(a-h,o-z),integer*4(i-n) 
      arg=(1.d0+y)*pi/2.d0
      res = -dcos(arg)
      return
      end       
c--------------------------------------------------------
      subroutine gauss(n,dd)
      implicit real(16) (a-h,o-z)
      real(8)  dd
      dimension d(n+1),e(n+1),a(n+1,n+1), dd(n)
      
      np=n+1
      
      do i=1, n
       di=float(i)
       d(i)=0.0_16
       e(i+1)=0.5_16*di/sqrt((di-0.5_16)*(di+0.5_16))
      end do
      
      call tqli(n,np,a,d,e)
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
      subroutine tqli(n,np,a,d,e)
      implicit real(16) (a-h,o-z)
      dimension a(np,np),d(np),e(np)
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
                f=a(k,i+1)
                a(k,i+1)=s*a(k,i)+c*f
                a(k,i)=c*a(k,i)-s*f
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
