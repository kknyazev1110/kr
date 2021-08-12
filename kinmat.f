      program kinmat
      
      include 'kinmat.h'
      
      pi=4.d0*datan(1.d0) 
      n=10
      m=10

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

      lr10=nxy*4      
      open(unit=10,file='kinmat.mat',form='unformatted',
     &    access='direct',recl=lr10,status='replace')

      tcycle=2.d0*pi/o_au
      dt=tcycle/ndt
   
      m=10

      allocate(a(0:nx+1,0:m))
      do i=1,nx
       a(nx-i+1,0)=dcos(pi*(4*i+1)/(4*nx+6))
      end do
      do i=1,nx
       do k=1,m 
         call Lezh(nx,a(i,k-1),res)
         ylezh=res
         call Lezh(nx+1,a(i,k-1),res) 
         wlezh=res
         call derLezh(nx+1,a(i,k-1),ylezh,wlezh,res)
         prlezh=res
         znam=-(nx+1)*(nx+2)*wlezh
         a(i,k)=a(i,k-1)-(1-a(i,k-1)**2)*prlezh/znam
       end do  ! do k=1,m 
      end do  ! do i=1,n
      
      x(0)=-1.d0
      x(nx+1)=1.d0
      do i=1, nx
      x(i)=a(i,m)
      end do
      deallocate(a)

      do i=0,nx+1
        call rad(Rm,Rb,x(i),res)
        r(i)=res
        call f(r(i),Rb,R0,pi,res)
        res1=res
c        write(*,*)res1
       write(*,*) x(i), r(i)
      end do



      allocate(a(0:ny+1,0:m))
      do i=1,ny
       a(ny-i+1,0)=dcos(pi*(4*i)/(4*ny+2))
      end do
      do i=1,ny
       do k=1,m 
         call Lezh(ny-1,a(i,k-1),res)
         ylezh=res
         call Lezh(ny,a(i,k-1),res) 
         wlezh=res
         call derLezh(ny,a(i,k-1),ylezh,wlezh,res)
         prlezh=res
         a(i,k)=a(i,k-1)-wlezh/prlezh
       end do  ! do k=1,m 
      end do  ! do i=1,n
      do j=1,ny
       y(j)=a(j,m) 	
       tetta(j)=dacos(y(j))
c       write (*,*) y(j), tetta(j)
      end do 
      deallocate(a)

      do i=0,nx+1
       call derrad(Rm,Rb,x(i),res)
       rpr(i)=res
c       write(*,*) rpr(i)
      end do      

      write(*,*)'a'
         
c--------TX----------------------------
      write(*,*)'T_x' 
      do i1=1,nx
       do i2=1,nx
        do n=0,nx+1
         call dx(x(n),x(i1),res)
         dx1=res
         call dx(x(n),x(i2),res)
         dx2=res
         tii(i1,i2)=tii(i1,i2)+dx1*dx2/rpr(n)
        end do
        tii(i1,i2)=tii(i1,i2)/(dsqrt(rpr(i1)*rpr(i2)))/2.d0        
       end do        
      end do

c--------TY---------------------------
      write(*,*)'T_y' 
      do j1=1,ny
       do j2=1,ny
        do n=1,ny
         call dy(y(n),y(j1),res)
         dy1=res
         call dy(y(n),y(j2),res)
         dy2=res
         tjj0(j1,j2)=tjj0(j1,j2)+dy1*dy2
         tjj1(j1,j2)=tjj1(j1,j2)+dy1*dy2*(1.d0-y(n)*y(n))
        end do
        tjj0(j1,j2)=tjj0(j1,j2)*dsqrt(1.d0-y(j1)*y(j1))*
     &              dsqrt(1.d0-y(j2)*y(j2))/2.d0
        tjj1(j1,j2)=tjj1(j1,j2)/2.d0
       end do
      end do

c-------H_0 calculation and diagonalization------------------
      open(8,file='energy.txt',status='new')
      open(6,file='cpsi0.txt',status='new')
      do 30 m=0, (nphi-1)/2 !m
      write(*,*)m            
      h0=0.d0
      if (mod(m,2).eq.0) then   !even angular momentum
      j=0
      do j2=1,ny
      do i2=1,nx
      j=j+1
      i=0 
       do j1=1,ny
       do i1=1,nx
        i=i+1
c         i=(i1-1)*ny+j1
c         j=(i2-1)*ny+j2
         if ((j1 .eq. j2)) then
         h0(i,j)=h0(i,j)+tii(i1,i2)
         end if
         if ((i1.eq.i2)) then
         h0(i,j)=h0(i,j)+tjj0(j1,j2)/(r(i1)**2)
         end if 
         if ((i1.eq.i2).and.(j1.eq.j2)) then
         h0(i,j)=h0(i,j)- 1.d0/(r(i1))+
     &           0.5d0*m*m/(r(i1)**2)/(1-y(j1)**2)
         end if
        enddo
        enddo
        enddo
        enddo
      else                     !odd angular momentum 
      j=0
      do j2=1,ny
      do i2=1,nx
      j=j+1
      i=0 
       do j1=1,ny
       do i1=1,nx
        i=i+1
c         i=(i1-1)*ny+j1
c         j=(i2-1)*ny+j2
         if ((j1 .eq. j2)) then
         h0(i,j)=h0(i,j)+tii(i1,i2)
         end if
         if ((i1.eq.i2)) then
         h0(i,j)=h0(i,j)+tjj1(j1,j2)/(r(i1)**2)
         end if 
         if ((i1.eq.i2).and.(j1.eq.j2)) then
         h0(i,j)=h0(i,j)- 1.d0/(r(i1))+ 1.d0/(r(i1)**2)+
     &           0.5d0*(m*m-1.d0)/(r(i1)**2)/(1-y(j1)**2)
         end if
        enddo
        enddo
        enddo
        enddo
      endif

      lwmax=100000
      lda=nxy
c      lwork=-1
      lwork=lwmax
      allocate(work(lwmax))
      write(*,*)'start dsyev', m
      call dsyev('V', 'U', nxy, h0, lda, val, work, lwork, info)
c      lwork=min(lwmax,int(work(1)))
      deallocate(work)

      if (info.ne.0) then
       write(*,*)info
      endif
      write(8,*)'E= ',val(1), -1.d0/2.d0/((m+1)**2)
      if(m.eq.0) then
       e0=val(1)
c       write(*,*)e0
       do i=1,nxy
        cpsi0(i)=dcmplx(h0(i,1),0.d0)
        write(6,*)dreal(cpsi0(i))
       end do
       write(10,rec=2) (cpsi0(k),k=1,nxy)	! GROUND STATE
      end if
      
      allocate(cvec(nxy,nxy))
      do i=1,nxy
       do j=1,nxy
        cvec(i,j)=dcmplx(h0(i,j),0.d0)
       end do
      end do
c      write(*,*)cvec(1,8),cvec(8,1)
      allocate(ce(nxy,nxy))
     
      do i=1,nxy
       ce(i,i)=dcmplx(dcos(dt*val(i)/2.d0),-dsin(dt*val(i)/2.d0))
      end do
c      write(*,*) ce(1,1)

c------Propagator matrix---------------------------------


      allocate(cpsim(nxy,nxy))
      allocate(cprop(nxy,nxy))
      allocate(cp(nxy,nxy))
     
      calpha=dcmplx(1.d0,0.d0)
      cbeta=dcmplx(0.d0,0.d0)
      call zgemm('n','n',nxy,nxy,nxy,calpha,cvec,nxy,
     &ce,nxy,cbeta,cpsim,nxy)
c      write(*,*) cpsim(2,2), cpsim(10,20)
      write(*,*)'1st zgemm'
      call zgemm('n','c',nxy,nxy,nxy,calpha,cpsim,nxy,
     &cvec,nxy,cbeta,cprop,nxy)
c      write(*,*) cprop(2,2), cprop(10,20),cprop(20,10)
      write(*,*)'2nd zgemm'
      call zgemm('n','c',nxy,nxy,nxy,calpha,cprop,nxy,
     & cprop,nxy,cbeta,cp,nxy)
      write(*,*) cp(2,2),cp(10,10), cp(10,20)

c      do i=1,nxy
c       do j=1, nxy
c        cpropall(m,i,j)=cprop(i,j)
c       end do
c      end do

      deallocate(ce)
      deallocate(cp)
      deallocate(cpsim)
c       do k1=1,nxy
        do k2=1,nxy
         write(10,rec=2+m*nxy+k2) (cprop(k1,k2),k1=1,nxy)
        end do
c       end do
      deallocate(cprop)
      deallocate(cvec)

30    end do    !m
      write(10,rec=1) lr10,Rm,Rb,R0,ndt,o_au,e0
      close(unit=10)

      end


c---------------------------------------------------------      
      subroutine rad(Rm,Rb,x,c)
      implicit real*8(a-h,o-z),integer*4(i-n)
      c=Rm*(1.d0+x)/(1.d0-x+(2.d0*Rm)/(Rb))
      return
      end
c---------------------------------------------------------
      subroutine derrad(Rm,Rb,x,c)
      implicit real*8(a-h,o-z),integer*4(i-n)
      c=2.d0*Rm*(1.d0+Rm/Rb)/((1.d0-x+2.d0*Rm/Rb)**2)
      return
      end
c---------------------------------------------------------
      subroutine derLezh(n,x,z1,z2,c)
      implicit real*8(a-h,o-z),integer*4(i-n)
      c=dble(n)/(1.d0-x**2)*(z1-x*z2)
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
      subroutine dx(x1,x2,c)
      implicit real*8(a-h,o-z),integer*4(i-n)
      if (x1 .eq. x2) then 
       c=0.d0
      else 
      c=1.d0/(x1-x2)
      end if
      return
      end
c---------------------------------------------------------
      subroutine dy(y1,y2,c)
      implicit real*8(a-h,o-z),integer*4(i-n)
      if (y1 .eq. y2) then 
      c=y1/(1.d0-y1**2)
      else 
      c=1.d0/(y1-y2)
      end if
      return
      end
c---------------------------------------------------------
      subroutine f(x,Rb,R0,pi,c)
      implicit real*8(a-h,o-z),integer*4(i-n)
      if (x<R0) then
      c=1
      else
      c=(dcos(((x-R0)/(Rb-R0))*pi/2.d0))**(0.25d0)
      end if
      return
      end



