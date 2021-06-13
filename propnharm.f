      program propnharm
      use omp_lib
      
      include 'propnharm.h'
      nxy=nyx

      
      open(97, file = 'moment.txt', status = 'unknown')
      open(98, file = 'acceleration.txt', status = 'unknown')
      open(87, file = 'hhg.txt', status = 'unknown')
      open(88, file = 'vn.txt', status = 'unknown')
c      open(77, file = 'coord0.txt', status = 'new')

      
      lr10=nxyphi*4     
      open(10,file='kinmat.mat',form='unformatted',
     >     access='direct',recl=lr10,status='old')

      read(10,rec=1) lr10,Rm,Rb,R0,ndt,o_au,e0
      close(10) 

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
      close(4)
      write(*,*) e0,o_au

      pi=4.d0*datan(1.d0)

c-------------------------------------

      tcycle=2.0d0*pi/o_au
      dt=tcycle/dble(ndt)                     !time step
      dw=o_au/20.d0                       !frequency step
      per=24.d0*pi/o_au

c------------------------------------
      
      w(1)=o_au
      npuls(1)=3
      tau(1)=2.d0*pi*npuls(1)/o_au
      etta(1)=1.d0
      t_pau(1)=2.0d0*pi/w(1)      
      w(2)=2.d0*o_au
      npuls(2)=3
      tau(2)=2.d0*pi*npuls(2)/o_au
      etta(2)=-1.d0
      t_pau(2)=2.0d0*pi/w(2)

c-----Time delay--------------------

      td=0.d0
      
c-----Grid------------------------------
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
c       write(*,*) x(i), r(i)
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

      do k=1, nphi
       phi(k)=2.d0*pi*(k-1)/dble(nphi)
      end do
      write(*,*)'a'
c-------------------------------

      open(10,file='kinmat.mat',form='unformatted',
     &     access='direct',recl=lr10,status='old')
      read(10,rec=2) (cpsi0(k),k=1,nxy)
      vn=dot_product(cpsi0,cpsi0)
      write(*,*)vn

c------Initial wave function---------------------------------

      dnphi=dble(nphi)
      do k=1, nphi
       do j=1, ny
        do i=1, nx
         ij=(j-1)*nx+i
         ijk=(k-1)*nxy+(j-1)*nx+i   !index
         cpsi(ijk)=cpsi0(ij)/dsqrt(dnphi)
        end do
       end do
      end do
      vn=dot_product(cpsi,cpsi)
      write(*,*)vn
      

c------Reading propagation matrix from file-------------------

      do m=0,(nphi-1)/2
       do k2=1,nxy
        read(10,rec=2+m*nxy+k2) (cpropall(m,k1,k2),k1=1,nxy)
       end do
      end do
      close(10)

c------Dipole moment and acceleration at start----------------

      time=-6*t_pau(1)
      call f_x(f_eau,w(1),time,tau(1),res)
      f1_x=res
      call f_y(f_eau,w(1),time,tau(1),etta(1),res)
      f1_y=res
      call f_x(f_eau,w(2),time-td,tau(2),res)
      f2_x=res
      call f_y(f_eau,w(2),time-td,tau(2),etta(2),res)
      f2_y=res
      write(*,*)'first step start'

      dip1=0.d0
      dip2=0.d0
      dip3=0.d0
      bz1=0.d0
      bz2=0.d0
      bz3=0.d0
c!$omp parallel do reduction(+:dip1,dip2,dip3,bz1,bz2,bz3)  
      do k=1,nphi
       do j=1,ny
        do i=1,nx
         m=(k-1)*nxy+(j-1)*nx+i   !index 
         call f(r(i),Rb,R0,pi,res)
         fun=res
         cpsi(m)=fun*cpsi(m)
         u0=cpsi(m)*dconjg(cpsi(m))
         xk=r(i)*dsin(tetta(j))*dcos(phi(k))
         yk=r(i)*dsin(tetta(j))*dsin(phi(k))
         zk=r(i)*dcos(tetta(j))
         dip1=dip1+xk*u0
         dip2=dip2+yk*u0
         dip3=dip3+zk*u0
         bz1=bz1-(xk/(r(i)**3))*u0!+f1_x+f2_x)*u0
         bz2=bz2-(yk/(r(i)**3))*u0!+f1_y+f2_y)*u0
         bz3=bz3-(zk/(r(i)**3))*u0
         write(77,*) xk/r(i),yk/r(i),zk/r(i)
         xk=0.d0
         yk=0.d0
         zk=0.d0
         u0=0.d0
         end do
        end do
       end do
c       close(77)
c!$end omp parallel do
      vn=dot_product(cpsi,cpsi)
      dip(1,0)=dip1
      dip(2,0)=dip2
      dip(3,0)=dip3
      bz(1,0)=bz1-(f1_x+f2_x)!*vn
      bz(2,0)=bz2-(f1_y+f2_y)!*vn
      bz(3,0)=bz3
      write(98,*) dip(1,0),dip(2,0),dip(3,0)
      write(97,*) bz(1,0),bz(2,0),bz(3,0)
      write(*,*)vn
      write(88,*)time/t_pau(1), vn
      vn=0.d0      

      write(*,*) 'end of first step'
     
c------Propagation--------------------------------------------

      do 10 n=1,ndt*12

       call halfv
      
        cpot=dcmplx(0.d0,0.d0)

       call f_x(f_eau,w(1),time+dt/2.d0,tau(1),res)
        f1_x=res
       call f_y(f_eau,w(1),time+dt/2.d0,tau(1),etta(1),res)
        f1_y=res
       call f_x(f_eau,w(2),time-td+dt/2.d0,tau(2),res)
        f2_x=res
       call f_y(f_eau,w(2),time-td+dt/2.d0,tau(2),etta(2),res)
        f2_y=res
        do k=1, nphi
         do j=1, ny
          do i=1, nx
          cpot=dcmplx(0.d0,0.d0)
          m=(k-1)*nxy+(j-1)*nx+i   !index
          xk=r(i)*dsin(tetta(j))*dcos(phi(k))
          yk=r(i)*dsin(tetta(j))*dsin(phi(k))
          pot=(f1_x+f2_x)*xk+(f1_y+f2_y)*yk
          cpot=dcmplx(dcos(dt*pot),-dsin(dt*pot)) 
          cpsi(m)=cpsi(m)*cpot 
          xk=0.d0
          yk=0.d0    
          end do
         end do
        end do

       call halfv

       time=time+dt

       call f_x(f_eau,w(1),time,tau(1),res)
       f1_x=res
       call f_y(f_eau,w(1),time,tau(1),etta(1),res)
       f1_y=res
       call f_x(f_eau,w(2),time-td,tau(2),res)
       f2_x=res
       call f_y(f_eau,w(2),time-td,tau(2),etta(2),res)
       f2_y=res

      dip1=0.d0
      dip2=0.d0
      dip3=0.d0
      bz1=0.d0
      bz2=0.d0
      bz3=0.d0

c!$omp parallel do reduction(+:dip1,dip2,dip3,bz1,bz2,bz3)  
      do k=1,nphi
       do j=1,ny
        do i=1,nx
         m=(k-1)*nxy+(j-1)*nx+i   !index 
         call f(r(i),Rb,R0,pi,res)
         fun=res
         cpsi(m)=fun*cpsi(m)
         u0=cpsi(m)*dconjg(cpsi(m))
         xk=r(i)*dsin(tetta(j))*dcos(phi(k))
         yk=r(i)*dsin(tetta(j))*dsin(phi(k))
         zk=r(i)*dcos(tetta(j))
         dip1=dip1+xk*u0
         dip2=dip2+yk*u0
         dip3=dip3+zk*u0
         bz1=bz1-(xk/((r(i))**3))*u0!+f1_x+f2_x)*u0
         bz2=bz2-(yk/((r(i))**3))*u0!+f1_y+f2_y)*u0
         bz3=bz3-(zk/((r(i))**3))*u0
         xk=0.d0
         yk=0.d0
         zk=0.d0
         u0=0.d0
        end do
       end do
      end do
c!$end omp parallel do
       vn=dot_product(cpsi,cpsi)
      dip(1,n)=dip1
      dip(2,n)=dip2
      dip(3,n)=dip3
      bz(1,n)=bz1-(f1_x+f2_x)!*vn
      bz(2,n)=bz2-(f1_y+f2_y)!*vn
      bz(3,n)=bz3	
      write(88,*)time/t_pau(1), vn
      write(98,*) dip(1,n),dip(2,n),dip(3,n)
      write(97,*) bz(1,n),bz(2,n),bz(3,n)
      vn=0.d0
      n0=mod(n,4096)
      if (n0.eq.0) then
       write(*,*) n/4096
      endif

10    end do

c------Harmonics calculation----------------------------------

      w0=0.d0
      do k1=0, 4000
       time=0.d0
       s1=0.d0       
       s2=0.d0
       cd=(0.d0,0.d0)
       cbz=(0.d0,0.d0)
       an=0.d0
       ann=0.d0
       do k2=0, 12*4096
        cres=dcmplx(dcos(w0*time),dsin(w0*time))
        call tt(o_au,time,per,res)
        fun1=res
         do i=1,3
         cd(i)=cd(i)+cres*fun1*dip(i,k2)
         cbz(i)=cbz(i)+cres*fun1*bz(i,k2)
         end do
        time=time+dt
       end do
       van=dot_product(cd,cd)
       vann=dot_product(cbz,cbz)
       s1=2.d0*(w0**4)*van/(3.d0*pi*((s_l*alpha)**3))     
       s2=2.d0*vann/(3.d0*pi*((s_l*alpha)**3))
c       s4=2.d0*annn/(3.d0*pi*((s_l*alpha)**3))
       write(87,*) w0/o_au ,s1, s2
c       write(97,*) w/o_au ,s2, s4
       w0=w0+dw
      end do






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
      subroutine f(x,Rb,R0,pi,c)
      implicit real*8(a-h,o-z),integer*4(i-n)
      if (x<R0) then
      c=1
      else
      c=(dcos(((x-R0)/(Rb-R0))*pi/2.d0))**(1/4)
      end if
      return
      end
c---------------------------------------------------------
      subroutine f_x(f_eau,w,t,tau,res)
      implicit real*8(a-h,o-z),integer*4(i-n)
      res1=2.d0*t/(tau**2)*dcos(w*t)*(-2.d0*dlog(2.d0))-w*dsin(w*t)
      res2=-f_eau/(w)*dexp(-2.d0*dlog(2.d0)*
     & (t**2)/(tau**2))
      res=res1*res2
      return
      end
c---------------------------------------------------------
      subroutine f_y(f_eau,w,t,tau,etta,res)
      implicit real*8(a-h,o-z),integer*4(i-n)
      res1=2.d0*t/(tau**2)*dsin(w*t)*(-2.d0*dlog(2.d0))+w*dcos(w*t)
      res2=-f_eau/(w)*dexp(-2.d0*dlog(2.d0)*
     & (t**2)/(tau**2))*etta
      res=res1*res2
      return
      end
c---------------------------------------------------------
      subroutine tt(o_au,time,per,res)
      implicit real*8(a-b,d-h,o-z),integer*4(i-n),complex*16(c)
      if ((time.ge.0).and.(time.lt.per/10.d0)) then
      res=(dsin(o_au*time/8.d0))**2
      else if ((time.ge.per/10.d0).and.(time.lt.per*9.d0/10.d0)) then
      res=1.d0
      else if ((time.ge.per*9.d0/10.d0).and.(time.le.per)) then
      res=(dsin(o_au*(per-time)/8.d0))**2
      end if
      return
      end
c--------------------------------------------------------
      subroutine dipole(ip,np,npp)
      include 'propnharm.h'
      nxy=nyx

      iip=nx/np
      iis=ip*iip+1
      iif=iis+iip-1
      if(ip.eq.np-1)then
        iif=nx
        npp=np
      endif
      ddx(ip)=0.d0
      ddy(ip)=0.d0
      ddz(ip)=0.d0
      dax(ip)=0.d0
      day(ip)=0.d0
      daz(ip)=0.d0
      do  k=1,nphi
        do  j=1,ny
          do  i=iis,iif
            ij=(j-1)*nx+i
            m=(k-1)*nxy+(j-1)*nx+i
            u0=cpsi(m)*dconjg(cpsi(m))
            xk=r(i)*dsin(tetta(j))*dcos(phi(k))
            yk=r(i)*dsin(tetta(j))*dsin(phi(k))
            zk=r(i)*dcos(tetta(j))
            ddx(ip)=ddx(ip)+xk*u0
            ddy(ip)=ddy(ip)+yk*u0
            ddz(ip)=ddz(ip)+zk*u0
            dax(ip)=dax(ip)-(xk/((r(i))**3))*u0
            day(ip)=day(ip)-(yk/((r(i))**3))*u0
            daz(ip)=daz(ip)-(zk/((r(i))**3))*u0
         write(77,*)xk,yk,zk
          write(*,*)'done'
          end do
        end do
      end do
      return
      end
c---------------------------------------------------------
      subroutine psi_to_m(ip,np,nxy,m,cnpsim,cua)
      include 'propnharm.h'
      nxy=nyx
      kp=nxy/np
      ks=ip*kp+1
      kf=ks+kp-1
      if(ip.eq.np-1) then
       kf=nxy
      endif
      do k=ks,kf
       cua(k,1)=cnpsim(k,-m)
       cua(k,2)=cnpsim(k,m)
      end do
      return
      end
c---------------------------------------------------------
      subroutine m_to_psi(ip,np,nxy,m,cnpsim,cub)
      include 'propnharm.h'
      nxy=nyx
      kp=nxy/np
      ks=ip*kp+1
      kf=ks+kp-1
      if(ip.eq.np-1) then
       kf=nxy
      endif
      do k=ks,kf
       cnpsim(k,-m)=cub(k,1)
       cnpsim(k,m)=cub(k,2)
      end do
      return
      end
c-------------------------------------------------------
      subroutine halfv
      use omp_lib
      use mkl_dfti
      include 'propnharm.h'
      type(dfti_descriptor),pointer::my_fft
      nxy=nyx
      dnphi=dble(nphi)      
      ist=DftiCreateDescriptor(my_fft,dfti_double,dfti_complex,1,nphi)
      ist=DftiCommitDescriptor(my_fft)
       do j=1,ny
        do i=1,nx
         ij=(j-1)*nx+i
         do k=1,nphi
          ijk=(k-1)*nxy+(j-1)*nx+i   !index 
          cnpsi(k)=cpsi(ijk)
         end do       
         ist=DftiComputeForward(my_fft,cnpsi)      
         do m=-(nphi-1)/2,(nphi-1)/2
          if(m.lt.0) then
           ii=nphi+m+1
          else
           ii=m+1
          endif
          cnpsim(ij,m)=cnpsi(ii)/dnphi
         end do
         cnpsi=0.d0
        end do
       end do
       do m=0,(nphi-1)/2
        ml=m
!$omp parallel default(private) shared(ml)
       np=omp_get_num_threads()
       ip=omp_get_thread_num()
       call psi_to_m(ip,np,nxy,m,cnpsim,cua)
!$omp end parallel
       nu=2
       allocate(cprop(nxy,nxy))
       do i=1,nxy
        do j=1,nxy
         cprop(i,j)=cpropall(m,i,j)
        end do
       end do
       call zgemm('n','n',nxy,nu,nxy,calpha,cprop,
     &      nxy,cua,nxy,cbeta,cub,nxy)
       deallocate(cprop)
!$omp parallel default(private) shared(ml)
       np=omp_get_num_threads()
       ip=omp_get_thread_num()
       call m_to_psi(ip,np,nxy,m,cnpsim,cub)
!$omp end parallel
      end do
      do j=1,ny
       do i=1,nx
        ij=(j-1)*nx+i
        do m=-(nphi-1)/2,(nphi-1)/2
         if(m.lt.0) then
          ii=nphi+m+1
         else
          ii=m+1
         endif
         cnpsi(ii)=cnpsim(ij,m)
        end do
        ist=DftiComputeBackward(my_fft,cnpsi)
        do k=1,nphi
         ijk=(k-1)*nxy+(j-1)*nx+i   !index 
         cpsi(ijk)=cnpsi(k)
        end do
       end do
      end do
      ist=DftiFreeDescriptor(my_fft)
      return
      end

      
      
      
      
      
      
      
      
      


