      program field
      implicit real*8(a-b,d-h,o-z),integer*4(i-n),complex*16(c)
      parameter(ndt=4096)
      parameter(ndw=20)

      pi=4.d0*datan(1.d0)
      q_m=9.109383701528282d-31
      e_h=4.3597447222071d-18
      h=6.62607015d-34
      hbar=h/(2.d0*pi)
      s_l=2.99792458d8
      q=1.602176634d-19
      qc=1.112650055451717d-10 
      alpha=q**2.d0/(qc*hbar*s_l) 
      a0=5.291772109038080d-11
c------------------------------------------------------
      open(4, file = 'fieldparams.dat', status = 'new')
c------------------------------------------------------

      fw_lenght=8d-7                      !wave length m
      f_int=1.d18                        !field intencity W/m^2
      omega=2.d0*pi*s_l/fw_lenght        !frequency s^-1
      o_au=omega*a0/(alpha*s_l)           !frequency a.u.
      f_e=dsqrt(8.d0*pi*qc*f_int/s_l)
      f_eau=f_e*a0*a0/q
      t_p=2.d0*pi/omega
      ttt=2.d0*pi/o_au
      dt=ttt/dble(ndt)                    !time step
      dw=o_au/dble(ndw)                      !frequency step
      per=24.d0*pi/o_au
c      write(*,*)t_p, t_pau, dt
      write(*,*) o_au-45.56335235d0/800.d0, 
     & f_eau-0.53380271076d-08*1d7
      write(4,*)ndt
      write(4,*)ndw
      write(4,*)o_au
      write(4,*)f_eau
      write(4,*)s_l
      write(4,*)alpha
      write(4,*)30.d0
      write(4,*)60.d0
      write(4,*)40.d0








      end      
