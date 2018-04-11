! this is 1DMC slit pore adsorption DFT program;
! can be used to calculate the gas adsorption in slit pore;
! you need input the gas type, gas number, adsorbent materials materials;
! we can output the adsorption isotherm, amount of adsorption, density
!   profile at certain pressure, grand free energy;
!
! developed by Hongguan Wu, Feb, 2018.
program main
  use global_module_energy
  implicit none
  real*8 err, err1, err_adsop, err_density, sum0, sum1
  real*8 adsp,adspp, VDP, SDP, VDS
  real*8 dllw, cpid, aa, bb, cc, dd, F
  real*8 time(0:1)
  integer i, number, m, j, k

  call cpu_time (time(0))
  write(*,*) "input the objects' parameters"
  call system_def
  !write(*,*) '1 number_of_comp is:' , number_of_comp

  write(*,*) 'calculate external potential'
  call comp_external_pot
  !write(*,*) '2 number_of_comp is:' , number_of_comp

  write(*,*) 'calculate bulk chemical potential, rho0, Miu'
  call comp_bulk_pot

  write(*,*) 'calculate miu_att'
  call miu_att
  ! initial the density of fluids in all position
  !call initial_density(P)

  open(31, file="out_iteration_process.dat")
  open(39, file="out_calculation_process.dat")
  open(32, file="out_density_profile.dat")
  open(33, file="out_thermodynamic_properties.dat")

  write(*,*) "Start the main process!"
  do i = 1, num_press ! presusre
    ! need to define the initial density distribution of components
    write(*,*) 'Now the pressure is pressure is', pres(i)
    F = pica_factor(i)
    if ( i > 1 ) then
      open(88, file = "density_exchange_file.dat")
      do j = 1, number_of_comp
        do k = 0, NL
          read(88,*) aa, bb, cc, P(j,k), dd ! i, j, k, P, PP
          P(j,k) = P(j,k) + rho0(i,j) - rho0(i-1,j)
          if(P(j,k) < 0.0)  P(j,k) = 0.0
        end do
      end do
      close(88)
    else
      do j =1, number_of_comp
        do k =0, NL
          P(j,k) = rho0(i,j)/2.0
        end do
      end do
    end if

    !write(*,*) 'start the',i,'th pressure calculate'
    do j = 1, number_of_comp ! components
      number = 0
      err = 1.0
      !write(*,*) 'start the', j, 'component'
      do while ( (number<=8) .or. ((number<Nmax).and.(err>errmax)) )
        number=number+1
        !write(*,*) 'start the hs calculation'
        ! calculate HS chemical and HS free energy
        call com_hs_pot
        !write(*,*) 'start the att calculation'
        ! calculate Att chemical potetial and ATT free energy
        call com_att_pot


        err1=0.0
        sum0 = xnf*P(j,0)
        PP(j,0) = rho0(i,j) * exp( Miu(i,j) - cp_hs(j,0) - cp_att(j,0) - cp_ext(j,0) )
        sum1 = xnf*PP(j,0)

        do k = 1, NL
          PP(j,k) = rho0(i,j) * exp( Miu(i,j) - cp_hs(j,k) - cp_att(j,k) - cp_ext(j,k) )
          if ( PP(j,k)>1.0e-10 ) then
            err_density = abs(PP(j,k)-P(j,k))/max(PP(j,k), P(j,k))
          else
            err_density = abs(PP(j,k)-P(j,k))
          end if
          err1=max(err1, err_density)
          sum0 = sum0+P(j,k)*deltx
          sum1 = sum1+PP(j,k)*deltx
          P(j,k) = P(j,k)*(1.0-F) + PP(j,k)*F
          !write(*,*) 'now we are in position:', k
        end do

        err_adsop = abs(sum1-sum0)/max(abs(sum1), abs(sum0))
        err=max(err1, err_adsop)
        write(*,*) "I'm in the picard iteration", 'err1:',err1, 'err_adsop:',err_adsop, 'err',err
        ! output the iteration process
        if (ITERA(i).eq.1) then ! ITERA not define
          write(*,111) i, j, number, F, err1, err_adsop, err
          write(31, 111) i, j,number, F, err1, err_adsop, err
        end if
        111 format(3(1X, I5), 1X, F6.4, 3(1X, E13.6))
      end do ! end while loop
      write(*,*) 'end the while loop'

      ! 39 out put the calculate process
      write(*,999) i,j, pres(i), rho0(i,j), number, err
      write(39, 999) i,j, pres(i), rho0(i,j), number, err
      999 format(2(1X, I5), 2(1X, E13.6), 1X, I5, 1X, E13.6)
      write(*,*) 'end the output of out_calculation_process.dat'

      ! 38 is a exchange file, gives initial density value for next rel_prees
      open (38, file="density_exchange_file.dat")
      do m = 0, NL
        write(38, 888) i, j, m, P(j,m), PP(j,m)
      end do
      888 format(3(1X, I3), 2(1X, E13.6))
      close(38)

      ! 32 out put density profile
      if ( IFDEN(i)==1 ) then
        do m = 0, NL
          write(32, 222) i, j, pres(i), rho0(i,j), xnf+m*deltx, P(j,m), PP(j,m)
        end do
      end if
      222 format(2(1X, I5), 5(1X, E13.6))

      write(*,*) 'start calculate amount of adsorption'
      ! calculate adsorption amount
      adsp=P(j,0)*xnf
      do m=0, NL
        if ( m==0 .or. m==NL ) then
          adsp=adsp + 0.5*P(j,m)*deltx
        else
          adsp=adsp + P(j,m)*deltx
        end if
      end do
      write(*,*) "total adsp amount is:", adsp
      SDP=2.0
      VDP=length_of_stru!-2.0*RS
      VDS=2*width_of_wall
      adspp=rho0(i,j)*VDP ! bulk
      ADSP0(j,i)=adspp/(VNMOL*SDP*d_min**2) !per pore area
      ADSP1(j,i)=2*adsp/(VNMOL*SDP*d_min**2) ! per pore area
      ADSP2(j,i)=2*adsp/(VNMOL*VDP*d_min**3) ! per pore volume
      ADSP3(j,i)=2*adsp/(VNMOL*VDS*d_min**3) ! per adsorbent volume


      ! calculate chemical potential, weight_of_comp unit should be ISO
      dllw = sqrt(plank**2/(2.0*pi*weight_of_comp(j)*kB*temperature))
      cpid = kB*temperature*log(rho0(i,j)*(dllw/d_min)**3)
      CP(i,j)= Miu(i,j) + cpid

      ! calculate the grand potential
      GRP(i,j) = 0.0!comp_grand_free_energy(P)

      ! 33 output the isothermodynamic properties and adsorption amount
      write(33, 333) i, j, pres(i), rho0(i,j), CP(i,j), GRP(i,j), ADSP0(j,i), ADSP1(j,i),&
      ADSP2(j,i), ADSP3(j,i)
      333  format(2(1X, I5), 3(1X, E13.6), 1X, E19.12, 4(1X, E15.6))

    end do ! end j components
  end do ! end i pressure
  close(31)
  close(32)
  close(33)
  close(39)
  call cpu_time (time(1))
  !time_diff = time(1)- time(0)
  write(*,*) 'the total time cost:', time(1)-time(0), 'seconds'
end program
