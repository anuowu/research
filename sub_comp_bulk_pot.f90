! this subroutine calculate the excess chemical potential of LJ fluids
! based on BMWR EOS
! aim to get miu(i) and rho0(i)
subroutine comp_bulk_pot
  use global_module_energy
  implicit none
  real*8 re_pres, re_Ar, re_Ur
  real*8 si, ep, errps, Ts, Ps, rhos0, Ps0, rhos, Ar, Ur
  real*8 sigxrho(max_num_press),epsixrho(max_num_press), Arxrho(max_num_press)
  integer nu, i, j, k
  external re_pres, re_Ar, re_Ur

  !write(*,*)  "get the sigma_x and eplison_x for mixture in MBWR"
  ! get the sigma_x and eplison_x for mixture in MBWR
  si = 0.0
  ep = 0.0
  !write(*,*) '3 number_of_comp', number_of_comp
  do i = 1, number_of_comp
    do j= 1, number_of_comp
      si = si +X(i)*X(j)*((sigma_gas(i)+sigma_gas(j))/2.0)**3
      ep = ep +X(i)*X(j)*((sigma_gas(i)+sigma_gas(j))/2.0)**3 &
      *sqrt(epsilon_gas(i)*epsilon_gas(j))
    end do
  end do
  write(*,*) 'here si is:', si
  ep=ep/si

  !write(*,*) 'open a bulk file'
  open(22, file = "out_bulk_property.dat")
  write(22,"2(A15,5X),3(A20,5X),2(A10,5X)" ) "number","pressure", "density", "real density", "bulk chemical potential", "iteration times:nu", "error"
  !! calculate the different components bulk density
  Ts = temperature/ep ! because here ep unit is K
  do i = 1, num_press

    Ps = pres(i)*(si*1.0E-30)/(ep*kB) !because Pres is ISO，so si and ep is transfered to ISO
    rhos0 = pres(i)/temperature/kB*1.0E-30*si*d_min**3 !here unit is rho*sigma^3, it is reduced
    Ps0 = re_pres(rhos0, Ts)
    errps = (Ps - Ps0)
    write(*,*) 'i is:', i, 'Ps is:', Ps, "Ps0 is:", Ps0
    write(*,*) 'the reduced bulk density is', rhos0, 'errps is:', errps
    nu=0
    do while (abs(errps) > 0.001 .and. nu < 2000)
      nu=nu+1
      rhos = rhos0 * (1+errps)
      Ps0 = re_pres(rhos, Ts)
      errps = (Ps - Ps0)
      rhos0 = rhos
      write(*,*) 'Ps0 is:', Ps0
      write(*,*) 'errps is:', errps, 'rhos0 is:', rhos0
    end do

    ! get pre components' bulk density
    do j = 1, number_of_comp
      rho0(i,j) = X(j) * (rhos0/(si*d_min**3)) ! rho0, X(j) are global varies.rho0 here is reduced density
      rho_real(i,j) = rho0(i,j)*1.0E30/VNMOL ! here unit is mol/m^3
      write(*,*) 'rhos0 is:', rhos0, 'rho0 is:', rho0(i,j), 'real density(mol/m^3):', rho_real(i,j)
    end do ! because si here is a reduced parameter
    ! here si is 1, becaust reduced sigma_gas is 1, we make rho0 become reduced density

    Ar = re_Ar(rhos0, Ts) ! Ar is the reduced excess Helmholtz free energy
    Ur = re_Ur(rhos0, Ts) ! Ur is the excess internal energy

    do j= 1, number_of_comp
      sigxrho(j) = 0.0
      epsixrho(j) = 0.0

      do k = 1, number_of_comp
        epsixrho(j) = epsixrho(j)+X(k)*sqrt(epsilon_gas(j)*epsilon_gas(k)) &
        *(0.5*(sigma_gas(j)+sigma_gas(k)))**3
        sigxrho(j) = sigxrho(j) + X(k)*((sigma_gas(j)+sigma_gas(k))/2.0)**3
      end do
      sigxrho(j) = (sigxrho(j) - si) * 2 /(rhos0/si)
      epsixrho(j) = (epsixrho(j)/si-ep)*2/(rhos0/si) - sigxrho(j)*ep/si
      Arxrho(j) = (Ps/rhos0**2-Ts/rhos0)*(si+(rhos0/si)*sigxrho(j)) - (Ar-Ur)*epsixrho(j)/ep
      Miu(i,j) = (Ar*ep + Ar*rhos0/si*epsixrho(j) + rhos0/si*ep*Arxrho(j))/temperature ! unit of Miu(i,j) depends on the epsilon_gas,
      !now it is unitless after we divided Temperature
    end do

    ! out put bulk density
    write(22,229) i, pres(i), (rho0(i,j), j=1,number_of_comp), (rho_real(i,j),j=1,number_of_comp), (Miu(i,j), j=1, number_of_comp), nu, errps

  end do
  close(22)
  229 format(I5,4X,F13.6,4X,3(F15.6,4X),I5,4X,F13.6)
end subroutine

! function re_Ur get the reduced internal energy
function re_Ur(rhos, Ts)
  use local_module_parameters
  implicit none
  real*8 rhos, Ts, re_Ur
  real*8 c(8), d(6), G(6), F
  integer i

  c(1) = x2*sqrt(Ts)/2 + x3 + 2*x4/Ts + 3*x5/Ts**2
  c(2) = x7 + 2*x8/Ts + 3*x9/Ts**2
  c(3) = x11 + 2*x12/Ts
  c(4) = x13
  c(5) = 2*x14/Ts + 3*x15/Ts**2
  c(6) = 2*x16/Ts
  c(7) = 2*x17/Ts + 3*x18/Ts**2
  c(8) = 3*x19/Ts**2
  d(1) = 3*x20/Ts**2 + 4*x21/Ts**3
  d(2) = 3*x22/Ts**2 + 5*x23/Ts**4
  d(3) = 3*x24/Ts**2 + 4*x25/Ts**3
  d(4) = 3*x26/Ts**2 + 5*x27/Ts**4
  d(5) = 3*x28/Ts**2 + 4*x29/Ts**3
  d(1) = 3*x30/Ts**2 + 4*x31/Ts**3 + 5*x32/Ts**4

  F=exp(-3*rhos**2)
  G(1) = (1-F)/(2*3)
  G(2) = -(F*rhos**2 - 2*G(1))/(2*3)
  G(3) = -(F*rhos**4 - 4*G(2))/(2*3)
  G(4) = -(F*rhos**6 - 6*G(3))/(2*3)
  G(5) = -(F*rhos**8 - 8*G(4))/(2*3)
  G(6) = -(F*rhos**10 - 10*G(5))/(2*3)

  re_Ur = 0.0
  do i = 1, 8
    re_Ur = re_Ur + c(i)*rhos**i/i
  end do
  do i =1, 6
    re_Ur = re_Ur + d(i)*G(i)
  end do

end function

! re_pres is reduced pressure
function re_pres(rhos, Ts)
  use local_module_parameters
  implicit none
  real*8 rhos, Ts, re_pres
  real*8 a(8), b(6), F
  integer i

  F=exp(-3*rhos**2)
  a(1) = x1*Ts + x2*sqrt(Ts) + x3 + x4/Ts + x5/Ts**2
  a(2) = x6*Ts + x7 + x8/Ts + x9/Ts**2
  a(3) = x10*Ts + x11 + x12/Ts
  a(4) = x13
  a(5) = x14/Ts + x15/Ts**2
  a(6) = x16/Ts
  a(7) = x17/Ts + x18/Ts**2
  a(8) = x19/Ts**2
  b(1) = x20/Ts**2 + x21/Ts**3
  b(2) = x22/Ts**2 + x23/Ts**4
  b(3) = x24/Ts**2 + x25/Ts**3
  b(4) = x26/Ts**2 + x27/Ts**4
  b(5) = x28/Ts**2 + x29/Ts**3
  b(6) = x30/Ts**2 + x31/Ts**3 + x32/Ts**4

  re_pres = 0.0
  do i=1,8
    re_pres = re_pres+a(i)*rhos**(i+1)
  end do
  do i=1,6
    re_pres = re_pres + F*b(i)*rhos**(2*i+1)
  end do
   re_pres = re_pres + rhos*Ts
   !re_pres = re_pres*kB*1.0E30 !transfer to international unit
end function

! re_Ar is reduced excess Helmholtz free energy
function re_Ar(rhos, Ts)
  use local_module_parameters
  implicit none
  real*8 rhos, Ts
  real*8 a(8), b(6), G(6), F, re_Ar
  integer i

  a(1) = x1*Ts + x2*sqrt(Ts) + x3 + x4/Ts + x5/Ts**2
  a(2) = x6*Ts + x7 + x8/Ts + x9/Ts**2
  a(3) = x10*Ts + x11 + x12/Ts
  a(4) = x13
  a(5) = x14/Ts + x15/Ts**2
  a(6) = x16/Ts
  a(7) = x17/Ts + x18/Ts**2
  a(8) = x19/Ts**2
  b(1) = x20/Ts**2 + x21/Ts**3
  b(2) = x22/Ts**2 + x23/Ts**4
  b(3) = x24/Ts**2 + x25/Ts**3
  b(4) = x26/Ts**2 + x27/Ts**4
  b(5) = x28/Ts**2 + x29/Ts**3
  b(6) = x30/Ts**2 + x31/Ts**3 + x32/Ts**4

  F=exp(-3*rhos**2)
  G(1) = (1-F)/(2*3)
  G(2) = -(F*rhos**2 - 2*G(1))/(2*3)
  G(3) = -(F*rhos**4 - 4*G(2))/(2*3)
  G(4) = -(F*rhos**6 - 6*G(3))/(2*3)
  G(5) = -(F*rhos**8 - 8*G(4))/(2*3)
  G(6) = -(F*rhos**10 - 10*G(5))/(2*3)

  re_Ar = 0.0
  do i = 1, 8
    re_Ar = re_Ar + a(i)*rhos**i/i
  end do
  do i = 1, 6
    re_Ar = re_Ar + b(i)*G(i)
  end do
end function
