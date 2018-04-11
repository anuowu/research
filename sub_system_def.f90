! this subroutine gives the gas and solid properity, interaction parameters
! and it defines the structure of slit pore
! it defines the adsorption conditions
subroutine system_def
  use global_module_energy
  implicit none

  real*8 ls, lt
  integer i, iflag,ii

  ! input the parameters
  open(42, file="input_parameters.dat")
  read(42,*)
  read(42,*)
  read(42,*) name_of_stru, dimension, length_of_stru, temperature
  read(42,*)
  read(42,*) type_of_wall, width_of_wall, number_density_wall, name_of_wall, delta_wall !number_density_wall should be the origin density,do not consider relative dimension
  read(42,*)
  read(42,*) type_of_fluids, number_of_comp ! number_of_comp < compmax=5
  ! 此处可以加个判断组分是否超过设置的最大组分
  read(42,*)
  do i = 1,number_of_comp
    read(42,*) name_of_comp(i), weight_of_comp(i), sigma_gas(i), epsilon_gas(i), dhs_gas(i), sigma_gw(i), epsilon_gw(i),X(i)
  end do
  read(42,*)
  read(42,*)
  read(42,*) appr_of_hs, appr_of_attr, cut_off_gg, cut_off_wg
  read(42,*)
  read(42,*) Nmax, errmax, num_deltx
  close(42)

  write(*,*) "this is ", number_of_comp, "components adsorption happened at ", name_of_stru
  write(*,*) "the length of stru is:", length_of_stru
  do ii = 1, number_of_comp
    write(*,*) "the gas are: ", name_of_comp(ii)
    write(*,*) "the gas parameter is:", sigma_gas(ii), epsilon_gas(ii), X(ii)
  end do
  write(*,*) "the approximation of attraction between fluids is: ", appr_of_attr
  write(*,*) "treat of wall: ", type_of_wall
  write(*,*) 'cut_off_gg', cut_off_gg, 'cut_off_wg', cut_off_wg

  open(43, file="input_pressure.dat")
  num_press=0
  iflag=1
  read(43,*)
  do while(iflag==1)
    num_press=num_press+1
    read(43,*) nb(num_press), pres(num_press), c, pica_factor(num_press), IFDEN(num_press), ITERA(num_press)
    if ( nb(num_press)==0 ) then
      iflag = 0
      num_press=num_press-1
    end if
  end do
  close(43)

  write(*,*) "adsorption isotherm has:", num_press, " data points"

  sigmamax = sigma_gas(1)
  do i =1, number_of_comp
    if(sigma_gas(i) > sigmamax) then
      sigmamax = sigma_gas(i)
    end if
  end do

  ! give a relative dimension
  d_min = dhs_gas(1)
  do i = 2, number_of_comp
    if ( dhs_gas(i) < d_min ) then
      d_min = dhs_gas(i)
    end if
  end do
  do i = 1, number_of_comp
    dhs_gas(i) = dhs_gas(i)/d_min
    sigma_gas(i) = sigma_gas(i)/d_min
    sigma_gw(i) = sigma_gw(i)/d_min
  end do
  length_of_stru = length_of_stru/d_min
  width_of_wall = width_of_wall/d_min
  delta_wall = delta_wall/d_min
  write(*,*) 'sigmamax is:', sigmamax, 'd_min is:', d_min
  sigmamax = sigmamax/d_min
  write(*,*) 'sigmamax is:', sigmamax
  deltx = 1.0 / num_deltx
  write(*,*) 'deltx=', deltx
  lt = length_of_stru/2.0+width_of_wall
  ls = length_of_stru/2.0
  if ( ls>30 ) then
    xnf = ls - 30.0
  else
    xnf = 0.0
  end if
  NBt = int((ls - xnf)/deltx)
  deltx = (ls-xnf)/NBt
  NL = int((lt-xnf)/deltx)

  write(*,*) 'here the step length and total number are:', deltx, NL

  ! calculate the saturation properity: P0, rho_bulk

end subroutine
