module global_module_parameters
  use global_module_constants

  ! defination of the structure
  character(10), save:: name_of_stru, name_of_wall, type_of_wall, dimension
  real*8,save:: length_of_stru, temperature, number_density_wall, width_of_wall, delta_wall
  !common name_of_stru, name_of_wall, type_of_wall, dimension
  !common length_of_stru, temperature, number_density_wall, width_of_wall, delta_wall

  ! defination of fluids
  ! integer, parameter::compmax=5
  character(10), save:: type_of_fluids, name_of_comp(compmax), appr_of_hs, appr_of_attr
  integer, save:: number_of_comp, cut_off_gg, cut_off_wg ! number_of_comp <= compmax
  real*8, save:: weight_of_comp(compmax), sigma_gas(compmax), epsilon_gas(compmax), dhs_gas(compmax)
  real*8,save:: sigma_gw(compmax),epsilon_gw(compmax)
  !common type_of_fluids, name_of_comp, appr_of_hs, appr_of_attr
  !common number_of_comp
  !common weight_of_comp, sigma_gas, epsilon_gas, dhs_gas, sigma_gw, epsilon_gw

  ! defination of iteration
  integer,save:: Nmax, num_deltx
  real*8,save:: errmax
  !common Nmax, errmax, num_deltx

  ! definition of the system step length
  integer, save:: NL, NBt
  real*8,save:: xnf, deltx

  ! define the pressure
  ! integer, parameter:: max_num_press = 200
  integer, save:: num_press, nb(max_num_press), c, IFDEN(max_num_press),ITERA(max_num_press)
  real*8, save:: pres(max_num_press), pica_factor(max_num_press)

  real*8, save:: d_min, sigmamax
  !common d_min, sigmamax
end module global_module_parameters
