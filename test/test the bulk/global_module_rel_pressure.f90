module global_module_rel_pressure
  ! this module define varies in relative pressure
  integer,parameter:: num_press_max=200
  integer num_press, nb(num_press_max), c, IFDEN(num_press_max)
  real*8 rel_prees(num_press_max), press(num_press_max), pica_factor(num_press_max)
  real*8 P0, rho_bulk(num_press)
  common num_press, nb, c, IFDEN
  common rel_prees, press, pica_factor, P0, rho_bulk

end module global_module_rel_pressure
