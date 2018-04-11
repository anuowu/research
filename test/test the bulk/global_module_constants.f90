module global_module_constants
  integer, parameter::pi = 3.1415926
  integer, parameter::max_NL = 9999
  integer, parameter::compmax = 5
  integer, parameter::max_num_press = 200

  ! define energy constants
  real*8, parameter::kB = 1.38064E-23 ! unit: J/K
  real*8, parameter::R = 8.3144598 ! unit: J/(K*MOL)
  real*8, parameter::plank = 6.6260693E-34 ! iso unit
  real*8, parameter::VNMOL = 6.0221367E23
end module
