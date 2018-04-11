module global_module_energy
  use global_module_parameters

  ! define external varies
  !integer,parameter::max_NL = 10000
  real*8, save:: cp_ext(compmax, 0:max_NL)
  !common cp_ext

  ! define bulk property
  real*8,save:: Miu(max_num_press, compmax), rho0(max_num_press, compmax), X(compmax)
  !common Miu, rho0, X

  !define density
  real*8,save:: P(compmax, 0:max_NL), PP(compmax, 0:max_NL)
  !common P, PP

  ! define hs terms
  real*8, save:: fi0(max_NL), fi1(max_NL), fi2(max_NL), fi3(max_NL), fiv1(max_NL), fiv2(max_NL) ! varies for MFMT
  real*8, save:: cp_hs(compmax, 0:max_NL), phi(0:max_NL)
  !common fi0, fi1, fi2, fi3, fiv1, fiv2
  !common cp_hs, phi

  ! define att terms
  real*8,save:: miu_attz(compmax, compmax, 0:max_NL), cp_att(compmax, 0:max_NL)
  real*8,save:: energy_att
  !common miu_attz, cp_att, energy_att

  !define the adsorption amount
  real*8,save:: ADSP0(compmax, max_num_press), ADSP1(compmax, max_num_press)
  real*8,save:: ADSP2(compmax, max_num_press), ADSP3(compmax, max_num_press)
  real*8,save:: CP(compmax, max_num_press), GRP(compmax, max_num_press)
  !common ADSP0, ADSP1, ADSP2, ADSP3, CP, GRP

end module global_module_energy
