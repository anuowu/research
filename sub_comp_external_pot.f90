! Attention:
! + here the unit of cp_ext depends on epsilon_gw, and used in exp() must be unitless
! + we multiply d_min, because delta_wall and sigma_gw are reduced by d_min, number_density_wall is
!   original unit.
subroutine comp_external_pot
  use global_module_energy
  implicit none
  real*8 xi, re_position, position
  integer i, j

  ! need the definition of the system
  ! need parameters, if it is a LJ wall
  open(41,file='external_potential.dat')
  do j=1, number_of_comp
    write(*,*) 'NL=', NL
    do i=0,NL ! start from 0 or 1
      xi = xnf+i*deltx
      if ( xi >= length_of_stru/2.0 ) then
        cp_ext(j,i) = 1.0E30 ! here take care of the unit
      else
        position = length_of_stru/2.0 - xi
        re_position = position/sigma_gw(j)
        if ( re_position > cut_off_wg ) then
          cp_ext(j,i) = 0.0
        else !!cp_ext(j,i) unit depends on epsilon_gw
          cp_ext(j,i) = 2.0*pi*number_density_wall*delta_wall*sigma_gw(j)**2*epsilon_gw(j)
          cp_ext(j,i) = cp_ext(j,i)*(0.4*re_position**(-10) - re_position**(-4) &
          - sigma_gw(j)**4/(3.0*delta_wall*(position+0.61*delta_wall)**3))
        end if

        position = length_of_stru/2.0+xi
        re_position = position/sigma_gw(j)
        if ( re_position > cut_off_wg ) then
          cp_ext(j,i) = cp_ext(j,i)+0.0
        else
          cp_ext(j,i) = cp_ext(j,i)+2.0*pi*number_density_wall*delta_wall*sigma_gw(j)**2 &
          *epsilon_gw(j)*(0.4*re_position**(-10) - re_position**(-4) &
          - sigma_gw(j)**4/(3.0*delta_wall*(position+0.61*delta_wall)**3))
        end if
        cp_ext(j,i)=cp_ext(j,i)*d_min**3/temperature !! becaust here the sigma_gw is reduced diameter
        ! here cp_ext(j,i) is unitless, canbe used in iteration equation in main.f90
        write(41,*) j, xi, cp_ext(j,i)
      end if
      !write(*,*) '5 number_of_comp', number_of_comp
    end do
    !write(*,*) '6 number_of_comp', number_of_comp
  end do
  !write(*,*) '7 number_of_comp', number_of_comp
  close(41)
end subroutine
