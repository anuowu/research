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
  open(41,file='out_external_potential.dat')
  write(41, "(A20,4X,A12,4X,A20)") "type of components", "position", "unitless potential"
  do j=1, number_of_comp
    write(*,*) 'NL=', NL
    do i=0,NL ! start from 0 or 1
      xi = xnf+i*deltx
      if ( xi >= length_of_stru/2.0 ) then
        cp_ext(j,i) = 200 ! here take care of the unit is kBT
      else
        position = length_of_stru/2.0 - xi
        re_position = position/sigma_gw(j)
        write(*,*) "xi is:", xi, "re_position is:", re_position, 'cut_off_wg:', cut_off_wg
        if ( re_position > cut_off_wg ) then
          cp_ext(j,i) = 0.0
        else !!cp_ext(j,i) unit depends on epsilon_gw
          cp_ext(j,i) = 2.0*pi*number_density_wall*delta_wall*sigma_gw(j)**2*d_min**3*epsilon_gw(j)/temperature
          cp_ext(j,i) = cp_ext(j,i)*(0.4*re_position**(-10) - re_position**(-4) &
          - sigma_gw(j)**4/(3.0*delta_wall*(position+0.61*delta_wall)**3))
        end if

        position = length_of_stru/2.0+xi
        re_position = position/sigma_gw(j)
        if ( re_position > cut_off_wg ) then
          cp_ext(j,i) = cp_ext(j,i)+0.0
        else
          cp_ext(j,i) = cp_ext(j,i)+2.0*pi*number_density_wall*delta_wall*sigma_gw(j)**2*d_min**3 &
          *(epsilon_gw(j)/temperature)*(0.4*re_position**(-10) - re_position**(-4) &
          - sigma_gw(j)**4/(3.0*delta_wall*(position+0.61*delta_wall)**3))
        end if
        !! becaust here the sigma_gw is reduced diameter
        ! here cp_ext(j,i) is kBT, canbe used in iteration equation in main.f90
      end if
      write(41,419) j, xi, cp_ext(j,i)
      write(*,*) j, xi, cp_ext(j,i)
    end do
    !write(*,*) '6 number_of_comp', number_of_comp
  end do
  !write(*,*) '7 number_of_comp', number_of_comp
  close(41)
  419 format(I5, 5X, F13.6, 5X, E12.6)
end subroutine
