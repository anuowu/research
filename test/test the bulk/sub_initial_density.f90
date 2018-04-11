subroutine initial_density(P)
  use global_module_energy
  implicit none
  integer j,i

  do j = 1, compmax
    do i =0, max_NL
      P(j,i) = 0.0
    end do
  end do

end subroutine
