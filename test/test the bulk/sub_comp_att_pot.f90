subroutine com_att_pot
  use global_module_energy
  implicit none
  real*8 x_cut, ddx, xz, xzz, kz, r_cut
  integer j, k, m, xc, ikz, iz, mj

  r_cut = cut_off_gg
  x_cut = sigmamax*r_cut
  write(*,*) 'r_cut =', r_cut, 'x_cut=', x_cut, 'cut_off_gg', cut_off_gg
  xc = int(x_cut/(0.999*deltx))
  write(*,*) 'xc=', xc, 'deltx=', deltx, 'x_cut=', x_cut
  ddx= x_cut/xc
  write(*,*) 'ddx =', ddx
  cp_att(:,:) = 0.0
  do j=1, number_of_comp
    !write(*,*) 'piosition data point is:', NL
    do k=0, NL
      !write(*,*) 'start the', k, 'position'
      xz = xnf + k*deltx
      !cp_att(j,k) = 0.0
      !write(*,*) 'xc=', xc
      do m = -xc, xc
        xzz = xz + m*ddx
        kz = m*ddx
        ikz = int(abs(kz)/ddx)+1
        iz = int(xzz/deltx)+1
        do mj = 1, number_of_comp
          cp_att(j,k) = cp_att(j,k) + 0.5*P(mj, iz)*miu_attz(j,mj,ikz)*deltx
        end do
      end do

    end do
  end do
end subroutine
