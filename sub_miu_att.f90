subroutine miu_att
  use global_module_energy
  implicit none
  integer xc,j,k,m,i
  real*8 x_cut, ddx, sigma_jk, eplison_jk, r1, uffc, xi, xm, rr, rs, uff,r_cut

  miu_attz(:,:,:) =0.0
  r_cut = cut_off_gg
  x_cut = sigmamax*r_cut
  xc = int(x_cut/(0.999*deltx))
  ddx= x_cut/xc
  write(*,*) 'x_cut=', x_cut, 'xc=', xc, 'ddx=', ddx, 'r_cut is', r_cut
  do j = 1, number_of_comp
    do k = 1, number_of_comp
      write(*,*) '****'
      sigma_jk = (sigma_gas(j)+sigma_gas(k))/2.0
      eplison_jk = sqrt(epsilon_gas(j)*epsilon_gas(k))
      r1 = 2.**(1./6.)*sigma_jk
      uffc = 4.0*eplison_jk*((sigma_jk/x_cut)**12-(sigma_jk/x_cut)**6)
      write(*,*) 'start i'
      do i = 0, xc
        xi = i*ddx
        !write(*,*) 'start m'
        do m = 0,xc
          xm = m*ddx
          rr = sqrt(xi**2+xm**2)
          rs = rr/sigma_jk
          if ( rr < r1 ) then
            uff = -eplison_jk -uffc
          else if(rr < x_cut) then
            uff = 4.0*eplison_jk*((1.0/rs)**12 - (1.0/rs)**6)-uffc
          else
            uff = 0.0
          end if
          miu_attz(j,k,i) = miu_attz(j,k,i) + uff*xm
        end do ! integral of miu_att as a function of z
        miu_attz(j,k,i) = miu_attz(j,k,i)*2*pi*ddx
        !write(*,*) 'miu_attz(j,k,i) is',miu_attz(j,k,i)
      end do !end i
    end do !end component k
  end do ! end component j

end subroutine
