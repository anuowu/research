! subroutine com_hs_pot calculate the all position's chemical potential of HS
! terms at certain pressure, cp_hs(j,i)
subroutine com_hs_pot
  use global_module_energy
  use local_module_mfmt
  implicit none

  do i = 0, NL
    xi = xnf+i*deltx
    xc = 100

    n_0 = 0.0
    n_1=0.0
    n_2=0.0
    n_3=0.0
    n_v1=0.0
    n_v2=0.0
    do j=1, number_of_comp

      n2(j) = 0.0
      n3(j) = 0.0
      nv2(j) =0.0
      ddx = sigma_gas(j)/(2*xc)
      do k = -xc, xc
        xxi= xi + k*ddx
        xz = k*ddx
        call qrho(rho, abs(xxi), j) ! calculate the density of
        n2(j) = n2(j) + ddx*rho
        n3(j) = n3(j) + ddx*rho*( (sigma_gas(j)/2.0)**2 - xz**2 )
        nv2(j) = nv2(j) + ddx*rho*(xz)
      end do ! z' relative position
      n2(j) = pi*sigma_gas(j)*n2(j)
      n3(j) = pi *n3(j)
      nv2(j) = 2*pi*nv2(j)
      n0(j) = n2(j) /(pi*sigma_gas(j)**2)
      n1(j) = n2(j) /(2*pi*sigma_gas(j))
      nv1(j) = nv2(j) /(2*pi*sigma_gas(j))

      n_2 = n_2 + n2(j)
      n_3 = n_3 + n3(j)
      n_v2 = n_v2 + nv2(j)
      n_0 = n_0 + n0(j)
      n_1 = n_1 + n1(j)
      n_v1 = n_v1 +nv1(j)

    end do ! components
    fi0(i) = -log(1-n_3) !delta_phi/delta_n0
    fi1(i) = n_2/(1-n_3) !delta_phi/delta_n1
    fi2(i) = n_1/(1-n_3) + (n_2**2-n_v2**2)/(12*pi*n_3)*(log(1-n_3)/n_3 + 1/(1-n_3)**2) !delta_phi/delta_n2
    fi3(i) = n_0/(1-n_3) + (n_1*n_2 - n_v1*n_v2)/(1-n_3)**2 - (n_2**2-3*n_2*n_v2**2) &
          *(log(1-n_3)/(18*pi*n_3**3) + (1-3*n_3+(1-n_3)**2)/(36*pi*n_3**2*(1-n_3)**3)) !delta_phi/delta_n3
    fiv1(i) = - n_v2/(1-n_3) !delta_phi/delta_nv1
    fiv2(i) = - n_v1/(1-n_3) - (log(1-n_3)/n_3 + 1/(1-n_3)**2)*(n_2*n_v2)/(6*pi*n_3) !delta_phi/delta_nv2
    ! phi(i) is the energy density
    phi(i) = -n_0*log(1-n_3) + n_2*n_1/(1-n_3) + (log(1-n_3)/n_3+1/(1-n_3)**2)/(36*pi*n_3)*n_2**2 &
           -n_v1*n_v2/(1-n_3) - (n_2*n_v2*n_v2)*(log(1-n_3)/n_3 + 1/(1-n_3)**2)/(12*pi*n_3)

  end do ! z position

  do i = 0, NL
    xi = xnf + i*deltx
    xc = 100

    do j = 1, number_of_comp
      ddx = sigma_gas(j) / (2*xc)
      dfn0 = 0.0
      dfn1 = 0.0
      dfn2 = 0.0
      dfn3 = 0.0
      dfnv1 = 0.0
      dfnv2 = 0.0
      do k = -xc, xc
        xxi = xi + k*ddx
        xz = k*ddx
        call qfi(fi, abs(xxi))
        dfn0 = dfn0 + fi(0) *ddx
        dfn1 = dfn1 + fi(1) *ddx
        dfn2 = dfn2 + fi(2) *ddx
        dfn3 = dfn3 + fi(3) *ddx*((sigma_gas(j)/2)**2-(xz)**2)
        dfnv1 = dfnv1 - fi(4) * ddx * (xz)
        dfnv2 = dfnv2 - fi(5) * ddx * (xz)
      end do
      dfn0 = dfn0 /sigma_gas(j)
      dfn1 = dfn1 /2
      dfn2 = dfn2 * pi *sigma_gas(j)
      dfn3 = dfn3 * pi
      dfnv1 = dfnv1/sigma_gas(j)
      dfnv2 = dfnv2*2*pi
      cp_hs(j,i) = dfn0+dfn1+dfn2+dfn3+dfnv1+dfnv2 !here cp_hs is unitless

    end do ! end components

  end do ! end position
end subroutine

subroutine qrho(rho, xxi, j)
  use global_module_energy ! include P(J,I), xnf, deltx,NL
  implicit none
  integer j, L
  real*8 rho, xxi, xi, x1, x2, xx1, xx2

  if ( xxi <= xnf ) then
    rho = P(j,0)
  else if (xxi>=(xnf+NL*deltx) ) then
    rho = P(j, NL)
  else
    L = int((xxi - xnf)/deltx)
    xi = xxi - xnf
    x1 = L*deltx
    x2 = (L+1)*deltx
    xx1 = (x2 - xi)/deltx
    xx2 = (xi-x1)/deltx
    rho = xx1*P(j,L) + xx2*P(j,L+1)
  end if

end subroutine

subroutine qfi(fi, xxi)
  use global_module_energy ! including fi0, fi1,fi2, fi3,fiv1,fiv2,xnf, NL，deltx
  implicit none
  real*8 fi(0:5), xxi, xi, x1, x2, xx1, xx2
  integer L

  if ( xxi <= xnf ) then
    fi(0) = fi0(0)
    fi(1) = fi1(0)
    fi(2) = fi2(0)
    fi(3) = fi3(0)
    fi(4) = fiv1(0)
    fi(5) = fiv2(0)
  else if(xxi >= (xnf+NL*deltx)) then
    fi(0) = fi0(NL)
    fi(1) = fi1(NL)
    fi(2) = fi2(NL)
    fi(3) = fi3(NL)
    fi(4) = fiv1(NL)
    fi(5) = fiv2(NL)
  else
    L = int((xxi -xnf)/deltx)
    xi = xxi - xnf
    x1 = L*deltx
    x2 = (L+1)*deltx
    xx1 = (x2-xi)/deltx
    xx2 = (xi-x1)/deltx
    fi(0) = fi0(L)*xx1 + fi0(L+1)*xx2
    fi(1) = fi1(L)*xx1 + fi1(L+1)*xx2
    fi(2) = fi2(L)*xx1+fi2(L+1)*xx2
    fi(3) = fi3(L)*xx1+fi3(L+1)*xx2
    fi(4) = fiv1(L)*xx1+fiv1(L+1)*xx2
    fi(5) = fiv2(L)*xx1+fiv2(L+1)*xx2
  end if

end subroutine
