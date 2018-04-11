program main
  implicit none
  real*8 pres(200), pica_factor(200),b,d
  integer num_press, iflag, nb(200), c, IFDEN(200), ITERA(200),a,e,f
  open(43, file="input_pressure.dat")
  num_press=0
  iflag=1
  read(43,*)
  do while(iflag==1)
    num_press=num_press+1
    read(43,*) a, b, c, d, e, f
    nb(num_press) = a
    pres(num_press) = b
    pica_factor(num_press) = d
    IFDEN(num_press) = e
    ITERA(num_press) = f
    if ( nb(num_press)==0 ) then
      iflag = 0
      num_press=num_press-1
    end if
  end do
  close(43)
end
