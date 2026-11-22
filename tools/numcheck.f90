! Driver: for each (real(8) x, integer n) record on stdin (stream), print the
! results of NJOY's a11, a10, sigfig, x**n and the e11.0 read-back of a11.
program numcheck
  use locale
  use endf, only: a11
  use util, only: a10, sigfig
  implicit none
  real(kr)::x,r
  integer::n,ios
  character(11)::h
  character(10)::t
  open(10,file='/dev/stdin',access='stream',form='unformatted')
  do
    read(10,iostat=ios) x,n
    if (ios.ne.0) exit
    h='           '
    call a11(x,h)
    call a10(x,t)
    read(h,'(e11.0)') r
    write(*,'(a,"|",a,"|",8(z16.16,"|"),1p,e12.5,"|",0p,f10.3)') h,t,&
      sigfig(x,7,0),sigfig(x,9,-1),x**n,x**2,x**3,x**5,x**13,r,x,x
  enddo
end program
