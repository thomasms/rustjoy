! Fuzz driver for the BROADR numeric kernels: prints, in hex, NJOY's
! erfc(x) and hnabb(n,a,b) for n=0..4, for each (x,a,b) triple read from
! standard input (little-endian real(8) triples, see numcheck_gen.py).
! Build against a copy of broadr.f90 whose module makes hnabb public.
program numcheck_broadr
  use locale
  use mathm, only: erfc
  use broadm, only: hnabb
  implicit none
  real(kr)::x,a,b
  integer::ios,n
  open(10,file='/dev/stdin',access='stream',form='unformatted')
  do
    read(10,iostat=ios) x,a,b
    if (ios.ne.0) exit
    write(*,'(z16.16,5(1x,z16.16))') erfc(x),(hnabb(n,a,b),n=0,4)
  enddo
end program
