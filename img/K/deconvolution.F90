program diffusion
implicit none

integer :: i, j, stat
real rc, e0, conv
real t

rc=0.82
open(1,file='14103104.txt')
open(2,file='14103104_deconvoluted.txt')
read(1, *) i, j, e0
do
   read(1, *, iostat=stat) i, j, conv
   if (stat /= 0) exit
   write(2, *) i, j, ((conv - e0*rc)/(1-rc))
   e0=conv
end do
close(1) 
close(2)

rc=0.85
open(1,file='14103107.txt')
open(2,file='14103107_deconvoluted.txt')
read(1, *) i, j, e0
do
   read(1, *, iostat=stat) i, j, conv
   if (stat /= 0) exit
   write(2, *) i, j, ((conv - e0*rc)/(1-rc))
   e0=conv
end do
close(1) 
close(2)


end program diffusion
