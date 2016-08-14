program diffusion
implicit none

integer :: i, j, stat
real rc, e0, conv
real t

rc=0.85
open(1,file='raw.txt')
open(2,file='raw_deconvoluted085.txt')
read(1, *) i, j, e0
do
   read(1, *, iostat=stat) i, j, conv
   if (stat /= 0) exit
   write(2, *) i, j, ((conv - e0*rc)/(1-rc))
   e0=conv
end do
close(1) 
close(2)

rc=0.75
open(1,file='raw.txt')
open(2,file='raw_deconvoluted075.txt')
read(1, *) i, j, e0
do
   read(1, *, iostat=stat) i, j, conv
   if (stat /= 0) exit
   write(2, *) i, j, ((conv - e0*rc)/(1-rc))
   e0=conv
end do
close(1) 
close(2)

rc=0.95
open(1,file='raw.txt')
open(2,file='raw_deconvoluted095.txt')
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
