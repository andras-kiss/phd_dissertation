program diffusion
implicit none

integer :: stat, counter
real t, rc, e0, conv

counter=0


rc=0.7
open(1,file='400ms.txt')
open(2,file='400ms_deconvoluted.txt')
read(1, *) t, e0
write(2, *) t, e0
do
   read(1, *, iostat=stat) t, conv
   if (stat /= 0) exit
   write(2, *) t, ((conv - e0*rc)/(1-rc))
   e0=conv
end do
close(1) 
close(2)

rc=0.6
open(1,file='1900ms.txt')
open(2,file='1900ms_deconvoluted.txt')
read(1, *) t, e0
write(2, *) t, e0
do
   read(1, *, iostat=stat) t, conv
   if (stat /= 0) exit
   write(2, *) t, ((conv - e0*rc)/(1-rc))
   e0=conv
end do
close(1) 
close(2)

rc=0.3
open(1,file='4900ms.txt')
open(2,file='4900ms_deconvoluted.txt')
read(1, *) t, e0
write(2, *) t, e0
do
   read(1, *, iostat=stat) t, conv
   if (stat /= 0) exit
   write(2, *) t, ((conv - e0*rc)/(1-rc))
   e0=conv
end do
close(1) 
close(2)

end program diffusion
