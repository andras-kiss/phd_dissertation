program diffusion
implicit none

integer :: stat, counter
real t, rc, e0, conv

counter=0


rc=0.87
open(1,file='500ms.txt')
open(2,file='500ms_deconvoluted.txt')
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

rc=0.81
open(1,file='1000ms.txt')
open(2,file='1000ms_deconvoluted.txt')
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

rc=0.39
open(1,file='5000ms.txt')
open(2,file='5000ms_deconvoluted.txt')
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
