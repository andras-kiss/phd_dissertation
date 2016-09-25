program diffusion
implicit none

integer :: stat
real t, tau1, tau2, rc, e0, conv

tau1=9
tau2=0.8

open(1,file='original_1sec_cut.txt')
open(2,file='deconvoluted.txt')
open(3,file='absolute_deconvoluted.txt')
read(1, *) t, e0
write(2, *) t, e0
write(3, *) t, e0
do
   read(1, *, iostat=stat) t, conv
   if (stat /= 0) exit
   write(2, *) t, (conv - e0*(exp(-t/tau1)+exp(-t/tau2))/2)/(1-(exp(-t/tau1)+exp(-t/tau2))/2)
   !e0=conv
   print *, t, e0
end do
close(1) 
close(2)
close(3)
end program diffusion
