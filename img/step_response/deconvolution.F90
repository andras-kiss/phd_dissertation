program diffusion
implicit none

integer :: stat, counter
real t, rc, e0, conv, sel, sey

counter=0
rc=0.7

open(1,file='transient3.txt')
open(2,file='transient3_deconvoluted07.txt')
!open(3,file='transient3.txt')
read(1, *) t, e0
do
   read(1, *, iostat=stat) t, conv
   if (stat /= 0) exit
   write(2, *) t, ((conv - e0*rc)/(1-rc))
   e0=conv
!   counter=counter+1
!   if (counter==21) then
!     write(3,*) t, conv
!     counter=0
!   endif
end do
close(1) 
close(2)
!close(3)


sel=0
sey=0

open(1,file='95.txt')
read(1, *) t, e0
do
   read(1, *, iostat=stat) conv
   if (stat /= 0) exit
   sel=sel+(conv+280)**2
end do
close(1) 


print *,sel/107
end program diffusion
