program diffusion
implicit none

integer :: i,j,x,y
real, dimension(0:44,0:20) :: sim, deconvoluted, eq
real rc, conv, e0

! mérés beolvasása
open(1,file='1.txt')
do x=1, 945
	read(1, *) i, j, conv
	sim(-i/50, -j/50)=conv 
end do
close(1)

! mérés beolvasása
open(1,file='1_check.txt')
do y=0,20
	do x=0,44
		write(1, *) x*50, y*50, sim(x,y)
	end do
end do
close(1)

rc=0.8
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='1_deconvoluted.txt')
do y=0,20 
	e0=sim(0,y)
	do x=0, 44
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x*50, y*50, deconvoluted(x,y)
	end do
end do
close(1)

end program diffusion
