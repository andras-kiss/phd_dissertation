program diffusion
implicit none

integer :: i,j,x,y,direction
real, dimension(0:100,0:100) :: sim, deconvoluted, eq
real rc, conv, e0

! mérés beolvasása
open(1,file='1.txt')
do x=1, 945
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(-i/50, -j/50)=conv 
end do
close(1)

rc=0.9

! MEANDER deconvolution
e0=sim(0,0)
direction=-1
open(1,file='1_deconvoluted.txt')
do y=0, 20, 1
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 44, 1
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			write(1, *) -x*50, -y*50, deconvoluted(x,y)
		end do
	else
		do x=44, 0, -1
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			write(1, *) -x*50, -y*50, deconvoluted(x,y)
		end do
	endif
end do
close(1)

end program diffusion
