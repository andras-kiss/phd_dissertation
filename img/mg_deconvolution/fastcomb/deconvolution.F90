program diffusion
implicit none

integer :: i,j,x,y,direction
real, dimension(0:100,0:100) :: sim, deconvoluted, eq
real, dimension(0:200,0:200) :: sim4000, deconvoluted4000
real rc, conv, e0, error_b, error_a
real, dimension(1:341,1:3) :: arcpol, arcpoldec
open(4,file='error.txt')

! mérés beolvasása
open(1,file='14111805.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

rc=0.8
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='14111805_deconvoluted.txt')
open(2,file='14111805_dd.txt')
do y=0, 100, 5
	e0=sim(0,y)
	do x=0, 100, 5
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x/5, y/5, deconvoluted(x,y)
		write(2, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
	end do
end do
close(1)
close(2)

! mérés beolvasása
open(1,file='14111806.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

rc=0.75
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='14111806_deconvoluted.txt')
open(2,file='14111806_dd.txt')
do y=0, 100, 5
	e0=sim(0,y)
	do x=0, 100, 5
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x/5, y/5, deconvoluted(x,y)
		write(2, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
	end do
end do
close(1)
close(2)

! mérés beolvasása
open(1,file='14111807.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

rc=0.65
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='14111807_deconvoluted.txt')
open(2,file='14111807_dd.txt')
do y=0, 100, 5
	e0=sim(0,y)
	do x=0, 100, 5
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x/5, y/5, deconvoluted(x,y)
		write(2, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
	end do
end do
close(1)
close(2)

! mérés beolvasása
open(1,file='14111808.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

rc=0.3
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='14111808_deconvoluted.txt')
open(2,file='14111808_dd.txt')
do y=0, 100, 5
	e0=sim(0,y)
	do x=0, 100, 5
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x/5, y/5, deconvoluted(x,y)
		write(2, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
	end do
end do
close(1)
close(2)

rc=0.3

! mérés beolvasása
open(1,file='14070706.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv
	eq(i*5, j*5)=conv 
end do
close(1)


! MEANDER deconvolution
error_a=0
error_b=0
e0=sim(0,0)
direction=-1
open(1,file='14070706_deconvoluted.txt')
open(2,file='14070706_diffo.txt')
open(3,file='14070706_diff.txt')
open(5,file='14070706_dd.txt')
do y=0, 100, 5
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 100, 5
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	else
		do x=100, 0, -5
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	endif
end do
close(1)
close(2)
close(3)
close(5)
write(4, *) "5s", error_b/441, error_a/441

rc=0.5

! mérés beolvasása
open(1,file='14070707.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

! MEANDER deconvolution
error_a=0
error_b=0
e0=sim(0,0)
direction=-1
open(1,file='14070707_deconvoluted.txt')
open(2,file='14070707_diffo.txt')
open(3,file='14070707_diff.txt')
open(5,file='14070707_dd.txt')
do y=0, 100, 5
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 100, 5
			if (e0<=sim(x,y)) then
				rc=0.5
			else
				rc=0.5
			endif
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	else
		do x=100, 0, -5
			if (e0<=sim(x,y)) then
				rc=0.5
			else
				rc=0.5
			endif
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	endif
end do
close(1)
close(2)
close(3)
close(5)
write(4, *) "2s", error_b/441, error_a/441

rc=0.65

! mérés beolvasása
open(1,file='14070708.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

! MEANDER deconvolution
error_a=0
error_b=0
e0=sim(0,0)
direction=-1
open(1,file='14070708_deconvoluted.txt')
open(2,file='14070708_diffo.txt')
open(3,file='14070708_diff.txt')
open(5,file='14070708_dd.txt')
do y=0, 100, 5
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 100, 5
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	else
		do x=100, 0, -5
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	endif
end do
close(1)
close(2)
close(3)
close(5)
write(4, *) "1s", error_b/441, error_a/441


rc=0.8

! mérés beolvasása
open(1,file='14070709.txt')
do x=1, 441
	read(1, *) i, j, conv
	!write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

! MEANDER deconvolution
error_a=0
error_b=0
e0=sim(0,0)
direction=-1
open(1,file='14070709_deconvoluted.txt')
open(2,file='14070709_diffo.txt')
open(3,file='14070709_diff.txt')
open(5,file='14070709_dd.txt')
do y=0, 100, 5
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 100, 5
			if (e0<=sim(x,y)) then
				rc=0.8
			else
				rc=0.8
			endif
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	else
		do x=100, 0, -5
			if (e0<=sim(x,y)) then
				rc=0.8
			else
				rc=0.8
			endif
			deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
			e0=sim(x,y)
			error_a=error_a+(deconvoluted(x,y)-eq(x,y))**2
			error_b=error_b+(sim(x,y)-eq(x,y))**2
			write(1, *) x/5, y/5, deconvoluted(x,y)
			write(2, *) x/5, y/5, deconvoluted(x,y)-eq(x,y)
			write(3, *) x/5, y/5, sim(x,y)-eq(x,y)
			write(5, *) x/5, y/5, -deconvoluted(x,y)+sim(x,y)
		end do
	endif
end do
close(1)
close(2)
close(3)
close(5)
write(4, *) "0.5s", error_b/441, error_a/441
close(4)

end program diffusion
