program diffusion
implicit none

integer :: i,j,x,y
real, dimension(0:40,0:30) :: sim, deconvoluted
real, dimension(0:40,0:20) :: sim_mg, deconvoluted_mg
real rc, conv, e0

! mérés beolvasása
open(1,file='solid_coupled.txt')
do x=1, 1271
	read(1, *) i, j, conv
	sim(i/25, j/50)=conv 
end do
close(1)

rc=0.5
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='solid_coupled_deconvoluted.txt')
do y=0, 30
	e0=sim(0,y)
	do x=0, 40
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x*25, y*50, deconvoluted(x,y)
	end do
end do
close(1)
close(2)

! mérés beolvasása
open(1,file='solid_uncoupled.txt')
do x=1, 1271
	read(1, *) i, j, conv
	sim(i/25, j/50)=conv 
end do
close(1)

rc=0.8
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='solid_uncoupled_deconvoluted.txt')
do y=0, 30
	e0=sim(0,y)
	do x=0, 40
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x*25, y*50, deconvoluted(x,y)
	end do
end do
close(1)
close(2)


! LIQUID CONTACT
! mérés beolvasása
open(1,file='liquid_coupled.txt')
do x=1, 1271
	read(1, *) i, j, conv
	sim(i/25, j/50)=conv 
end do
close(1)

rc=0.9
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='liquid_coupled_deconvoluted.txt')
do y=0, 30
	e0=sim(0,y)
	do x=0, 40
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x*25, y*50, deconvoluted(x,y)
	end do
end do
close(1)
close(2)

! mérés beolvasása
open(1,file='liquid_uncoupled.txt')
do x=1, 1271
	read(1, *) i, j, conv
	sim(i/25, j/50)=conv 
end do
close(1)

rc=0.9
! FAST COMB deconvolution
e0=sim(0,0)
open(1,file='liquid_uncoupled_deconvoluted.txt')
do y=0, 30
	e0=sim(0,y)
	do x=0, 40
		deconvoluted(x,y) = (sim(x,y) - e0*rc)/(1-rc)
		e0=sim(x,y)
		write(1, *) x*25, y*50, deconvoluted(x,y)
	end do
end do
close(1)
close(2)

! Mg
! mérés beolvasása
open(1,file='liquid_Mg.txt')
do x=1, 861
	read(1, *) i, j, conv
	sim_mg(i/25, j/50)=conv 
end do
close(1)

rc=0.8
! FAST COMB deconvolution
e0=sim_mg(0,0)
open(1,file='liquid_Mg_deconvoluted.txt')
do y=0, 20
	e0=sim_mg(0,y)
	do x=0, 40
		deconvoluted_mg(x,y) = (sim_mg(x,y) - e0*rc)/(1-rc)
		e0=sim_mg(x,y)
		write(1, *) x*25, y*50, deconvoluted_mg(x,y)
	end do
end do
close(1)
close(2)

! mérés beolvasása
open(1,file='solid_Mg.txt')
do x=1, 861
	read(1, *) i, j, conv
	sim_mg(i/25, j/50)=conv 
end do
close(1)

rc=0.8
! FAST COMB deconvolution
e0=sim_mg(0,0)
open(1,file='solid_Mg_deconvoluted.txt')
do y=0, 20
	e0=sim_mg(0,y)
	do x=0, 40
		deconvoluted_mg(x,y) = (sim_mg(x,y) - e0*rc)/(1-rc)
		e0=sim_mg(x,y)
		write(1, *) x*25, y*50, deconvoluted_mg(x,y)
	end do
end do
close(1)
close(2)



end program diffusion
