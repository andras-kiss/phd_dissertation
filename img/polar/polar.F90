program diffusion

! input file, mask, D, h, time, output lower left, output top right, output height, output file

implicit none

integer :: x_size, y_size, z_size
integer :: res, height
integer :: time_res

! 0.000703 mm^2/s = 703um^2/s for resolution of 1umx1um and 1s
! = 0.00000703 cm^2/s for resolution of 1umx1um and 1s
! CRC Handbook of chemistry and physics 87th 2006-2007.pdf page 848 (acrobat reader p#)
! 1/2 Zn2+ -ra vonatkozik.
real :: const=0.017575 !for resolution of 20umx20um and 0.1s
! real :: const=0.02812 !for resolution of 50umx50um and 0.1s

integer :: h,i,j,k,x,y,m,switch,cells, intercounter
real, dimension(0:101,0:101,0:101) :: a, b
real, dimension(1:100,1:100,1:100) :: flux
integer, dimension(1:100,1:100,1:100) :: mask
real :: pi=3.1415926535897932384626433832795
real maximum
real e0, x_real, y_real
real, dimension(0:100,0:100) :: sim, deconvoluted
integer direction, divisions, alpha
real circumference, r_real
real conv
integer r, n
real rc
real meander_delta, comb_delta, fastcomb_delta, web_delta, arc_delta

!interpoláláshoz:
!integer :: nd =352, ni=441
integer nd
integer ni
real xd(1:352)
real xi(1:441)
real yd(1:352)
real yi(1:441)
real zd(1:352)
real zi(1:441)
real p
nd = 352
ni = 441
p=10.0000


rc=0.8
direction=-1
maximum=0.00000000000000
height=50


a=0.
b=0.
mask=1

! Tömb feltöltése. Ha azt akarod, hogy egy cella ne változzon (Dirichlet), akkor azon a cellán mask=0.
! Fluxus tömb: minden körben hozzáadódik a fluxusban tárolt érték a c-tömbhöz.
! Térbeli felbontást a D (const) határozza meg.
do i=1, 100
	do j=1, 100
		if ( ((i-50)**2+(j-50)**2) < 17**2 ) then !20 itt a source sugara = 400um
			flux(i,j,1)=0.1
		endif
	end do
end do

!open(1,file='flux at 1um.txt')
!do i=0, 100
!	do j=0, 100
!		write(1, *) i*20, j*20, flux(i,j,1)  
!	end do
!end do
close(1)

!open(1,file='flux at 1um_final res.txt')
!do i=0, 100, 5
!	do j=0, 100, 5
!		write(1, *) i*20, j*20, flux(i,j,1)  
!	end do
!end do
!close(1)

! Dirichlet feltétel
! front és hátsó lap
do k=1, 100
	do i=1, 100
		mask(i, 1, k)=0
		mask(i, 100, k)=0
	end do
end do
! jobb és bal lap
do k=1, 100
	do j=1, 100
		mask(1, j, k)=0
		mask(100, j, k)=0
	end do
end do
! fedőlap
do j=1, 100
	do i=1, 100
		mask(i, j, 100)=0
	end do
end do

! Számolás. Az időbeli felbontást D (const) határozza meg.
! b és a tömbök váltakoznak, így nem kell az egész tömböt másolgatni.
b=a
switch=0
! x=i, y=j, z=k, h=time
do h=1, 200 ! M A I N  L O O P

	! All real cells computed. Not cycled: borders, which are all zeros.
	do k=1, 100
		do j=1, 100
			do i=1, 100
				if (mask(i,j,k)==1) then
					cells=6
					if ((k==1) .or. (k==100)) then 
						cells=cells-1
					endif
					if ((j==1) .or. (j==100)) then 
						cells=cells-1
					endif
					if ((i==1) .or. (i==100)) then 
						cells=cells-1
					endif
					b(i,j,k)=a(i,j,k)+const*(a(i,j+1,k)+a(i-1,j,k)+a(i+1,j,k) &
					+a(i,j-1,k)+a(i,j,k-1)+a(i,j,k+1)-cells*a(i,j,k))
				endif
			end do
		end do
	end do
	
	! Fluxusok hozzáadása.

	do j=1, 100
		do i=1, 100
			b(i,j,1)=b(i,j,1)+flux(i,j,1)  
		end do
	end do	
	
	! a és b tömbök cseréje
	a=b
	
	print *,h
end do
	
do i=0, 100
	do j=0, 100
		if (a(i,j,height)>maximum) then
			maximum=a(i,j,height)
		endif
	end do
end do	
print *,maximum
do i=0, 100
	do j=0, 100
		a(i,j,height)=a(i,j,height)/maximum
	end do
end do	
					
	
! O U T P U T	
!open(1,file='real_100um_fullres.txt')
!do i=0, 100
!	do j=0, 100
!		write(1, *) i*20, j*20, a(i,j,height)
!	end do
!end do		
!close(1)

open(1,file='real_sim.txt')
do i=0, 100, 5
	do j=0, 100, 5
		write(1, *) i*20, j*20, a(i,j,height)
	end do
end do		
close(1)


! SECM scanning simulation
! FAST COMB
sim=0
e0=a(0,0,height)
fastcomb_delta=0
open(1,file='fastcomb_sim.txt')
open(2,file='fastcomb_delta.txt')
do y=0, 100, 5
	do x=0, 100, 5
		sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
		e0=sim(x,y)
		write(1, *) x*20, y*20, sim(x,y)
		write(2, *) x*20, y*20, sim(x,y)-a(x,y,height)
		fastcomb_delta=fastcomb_delta+(sim(x,y)-a(x,y,height))**2
	end do
	do x=100, 0, -5
 		e0 = a(x,y,height) + (e0-a(x,y,height))*0.99
	end do
end do
close(1)
close(2)

! MEANDER
sim=0
e0=a(0,0,height)
open(1,file='meander_sim.txt')
open(2,file='meander_delta.txt')
open(3,file='meander_real.txt')
meander_delta=0
do y=0, 100, 5
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 100, 5
			sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
			e0=sim(x,y)
			write(1, *) x*20, y*20, sim(x,y)
			write(2, *) x*20, y*20, sim(x,y)-a(x,y,height)
			write(3, *) x*20, y*20, a(x, y, height)
			meander_delta=meander_delta+(sim(x,y)-a(x,y,height))**2
		end do
	else
		do x=100, 0, -5
			sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
			e0=sim(x,y)
			write(1, *) x*20, y*20, sim(x,y)
			write(2, *) x*20, y*20, sim(x,y)-a(x,y,height)
			write(3, *) x*20, y*20, a(x, y, height)
			meander_delta=meander_delta+(sim(x,y)-a(x,y,height))**2
		end do
	endif
end do
close(1)
close(2)
close(3)

open(1,file='meander.txt')
do x=1, 441
	read(1, *) i, j, conv
	write(2, *) i, j, conv
	sim(i*5, j*5)=conv 
end do
close(1)

! MEANDER deconvolution
e0=sim(0,0)
direction=-1
open(1,file='meander_deconvoluted.txt')
do y=0, 100, 5
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 100, 5
			deconvoluted(x,y) = (sim(x,y) - e0*0.75)/0.25
			e0=sim(x,y)
			write(1, *) x*20, y*20, deconvoluted(x,y)
		end do
	else
		do x=100, 0, -5
			deconvoluted(x,y) = (sim(x,y) - e0*0.75)/0.25
			e0=sim(x,y)
			write(1, *) x*20, y*20, deconvoluted(x,y)
		end do
	endif
end do
close(1)

! MEANDER sim deconvolution
e0=sim(0,0)
direction=-1
open(1,file='meander_sim_deconvoluted.txt')
do y=0, 100, 5
	direction=direction*(-1)
	if (direction==1) then
		do x=0, 100, 5
			deconvoluted(x,y) = (sim(x,y) - e0*0.8)/0.2
			e0=sim(x,y)
			write(1, *) x*20, y*20, deconvoluted(x,y)
		end do
	else
		do x=100, 0, -5
			deconvoluted(x,y) = (sim(x,y) - e0*0.8)/0.2
			e0=sim(x,y)
			write(1, *) x*20, y*20, deconvoluted(x,y)
		end do
	endif
end do
close(1)


! COMB
sim=0
e0=a(0,0,height)
open(1,file='comb_sim.txt')
open(2,file='comb_pattern.txt')
open(3,file='comb_delta.txt')
do y=0, 100, 5
	do x=0, 100, 5
		sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
		e0=sim(x,y)
		!write(1, *) x*20, y*20, sim(x,y)
		write(2, *) x*20, y*20, sim(x,y)
	end do
	do x=100, 0, -5
		sim(x,y) = (sim(x,y) + a(x,y,height) + (e0-a(x,y,height))*rc)/2
		e0=a(x,y,height) + (e0-a(x,y,height))*rc
		write(1, *) x*20, y*20, sim(x,y)
		write(2, *) x*20, y*20, sim(x,y)
		write(3, *) x*20, y*20, sim(x,y)-a(x,y,height)
		comb_delta=comb_delta+(sim(x,y)-a(x,y,height))**2
	end do
end do
close(1)
close(2)

! SPIDER-NET integer loop
x_real=0
y_real=0
e0=a(50,50,height)
web_delta=0
open(1,file='web_sim.txt')
open(2,file='web_delta.txt')
do r=0, 10
	do alpha=0, 9
		x_real=100*r*cos(2*pi*alpha/10)
		y_real=100*r*sin(2*pi*alpha/10)
		x=nint((x_real+1000)/20)
		y=nint((y_real+1000)/20)
		!sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
		e0=a(x,y,height)+(e0-a(x,y,height))*rc
		write(1, *) nint(x_real), nint(y_real), e0
		write(2, *) nint(x_real), nint(y_real), e0-a(x,y,height)
		web_delta=web_delta+(e0-a(x,y,height))**2
		print *,x_real,y_real
	end do
end do
close(1)

! SPIDER-NET
!x_real=0
!y_real=0
!e0=a(50,50,height)
!open(1,file='spider_net.txt')
!do r=0, 10
!	do alpha=0, 2*pi-pi/10, 2*pi/10
!		x_real=100*r*cos(alpha)
!		y_real=100*r*sin(alpha)
!		x=nint((x_real+1000)/20)
!		y=nint((y_real+1000)/20)
!		!sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
!		e0=a(x,y,height)+(e0-a(x,y,height))*rc
!		write(1, *) nint(x_real), nint(y_real), e0
!		print *,x_real,y_real
!	end do
!end do
!close(1)

! STONEHENGE
!x_real=0
!y_real=0
!e0=a(50,50,height)
!open(1,file='stonehenge.txt')
!write(1, *) 0, 0, e0
!do r=0, 10
!	circumference=2*r*100*pi
!	divisions=circumference/100
!	do alpha=0, 2*pi-2*pi/divisions, 2*pi/divisions
!		x_real=100*r*cos(alpha)
!		y_real=100*r*sin(alpha)
!		x=nint((x_real+1000)/20)
!		y=nint((y_real+1000)/20)
!		!sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
!		e0=a(x,y,height)+(e0-a(x,y,height))*rc
!		write(1, *) nint(x_real), nint(y_real), e0
!		print *,x_real,y_real
!	end do
!end do
close(1)

! STONEHENGE integer

x_real=0
y_real=0
e0=a(50,50,height)
intercounter=1
arc_delta=0
open(1,file='arc_sim.txt')
open(2,file='arc_delta.txt')
write(1, *) 0, 0, e0
do r=0, 10
	circumference=2*r*100*pi
	divisions=circumference/100
	do alpha=0, divisions-1
		x_real=100*r*cos(2*pi*alpha/divisions)
		y_real=100*r*sin(2*pi*alpha/divisions)
		x=nint((x_real+1000)/20)
		y=nint((y_real+1000)/20)
		!sim(x,y) = a(x,y,height) + (e0-a(x,y,height))*rc
		e0=a(x,y,height)+(e0-a(x,y,height))*rc
		xd(intercounter)=nint(x_real)
		yd(intercounter)=nint(y_real)
		zd(intercounter)=e0
		intercounter=intercounter+1
		write(1, *) nint(x_real), nint(y_real), e0
		write(2, *) nint(x_real), nint(y_real), e0-a(x,y,height)
		arc_delta=arc_delta+(e0-a(x,y,height))**2
		print *,x_real,y_real
	end do
end do
close(1)

intercounter=1
do x=-1000, 1000, 100
	do y=-1000, 1000, 100
		xi(intercounter)=x
		yi(intercounter)=y
		intercounter=intercounter+1
	end do
end do	

call shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi, zi)


open(1,file='arc_sim_interpolated.txt')
do x=1, ni
 write(1, *) xi(x), yi(x), zi(x)
end do
close(1)

open(1,file='deltas.txt')
write(1, *) meander_delta/441, comb_delta/441, fastcomb_delta/441, web_delta/110, arc_delta/341
close(1)

end program diffusion

subroutine shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi, zi )

!*****************************************************************************80
!
!! SHEPARD_INTERP_2D evaluates a 2D Shepard interpolant.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), YD(ND), the data points.
!
!    Input, real ( kind = 8 ) ZD(ND), the data values.
!
!    Input, real ( kind = 8 ) P, the power.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), YI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) ZI(NI), the interpolated values.
!
  implicit none

  integer nd
  integer ni

  integer i
  integer j
  real p
  real s
  real w(nd)
  real xd(nd)
  real xi(ni)
  real yd(nd)
  real yi(ni)
  integer z
  real zd(nd)
  real zi(ni)

  do i = 1, ni

    z = -1
    do j = 1, nd
      w(j) = sqrt ( ( xi(i) - xd(j) ) ** 2 + ( yi(i) - yd(j) ) ** 2 )
      if ( w(j) == 0.0) then
        z = j
        exit
      end if
    end do
    if ( z /= -1 ) then
      w(1:nd) = 0.0
      w(z) = 1.0
    else
      w(1:nd) = 1.0 / w(1:nd) ** p
      s = sum ( w )
      w(1:nd) = w(1:nd) / s
    end if

    zi(i) = dot_product ( w, zd )

  end do

  return
end
