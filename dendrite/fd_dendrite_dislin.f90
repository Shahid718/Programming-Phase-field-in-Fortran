!-------------------------------------------------------------------------------
!      
!   Finite Difference Phase Field Code of Dendrite Solidification.
!
!   The program uses Dislin plotting library.
!
!   Author  :
!              Shahid Maqbool
! 
!   Modified   :
!                    22 December 2022
!
!   To compile and run :
!                            check ReadMe file
!
!-------------------------------------------------------------------------------

program fd_Kobayashi_model_test
  use Dislin
  implicit none

  !--- simulation cell parameters

  integer ( kind = 4 ), parameter :: Nx = 300
  integer ( kind = 4 ), parameter :: Ny = 300
  real ( kind = 8 )               :: dx = 0.03
  real ( kind = 8 )               :: dy = 0.03

  !--- time integeration parameters

  integer (kind = 4 ) :: nsteps = 2000
  integer (kind = 4 ) :: nprint = 100
  integer (kind = 4 ) :: tsteps 
  real ( kind = 8 )   :: dtime  = 1.0e-4
  real ( kind = 8 )   :: start, finish

  !--- material specific parameters

  real ( kind = 8 )   :: tau   = 0.0003
  real ( kind = 8 )   :: epsilonb = 0.01
  real ( kind = 8 )   :: mu    = 1.0
  real ( kind = 8 )   :: kappa = 1.8
  real ( kind = 8 )   :: delta = 0.02
  real ( kind = 8 )   :: aniso = 6.0
  real ( kind = 8 )   :: alpha = 0.9
  real ( kind = 8 )   :: gama  = 10.0
  real ( kind = 8 )   :: teq   = 1.0
  real ( kind = 8 )   :: theta0= 0.2 
  real ( kind = 8 )   :: seed  = 5.0

  real ( kind = 8 )   :: pix   = 4.0*atan(1.0)

  !--- initial nuclei and evolution parameters

  real ( kind = 8 ) , dimension( Nx, Ny ) :: phi, tempr
  real ( kind = 8 ) , dimension( Nx, Ny ) :: lap_phi, lap_tempr
  real ( kind = 8 ) , dimension( Nx, Ny ) :: phidx, phidy
  real ( kind = 8 ) , dimension( Nx, Ny ) :: epsil, epsilon_deriv
  real ( kind = 8 )                       :: phi_old, term1, term2
  real ( kind = 8 )                       :: theta, m
  integer ( kind = 4 )                    :: i, j, istep, ip, im, jp, jm

  open ( 1, file = 'phi.dat' )
  open ( 2, file = 'temperature.dat')

  call cpu_time ( start )

  !--- initialize and introduce initial nuclei

  phi = 0.0
  tempr = 0.0

  do i = 1, Nx
     do j = 1, Ny
        if ( (i - Nx/2.0)*(i - Nx/2.0) + (j - Ny/2.0)*(j - Ny/2.0)&
             & < seed ) then
           phi(i,j) = 1.0
        end if
     end do
  end do

  !--- start microstructure evolution

  time_loop: do tsteps = 1, nsteps

     do i = 1, Nx
        do j = 1, Ny

           jp = j + 1
           jm = j - 1

           ip = i + 1
           im = i - 1

           if ( im == 0 ) im = Nx
           if ( ip == ( Nx + 1) ) ip = 1
           if ( jm == 0 ) jm = Ny
           if ( jp == ( Ny + 1) ) jp = 1

           !--- laplacian

           lap_phi(i,j) = ( phi(ip,j) + phi(im,j) + phi(i,jm) + phi(i,jp)&
                & - 4.0*phi(i,j)) / ( dx*dy )
           lap_tempr(i,j) = ( tempr(ip,j) + tempr(im,j) + tempr(i,jm) + &
                & tempr(i,jp) - 4.0*tempr(i,j)) / ( dx*dy )

           !--- gradients

           phidx(i,j) = ( phi(ip,j) - phi(im,j) ) / dx
           phidy(i,j) = ( phi(i,jp) - phi(i,jm) ) / dy

           !--- angle

           theta  = atan2( phidy(i,j),phidx(i,j) )

           !--- epsilon and its derivative

           epsil(i,j) = epsilonb*( 1.0 + delta*cos(aniso*&
                & ( theta - theta0 ) ) )
           epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin&
                & ( aniso*( theta - theta0 ) )

        end do
     end do

     do i = 1, Nx
        do j = 1, Ny

           jp = j + 1
           jm = j - 1

           ip = i + 1
           im = i - 1

           if ( im == 0 ) im = Nx
           if ( ip == ( Nx + 1) ) ip = 1
           if ( jm == 0 ) jm = Ny
           if ( jp == ( Ny + 1) ) jp = 1

           phi_old = phi(i,j)

           !--- term1 and term2

           term1 = ( epsil(i,jp)*epsilon_deriv(i,jp)*phidx(i,jp)&
                & - epsil(i,jm)*epsilon_deriv(i,jm)*phidx(i,jm) ) / dy
           term2 = -( epsil(ip,j)*epsilon_deriv(ip,j)*phidy(ip,j)&
                & - epsil(im,j)*epsilon_deriv(im,j)*phidy(im,j) ) / dx

           !--- factor m

           m = alpha/pix*atan( gama*( teq - tempr(i,j) ) )

           !--- time integration

           phi(i,j) = phi(i,j) + ( dtime/tau )*( term1 + term2 +&
                & epsil(i,j)**2*lap_phi(i,j) ) + &
                & phi_old*( 1.0 - phi_old )*( phi_old -0.5 + m )
           tempr(i,j) = tempr(i,j) + dtime*lap_tempr(i,j) &
                & + kappa*( phi(i,j) - phi_old )

        end do
     end do

     !--- print steps

     if ( mod( tsteps, nprint ) .eq. 0 ) print *, 'Done steps  =  ', tsteps

     !--- end microstructure evolution

  end do time_loop

  call cpu_time ( finish )

  !--- write phi and temperature on the files and closes it

  do i = 1, Nx
     write( 1, * ) ( phi(i,j),j = 1, Ny )
     write( 2, * ) ( tempr(i,j),j = 1, Ny )
  end do

  close( 1 )
  close( 2 )

  !--- prints computed time on the screen

  print*,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start

  !--- dislin multiplot

  call metafl ( 'cons' )
  call scrmod ( 'revers' )
  call disini ( )

  call complx ( )

  call axspos ( 400, 1500 )
  call ax3len ( 600, 600, 600 )
  call graf3 ( 0.d0, 300.d0, 0.d0, 100.d0, 0.d0, 300.d0,&
       & 0.d0, 100.d0, 0.0d0, 1.2d0, 0.0d0, 0.2d0 )
  call crvmat ( phi, Nx, Ny, 1, 1 )
  call endgrf

  call axspos ( 1700, 1500 )
  call ax3len ( 600, 600, 600 )
  call graf3 ( 0.d0, 300.d0, 0.d0, 100.d0, 0.d0, 300.d0,&
       & 0.d0, 100.d0, 0.0d0, 1.2d0, 0.0d0, 0.2d0 )
  call crvmat ( tempr, Nx, Ny, 1, 1 )
  call endgrf

  call disfin ( )

end program fd_Kobayashi_model_test
