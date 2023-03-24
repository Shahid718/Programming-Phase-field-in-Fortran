!-------------------------------------------------------------------------------
!      
!   Finite Difference Phase Field Code of Allen-Cahn Equation.
!
!   The program  uses Dislin plotting library
!
!   Author :
!               Shahid Maqbool
!
!   Modified   :
!                    20 December 2022
!
!   To compile and run :
!                            check ReadMe file
!
!-------------------------------------------------------------------------------

program fd_ac_test
  use Dislin
  implicit none

  !--- simulation cell parameters

  integer ( kind = 4 ), parameter :: Nx = 128
  integer ( kind = 4 ), parameter :: Ny = 128
  integer ( kind = 4 ), parameter :: dx = 2
  integer ( kind = 4 ), parameter :: dy = 2

  !--- time integration parameters

  integer ( kind = 4 ), parameter :: nsteps = 1500
  integer ( kind = 4)             :: nprint = 100
  integer (kind = 4 )             :: tsteps 
  real ( kind = 8 )   , parameter :: dt = 0.01
  real ( kind = 8 )               :: start, finish

  !--- material specific parameters

  real ( kind = 8 )   , parameter :: phi_0 = 0.5
  real ( kind = 8 )   , parameter :: mobility = 1.0
  real ( kind = 8 )   , parameter :: grad_coef = 1.0

  !--- microstructure parameters

  real ( kind = 8 )   , parameter :: noise = 0.02
  real ( kind = 8 )   , parameter :: A  = 1.0
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: r, phi, dfdphi
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: lap_phi, dummy_phi
  integer ( kind = 4 )            :: i, j, jp, jm, ip, im

  call cpu_time ( start )

  !--- initial microstructure

  do i = 1 , Nx
     do j = 1, Ny
 
        call random_number ( r (i,j) )
        phi(i,j) = phi_0 + noise*( 0.5 - r(i,j) )
        
     end do
  end do

  !--- start microstructure evolution

  time_loop: do tsteps = 1, nsteps

     row: do i = 1, Nx
        column: do j = 1, Ny

           !--- free energy derivative

           dfdphi(i,j) = A*( 2.0*phi(i,j)*( 1.0 - phi(i,j) )**2 &
                *( 1.0 - 2*phi(i,j) ) )

           !--- laplace evaluation

           jp = j + 1
           jm = j - 1

           ip = i + 1
           im = i - 1

           if ( im == 0 ) im = Nx
           if ( ip == ( Nx + 1 ) ) ip = 1
           if ( jm == 0 ) jm = Ny
           if ( jp == ( Ny + 1 ) ) jp = 1

           lap_phi(i,j) = ( phi(ip,j) + phi(im,j) + phi(i,jm) + &
                phi(i,jp) - 4.0*phi(i,j)) /( dx*dy )              

           dummy_phi(i,j) = dfdphi(i,j) - grad_coef*lap_phi(i,j)

           !--- time integration

           phi(i,j) = phi(i,j) - dt*mobility*dummy_phi(i,j)

           !--- for small deviations

           if ( phi(i,j) >= 0.99999 ) phi(i,j) = 0.99999
           if ( phi(i,j) < 0.00001 )  phi(i,j) = 0.00001

        end do column
     end do row

     !--- print steps

     if ( mod( tsteps, nprint ) .eq. 0 )  print *, 'Done steps  =  ', tsteps

     !--- end microstructure evolution

  end do time_loop

  call cpu_time ( finish )

  !--- write phi on the file and closes it

  open ( 1, file = "ac.dat" )
  
  do i = 1, Nx
     write( 1, * ) ( phi(i,j),j = 1, Ny )
  end do

  close( 1 )

  !--- prints computed time on the screen

  print*,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start

  !--- dislin plot

  call qplclr ( phi, Nx, Ny )

end program fd_ac_test
