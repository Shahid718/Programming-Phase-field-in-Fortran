!-------------------------------------------------------------------------------
!      
!   Finite Difference Phase Field Code of Cahn-Hilliard Equation.
!
!   The program  uses Dislin plotting library
!
!   Author  :
!               Shahid Maqbool
!
!   Modified   :
!                    20 December 2022
!
!   To compile and run :
!                            check ReadMe file
!
!-------------------------------------------------------------------------------

program fd_ch_test
  use Dislin
  implicit none

  !--- simulation cell parameters

  integer ( kind = 4 ), parameter :: Nx = 64
  integer ( kind = 4 ), parameter :: Ny = 64
  integer ( kind = 4 ), parameter :: dx = 1
  integer ( kind = 4 ), parameter :: dy = 1

  !--- time integration parameters

  integer ( kind = 4 ), parameter :: nsteps = 10000
  integer ( kind = 4)             :: nprint = 1000
  integer (kind = 4 )             :: tsteps 
  real ( kind = 8 )   , parameter :: dt = 0.01
  real ( kind = 8 )               :: start, finish

  !--- material specific parameters

  real ( kind = 8 )   , parameter :: c0 = 0.4
  real ( kind = 8 )   , parameter :: mobility = 1.0
  real ( kind = 8 )   , parameter :: grad_coef = 0.5

  !--- microstructure parameters

  real ( kind = 8 )   , parameter :: noise = 0.02
  real ( kind = 8 )   , parameter :: A  = 1.0
  integer ( kind = 4 )            :: i, j, jp, jm, ip, im
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: r, con, lap_con, dfdcon
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: dummy_con, lap_dummy

 
  call cpu_time ( start )
  
  !--- initial microstructure
  
   do i = 1 , Nx
     do j = 1, Ny

        call random_number ( r (i,j) )

        con(i,j) = c0 + noise*( 0.5 - r(i,j) )

     end do
  end do

  !--- start microstructure evolution

  time_loop: do tsteps = 1, nsteps

     row:do i = 1, Nx
        column:do j = 1, Ny

           !--- free energy derivative

           dfdcon(i,j) = A*( 2.0*con(i,j)*( 1.0 - con(i,j) )**2 &
                - 2.0*con(i,j)**2*( 1.0 - con(i,j) ) )

           !--- laplace evaluation

           jp = j + 1
           jm = j - 1

           ip = i + 1
           im = i - 1

           if ( im == 0 ) im = Nx
           if ( ip == ( Nx + 1) ) ip = 1
           if ( jm == 0 ) jm = Ny
           if ( jp == ( Ny + 1) ) jp = 1

           lap_con(i,j)   = ( con(ip,j) + con(im,j) + con(i,jm) + con(i,jp) - &
                4.0*con(i,j) ) /( dx*dy )

           dummy_con(i,j) = dfdcon(i,j) - grad_coef*lap_con(i,j)

           lap_dummy(i,j) = ( dummy_con(ip,j) + dummy_con(im,j) + dummy_con(i,jm) &
                + dummy_con(i,jp) - 4.0*dummy_con(i,j) ) / ( dx*dy )

           !--- time integration

           con(i,j) =  con(i,j) + dt*mobility*lap_dummy(i,j)

           !--- for small deviations

           if ( con(i,j) >= 0.99999 )  con(i,j) = 0.99999
           if ( con(i,j) < 0.00001 )   con(i,j) = 0.00001

        end do column
     end do row

     !--- print steps

     if ( mod( tsteps, nprint ) .eq. 0 )  print *, 'Done steps  =  ', tsteps

     !--- end microstructure evolution

  end do time_loop

  call cpu_time ( finish )

  !--- Open, write concentration on the file and closes it
  
  open ( 1, file = "ch.dat" )
   
  do i = 1, Nx
     write( 1, * ) ( con(i,j),j = 1, Ny )
  end do

  close( 1 )

  !--- prints computed time on the screen
  
  print*,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start

  !--- dislin plot

  call qplclr ( con, Nx, Ny )

end program fd_ch_test
