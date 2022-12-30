!-------------------------------------------------------------------------------
!      
!   Finite Difference Phase Field Code of Fe-Cr phase separation.
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

program fd_feCr_test
  use Dislin
  implicit none

  !--- simulation cell parameters

  integer ( kind = 4 ), parameter :: Nx = 128
  integer ( kind = 4 ), parameter :: Ny = 128
  integer ( kind = 4 ), parameter :: dx = 1
  integer ( kind = 4 ), parameter :: dy = 1
  integer ( kind = 4 )            :: dxdy 

  !--- time integration parameters

  integer ( kind = 4 ), parameter :: nsteps = 30000
  integer ( kind = 4 )            :: nprints = 1000
  integer (kind = 4 )             :: tsteps 
  real ( kind = 8 )   , parameter :: dt = 0.01
  real ( kind = 8 )               :: start, finish

  !--- material specific parameters

  real ( kind = 8 )   , parameter :: c0 = 0.2
  real ( kind = 8 )   , parameter :: mobility = 0.50
  real ( kind = 8 )   , parameter :: grad_coef = 2.0
  real ( kind = 8 )   , parameter :: temperature = 535.0 
  real ( kind = 8 )   , parameter :: gas_constant = 8.314462
  real ( kind = 8 )               :: RT

  !--- microstructure parameters

  real ( kind = 8 )   , parameter :: noise = 0.02
  real ( kind = 8 )   , parameter :: A  = 1.0
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: r, cr, lap_cr
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: dfdcr, dummy_cr, lap_dummy
  integer ( kind = 4 )            :: i, j, jp, jm, ip, im

  open ( 1, file = "feCr.dat" )

  !--- initial microstructure

  call cpu_time ( start )

  call random_number ( r )

  cr = c0 + noise*( 0.5 - r )

  !--- start microstructure evolution

  RT = gas_constant*temperature
  dxdy = dx*dy

  temporal_loop: do tsteps = 1, nsteps

     spatial_loop: do concurrent ( i = 1 : Nx, j = 1 : Ny )

        !--- free energy derivative

        dfdcr(i,j) = ( -cr(i,j)*( 20500.0 - 9.68*temperature ) + &
             & ( 1.0 - cr(i,j) )* ( 20500.0 - 9.68*temperature ) + &
             & ( log(cr(i,j)) - log(1.0 - cr(i,j)) ) *RT ) / RT

        !--- laplace evaluation

        jp = j + 1
        jm = j - 1

        ip = i +1
        im = i-1

        if ( im == 0 ) im = Nx
        if ( ip == ( Nx + 1 ) ) ip = 1
        if ( jm == 0 ) jm = Ny
        if ( jp == ( Ny + 1 ) ) jp = 1

        lap_cr(i,j) = ( cr(ip,j) + cr(im,j) + cr(i,jm) + &
             & cr(i,jp) - 4.0*cr(i,j) ) / ( dxdy )

        dummy_cr(i,j) = dfdcr(i,j) - grad_coef*lap_cr(i,j)

        lap_dummy(i,j) = (dummy_cr(ip,j) + dummy_cr(im,j) &
             & + dummy_cr(i,jm) + dummy_cr(i,jp) &
             & - 4.0*dummy_cr(i,j) ) / ( dxdy )

        !--- time integration

        cr(i,j) =  cr(i,j) + dt*mobility*lap_dummy(i,j)

        !--- for small deviations

        if ( cr(i,j) >= 0.99999 ) cr(i,j) = 0.99999
        if ( cr(i,j) < 0.00001)  cr(i,j) = 0.00001

     end do spatial_loop

     !--- print done steps

     if ( mod( tsteps, nprints ).eq. 0 ) Print*, 'Done steps  =  ', tsteps

     !--- end microstructure evolution

  end do temporal_loop

  call cpu_time ( finish )

  !--- write concentration on the file and closes it

  do i = 1, Nx
     write( 1, * ) ( cr(i,j), j = 1, Ny )
  end do

  close( 1 )

  !--- prints computed time on the screen

  print*,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start

  !--- quick Dislin color plot

  call qplclr ( cr, Nx, Ny )

end program fd_feCr_test
