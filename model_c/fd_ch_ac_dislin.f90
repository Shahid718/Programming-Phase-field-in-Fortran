!                                                                             
!   FINITE DIFFERENCE PHASE-FIELD CODE FOR SOLVING MODEL C
!
!   The program uses Dislin graphical library
!                                                                             
!
!   Author  :
!               Shahid Maqbool
!
!   Modified   :
!                    25 May 2023
!
!   To compile and run :
!                            check ReadMe file
!
!------------------------------------------------------------------------------


program fd_ch_ac_test
  use Dislin
  implicit none


  !--- simulation cell parameters

  integer ( kind = 4 ), parameter :: Nx = 128
  integer ( kind = 4 ), parameter :: Ny = 128
  integer ( kind = 4 ), parameter :: dx = 1
  integer ( kind = 4 ), parameter :: dy = 1


  !--- time integration parameters

  integer ( kind = 4 ), parameter :: nsteps = 10000
  integer ( kind = 4)             :: nprint = 1000
  integer (kind = 4 )             :: istep 
  real ( kind = 8 )   , parameter :: dt = 0.03
  real ( kind = 8 )               :: start, finish


  !--- material specific parameters

  real ( kind = 8 )   , parameter :: A = 1.0
  real ( kind = 8 )   , parameter :: B = 1.0
  real ( kind = 8 )   , parameter :: D = 1.0
  real ( kind = 8 )   , parameter :: mobility_con = 0.5
  real ( kind = 8 )   , parameter :: mobility_phi = 0.5
  real ( kind = 8 )   , parameter :: grad_coef_con = 1.5
  real ( kind = 8 )   , parameter :: grad_coef_phi = 1.5
  real ( kind = 8 )   , parameter :: radius = 10.0


  !--- microstructure parameters

  integer ( kind = 4 )                       :: i, j, jp, jm, ip, im
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: con, phi, dfdcon, dfdphi
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: lap_con, lap_phi
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: dummy_con, lap_dummy, phi_dummy
  real ( kind = 8 )   , dimension ( Nx )     :: con_profile, phi_profile, x


  call cpu_time ( start )


  !--- initial microstructure


  do i = 1, Nx
     do j = 1, Ny


        if ( (i - Nx/2)*(i - Nx/2) + (j - Ny/2)*(j - Ny/2) < radius**2 ) then
           con(i,j) = 1.0
           phi(i,j) = 1.0
        else
           con(i,j) = 0.02
           phi(i,j) = 0.0
        endif


     end do
  end do


  !---- setting initial dislin routines for multiplot 

  call metafl ( 'cons' )
  call scrmod ( 'revers' )
  call disini ( )


  !--- start microstructure evolution


  time_loop: do istep = 1, nsteps


     do i = 1, Nx
        do j = 1, Ny


           !--- derivative wrt concentration and phi

           dfdcon(i,j) = 2*A*con(i,j)*(1-( phi(i,j)**3*( 10 - 15*phi(i,j) + &
                6*phi(i,j)**2 ) )) - 2*B*(1 - con(i,j))* &
                ( phi(i,j)**3*( 10 - 15*phi(i,j) + 6*phi(i,j)**2 ) )

           dfdphi(i,j) = -A*con(i,j)*con(i,j)*( 3*phi(i,j)**2*( 10 - &
                15*phi(i,j) + 6*phi(i,j)**2 ) + phi(i,j)**3* &
                ( 12*phi(i,j) - 15 )) + 2*B*(1 - con(i,j))* &
                (1 - con(i,j))*( 3*phi(i,j)**2*( 10 - 15*phi(i,j) + &
                6*phi(i,j)**2 ) + phi(i,j)**3*( 12*phi(i,j) - 15 )) + &
                2*D*phi(i,j)*(1 - phi(i,j))*(1 - 2*phi(i,j) )


           !--- laplace evaluation

           jp = j + 1
           jm = j - 1

           ip = i + 1
           im = i - 1

           if ( im == 0 ) im = Nx
           if ( ip == ( Nx + 1) ) ip = 1
           if ( jm == 0 ) jm = Ny
           if ( jp == ( Ny + 1) ) jp = 1


           !--- concentration

           lap_con(i,j)   = ( con(ip,j) + con(im,j) + con(i,jm) + &
                con(i,jp) - 4.0*con(i,j) ) / ( dx*dy )
           dummy_con(i,j) = dfdcon(i,j) - grad_coef_con*lap_con(i,j)


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


           lap_dummy(i,j) = ( dummy_con(ip,j) + dummy_con(im,j) + &
                dummy_con(i,jm) + dummy_con(i,jp) - 4.0*dummy_con(i,j) ) &
                / ( dx*dy )


           !--- phi

           lap_phi(i,j) = ( phi(ip,j) + phi(im,j) + phi(i,jm) + &
                phi(i,jp) - 4.0*phi(i,j) ) / ( dx*dy )
           phi_dummy(i,j) = dfdphi(i,j) - grad_coef_phi*lap_phi(i,j)


           !--- time integration

           con(i,j) = con(i,j) + dt*mobility_con*lap_dummy(i,j)
           phi(i,j) = phi(i,j) - dt*mobility_phi*phi_dummy(i,j)


        end do
     end do


     !--- for small deviations

     do i = 1, Nx
        do j = 1, Ny

           if ( phi(i,j) >= 0.99999 )  phi(i,j) = 0.99999
           if ( phi(i,j) < 0.00001 )   phi(i,j) = 0.00001

        end do
     end do


     !--- print steps

     if ( mod( istep, nprint ) .eq. 0 )  print *, 'Done steps  =  ', istep


     !--- dislin multiplot

     do i = 1, Nx
        x(i)  = i
     end do


     con_profile = con( : , Nx/2 )
     phi_profile = phi( : , Nx/2 )


     if ( istep == 1 ) then

        call complx ( )
        call axspos ( 500, 800 )
        call axslen ( 600, 600 )
        call ticks( 0,'XY' )
        call titlin ( 'concentration', 4 )

        call frame ( 3 )                 
        call height ( 50 )
        call hname ( 50 )                
        call name ( 'Nx', 'x' )
        call name ( 'con', 'y' )

        call labdig ( -1, 'x' )  
        call graf ( 0.d0, 128.d0, 0.d0, 32.d0, 0.0d0, 1.3d0, 0.0d0, 0.2d0 )

        call thkcrv ( 3 )                
        call color ( 'blue' )
        call curve ( x, con_profile, Nx )
        call color ( 'fore' )
        call title
        call endgrf

     end if


     if ( istep == 1 ) then

        call autres ( Nx, Ny )
        call labdig ( 2, 'Z' )
        call labdig ( -1, 'y' ) 
        call ticks( 0,'XY' )
        call titlin ( 'phi', 4 )
        call complx ( )

        call axspos ( 1700, 800 )
        call ax3len ( 600, 600, 600 )
        call name ( 'Ny', 'y' )
        call graf3 ( 0.d0, 128.d0, 0.d0, 32.d0, 0.d0, 128.d0,&
             & 0.d0, 32.d0, 0.0d0, 1.0d0, 0.0d0, 0.2d0 )
        call crvmat ( phi, Nx, Ny, 1, 1 )
        call title ( )
        call endgrf ( )


     end if

     if ( istep == 10000 ) then
        call complx ( )
        call axspos ( 500, 1700 )
        call axslen ( 600, 600 )

        call frame ( 3 )                 
        call height ( 50 )
        call hname ( 50 )                
        call name ( 'Nx', 'x' )
        call name ( 'con', 'y' )
        call labdig ( -1, 'x' )
        call labdig ( 1, 'y' )
        call graf ( 0.d0, 128.d0, 0.d0, 32.d0, 0.0d0, 1.3d0, 0.d0, 0.2d0 )

        call thkcrv ( 3 )                 
        call color ( 'red' )
        call curve ( x, con_profile, Nx )
        call color ( 'fore' )
        call endgrf


     end if


     if ( istep == 10000 ) then

        call autres ( Nx, Ny )
        call labdig ( 2, 'Z' )
        call labdig ( -1, 'y' ) 
        call ticks( 0,'XY' )
        call complx ( )

        call axspos ( 1700, 1700 )
        call ax3len ( 600, 600, 600 )
        call name ( 'Ny', 'y' )
        call graf3 ( 0.d0, 128.d0, 0.d0, 32.d0, 0.d0, 128.d0,&
             & 0.d0, 32.d0, 0.0d0, 1.0d0, 0.0d0, 0.2d0 )
        call crvmat ( phi, Nx, Ny, 1, 1 )
        call endgrf

        call disfin ( )

     end if


     !--- write the values on files

     if ( istep == 1 ) then
        open ( 1, file = 'phi_1.dat' )
        do i = 1, Nx    
           write ( 1, * ) ( phi( i, j), j = 1, Ny )  
        end do
        close ( 1 )
     end if

     if ( istep == 10000 ) then      
        open ( 2, file = 'phi_10000.dat' )
        do i = 1, Nx            
           write ( 2, * ) ( phi( i, j), j = 1, Ny )
        end do
        close ( 2 )
     end if


     !--- end microstructure evolution

  end do time_loop


  call cpu_time ( finish )


  !--- prints computed time on the screen

  print*,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start


end program fd_ch_ac_test
