!-------------------------------------------------------------------------------
!      
!   Fortran Finite Difference Phase Field Code for 2 Grains.
!
!   Author :
!               Shahid Maqbool
!
!   Modified   :
!                    23 December 2022
!
!   To compile and run :
!                            check ReadMe file
!
!-------------------------------------------------------------------------------

program fd_grain_2_test
  implicit none

  !-- simulation cell parameters

  integer( kind = 4 ), parameter :: Nx = 64
  integer( kind = 4 ), parameter :: Ny = 64
  integer( kind = 4 ), parameter :: NxNy = Nx*Ny
  real   ( kind = 8 ), parameter :: dx = 0.5
  real   ( kind = 8 ), parameter :: dy = 0.5

  !--- time integration parameters

  integer( kind = 4 ), parameter :: nstep  = 5000
  integer( kind = 4 ), parameter :: nprint = 500
  real   ( kind = 8 ), parameter :: dt  = 0.005
  integer ( kind = 4 )           :: istep

  !--- material parameters

  real   ( kind = 8 ), parameter :: mobility  = 5.0
  real   ( kind = 8 ), parameter :: grad_coef = 0.1

  !--- initial grain structure parameters

  integer( kind = 4 ), parameter :: ngrain = 2
  real   ( kind = 8 ), parameter :: radius = 14.0
  real   ( kind = 8 ), dimension( Nx,Ny,ngrain ) :: etas
  real   ( kind = 8 ) :: xlength, x0, y0 

  !--- evolution and free energy parameters

  real   ( kind = 8 ), dimension( ngrain ) :: glist
  real   ( kind = 8 ), dimension( Nx,Ny )  :: eta, dfdeta, lap_eta
  real   ( kind = 8 ), parameter :: A = 1.0, B = 1.0
  integer( kind = 4 ) :: i, j, jp, jm, ip, im, igrain, jgrain
  real   ( kind = 8 ) :: summ, grain_sum

  !--- inital microstructure

  x0 = Nx/2
  y0 = Ny/2

  do i = 1, Nx
     do j = 1, Ny

        etas(i,j,1) = 1.0
        etas(i,j,2) = 0.0

        xlength = sqrt( ( i - x0 )**2 + ( j - y0 )**2 )

        if ( xlength <= radius ) then
           etas(i,j,1) = 0.0
           etas(i,j,2) = 1.0
        end if

     end do
  end do

  !--- initialize glist

  do igrain = 1, ngrain
     glist(igrain) = 1.0
  end do

  !--- starts microstructure evolution

  time_loop: do istep = 1, nstep

     do igrain = 1, ngrain
        if ( glist(igrain) == 1 ) then
           do i = 1, Nx
              do j = 1, Ny
                 eta(i,j) = etas(i,j, igrain)
              end do
           end do

           do i = 1, Nx
              do j = 1, Ny

                 jp = j + 1
                 jm = j - 1

                 ip = i + 1
                 im = i - 1

                 if ( im == 0 ) im = Nx
                 if ( ip == ( Nx +1 ) ) ip = 1
                 if ( jm == 0) jm = Ny
                 if ( jp == ( Ny + 1) ) jp = 1

                 !--- laplace evaluation

                 lap_eta(i,j) = ( eta(ip,j) + eta(im,j) + eta(i,jm) + &
                      eta(i,jp) - 4.0*eta(i,j) ) / ( dx*dy )

                 !--- derivative of free energy

                 summ = 0.0

                 do jgrain = 1,ngrain
                    if ( jgrain /= igrain ) then
                       summ = summ + etas(i,j,jgrain)**2
                    end if
                 end do

                 dfdeta(i,j) = A*( 2.0*B* eta(i,j)*summ  + eta(i,j)**3 &
                      - eta(i,j) )

                 !--- time integration

                 eta(i,j) = eta(i,j) - dt*mobility*( dfdeta(i,j) &
                      - grad_coef*lap_eta(i,j) )

                 !-- for small deviations

                 if ( eta(i,j) >= 0.9999 ) eta(i,j) = 0.9999
                 if ( eta(i,j) < 0.00001 ) eta(i,j) = 0.00001

              end do
           end do

           grain_sum = 0.0

           do i = 1, Nx
              do j = 1, Ny
                 etas(i,j,igrain) = eta(i,j)
                 grain_sum = grain_sum + eta(i,j)
              end do
           end do

           !--- check volume fraction of current grain

           grain_sum = grain_sum / NxNy

           if ( grain_sum <= 0.001 ) then
              glist(igrain) = 0
           end if

        end if
     end do

     !--- print done steps and write the values on files

     writing_files:block
       character ( len = 80 ) :: filename 
       integer ( kind = 4 )   :: lun

       if ( mod ( istep, nprint ).eq.0 ) then
          print*, 'done step =', istep
          write ( filename, '( "grain_", i0.3, ".dat" )') istep
          open ( newunit = lun, file = filename, status = 'replace' )
          do i = 1, Nx    
             write ( lun, * ) ( eta( i, j), j = 1, Ny )
          end do
          close ( unit = lun )
       endif

     end block writing_files

  end do time_loop

end program fd_grain_2_test
