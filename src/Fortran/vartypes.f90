module vartypes
!
! ==================================================================
!
! This module defines the data type used to represent floating point
! double precision numbers in the package.
!
! ==================================================================


    implicit none

    integer, parameter :: dp = kind(1.d0) ! Double precision
    integer, parameter :: qp = kind(1.q0) ! Quadruple precision

end module vartypes
