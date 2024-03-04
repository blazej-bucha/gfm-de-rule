function dlog1p(x) result(log1p)
!
! ==================================================================
!
!   Double precision function program
!   to compute a special exponential function,
!       log1p(x) = log(1+x),
!   accurately by minimax rational approximation
!
!   "dlog1p" keeps 16 digit accuracy for arbitrary argument
!   while runs 20 % faster than dlog(x) itself
!
!   Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!   Date: 2019/10/30
!
!
!
!   Minor (formal) modifications by Blazej Bucha.
!
! ==================================================================

    ! Import modules
    ! ===============================================
    use vartypes  ! Module defining floating point double and quadruple
                  ! precision numbers
    use constants ! Module defining constants (pi, etc.)
    ! ===============================================


    ! Declaration of variables
    ! ===============================================
    real(dp), intent(in) :: x
    real(dp) :: x1, x2, x4
    real(dp) :: a1, a2, a3, a4, a5, a6, a7
    real(dp) :: b1, b2, b3, b4, b5, b6
    real(dp) :: log1p
    ! ===============================================


    ! Main computational part of the routine
    ! ===============================================
    parameter (a1 = 0.99999999999999999405_dp, a2 = 2.45235728562912886048_dp)
    parameter (a3 = 2.17053627298972253249_dp, a4 = 0.83928994566440838378_dp)
    parameter (a5 = 0.13520496594993836479_dp, a6 = 0.00682631751459270270_dp)
    parameter (a7 = 0.00002291289324181940_dp)
    parameter (b1 = 2.95235728562912599232_dp, b2 = 3.31338158247117791600_dp)
    parameter (b3 = 1.76186164168333482938_dp, b4 = 0.44976458082070468584_dp)
    parameter (b5 = 0.04896199808811261680_dp, b6 = 0.00157389087429218809_dp)

    if(x .lt. -0.5_dp .or. x .gt. 1.0_dp) then
        log1p = log(1.0_dp + x)
    elseif(x .ge. 0.0_dp) then
        x1 = x
        x2 = x1 * x1
        x4 = x2 * x2

        log1p = x1 * (((a1 + x1 * a2) + x2 * (a3 + x1 * a4)) &
                 + x4 * ((a5 + x1 * a6) + x2 * a7)) &
                 / (((1.0_dp + x1 * b1) + x2 * (b2 + x1 * b3)) &
                 + x4 * ((b4 + x1 * b5) + x2 * b6))
    else
        x1 = -x / (1.0_dp + x)
        x2 = x1 * x1
        x4 = x2 * x2

        log1p = -x1 * (((a1 + x1 * a2) + x2 * (a3 + x1 * a4)) &
                 + x4 * ((a5 + x1 * a6) + x2 * a7)) &
                 / (((1.0_dp + x1 * b1) + x2 * (b2 + x1 * b3)) &
                 + x4 * ((b4 + x1 * b5) + x2 * b6))
    endif
    ! ===============================================


    return

end function dlog1p

