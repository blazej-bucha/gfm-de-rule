function qlog1p(x) result(log1p)
!
! ==================================================================
!
!   Quadruple precision function program
!   to compute a special exponential function,
!       log1p(x) = log(1+x),
!   accurately by minimax rational approximation
!
!   "qlog1p" keeps 33 digit accuracy for arbitrary argument
!   while runs 9 times slower than qlog(x) itself
!
!   Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>
!   Date: 2019/10/30
!
!
!
!   Minor (formal) modifications by Blazej Bucha
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
    real(qp) :: x, x1, x2, x4, x8
    real(qp) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13
    real(qp) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13
    real(qp) :: log1p
    ! ===============================================


    ! Main computational part of the routine
    ! ===============================================
    parameter (a1 = 0.999999999999999999999999999999999989679405409_qp)
    parameter (a2 = 5.43952695470870990343440472736508508440853974_qp)
    parameter (a3 = 12.8687954156331531486576454765577285904854589_qp)
    parameter (a4 = 17.3813399893763161892491716455513072747241961_qp)
    parameter (a5 = 14.7937221247975315812946713068099324296210500_qp)
    parameter (a6 = 8.26554125834405953536667781138761190677622028_qp)
    parameter (a7 = 3.06422110888383600247671821221533143683683084_qp)
    parameter (a8 = 0.745288118102087966956230766831419822309727243_qp)
    parameter (a9 = 0.115047592893588673011032562622446366266055797_qp)
    parameter (a10 = 0.0105968202344748100801276898051500664263281933_qp)
    parameter (a11 = 0.000522466799810767549674343323841954134608619668_qp)
    parameter (a12 = 0.0000112217383140869479740242371068976645949743233_qp)
    parameter (a13 = 6.39672301036246645989908315040681256718633817e-8_qp)
    parameter (b1 = 5.93952695470870990343440472736506510635091212_qp)
    parameter (b2 = 15.5052255596541747670415145069133006422898887_qp)
    parameter (b3 = 23.4041104509671669382917939890824601741965693_qp)
    parameter (b4 = 22.6122505690735676039519980345629102098426624_qp)
    parameter (b5 = 14.6253640581969227356522409769061106689586552_qp)
    parameter (b6 = 6.43653279869637374357975617245086508089670491_qp)
    parameter (b7 = 1.92137412606261539134961029284884346128770878_qp)
    parameter (b8 = 0.381097234341555148065443000297505125713856247_qp)
    parameter (b9 = 0.0482176170045436717981339396150227533590537972_qp)
    parameter (b10 = 0.00362957310087890393686248940312126067485435786_qp)
    parameter (b11 = 0.000144194286147019122572434476322438425343273856_qp)
    parameter (b12 = 2.40586135243166027564864363686803405251652763e-6_qp)
    parameter (b13 = 9.53384644636009246025258296194840899175886989e-9_qp)
    !
    !   bit loss occurs when -1/2 <= x <= 1
    !
    if(x .lt. -0.5_qp .or. x .gt. 1.0_qp) then
        log1p = log(1.0_qp + x)
    !
    !   1/2 <= x < 0: log1p(x)=-log1p(x/(1+x))
    !
    elseif(x .lt. 0.0_qp) then
        x1 = -x / (1.0_qp + x)
        x2 = x1 * x1
        x4 = x2 * x2
        x8 = x4 * x4
        log1p = -x1 * ((((a1 + x1 * a2) + x2 * (a3 + x1 * a4)) &
                + x4 * ((a5 + x1 * a6) + x2 * (a7 + x1 * a8))) &
                + x8 * (((a9 + x1 * a10) + x2 * (a11 + x1 * a12)) + x4 * a13)) &
                / ((((1.0_qp + x1 * b1) + x2 * (b2 + x1 * b3)) &
                + x4 * ((b4 + x1 * b5) + x2 * (b6 + x1 * b7))) &
                + x8 * (((b8 + x1 * b9) + x2 * (b10 + x1 * b11)) + x4 &
                * (b12 + x1 * b13)))
    !
    !   0 <= x <= 1
    !
    else
        x1 = x
        x2 = x1 * x1
        x4 = x2 * x2
        x8 = x4 * x4
        log1p = x1 * ((((a1 + x1 * a2) + x2 * (a3 + x1 * a4)) &
                + x4 * ((a5 + x1 * a6) + x2 * (a7 + x1 * a8))) &
                + x8 * (((a9 + x1 * a10) + x2 * (a11 + x1 * a12)) + x4 * a13)) &
                / ((((1.0_qp + x1 * b1) + x2 * (b2 + x1 * b3)) &
                + x4 * ((b4 + x1 * b5) + x2 * (b6 + x1 * b7))) &
                + x8 * (((b8 + x1 * b9) + x2 * (b10 + x1 * b11)) &
                + x4 * (b12 + x1 * b13)))
    endif
    ! ===============================================


    return

end function qlog1p
