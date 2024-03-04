subroutine qanm_bnm(lat, nmax, cs, anm, bnm, lcanm, lcbnm)
!
! ==================================================================
!
! DESCRIPTION: This subroutine computes lumped coefficients for a single
!              latitude. The computation is performed in quadruple precision.
!
!              The algorithm used compute Legendre functions in this subroutine
!              is numerically stable up to degree 21600 at most. Beyond this
!              degree, one can use, for instance, the extended-range artihmetic
!              approach by Fukushima (2012).
!
!
! INPUTS: "lat"  -- Latitude of the computation point (a scalar).
!
!         "nmax" -- Maximum harmonic degree up to which the lumped coefficients
!                   are computed.
!
!         "cs"   -- Matrix of dimensions (nmax + 1, nnmax + 1) with
!                   the input spherical harmonic coefficients. The structure
!                   of the matrix is as follows:
!
!                   [C(0, 0)    S(1, 1) S(2, 1)    ...        S(nmax, 1)]
!                   [C(1, 0)    C(1, 1) S(2, 2)    ...        S(nmax, 2)]
!                   [C(2, 0)    C(2, 1) C(2, 2)    ...        S(nmax, 3)]
!                   [  .          .      .         .             .      ]
!                   [  .          .      .          .            .      ]
!                   [  .          .      .           .           .      ]
!                   [C(nmax, 0)                    ...     C(nmax, nmax)]
!
!                   where C(n, m) is the coefficient "C" of degree "n" and
!                   order "m", etc.
!
!          "anm" -- Coefficients for recurrence relations to compute Legendre
!                   functions (e.g., Fukushima 2012). The algorithm used
!                   compute Legendre functions is numerically stable up to
!                   degree 21600 at most (see Fukushima 2012). Beyond this
!                   degree, one can use, for instance, the extended-range
!                   arithmetic approach by Fukushima (2012).
!
!                   "anm" is a matrix of dimmensions (nmax + 1, nmax + 1) and
!                   has the following form
!
!                         [a(0, 0)    0.0_qp     0.0_qp   ...           0.0_qp]
!                   anm = [a(1, 0)    a(1, 1)    0.0_qp   ...           0.0_qp]
!                         [a(2, 0)    a(2, 1)    a(2, 2)  ...           0.0_qp]
!                         [   .          .          .     .              .    ]
!                         [   .          .          .      .             .    ]
!                         [   .          .          .       .            .    ]
!                         [a(nmax, 0) a(nmax, 1)          ...    a(nmax, nmax)]
!
!                   where "a(n, m)" is the "anm "coefficient of degree "n" and
!                   order "m" (see, e.g., Fukushima 2012).
!
!          "bnm" -- Coefficients for recurrence relations to compute Legendre
!                   functions (e.g., Fukushima 2012). The algorithm used
!                   compute Legendre functions is numerically stable up to
!                   degree 21600 at most (see Fukushima 2012). Beyond this
!                   degree, one can use, for instance, the extended-range
!                   arithmetic approach by Fukushima (2012).
!
!                   "bnm" is a matrix of dimmensions (nmax + 1, nmax + 1) and
!                   has the following form
!
!                         [b(0, 0)    0.0_qp     0.0_qp   ...           0.0_qp]
!                   bnm = [b(1, 0)    b(1, 1)    0.0_qp   ...           0.0_qp]
!                         [b(2, 0)    b(2, 1)    b(2, 2)  ...           0.0_qp]
!                         [   .          .          .     .              .    ]
!                         [   .          .          .      .             .    ]
!                         [   .          .          .       .            .    ]
!                         [b(nmax, 0) b(nmax, 1)          ...    b(nmax, nmax)]
!
!                   where "b(n, m)" is the "bnm" coefficient of degree "n" and
!                   order "m" (see, e.g., Fukushima 2012).
!
!
! OUTPUTS: "lcanm" -- A vector of dimmension (nmax + 1) with lumped 
!                     coefficients for the cosine spherical harmonic 
!                     coefficients.
!
!          "lcbnm" -- A vector of dimmension (nmax + 1) with lumped 
!                     coefficients for the sine spherical harmonic 
!                     coefficients.
!
!
! REFERENCES: Fukushima, T., 2012. Numerical computation of spherical
!                 harmonics of arbitrary degree and order by extending exponent
!                 of floating point numbers. Journal of Geodesy 86:271--285.
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
    implicit none

    integer,  intent(in) :: nmax
    real(qp), intent(in) :: lat
    real(qp), intent(in) :: cs(nmax + 1, nmax + 1)
    real(qp), intent(in) :: anm(nmax + 1, nmax + 1), bnm(nmax + 1, nmax + 1)

    integer  :: n, n2, m
    real(qp) :: pnm0, pnm1, pnm2
    real(qp) :: u, um, t, sqrt3, prod, atmp, btmp

    real(qp), intent(out) :: lcanm(nmax + 1), lcbnm(nmax + 1)
    ! ===============================================


    ! Main computational part of the routine
    ! ===============================================
    u = cos(lat)
    t = sin(lat)
    sqrt3 = sqrt(3.0_qp)

    ! Loop over harmonic orders
    do m = 0, nmax

        if (m == 0) then

            um = u

            ! Zonal harmonics
            ! --------------------------------------------------------------
            ! P00
            pnm0 = 1.0_qp
            atmp = pnm0 * cs(1, m + 1)


            ! P10
            if (nmax >= 1) then
                pnm1 = sqrt3 * t
                atmp = atmp + pnm1 * cs(2, m + 1)
            end if


            ! P20, P30, ..., Pnmax,0
            if (nmax >= 2) then
                do n = 2, nmax
                    pnm2 = anm(n + 1, m + 1) * t * pnm1 &
                           - bnm(n + 1, m + 1) * pnm0

                    atmp = atmp + pnm2 * cs(n + 1, m + 1)

                    pnm0 = pnm1
                    pnm1 = pnm2
                end do
            end if

            btmp = 0.0_qp
            ! --------------------------------------------------------------

        else ! Non-zonal harmonics

            ! Sectorial harmonics
            ! --------------------------------------------------------------
            prod = 1.0_qp
            do n = 2, m
                n2 = 2 * n
                prod = prod * (real(n2, qp) + 1.0_qp) / real(n2, qp)
            end do

            ! Pmm
            pnm0 = um * sqrt3 * sqrt(prod)
            um = um * u

            atmp = pnm0 * cs(m + 1, m + 1)
            btmp = pnm0 * cs(m, m + 1)
            ! --------------------------------------------------------------


            ! Tesseral harmonics
            ! --------------------------------------------------------------
            if (m < nmax) then

                ! Pm+1,m
                pnm1 = anm(m + 2, m + 1) * t * pnm0

                atmp = atmp + pnm1 * cs(m + 2, m + 1)
                btmp = btmp + pnm1 * cs(m, m + 2)


                ! Pm+2,m, Pm+3,m, ..., Pnmax,m
                do n = (m + 2), nmax
                    pnm2 = anm(n + 1, m + 1) * t * pnm1 &
                           - bnm(n + 1, m + 1) * pnm0

                    atmp = atmp + pnm2 * cs(n + 1, m + 1)
                    btmp = btmp + pnm2 * cs(m, n + 1)

                    pnm0 = pnm1
                    pnm1 = pnm2
                end do

            end if
            ! --------------------------------------------------------------

        end if


        ! Lumped coefficients
        ! --------------------------------------------------------------
        lcanm(m + 1) = atmp
        lcbnm(m + 1) = btmp
        ! --------------------------------------------------------------

    end do
    ! ===============================================

end subroutine qanm_bnm
