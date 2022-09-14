subroutine qqde1(a, b, delta, lat, lon, r, cs, nmax, rb, r0, anm, bnm, s)
!
! ==================================================================
!
! DESCRIPTION: This subroutine computes one of the two integrals (depends on 
!              the specified limits) from Eq. (28) of Fukushima (2017) via the
!              double exponential rule. The computation is performed in
!              quadruple precision.
!
!              The subroutine is designed such that it can be used to
!              compute the gravitational potential implied by a body defined
!              via a finite surface spherical harmonic expansion or a
!              gravitational potential of the same body but residual to a ball
!              of a specified radius.
!
!              The routine was take from Fukushima (2017), while some necessary
!              modifications were introduced. Note that the routine can only
!              be used for bodies with a constant mass density.
!
!
! INPUTS: "a"     -- Lower integration limit, a scalar.
!
!         "b"     -- Upper integration limit, a scalar.
!
!         "delta" -- Relative error tolerance (see Fukushima 2017). For
!                    instance, "delta = 1e-32_qp" ensures about 32-digit
!                    accuracy of the output gravitational potential.
!
!         "lat"   -- Latitude of the evaluation point in radians, a scalar.
!
!         "lon"   -- Longitude of the evaluation points in radians, a scalar.
!
!         "r"     -- Spherical radius of the evaluation point in metres,
!                    a scalar.
!
!         "cs"    -- Matrix of dimensions (nmax + 1, nnmax + 1) with
!                    the input spherical harmonic coefficients defining the
!                    shape of the gravitating body. Each coefficients is given
!                    in metres. The structure of the matrix is as follows:
!
!                    [C(0, 0)    S(1, 1) S(2, 1)    ...        S(nmax, 1)]
!                    [C(1, 0)    C(1, 1) S(2, 2)    ...        S(nmax, 2)]
!                    [C(2, 0)    C(2, 1) C(2, 2)    ...        S(nmax, 3)]
!                    [  .          .      .         .             .      ]
!                    [  .          .      .          .            .      ]
!                    [  .          .      .           .           .      ]
!                    [C(nmax, 0)                    ...     C(nmax, nmax)]
!
!                    where C(n, m) is the coefficient "C" of degree "n" and
!                    order "m", etc.
!
!         "nmax"  -- Maximum harmonic degree of the surface spherical harmonic
!                    expansion of the body (must correspond with the "cs"
!                    matrix).
!
!         "rb"    -- Radius (in metres) of a ball, the gravitational potential
!                    of which is subtracted from the output potential.
!                    If "rb = 0.0_qp", the gravitational field of the entire
!                    body (as defined by its spherical harmonic coefficients)
!                    is computed.
!
!         "r0"    -- Reference radius in metres, a scalar. This value can
!                    efficiently be set to the mean radius of the gravitating
!                    body.
!
!         "anm"   -- Coefficients for recurrence relations to compute Legendre
!                    functions (e.g., Fukushima 2012). The algorithm used
!                    compute Legendre functions is numerically stable up to
!                    degree 21600 at most. Beyond this degree, one can use,
!                    for instance, the extended-range arithmetic approach by
!                    Fukushima (2012).
!
!                    "anm" is a matrix of dimmensions (nmax + 1, nmax + 1) and
!                    has the following form
!
!                         [a(0, 0)    0.0_qp     0.0_qp   ...           0.0_qp]
!                   anm = [a(1, 0)    a(1, 1)    0.0_qp   ...           0.0_qp]
!                         [a(2, 0)    a(2, 1)    a(2, 2)  ...           0.0_qp]
!                         [   .          .          .     .              .    ]
!                         [   .          .          .      .             .    ]
!                         [   .          .          .       .            .    ]
!                         [a(nmax, 0) a(nmax, 1)          ...    a(nmax, nmax)]
!
!                    where "a(n, m)" is the "anm "coefficient of degree "n" and
!                    order "m" (see, e.g., Fukushima 2012).
!
!         "bnm"   -- Coefficients for recurrence relations to compute Legendre
!                    functions (e.g., Fukushima 2012). The algorithm used
!                    compute Legendre functions is numerically stable up to
!                    degree 21600 at most. Beyond this degree, one can use,
!                    for instance, the extended-range arithmetic approach by
!                    Fukushima (2012).
!
!                    "bnm" is a matrix of dimmensions (nmax + 1, nmax + 1) and
!                    has the following form
!
!                         [b(0, 0)    0.0_qp     0.0_qp   ...           0.0_qp]
!                   bnm = [b(1, 0)    b(1, 1)    0.0_qp   ...           0.0_qp]
!                         [b(2, 0)    b(2, 1)    b(2, 2)  ...           0.0_qp]
!                         [   .          .          .     .              .    ]
!                         [   .          .          .      .             .    ]
!                         [   .          .          .       .            .    ]
!                         [b(nmax, 0) b(nmax, 1)          ...    b(nmax, nmax)]
!
!                    where "b(n, m)" is the "bnm" coefficient of degree "n" and
!                    order "m" (see, e.g., Fukushima 2012).
!
!
! OUTPUTS: "s" -- Output value of the integral in m**2 * s**-2.
!
!
! REFERENCES: Fukushima, T., 2012. Numerical computation of spherical
!                 harmonics of arbitrary degree and order by extending exponent
!                 of floating point numbers. Journal of Geodesy 86:271--285.
!
!             Fukushima, T., 2017. Precise and fast computation of the
!                 gravitational field of a general finite body and its
!                 application to the gravitational study of asteroid Eros. The
!                 Astronomical Journal 154(145):15pp,
!                 doi: 10.3847/1538-3881/aa88b8
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
    real(qp), intent(in) :: a, b, delta
    real(qp), intent(in) :: lat, lon, r
    real(qp), intent(in) :: rb, r0
    real(qp), intent(in) :: cs(nmax + 1, nmax + 1)
    real(dp), intent(in) :: anm(nmax + 1, nmax + 1), bnm(nmax + 1, nmax + 1)

    integer  :: MMAX, m
    real(qp) :: TWOPI, SAFETY, hmax, sdelta, factor
    real(qp) :: deltah, h0, eph, emh, deltat
    real(qp) :: apb, bma, sr, h, sprev, srprev, t, et, ep, em, eppem, xw, xa
    real(qp) :: wg, fa, fb, fapfb, err, errt, errh, errd
    real(qp) :: fatmp, fbtmp

    real(qp), intent(out) :: s
    ! ===============================================


    ! Main computational part of the routine
    ! ===============================================
    parameter (TWOPI = 2.0_qp * pi_qp, MMAX = 4092, SAFETY = 10.0_qp)

    if(delta .gt. 1e-4_qp) then
        hmax = 2.75_qp + (-log10(delta)) * 0.75_qp
    elseif(delta .gt. 1e-6_qp) then
        hmax = 3.75_qp + (-log10(delta)) * 0.5_qp
    elseif(delta .gt. 1e-10_qp) then
        hmax = 5.25_qp + (-log10(delta)) * 0.25_qp
    else
        hmax = 6.5_qp + (-log10(delta)) * 0.125_qp
    endif

    sdelta = SAFETY * delta
    factor = 1.0_qp - log(sdelta)
    deltah = sqrt(sdelta)
    h0 = hmax / factor
    eph = exp(h0)
    emh = 1.0_qp / eph
    deltat = exp(-emh * factor)
    apb = a + b
    bma = b - a
    m = 1

    call qqde2(0.0_qp, TWOPI, delta, lat, lon, r, apb * 0.5_qp, &
               cs, nmax, rb, r0, anm, bnm, &
               sr)

    sr = sr * bma * 0.25_qp
    s = sr * (2.0_qp * TWOPI)
    err = abs(s) * deltat
    h = 2.0_qp * h0

  1 continue

    sprev = s
    srprev = sr
    t = h * 0.5_qp

  2 continue

    et = exp(t)
    ep = TWOPI * et
    em = TWOPI / et

  3 continue

    xw = 1.0_qp / (1.0_qp + exp(ep - em))
    eppem = ep + em
    xa = bma * xw
    wg = xa * (1.0_qp - xw)

    call qqde2(0.0_qp, TWOPI, delta, lat, lon, r, a + xa, &
               cs, nmax, rb, r0, anm, bnm, &
               fatmp)
    fa = fatmp * wg

    call qqde2(0.0_qp, TWOPI, delta, lat, lon, r, b - xa, &
               cs, nmax, rb, r0, anm, bnm, &
               fbtmp)
    fb = fbtmp * wg

    fapfb = fa + fb
    sr = sr + fapfb
    s = s + fapfb * eppem
    errt = (abs(fa) + abs(fb)) * eppem

    if (m .eq. 1) then
        err = err + errt * deltat
    end if

    ep = ep * eph
    em = em * emh

    if (errt .gt. err .or. xw .gt. deltah) then
        goto 3
    end if

    t = t + h

    if (t .lt. h0) then
        goto 2
    end if

    if (m .eq. 1) then
        errh = (err / deltat) * deltah * h0
        errd = 1.0_qp + 2.0_qp * errh
    else
        errd = h * (abs(s - 2.0_qp * sprev) &
               + 4.0_qp * abs(sr - 2.0_qp * srprev))
    end if

    h = h * 0.5_qp
    m = m * 2

    if (errd .gt. errh .and. m .lt. MMAX) then
        goto 1
    end if

    s = s * h

    return
    ! ===============================================

end subroutine qqde1
