subroutine dqde2(a, b, delta, lat, lon, r, xi, cs, nmax, rb, r0, anm, bnm, s)
!
! ==================================================================
!
! DESCRIPTION: This subroutine computes the integral from Eq. (29) of Fukushima
!              (2017) (except for the "G * R_0^2 * rho_n" factor from that
!              equation) via the double exponential rule. The computation is
!              performed in double precision.
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
!                    instance, "delta = 1e-16_dp" ensures about 15- or 16-digit
!                    accuracy of the output gravitational potential.
!
!         "lat"   -- Latitude of the evaluation point in radians, a scalar.
!
!         "lon"   -- Longitude of the evaluation points in radians, a scalar.
!
!         "r"     -- Spherical radius of the evaluation point in metres,
!                    a scalar.
!
!         "xi"    -- Temporary latitude (dummy integration variable) in radians
!                    resulting from the "dqde1.f90" subroutine.
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
!                    If "rb = 0.0_dp", the gravitational field of the entire
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
!                    degree 1800 at most. Beyond this degree, one can use,
!                    for instance, the extended-range arithmetic approach by
!                    Fukushima (2012).
!
!                    "anm" is a matrix of dimmensions (nmax + 1, nmax + 1) and
!                    has the following form
!
!                         [a(0, 0)    0.0_dp     0.0_dp   ...           0.0_dp]
!                   anm = [a(1, 0)    a(1, 1)    0.0_dp   ...           0.0_dp]
!                         [a(2, 0)    a(2, 1)    a(2, 2)  ...           0.0_dp]
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
!                    degree 1800 at most. Beyond this degree, one can use,
!                    for instance, the extended-range arithmetic approach by
!                    Fukushima (2012).
!
!                    "bnm" is a matrix of dimmensions (nmax + 1, nmax + 1) and
!                    has the following form
!
!                         [b(0, 0)    0.0_dp     0.0_dp   ...           0.0_dp]
!                   bnm = [b(1, 0)    b(1, 1)    0.0_dp   ...           0.0_dp]
!                         [b(2, 0)    b(2, 1)    b(2, 2)  ...           0.0_dp]
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

    integer, intent(in) :: nmax
    real(dp), intent(in) :: a, b, delta
    real(dp), intent(in) :: lat, lon, r, xi
    real(dp), intent(in) :: cs(nmax + 1, nmax + 1)
    real(dp), intent(in) :: rb, r0
    real(dp), intent(in) :: anm(nmax + 1, nmax + 1), bnm(nmax + 1, nmax + 1)

    integer  :: MMAX, n, m
    real(dp) :: TWOPI, HALFPI, SAFETY, hmax, sdelta, factor
    real(dp) :: deltah, h0, eph, emh, deltat
    real(dp) :: apb, bma, sr, h, sprev, srprev, t, et, ep, em, eppem, xw, xa
    real(dp) :: wg, fa, fb, fapfb, err, errt, errh, errd
    real(dp) :: lattmp, xitmp, clattmp, sxi22, clattmp_clat, eta, eta1, eta2
    real(dp) :: alpha, alpha2, alpha3, alpha_pow2
    real(dp) :: lcanm(nmax + 1), lcbnm(nmax + 1)
    real(dp) :: r_surf, zeta_t, zeta_b, ac, bc, c0c, d0_tc, s_tc, s_bc
    real(dp) :: beta, e0, ta, tb, tt, ta_temp, tb_temp, tt_temp, k
    real(dp) :: etatmp1, etatmp2
    real(dp) :: dlog1p

    real(dp), intent(out) :: s
    ! ===============================================


    ! Main computational part of the routine
    ! ===============================================
    parameter (TWOPI = 2.0_dp * pi_dp, HALFPI = pi_dp / 2.0_dp, &
               MMAX = 1024, SAFETY = 10.0_dp)

    if(delta .gt. 1e-4_dp) then
        hmax = 2.75_dp + (-log10(delta)) * 0.75_dp
    elseif(delta .gt. 1e-6_dp) then
        hmax = 3.75_dp + (-log10(delta)) * 0.5_dp
    elseif(delta .gt. 1e-10_dp) then
        hmax = 5.25_dp + (-log10(delta)) * 0.25_dp
    else
        hmax = 6.5_dp + (-log10(delta)) * 0.125_dp
    endif

    sdelta = SAFETY * delta
    factor = 1.0_dp - log(sdelta)
    deltah = sqrt(sdelta)
    h0 = hmax / factor
    eph = exp(h0)
    emh = 1.0_dp / eph
    deltat = exp(-emh * factor)
    apb = a + b
    bma = b - a
    m = 1

    ! Computation of "sr"
    ! -----------------------------------------------
    lattmp = lat + xi
    clattmp = cos(lattmp)
    sxi22 = sin(xi / 2.0_dp)**2
    clattmp_clat = cos(lattmp) * cos(lat)
    eta = apb / 2.0_dp
    alpha = r / r0
    alpha2 = 2.0_dp * alpha
    alpha3 = 3.0_dp * alpha
    alpha_pow2 = alpha**2

    xitmp = xi + lat
    call danm_bnm(xitmp, nmax, cs, anm, bnm, lcanm, lcbnm)

    r_surf = 0.0_dp
    do n = 0, nmax
        r_surf = r_surf + cos(real(n, dp) * (eta + lon)) * lcanm(n + 1) &
                        + sin(real(n, dp) * (eta + lon)) * lcbnm(n + 1)
    end do

    zeta_t = (r_surf - r) / r0
    bc = alpha2 * (sxi22 + clattmp_clat * sin(eta / 2.0_dp)**2)
    ac = bc * (alpha2 - bc)
    c0c = alpha_pow2 - alpha3 * bc + 1.5_dp * bc**2
    d0_tc = (alpha2 - 1.5_dp * bc) + 0.5_dp * zeta_t
    s_tc = sqrt(ac + (bc + zeta_t)**2)
    zeta_b = (rb - r) / r0
    s_bc = sqrt(ac + (bc + zeta_b)**2)
    beta = zeta_t - zeta_b
    e0 = 0.5_dp
    tt = (zeta_t + zeta_b + 2.0_dp * bc) / (s_tc + s_bc)
    if (tt .le. -1.0_dp) then
        s = 0.0_dp ! This prevents possible numerical inaccuracies
        return
    end if
    if ((zeta_b + bc) .gt. 0.0_dp) then
        tt_temp = ((1.0_dp + tt) * beta) / (zeta_b + bc + s_bc)
    else
        tt_temp = ((1.0_dp + tt) * beta) * ((-zeta_b - bc + s_bc) / ac)
    end if
    k = c0c * dlog1p(tt_temp) + (d0_tc * tt + s_bc * e0) * beta

    sr = k * (bma * 0.25_dp)
    ! -----------------------------------------------

    s = sr * (2.0_dp * TWOPI)
    err = abs(s) * deltat
    h = 2.0_dp * h0

  1 continue

    sprev = s
    srprev = sr
    t = h * 0.5_dp

  2 continue

    et = exp(t)
    ep = TWOPI * et
    em = TWOPI / et

  3 continue

    xw = 1.0_dp / (1.0_dp + exp(ep - em))
    eppem = ep + em
    xa = bma * xw
    wg = xa * (1.0_dp - xw)

    ! Computation of "fa"
    ! -----------------------------------------------
    eta1 = a + xa
    etatmp1 = eta1 + lon

    r_surf = 0.0_dp
    do n = 0, nmax
        r_surf = r_surf + cos(real(n, dp) * etatmp1) * lcanm(n + 1) &
                        + sin(real(n, dp) * etatmp1) * lcbnm(n + 1)
    end do

    zeta_t = (r_surf - r) / r0
    bc = alpha2 * (sxi22 + clattmp_clat * sin(eta1 / 2.0_dp)**2)
    ac = bc * (alpha2 - bc)
    c0c = alpha_pow2 - alpha3 * bc + 1.5_dp * bc**2
    d0_tc = (alpha2 - 1.5_dp * bc) + 0.5_dp * zeta_t
    s_tc = sqrt(ac + (bc + zeta_t)**2)
    s_bc = sqrt(ac + (bc + zeta_b)**2)
    beta = zeta_t - zeta_b
    ta = (zeta_t + zeta_b + 2.0_dp * bc) / (s_tc + s_bc)
    if ((zeta_b + bc) .gt. 0.0_dp) then
        ta_temp = ((1.0_dp + ta) * beta) / (zeta_b + bc + s_bc)
    else
        ta_temp = ((1.0_dp + ta) * beta) * ((-zeta_b - bc + s_bc) / ac)
    end if
    k = c0c * dlog1p(ta_temp) + (d0_tc * ta + s_bc * e0) * beta

    fa = k * wg
    ! -----------------------------------------------

    ! Computation of "fb"
    ! -----------------------------------------------
    eta2 = b - xa
    etatmp2 = eta2 + lon

    r_surf = 0.0_dp
    do n = 0, nmax
        r_surf = r_surf + cos(real(n, dp) * etatmp2) * lcanm(n + 1) &
                        + sin(real(n, dp) * etatmp2) * lcbnm(n + 1)
    end do

    zeta_t = (r_surf - r) / r0
    bc = alpha2 * (sxi22 + clattmp_clat * sin(eta2 / 2.0_dp)**2)
    ac = bc * (alpha2 - bc)
    c0c = alpha_pow2 - alpha3 * bc + 1.5_dp * bc**2
    d0_tc = (alpha2 - 1.5_dp * bc) + 0.5_dp * zeta_t
    s_tc = sqrt(ac + (bc + zeta_t)**2)
    s_bc = sqrt(ac + (bc + zeta_b)**2)
    beta = zeta_t - zeta_b
    tb = (zeta_t + zeta_b + 2.0_dp * bc) / (s_tc + s_bc)
    if ((zeta_b + bc) .gt. 0.0_dp) then
        tb_temp = ((1.0_dp + tb) * beta) / (zeta_b + bc + s_bc)
    else
        tb_temp = ((1.0_dp + tb) * beta) * ((-zeta_b - bc + s_bc) / ac)
    end if
    k = c0c * dlog1p(tb_temp) + (d0_tc * tb + s_bc * e0) * beta

    fb = k * wg
    ! -----------------------------------------------

    if ((ta .le. -1.0_dp) .or. (tb .le. -1.0_dp)) then
        fa = 0.0_dp ! This prevents possible numerical inaccuracies
        fb = 0.0_dp ! This prevents possible numerical inaccuracies
    end if

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
        errd = 1.0_dp + 2.0_dp * errh
    else
        errd = h * (abs(s - 2.0_dp * sprev) &
               + 4.0_dp * abs(sr - 2.0_dp * srprev))
    end if

    h = h * 0.5_dp
    m = m * 2

    if (errd .gt. errh .and. m .lt. MMAX) then
        goto 1
    end if

    s = s * clattmp * h

    return
    ! ===============================================

end subroutine dqde2
