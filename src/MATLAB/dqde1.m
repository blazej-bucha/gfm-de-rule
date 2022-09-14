function s = dqde1(a, b, delta, lat, lon, r, cs, nmax, rB, R0, ...
                   anm, bnm)
%
% ==================================================================
%
% DESCRIPTION: This function computes one of the two integrals (depends on
%              the specified limits) from Eq. (28) of Fukushima (2017) via
%              the double exponential rule.
%
%              The function is designed such that it can be used to
%              compute the gravitational potential implied by a body
%              defined via a finite surface spherical harmonic expansion or
%              a gravitational potential of the same body but residual to
%              a ball of a specified radius.
%
%              The function is based on the Fortran routine by
%              Fukushima (2017), while some necessary modifications were
%              introduced. Note that the function can only be used for 
%              bodies with a constant mass density.
%
%
% INPUTS: "a"     -- Lower integration limit, a scalar.
%
%         "b"     -- Upper integration limit, a scalar.
%
%         "delta" -- Relative error tolerance (see Fukushima 2017). For
%                    instance, "delta = 1e-16" ensures about 15- or 
%                    16-digit accuracy of the output gravitational
%                    potential.
%
%         "lat"   -- Latitude of the evaluation point in radians, a scalar.
%
%         "lon"   -- Longitude of the evaluation points in radians, 
%                    a scalar.
%
%         "r"     -- Spherical radius of the evaluation point in metres,
%                    a scalar.
%
%         "cs"    -- Matrix of dimensions (nmax + 1, nnmax + 1) with
%                    the input spherical harmonic coefficients defining the
%                    shape of the gravitating body. Each coefficients is 
%                    given in metres. The structure of the matrix is as
%                    follows:
%
%                    [C(0, 0) S(1, 1) S(2, 1)    ...        S(nmax, 1)]
%                    [C(1, 0) C(1, 1) S(2, 2)    ...        S(nmax, 2)]
%                    [C(2, 0) C(2, 1) C(2, 2)    ...        S(nmax, 3)]
%                    [  .       .      .         .             .      ]
%                    [  .       .      .          .            .      ]
%                    [  .       .      .           .           .      ]
%                    [C(nmax, 0)                 ...     C(nmax, nmax)]
%
%                    where C(n, m) is the coefficient "C" of degree "n" and
%                    order "m", etc.
%
%         "nmax"  -- Maximum harmonic degree of the surface spherical 
%                    harmonic expansion of the body (must correspond with
%                    the "cs" matrix).
%
%         "rB"    -- Radius (in metres) of a ball, the gravitational 
%                    potential of which is subtracted from the output
%                    potential.
%                    If "rb = 0", the gravitational field of the
%                    entire body (as defined by its spherical harmonic
%                    coefficients) is computed.
%
%         "R0"    -- Reference radius in metres, a scalar. This value can
%                    efficiently be set to the mean radius of the 
%                    gravitating body.
%
%         "anm"   -- Coefficients for recurrence relations to compute
%                    Legendre functions (e.g., Fukushima 2012). The 
%                    algorithm used compute Legendre functions is 
%                    numerically stable up to degree 1800 at most. Beyond
%                    this degree, one can use, for instance, the 
%                    extended-range arithmetic approach by
%                    Fukushima (2012).
%
%                    "anm" is a matrix of dimmensions (nmax + 1, nmax + 1)
%                    and has the following form
%
%                     [a(0, 0)    0.0        0.0      ...           0.0   ]
%               anm = [a(1, 0)    a(1, 1)    0.0      ...           0.0   ]
%                     [a(2, 0)    a(2, 1)    a(2, 2)  ...           0.0   ]
%                     [   .          .          .     .              .    ]
%                     [   .          .          .      .             .    ]
%                     [   .          .          .       .            .    ]
%                     [a(nmax, 0) a(nmax, 1)          ...    a(nmax, nmax)]
%
%                    where "a(n, m)" is the "anm "coefficient of degree "n"
%                    and order "m" (see, e.g., Fukushima 2012).
%
%         "bnm"   -- Coefficients for recurrence relations to compute
%                    Legendre functions (e.g., Fukushima 2012). The 
%                    algorithm used compute Legendre functions is 
%                    numerically stable up to degree 1800 at most. Beyond
%                    this degree, one can use, for instance, the 
%                    extended-range arithmetic approach by
%                    Fukushima (2012).
%
%                    "bnm" is a matrix of dimmensions (nmax + 1, nmax + 1)
%                    and has the following form
%
%                     [b(0, 0)    0.0        0.0      ...           0.0   ]
%               bnm = [b(1, 0)    b(1, 1)    0.0      ...           0.0   ]
%                     [b(2, 0)    b(2, 1)    b(2, 2)  ...           0.0   ]
%                     [   .          .          .     .              .    ]
%                     [   .          .          .      .             .    ]
%                     [   .          .          .       .            .    ]
%                     [b(nmax, 0) b(nmax, 1)          ...    b(nmax, nmax)]
%
%                    where "b(n, m)" is the "bnm "coefficient of degree "n"
%                    and order "m" (see, e.g., Fukushima 2012).
%
%
% OUTPUTS: "s" -- Output value of the integral in m^2 * s^(-2).
%
%
% REFERENCES: Fukushima, T., 2012. Numerical computation of spherical
%                 harmonics of arbitrary degree and order by extending
%                 exponent of floating point numbers. Journal of 
%                 Geodesy 86:271--285.
%
%             Fukushima, T., 2017. Precise and fast computation of the
%                 gravitational field of a general finite body and its
%                 application to the gravitational study of asteroid Eros.
%                 The Astronomical Journal 154(145):15pp,
%                 doi: 10.3847/1538-3881/aa88b8.
%
% ==================================================================


twopi = 6.28318530717958648;
mmax = 1024;
safety = 10;

if delta > 1e-4
    hmax = 2.75 + (-log10(delta)) * 0.75;
elseif delta > 1e-6
    hmax = 3.75 + (-log10(delta)) * 0.5;
elseif delta > 1e-10
    hmax = 5.25 + (-log10(delta)) * 0.25;
else
    hmax = 6.5 + (-log10(delta)) * 0.125;
end

sdelta = safety * delta;
factor = 1 - log(sdelta);
deltah = sqrt(sdelta);
h0 = hmax / factor;
eph = exp(h0);
emh = 1 / eph;
deltat = exp(-emh * factor);
apb = a + b;
bma = b - a;
m = 1;
sr = dqde2(0, twopi, delta, lat, lon, r, apb * 0.5, cs, ...
           nmax, rB, R0, anm, bnm) * (bma * 0.25);
s = sr * (2 * twopi);
err = abs(s) * deltat;
h = 2 * h0;


con1 = 1;
while con1

    sprev = s;
    srprev = sr;
    t = h * 0.5;
    
    con2 = 1;
    while con2
        
        et = exp(t);
        ep = twopi * et;
        em = twopi / et;
        
        con3 = 1;
        while con3
            
            xw = 1 / (1 + exp(ep - em));
            eppem = ep + em;
            xa = bma * xw;
            wg = xa * (1 - xw);
            fa = dqde2(0, twopi, delta, lat, lon, r, a + xa, cs,...
                       nmax, rB, R0, anm, bnm) * wg;
            fb = dqde2(0, twopi, delta, lat, lon, r, b - xa, cs,...
                       nmax, rB, R0, anm, bnm) * wg;
            fapfb = fa + fb;
            sr = sr + fapfb;
            s = s + fapfb * eppem;
            errt = (abs(fa) + abs(fb)) * eppem;
            
            if m == 1
                err = err + errt * deltat;
            end
            
            ep = ep * eph;
            em = em * emh;
            
            if ~(errt > err || xw > deltah)
                con3 = 0;
            end
            
        end
        
        t = t + h;
        
        if ~(t < h0)
            con2 = 0;
        end

    end
    
    if m == 1
        errh = (err / deltat) * deltah * h0;
        errd = 1 + 2 * errh;
    else
        errd = h * (abs(s - 2 * sprev)) + 4 * abs(sr - 2 * srprev);
    end
    
    h = h * 0.5;
    m = m * 2;
    
    if ~(errd > errh && m < mmax)
        con1 = 0;
    end
end

s = s * h;
