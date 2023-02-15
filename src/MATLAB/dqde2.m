function s = dqde2(a, b, delta, lat, lon, r, xi, cs, nmax, ...
                   rB, R0, anm, bnm)
%
% ==================================================================
%
% DESCRIPTION: This subroutine computes the integral from Eq. (29) of Fukushima
%              (2017) (except for the "G * R_0^2 * rho_n" factor from that
%              equation) via the double exponential rule. The computation is
%              performed in double precision.
%
%              The function is designed such that it can be used to
%              compute the gravitational potential implied by a body
%              defined via a finite surface spherical harmonic expansion
%              or a gravitational potential of the same body but residual
%              to a ball of a specified radius.
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
%         "xi"    -- Temporary latitude (dummy integration variable) in 
%                    radians resulting from the "dqde1.m" function.
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

% Computation of "sr"
%--------------------------------------------------------------------------
lattmp = lat + xi;
clattmp = cos(lattmp);
sxi22 = sin(xi / 2)^2;
clattmp_clat = cos(lattmp) * cos(lat);
eta = apb / 2;
alpha = r / R0;
alpha2 = 2 * alpha;
alpha3 = 3 * alpha;
alpha_pow2 = alpha^2;

xitmp = xi + lat;
[Anm, Bnm] = danm_bnm(xitmp, nmax, cs, anm, bnm);
m_temp = [0:nmax]';
cosla = cos(m_temp * (eta + lon));
sinla = sin(m_temp * (eta + lon));
r_surf = Anm * cosla + Bnm * sinla;

zeta_T = (r_surf - r) / R0;
B = alpha2 * (sxi22 + clattmp_clat * sin(eta / 2)^2);
A = B * (alpha2 - B);
C0 = alpha_pow2 - alpha3 * B + 1.5 * B^2;
D0_T = (alpha2 - 1.5 * B) + 0.5 * zeta_T;
S_T = sqrt(A + (B + zeta_T)^2);
zeta_B = (rB - r) / R0;
S_B = sqrt(A + (B + zeta_B)^2);
beta = zeta_T - zeta_B;
E0 = 0.5;
T = (zeta_T + zeta_B + 2 * B) / (S_T + S_B);
if (T <= -1)
    s = 0; % This prevents possible numerical inaccuracies
    return
end
if (zeta_B + B) > 0
    T_temp = ((1 + T) * beta) / (zeta_B + B + S_B);
else
    T_temp = ((1 + T) * beta) * ((-zeta_B - B + S_B) / A);
end
K = C0 * log1p(T_temp) + (D0_T * T + S_B * E0) * beta;

sr = K * (bma * 0.25);
%--------------------------------------------------------------------------

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
            
            % Computation of "fa"
            %--------------------------------------------------------------
            eta1 = a + xa;
            eta_temp1 = eta1 + lon;
            eta2 = b - xa;
            eta_temp2 = eta2 + lon;
            
            m_temp_eta_temp1_eta_temp2 = m_temp * [eta_temp1 eta_temp2];
            cosla = cos(m_temp_eta_temp1_eta_temp2);
            sinla = sin(m_temp_eta_temp1_eta_temp2);
            r12 = Anm * cosla + Bnm * sinla;
            r_surf = r12(1);
            
            zeta_T = (r_surf - r) / R0;
            B = alpha2 * (sxi22 + clattmp_clat * sin(eta1 / 2)^2);
            A = B * (alpha2 - B);
            C0 = alpha_pow2 - alpha3 * B + 1.5 * B^2;
            D0_T = (alpha2 - 1.5 * B) + 0.5 * zeta_T;
            S_T = sqrt(A + (B + zeta_T)^2);
            S_B = sqrt(A + (B + zeta_B)^2);
            
            beta = zeta_T - zeta_B;
            Ta = (zeta_T + zeta_B + 2 * B) / (S_T + S_B);
            if (zeta_B + B) > 0
                Ta_temp = ((1 + Ta) * beta) / (zeta_B + B + S_B);
            else
                Ta_temp = ((1 + Ta) * beta) * ((-zeta_B - B + S_B) / A);
            end
            K = C0 * log1p(Ta_temp) + (D0_T * Ta + S_B * E0) * beta;

            fa = K * wg;
            %--------------------------------------------------------------

            
            %Computation of "fb"
            %--------------------------------------------------------------
            r_surf = r12(2);
            
            zeta_T = (r_surf - r) / R0;
            B = alpha2 * (sxi22 + clattmp_clat * sin(eta2 / 2)^2);
            A = B * (alpha2 - B);
            C0 = alpha_pow2 - alpha3 * B + 1.5 * B^2;
            D0_T = (alpha2 - 1.5 * B) + 0.5 * zeta_T;
            S_T = sqrt(A + (B + zeta_T)^2);
            S_B = sqrt(A + (B + zeta_B)^2);
            
            beta = zeta_T - zeta_B;
            Tb = (zeta_T + zeta_B + 2 * B) / (S_T + S_B);
            if (zeta_B + B) > 0
                Tb_temp = ((1 + Tb) * beta) / (zeta_B + B + S_B);
            else
                Tb_temp = ((1 + Tb) * beta) * ((-zeta_B - B + S_B) / A);
            end
            K = C0 * log1p(Tb_temp) + (D0_T * Tb + S_B * E0) * beta;
            
            fb = K * wg;
            %--------------------------------------------------------------
            
            if (Ta <= -1 || Tb <= -1)
                fa = 0; % This prevents possible numerical inaccuracies
                fb = 0; % This prevents possible numerical inaccuracies
            end
            
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
    
    if m==1
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

s = s * clattmp * h;
