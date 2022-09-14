function [Anm, Bnm] = danm_bnm(lat, nmax, cs, anm, bnm)
%
% ==================================================================
%
% DESCRIPTION: This function computes lumped coefficients for a single
%              latitude.
%
%              The algorithm used compute Legendre functions in this
%              subroutine is numerically stable up to degree 1800 at most.
%              Beyond this degree, one can use, for instance, the
%              extended-range artihmetic approach by Fukushima (2012).
%
%
% INPUTS: "lat"  -- Latitude of the computation point (a scalar).
%
%         "nmax" -- Maximum harmonic degree up to which the lumped
%                   coefficients are computed.
%
%         "cs"   -- Matrix of dimensions (nmax + 1, nnmax + 1) with
%                   the input spherical harmonic coefficients. The
%                   structure of the matrix is as follows:
%
%                   [C(0, 0)    S(1, 1) S(2, 1)    ...        S(nmax, 1)]
%                   [C(1, 0)    C(1, 1) S(2, 2)    ...        S(nmax, 2)]
%                   [C(2, 0)    C(2, 1) C(2, 2)    ...        S(nmax, 3)]
%                   [  .          .      .         .             .      ]
%                   [  .          .      .          .            .      ]
%                   [  .          .      .           .           .      ]
%                   [C(nmax, 0)                    ...     C(nmax, nmax)]
%
%                   where C(n, m) is the coefficient "C" of degree "n" and
%                   order "m", etc.
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
% OUTPUTS: "Anm" -- A vector of dimmension (nmax + 1) with lumped
%                   coefficients for the cosine spherical harmonic
%                   coefficients.
%
%          "Bnm" -- A vector of dimmension (nmax + 1) with lumped
%                   coefficients for the cosine spherical harmonic
%                   coefficients..
%
%
% REFERENCES: Fukushima, T., 2012. Numerical computation of spherical
%                 harmonics of arbitrary degree and order by extending
%                 exponent of floating point numbers.
%                 Journal of Geodesy 86:271--285.
%
% ==================================================================

% Initialization
Anm = zeros(1, nmax+1);
Bnm = zeros(1, nmax+1);

u = cos(lat);
t = sin(lat);
sqrt3 = sqrt(3);
um = u;

% Loop over harmonic orders
for m = 0:nmax
    
    if m == 0
               
        % -----------------------------------------------------------------
        % P00
        Pnm0 = 1;
        atmp = Pnm0 * cs(1, m + 1);
        
        
        % P10
        if nmax >= 1
            Pnm1 = sqrt3 * t;
            atmp = atmp + Pnm1 * cs(2, m + 1);
        end
         
        
        % P20, P30, ..., Pnmax,0
        if nmax >=2
            for n = 2:nmax
                Pnm2 = anm(n + 1, m + 1) * t * Pnm1 -...
                       bnm(n + 1, m + 1) * Pnm0;
                   
                atmp = atmp + Pnm2 * cs(n + 1, m + 1);
                   
                Pnm0 = Pnm1;
                Pnm1 = Pnm2;
            end            
        end
        
        btmp = 0;
        % -----------------------------------------------------------------
        
    else % Non-zonal harmonics
        
        % Sectorial harmonics
        % -----------------------------------------------------------------
        prod = 1;
        for n = 2:m
            n2 = 2 * n;
            prod = prod * (n2 + 1) / n2;
        end
        
        % Pmm
        Pnm0 = um * sqrt3 * sqrt(prod);
        um = um * u;
        
        atmp = Pnm0 * cs(m + 1, m + 1);
        btmp = Pnm0 * cs(m, m + 1);
        % -----------------------------------------------------------------
        
        
        % Tesseral harmonics
        % -----------------------------------------------------------------
        if m < nmax
            
            % Pm+1,m
            Pnm1 = anm(m + 2, m + 1) * t * Pnm0;
            
            atmp = atmp + Pnm1 * cs(m + 2, m + 1);
            btmp = btmp + Pnm1 * cs(m, m + 2);
            
            % Pm+2,m, Pm+3,m, ..., Pnmax,m
            for n = (m + 2):nmax
                Pnm2 = anm(n + 1, m + 1) * t * Pnm1 -...
                       bnm(n + 1, m + 1) * Pnm0;
                   
                atmp = atmp + Pnm2 * cs(n + 1, m + 1);
                btmp = btmp + Pnm2 * cs(m, n + 1);
                
                Pnm0 = Pnm1;
                Pnm1 = Pnm2;
            end
            
        end
        % -----------------------------------------------------------------
        
    end
    
    
    % ---------------------------------------------------------------------
    Anm(m + 1) = atmp;
    Bnm(m + 1) = btmp;
    % ---------------------------------------------------------------------
    
end
