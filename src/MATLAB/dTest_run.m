%
% ==========================================================================
%
% This script computes the gravitation potential of the Bennu asteroid.
% The shape of the asteroid is given by a surface spherical harmonic
% coefficients (in metres; the structure of the file is explained, for
% instance, in "./danm_bnm.f90")
%
% ../../data/Bennu_Shape_SHCs_to15.txt
%
% and the mass of the asteroid is assumed to be constant. Three evaluation
% points (one located inside and two located outside the body) are defined
% in
%
% ../../data/Bennu_Computing_points.txt
%
% The file contains three columns representing the spherical latitude,
% the spherical longitude (both in degrees) and the spherical radius
% (in metres) of the evaluation points.
%
% The loop over the evaluation points is paralellized using parfor. The
% parpool session is assumed to be open outside this script.
%
%
% REFERENCES: Fukushima, T., 2017. Precise and fast computation of the
%                 gravitational field of a general finite body and its
%                 application to the gravitational study of asteroid Eros.
%                 The Astronomical Journal 154(145):15pp,
%                 doi: 10.3847/1538-3881/aa88b8
%
% ==========================================================================

clear;
tic

% Input parameters
% -------------------------------------------------------------------------
nmax = 15;  % Maximum harmonic degree of the surface spherical
            % harmonic expansion of the body's surface (must
            % correspond to the imported "cs" matrix, see below)

rB = 0; % Radius (in metres) of a ball, the gravitational potential
        % of which will be subtracted from the output potential.
        % If "rb = 0", the gravitational field of the entire body
        % (as defined by its spherical harmonic coefficients) is computed.

G = 6.67384 * 10^(-11); % Newton's gravitational constant

rho = 1260; % Constant mass density

delta = 1e-16; % Relative error tolerance (see Fukushima 2017). For
               % instance, "delta = 1e-16" ensures about 15- or
               % 16-digit accuracy of the output gravitational potential.
% -------------------------------------------------------------------------


% Import spherical harmonic coefficients defining the gravitating body
% -------------------------------------------------------------------------
fprintf('Reading spherical harmonic coefficients...\n')

cs = load('../../data/Bennu_Shape_SHCs_to15.txt'); % Loads variable "cs" of
                                                % dimmensions
                                                % (nmax + 1, nmax + 1)
% -------------------------------------------------------------------------


% Import evaluation points
% -------------------------------------------------------------------------
fprintf('Reading evaluation points...\n')

points = load('../../data/Bennu_Computing_points.txt'); % Loads variables
                                                     %  "lat", "lon", "r"
                                                     % with spherical
                                                     % latitudes (degrees),
                                                     % longitudes (degrees)
                                                     % and radii (metres)
                                                     % of evaluation
                                                     % points, respectively

lat = points(:, 1);
lon = points(:, 2);
r = points(:, 3);
clear points
% -------------------------------------------------------------------------


% Transform latitudes and longitudes of the evaluation points to radians
% -------------------------------------------------------------------------
lat = lat * pi / 180;
lon = lon * pi / 180;
% -------------------------------------------------------------------------


% Set the reference radius to the mean sphere
% -------------------------------------------------------------------------
R0 = cs(1, 1);
% -------------------------------------------------------------------------


% Precompute "anm" and "bnm" coefficients for Legendre recurrences
% -------------------------------------------------------------------------
fprintf('Computing "anm" and "bnm" coefficients...\n')

anm = zeros(nmax + 1); % Initialization
bnm = anm;             % Initialization

for n = 0:nmax
    for m = 0:(n - 1)

        anm(n + 1, m + 1) = sqrt((2 * n + 1) * (2 * n -1) /...
                                 ((n + m) * (n - m)));

        bnm(n + 1, m + 1) = sqrt((2 * n + 1) * (n + m - 1) *...
                                 (n - m - 1) / ((n - m) * (n + m) *...
                                 (2 * n - 3)));
    end
end % End of the loop over harmonic orders
% -------------------------------------------------------------------------


% Computation of the gravitational potential
% -------------------------------------------------------------------------
fprintf('Computing the gravitational potential...\n')

npoints = length(lat); % Total number of evaluation points
V = zeros(npoints, 1); % Initialization

parfor i = 1:npoints

    lati = lat(i);
    loni = lon(i);
    ri = r(i);

    vi1 = dqde1(-pi / 2 - lati, 0, delta, lati, loni, ri, ...
                cs, nmax, rB, R0, anm, bnm);
    vi2 = dqde1(0, pi / 2 - lati, delta, lati, loni, ri, ...
                cs, nmax, rB, R0, anm, bnm);

    V(i) = vi1 + vi2;

end % End of the loop over evaluation points

V = G * R0^2 * rho * V;
% -------------------------------------------------------------------------


fprintf('\n')
fprintf('\n')
fprintf('Validation\n')
fprintf('=======================================\n')
fprintf('Latitudes of the evaluation points (in degrees)\n')
lat * 180 / pi
fprintf('Longitudes of the evaluation points (in degrees)\n')
lon * 180 / pi
fprintf('Radii of the evaluation points (in metres)\n')
r
fprintf('Computed values of the gravitation potential (in m^2 * s^-2)\n')
V
fprintf('Reference values of the gravitation potential (in m^2 * s^-2)\n')
Vref = [3.01907463537325338e-02; 1.7306053829097761e-02; ...
        5.23751595802091161e-03];
Vref
fprintf('=======================================\n')

toc
