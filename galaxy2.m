function [t,r] = galaxy2(tmax, level, mc, Nstars, gr_0, gv_0)
% galaxy: solves simple n-body problem of massless stars orbiting two
%         cores using an O(delta t^2) FDA.  Stars are given 
%         initially circular orbits in x-y plane around their respective 
%         cores.
%
% Inputs:
%   tmax:   (real scalar) Final solution time.
%   level:  (integer scalar) Discretization level.
%   mc:     (real vector) Vector of length cores x 1 containing core 
%           masses.
%   Nstars: (real vector) Vector of length Nstars x 1 containing number of
%           stars per core.
%   gr_0:   (real array) 2 x 3 array of initial core positions. 
%   gv_0:   (real array) 2 x 4 array of initial core velocities. 4th elt in
%           each row is for star rotation about core - if non-zero,
%           rotation is counter-clockwise when viewed from +z axis, else 
%           clockwise.

% Outputs:
%   t:    (real vector) Vector of length nt = 2^level + 1 containing time 
%         mesh.
%   r:    (real array) Array of size Nstars*length(mc)+length(mc) x 3 x nt
%         containing computed positions of all stars and cores at discrete
%         times t(n).
% /=======================================================================/

% Define number of time steps and create t, r, v, m, and a arrays of 
% appropriate size.
Nbodies = 2*Nstars + 2;

nt = 2^level + 1;
t = linspace(0.0, tmax, nt);
r = zeros(Nbodies, 3, nt);
v = zeros(Nbodies, 3, nt);
a2 = zeros(Nbodies, 3, nt);
m = zeros(Nbodies, 1);

% Determine discrete time step from t array
deltat = t(2) - t(1);

% Initialize core masses, positions
m(1) = mc(1);
m(Nstars+2) = mc(2);

% Initialize core positions
r(1, :, 1) = gr_0(1, :);
r(Nstars + 2, :, 1) = gr_0(2, :);

% Calculate and assign initial star positions s/t their orbits are circular
% about each core.
rmin = 10;
rmax = 30;

% Generate Nx1 array of random radii within the bounds [rmin, rmax]
r_rand = (rmax - rmin).*rand(Nstars, 1) + rmin;

% Generate Nx1 array of random angles within the bounds [0, 2*pi]
ang_rand = 2*pi.*rand(Nstars, 1);

% Assign initial x, y positions for stars.
% Compute values such that star orbits are circular
r(2:Nstars + 1, 1, 1) = r_rand.*cos(ang_rand) + r(1, 1, 1);  
r(2:Nstars + 1, 2, 1) = r_rand.*sin(ang_rand) + r(1, 2, 1); 

r(Nstars + 3:end, 1, 1) = r_rand.*cos(ang_rand) + r(Nstars + 2, 1, 1);  
r(Nstars + 3:end, 2, 1) = r_rand.*sin(ang_rand) + r(Nstars + 2, 2, 1);

% Adjust star z-coordinates to match their core
r(2:Nstars + 1, 3, 1) = r(1, 3, 1);
r(Nstars + 3:end, 3, 1) = r(Nstars + 2, 3, 1);

% Assign initial x, y velocities for stars about each core.
% Compute values such that star orbits are circular
vmag1 = sqrt(mc(1) ./ r_rand);
vmag2 = sqrt(mc(2) ./ r_rand);

% Direction of star rotation control:
if gv_0(1, 4) ~=0
    vmag1 = -vmag1;
end
if gv_0(2, 4) ~=0
    vmag2 = -vmag2;
end

v(2:Nstars + 1, 1, 1) = -vmag1.*sin(ang_rand);
v(2:Nstars + 1, 2, 1) = vmag1.*cos(ang_rand);


v(Nstars + 3:end, 1, 1) = -vmag2.*sin(ang_rand);
v(Nstars + 3:end, 2, 1) = vmag2.*cos(ang_rand);


% Add galactic velocities
v(1:Nstars + 1, :, 1) = v(1:Nstars + 1, :, 1) + gv_0(1);
v(Nstars + 2:end, :, 1) = v(Nstars + 2:end, :, 1) + gv_0(2);


% Initialize initial acceleration
a = nbodyaccn(m, r(:, :, 1));
a2(:, :, 1) = a;

% Initialize second position values
r(:, :, 2) = r(:, :, 1) + deltat.*v(:, :, 1)...
                + (0.5*deltat^2).*a2(:, :, 1);
            
% Begin computations
for n = 2 : nt - 1
    a = nbodyaccn(m, r(:, :, n));
    a2(:, :, n) = a;
    r(:, :, n+1) = 2.*r(:, :, n) - r(:, :, n-1) + (deltat^2).*a2(:, :, n);
end
end


