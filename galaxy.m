function [t,r] = galaxy(tmax, level, m, r_0, v_0, gv_0)
% galaxy: solves and simple n-body problem of massless stars orbiting
%         a core using an O(delta t^2) FDA.
%
% Inputs:
% % Inputs:
%   tmax:  (real scalar) Final solution time.
%   level: (integer scalar) Discretization level.
%   m:     (real vector) Vector of length Nstars + 1, first element is
%          core.
%   r_0:   (real array) Nstars+1 x 3 array of initial positions. First row
%          is core position (0, 0, 0).
%   v_0:   (real array) Nstars+1 x 3 array of initial velocities. First row
%          is core v0=(0, 0, 0).
%   gv_0:  (real vector) 1 x 4 vector of initial galaxy velocity, 4th elt 
%          is for star rotation about core - if non-zero, rotation is 
%          counter-clockwise when viewed from +z axis, else clockwise.

% Outputs:
%   t:    (real vector) Vector of length nt = 2^level + 1 containing time 
%         mesh.
%   r:    (real vector) Nstars+1 x 3 x nt array containing computed 
%         positions at  discrete times t(n).
% /=======================================================================/


% Define number of time steps and create t, r, v, and a arrays of 
% appropriate size.
N = length(m);
nt = 2^level + 1;
t = linspace(0.0, tmax, nt);
r = zeros(N, 3, nt);
v = zeros(N, 3, nt);
a2 = zeros(N, 3, nt);

% Determine discrete time step from t array
deltat = t(2) - t(1);

% Copy initial values 
r(:, 1, 1) = r_0(:, 1);
r(:, 2, 1) = r_0(:, 2);


% Control star rotation
if gv_0(4) ~= 0
    v(:, 1, 1) = -v_0(:, 1);
    v(:, 2, 1) = -v_0(:, 2);
else
    v(:, 1, 1) = v_0(:, 1);
    v(:, 2, 1) = v_0(:, 2);
end


% Add galactic velocity
v(:, :, 1) = v(:, :, 1) + gv_0(1:3);


% Initialize initial acceleration
a = nbodyaccn(m, r(:, :, 1));
a2(:, :, 1) = a;

% Initialize second position value
r(:, :, 2) = r(:, :, 1) + deltat.*v(:, :, 1)...
                + (0.5*deltat^2).*a2(:, :, 1);



% Begin computations
for n = 2 : nt - 1
    a = nbodyaccn(m, r(:, :, n));
    a2(:, :, n) = a;
    r(:, :, n+1) = 2.*r(:, :, n) - r(:, :, n-1) + (deltat^2).*a2(:, :, n);
end
end

