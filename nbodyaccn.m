function [a] = nbodyaccn(m, r)
% nbodyaccn: Compute accelerations for each particle at a discrete time

% Inputs:
%   m: Vector of length N containing particle masses
%   r: N x 3 array containing particle positions

% Output:
%   a: N x 3 array containing computed particle accelerations


N = length(m);

% Initialize output array of correct size
a = zeros(N, 3);


for i = 1 : N
    for j = 1 : N
        if m(j) ~=0 && i ~= j
            rmag3 = sum((r(j, :) - r(i, :)).^2)^(3/2);
            a(i, :) = a(i, :) + (m(j)/rmag3).*(r(j, :) - r(i, :));
        end
    end
end
end

