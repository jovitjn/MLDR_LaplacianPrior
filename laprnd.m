function r = laprnd(m, n, mu, b)
% LAPRND Generate laplacian random numbers.
%   R = LAPRND(M, N, MU, B) returns an M-by-N matrix of
%   Laplacian distributed random numbers with mean MU and scale
%   parameter B.
%
%   Default values for MU and B are 0 and 1, respectively.
%
%   Note: The Laplacian distribution is also known as the
%   double exponential distribution.

if nargin < 3
    mu = 0;
end
if nargin < 4
    b = 1;
end

U = rand(m, n) - 0.5;
r = mu - b * sign(U) .* log(1 - 2*abs(U));
