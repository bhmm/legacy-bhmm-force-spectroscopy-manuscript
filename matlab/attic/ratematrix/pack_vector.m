function x = pack_vector(Tij, tau)
% Convert a transition matrix into a packed version of Hummer's parameterization of the rate matrix.
% The transition matrix Tij must be determined from a proper rate matrix such that
%
% Tij = exp(Kij * tau)
%
% ARGUMENTS
%  Tij (MxM) - row-stochastic transition matrix that satsifes Tij = exp(Kij * tau) for rate matrix Kij with positive off-diagonal elements
%  tau - lag time relating Kij and Tij
%
% RETURN VALUES
%  x (M-1 + M(M-1)/2 vector) - packed form of Hummer's parameterization of a rate matrix satisfying detailed balance
%    stationary probabilities Pi(i) = x(i) for i = 1..(M-1)
%    upper triangle of rates K(i,j) are derived from x(M:end)   

% Determine number of states.
M = size(Tij,1);

% Determine number of unique elements for dimension of x.
Nunique = M-1 + M*(M-1)/2; % number of unique elements not determined by detailed balance

% Allocate storage for x.
x = zeros(Nunique, 1); % allocate storage for vector

% Determine rate matrix from Tij.
Kij = logm(Tij) / tau;

% Compute stationary probability.
%Pi = null(Kij');
options = optimset;
options.disp = 0;
[Pi,lambda] = eigs(Tij',1,'LM',options);
Pi = Pi / sum(Pi);

% Compute elementwise logarithms for working parameters.
logPi = log(Pi);
logKij = log(Kij - 2*diag(diag(Kij))); % take care to ensure diagonal of Kij is positive

% Store elements 1..(M-1) of logPi.
x(1:(M-1)) = logPi(1:(M-1));

% Store upper right triangle of logKij.
index = M;
for i = 1:M
  x(index:index+(M-i-1)) = logKij(i,i+1:end);
  index = index + M-i;
end

return
