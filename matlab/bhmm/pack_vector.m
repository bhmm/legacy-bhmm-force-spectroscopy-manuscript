function x = pack_vector(Tij, tau)
% Convert a transition matrix into a packed version of Hummer's parameterization of the rate matrix.
% The transition matrix Tij must be determined from a proper rate matrix such that
%
% Tij = exp(Kij * tau)
%
% ARGUMENTS
%  Tij (MxM) - row-stochastic transition matrix that satsifes Tij = exp(Kij * tau) for rate matrix Kij with positive off-diagonal elements
%    Note that Tij must come from the exponentiation of a rate matrix that satisfies detailed balance, or an error will be thrown.
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

% Compute stationary probability from Tij.
Pi = stationary_probability(Tij);

% Check to make sure Tij satisfies detailed balance.
D = diag(Pi);
TOLERANCE = 1.0e-3;
Tij_symmetrized = D^(+1/2)*Tij*D^(-1/2);
if norm(Tij_symmetrized - Tij_symmetrized') > TOLERANCE
  Tij
  Pi
  Tij_symmetrized
  disp('norm');
  norm(Tij_symmetrized - Tij_symmetrized')
  error('Tij does not satisfy detailed balance');
end

% Ensure all eigenvalues are positive
smallest_eigenvalue = eigs(Tij, 1, 'SR', struct('disp', 0));   
if (smallest_eigenvalue <= eps * M)
  smallest_eigenvalue  
  error('smallest eigenvalue is smaller than allowed tolerance to ensure real packed x');
end
    


% Determine rate matrix from Tij.
Kij = logm(Tij) / tau;

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

% DEBUG make sure x is strictly real
if ~isreal(x)
  Tij
  Kij
  logPi
  logKij
  x
  smallest_eigenvalue
  error('x is not strictly real');
end

return
