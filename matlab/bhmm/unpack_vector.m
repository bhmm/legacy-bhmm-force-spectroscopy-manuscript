function Tij = unpack_vector(x, tau)
% Convert packed form of Hummer's parameterization of rate matrix into transition matrix.

% Determine dimension of packed vector.
Nunique = length(x);

% Determine number of states.
M = round(-(1/2) + sqrt(1/4 + 2*(Nunique+1)));

% Allocate storage for Hummer's parameterization.
logPi = zeros(M,1);
logKij = zeros(M,M);

% Extract logarithms of stationary probabilities.
logPi(1:M-1) = x(1:M-1);
logPi(M) = log(1.0 - sum(exp(logPi(1:M-1))));

% Extract logarithms of unique rate matrix entries.
index = M;
for i = 1:M  
  logKij(i,i+1:end) = x(index:index+(M-i-1));
  index = index + M-i;
end

% Use detailed balance to determine lower triangle of rate matrix entries.
for i = 1:M
  for j = 1:i-1
    logKij(i,j) = logKij(j,i) + logPi(j) - logPi(i);
  end
end

% Compute rate matrix.
Kij = exp(logKij);
for i = 1:M
  Kij(i,i) = 0.0;
  Kij(i,i) = - sum(Kij(i,:));
end

% Check to make sure no state probabilities are too large.
if any(logPi(1:M-1) > 0)
  logPi
  error('Pi > 1');
end  

% Compute transition matrix.
Tij = expm(Kij * tau);

% Check normalization of transition matrix.
if any(sum(Tij,2) > 1.01) || any(sum(Tij,2) < 0.99)
  Kij
  Tij
  error('Tij no longer normalized');
end

return

