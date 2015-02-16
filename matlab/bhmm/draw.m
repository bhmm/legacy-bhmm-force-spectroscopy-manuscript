function outcome = draw(Pi)
% Draw a random outcome given vector of outcome probabilities.
%
% ARGUMENTS
%   Pi (N vector) - vector of outcome probabilities
%
% RETURNS
%   outcome - integer outcome \in 1..N

if any(Pi < 0)
  Pi  
end

% return random outcome
outcome = randsample(length(Pi), 1, true, Pi);
return

% Compute cumulative probability.
CMF = cumsum(Pi);

% Normalize just in case.
CMF = CMF / CMF(end);

% Draw an outcome.
r = rand;
[n,bin] = histc(r, CMF);
outcome = bin + 1;

return

