function Pi = stationary_probability(Tij)
% Compute stationary probability of a row-stochastic transition matrix.
%
% Pi = stationary_probability(Tij)
%
% ARGUMENTS
%  Tij (N x N matrix) - row-stochastic transition matrix
%
% RETURNS
%  Pi (N vector) - vector of stationary probabilities
%
% TODO
%  * Use a more numerically stable / efficient method here, such as null(Tij' - eye(N))?

% Compute stationary eigenvector.
options = struct();
options.disp = 0;
options.issym = false;
[Pi,lambda] = eigs(Tij',1,'LR',options);

% Normalize to obtain stationary probability.
Pi = Pi / sum(Pi); 

return
