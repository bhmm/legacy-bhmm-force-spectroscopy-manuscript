function [tvec, C_t] = compute_autocorrelation(o_t, tmax, options)
% Estimate normalized fluctuation autocorrelation function for observed timeseries data.
%
% [tvec, C_t] = compute_autocorrelation(o_t, tmax, options)
%

% Determine timeseries length.
T = length(o_t);

% Compute fluctuation timeseries.
o_mean = mean(o_t);
do_t = o_t - o_mean;

% Compute unnormalized correlation function.
C_t = zeros(1,tmax+1);
for t = 0:tmax
  C_t(t+1) = mean( do_t(1:T-t) .* do_t(1+t:T) );
end  

% Normalize correlation function.
C_t = C_t / C_t(1);

% Compute times.
tvec = (0:tmax) * options.tau;

return;



