function x = pack_vector2(logKij, logPi)

M = length(logPi);

% Store elements 1..(M-1) of logPi.
x(1:(M-1)) = logPi(1:(M-1));

% Store upper right triangle of logKij.
index = M;
for i = 1:M
  x(index:index+(M-i-1)) = logKij(i,i+1:end);
  index = index + M-i;
end

return
