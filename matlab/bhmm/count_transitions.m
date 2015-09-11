function Nij = count_transitions(trajectory, nstates)
% Count number of state-to-state transitions observed in a trajectory.
%
% ARGUMENTS
%   trajectory (T array of integers in 1..nstates) - trajectory of states visited (trajectory length T)
%   nstates - number of states
%
% RETURNS
%   Nij (nstate x nstates array of integers) - Nij(i,j) is number of times transition from i -> j was observed for lag time 1
%
% NOTES
%  The algorithm used here scales as O(nstates^2) and could potentially be improved, potentially using find() and histc().

% Allocate storage
Nij = zeros(nstates, nstates, 'int32');

% Count number of transitions observed in this trajectory.
for i = 1:nstates
  for j = 1:nstates
    indices = find((trajectory(1:end-1) == i) & (trajectory(2:end) == j));
    Nij(i,j) = length(indices);
  end
end

return

