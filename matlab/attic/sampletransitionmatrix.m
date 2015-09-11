
%************************************************************************************************************************************************************
%Arguments:
%counts: counts matrix with each cell corresponding the number of times a transition between each indexed state is observed
%************************************************************************************************************************************************************

function [transitionmatrix] = sampletransitionmatrix(counts)


%************************************************************************************************************************************************************
% PRIOR DISTRIBUTION FOR DIRICHLET DISTRIBUTION?
%************************************************************************************************************************************************************

  numstates = size(counts,2);

% MAXIMUM LIKELIHOOD MATRIX
% numsampled = sum(counts,2);
% normmatrix = numsampled * ones(1,numstates);
% transitionmatrix = counts ./ normmatrix;

  for i = 1:numstates,
    transitionmatrix(:,i) = dirichopt_rnd(counts(i,:)');
  end
  transitionmatrix = transitionmatrix';
