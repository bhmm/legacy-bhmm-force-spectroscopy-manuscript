% Test transition matrix updating scheme.

% Clear all data
clear;

% Define a test count matrix.
%
% This count matrix for the alanine dipeptide six-state decomposition comes from Eq. 27 of:
% Hinrichs NS and Pande VS. Calculation of the distribution of eigenvalues and eigenvectors in Markovian state models for molecular dynamics. JCP 126:244101, 2007.

Nij = [ 
    4380  153   15    2    0    0;
     211 4788    1    0    0    0;
     169    1 4604  226    3    0;
       3   13  158 4823    3    0;
       0    0    0    4 4978   18;
       7    5    0    0   62 4926;
    ]

tau = 0.1; % lag time (ps)

M = size(Nij, 1); % number of states

% Simple estimate of Tij that does not obey detailed balance.
Tij_reference = zeros(size(Nij));
for i = 1:M
  Tij_reference(i,:) = Nij(i,:) / sum(Nij(i,:));
end

% Generate an initial guess.
disp('Computing maximum-likelihood estimate of Tij...');
Tij = generate_transition_matrix_guess(Nij, tau);

% Compare the result.
disp('Tij estimate from counting (without detailed balance):')
Tij_reference
disp('Tij estimate from likelihood maximization (with detailed balance constraint)')
Tij

% Produce samples.
disp('Sampling with Hummer method...');
for i = 1:10
  Tij = hummer_rate_matrix_update(Tij, Nij, tau)
end

