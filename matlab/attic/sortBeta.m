function [sorted_mu sorted_sigma sorted_transition] = sortBeta(mu, sigma, transition)

%WILL SORT STATES IN ACCORDANCE TO RULES OF THE COMPARESTATES FUNCTION, ASCENDING
%RETURNS SORTED MU, SIGMA, AND TRANSITION MATRIX 

numstates = size(mu,2); 

states = 1:numstates;

sorted_mu = mu; 
sorted_sigma = sigma; 
sorted_transition = transition;


%MERGE SORT STATE DEFFINITIONS
new_indices = state_merge_sort(states, mu, sigma, transition);


%REORDER MU SIGMA AND TRANSITION MATRIX ACCORDING TO NEW STATE DEFFINITIONS 

[sorted_mu sorted_sigma sorted_transition] = order_states(new_indices, mu, sigma, transition);
