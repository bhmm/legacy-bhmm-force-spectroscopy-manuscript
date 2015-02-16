%GIVEN NEW DEFFINITIONS ORDERING OF STATES, RESORTS THE MU, SIGMA, AND TRANSITION MATRIX TO REFLECT. 


function [sorted_mu sorted_sigma sorted_transition] = order_states(new_indices, mu, sigma, transition)



numstates = size(new_indices,2); 

%sort mu values, sigma values, and columns of transition matrix 

for i = 1:numstates,
	new_index = new_indices(i); 
	sorted_mu(new_index) = mu(i); 
	sorted_sigma(new_index) = sigma(i); 
	sorted_transition(:,new_index) = transition(:,i); 
end

%sort rows of transition matrix

temp = sorted_transition; 

for i = 1:numstates,
	new_index = new_indices(i); 
	sorted_transition(new_index,:) = temp(i,:); 
end


