%FUNCTION TAKES A LIST OF INDEX VALUES AND SORTS THEM ACCORDING TO THE RULES IN COMPARE STATES IN ASCENDING ORDER

function sorted_states = state_merge_sort(state_indices, mu, sigma, transition)


if size(state_indices,2) <= 1,

	sorted_states = state_indices; 

else
    len = size(state_indices,2);
    middle = int32(len/2); 

    left = zeros(1,middle); 
    right = zeros(1, len-middle);

    for i = 1:middle,
    	left(i) = state_indices(i); 
    end

    for i = middle+1:len,
    	right(i-middle) = state_indices(i); 
    end

    left = state_merge_sort(left, mu, sigma, transition);
    right = state_merge_sort(right, mu, sigma, transition); 

    sorted_states = state_merge(left, right, mu, sigma, transition); 
end

