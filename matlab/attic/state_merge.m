%MERGES TWO LISTS OF STATE INDICES IN ASCENDING SORTED ORDER BASED ON COMPARE STATES SEE http://en.wikipedia.org/wiki/Merge_sort 


function merged_indices = state_merge(left, right, mu, sigma, transition)


    merged_indices = zeros(1); 
    merged_indices(1,:) = []; 


while size(left,2) > 0 && size(right, 2) > 0,
	
    
	if comparestates(mu(left(1)), sigma(left(1)), transition(left(1),left(1)), mu(right(1)), sigma(right(1)), transition(right(1),right(1)))>-1,
	
		merged_indices = [merged_indices left(1)] ;
		left(:,1) = []; 

	else
		merged_indices = [merged_indices right(1)]; 
		right(:,1) = []; 
	end

end

if size(left,2) > 0,
	merged_indices = [merged_indices, left];
end

if size(right,2) > 0,
	merged_indices = [merged_indices, right]; 
end	

