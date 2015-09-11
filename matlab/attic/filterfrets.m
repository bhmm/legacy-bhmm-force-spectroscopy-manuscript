
%***********************************************************************************************************************************************************
%Function will accept a data list and eliminate ludicrous fret values
%that is all fret data that is less or equal to zero or greater then or equal to one THAT IS NOT A DELIMITER VALUE
%************************************************************************************************************************************************************



%************************************************************************************************************************************************************%Arguments
%datalist: trajectory value
%************************************************************************************************************************************************************

function filtered_data = filterfrets(datalist)


	filtered_data = datalist;

	filtered_data(find(filtered_data >= 1)) = [];

	lower_indices = find(filtered_data<=0); 
	lower_indices(find(filtered_data(lower_indices) == -1000)) = []; 

	filtered_data(lower_indices) = []; 

