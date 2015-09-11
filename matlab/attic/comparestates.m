%FUNCTION WILL COMPARE TWO  STATES AND RETURN -1 IF THE FIRST (STATE 0) IS GREATER THEN THE SECOND (STATE 1) RETURN 1 IF THE SECOND STATE (STATE 1) IS GREATER THEN THE FIRST (STATE 0)
%RETURN 0 IF ALL COMPARABLE VALUES ARE WITHIN THRESHOLDS

function truth_value = comparestates(mu_0, sigma_0, self_transition_0, mu_1, sigma_1, self_transition_1) 

		      truth_value = 1; 

		      EPSILON = .08;  %Threshold value for comapring the difference in mu
		      EPSILON0 = .005; %Threshold value for comparing the differences in self transitions 



		      delta = abs(mu_0 - mu_1);
		      delta0 = abs(self_transition_0 - self_transition_1);


		      %IF THE TWO STATES ARE PROPERLY SEPERATED BY THEIR MU VALUES THEN COMPARE THEM AS SUCH
		      
		      if delta >= EPSILON,
		      	 
			 if mu_0 >= mu_1,

			 	truth_value = -1; 

			 end


		        %SHOULD THE MUS BE PACKED TOGETHER WE WILL COMPARE THE STATES BY THEIR STABILITY, THAT IS THEIR SELF TRANSITION VALUES. THE MORE STABLE, THE GREATER  

		    

			 elseif delta0 >= EPSILON0,
			
				if self_transition_0 >= self_transition_1,
				
					truth_value = -1; 

				end

			%SHOULD THE MUS AND SELF TRANSITIONS BE PACKED, WE FINALLY COMPARE SIGMA IN ASCENDING ORDER 

              elseif abs(sigma_0-sigma_1) >= EPSILON0,

				if sigma_0 >= sigma_1,
					
					truth_value = -1; 

				end

                else

				truth_value = 0; 
			
			end






