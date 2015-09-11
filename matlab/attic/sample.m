%************************************************************************************************************************************************************
%sample.m will run sampling on the specified data trajectory for a specified number of steps
%parameters are as follows
%sample_file_index: the index number of a trajectory file that you wish to sample from
%out_index: the index of the trajectory you wish to write results to
%num_samples: the number of samples you wish to take
%************************************************************************************************************************************************************


function sample(sample_file_index,out_index, num_samples,ref)

%************************************************************************************************************************************************************
%refreshlength specifies the interval of samples that program will write out
%************************************************************************************************************************************************************


refreshlength = num_samples; 


%************************************************************************************************************************************************************
%load trajectory to sample
%************************************************************************************************************************************************************

id = int2str(sample_file_index); 
inid = strcat(' trajectory_', id); 
outid = strcat(' trajectory_', int2str(out_index)); 
s = strcat('load', inid); 
eval(s);

%************************************************************************************************************************************************************
%allmunew,allsigmanew,allstatesnew,alltransitionsnew are temporary history arrays that are periodically appended to allmu,allsigma,allstates,alltransitions
%These histories are as long as the refreshlength
%************************************************************************************************************************************************************

allmunew = zeros(numstates,refreshlength); 
allsigmanew = zeros(numstates,refreshlength); 
alltransitionsnew = zeros(numstates,numstates,refreshlength); 



%************************************************************************************************************************************************************
%if we are using contrived data then load the actual model parameters
%************************************************************************************************************************************************************

if opt == 0, 
	actinid = strcat( ' actual_parameters_', int2str(origin_index));

	s = strcat('load', actinid); 

	eval(s); 

	allstatesnew = int8(zeros(datalength,refreshlength));
end

%************************************************************************************************************************************************************
%SAMPLE FOR THE NUMBER OF SAMPLES PROVIDED BY SAMPLE_NUM
%************************************************************************************************************************************************************


for samples = (sample_num+1):(num_samples+sample_num),
        samples

%************************************************************************************************************************************************************
% Generate sample of mu and sigma from state assignments and data
%************************************************************************************************************************************************************

	[muguess, sigmaguess] = samplemusigma(observations, stateguess,numstates);

%************************************************************************************************************************************************************
% Generate counts of transitions between states
% Initially with psuedocounts of .25
%************************************************************************************************************************************************************
  
  
  	counts = ones(numstates,numstates)./4; 

%************************************************************************************************************************************************************
%COUNT TRANSITIONS BETWEEN ALL OBSERVATIONS NOT INCLUDING ZERO DELIMETERS (SEE INITIALIZE.M)    	
%************************************************************************************************************************************************************

 
  
  for i = 2:datalength,
  	if stateguess(i) ~= 0 && stateguess(i-1) ~= 0,	
		counts(stateguess(i), stateguess(i-1)) = counts(stateguess(i), stateguess(i-1)) + 1;
  	end 
  end
  

%************************************************************************************************************************************************************
%Generate a Dirichlet sampled transition matrix from counts
%************************************************************************************************************************************************************
 
 transitionguess = sampletransitionmatrix(counts);

%************************************************************************************************************************************************************
%Sort states with mergesort
%************************************************************************************************************************************************************
 
 [muguess, sigmaguess, transitionguess] = sortBeta(muguess, sigmaguess, transitionguess);

%************************************************************************************************************************************************************
% Generate sample of state assignments from mu, sigma, fretdata, transitionmatri
%************************************************************************************************************************************************************
  
  if ~isempty(observation_sets(1).data),
  	stateguess = samplestateassignments(observation_sets(1).data, muguess, sigmaguess, transitionguess);
  end

  for i = 2:length(observation_sets),
  	if ~isempty(observation_sets(i).data),
  		stateguess = [stateguess 0  samplestateassignments(observation_sets(i).data, muguess, sigmaguess, transitionguess)];
  	else
		stateguess = [stateguess 0]; 
  	end
  end
  
  stateguess = int8(stateguess); 


  curind = mod(samples,refreshlength);  
  
  if curind == 0
  	curind = refreshlength; 
  end
  
  allmunew(:,curind) = muguess;
  allsigmanew(:,curind) = sigmaguess;
  alltransitionsnew(:,:,curind) = transitionguess;

 if opt == 0,
 	allstatesnew(:,curind) = stateguess;   
 end

%************************************************************************************************************************************************************
%if we are using contrived data
%get counts of correct and incorrect guesses for each state
%************************************************************************************************************************************************************ 
 if opt==0,
  for i = 1:numstates,
      for j = 1:size(stateguess,2),
	if stateguess(j) ~= 0, 
          if isequal(stateguess(j), i), %count positives
          	if isequal(realstates(j),i), %count true
          		guess_counts(1,i) = guess_counts(1,i) + 1; 
          	else
			correct = correct+1; 
          		guess_counts(2,i) = guess_counts(2,i) + 1; %count false
          	end
         elseif ~isequal(stateguess(j), i)
         	if ~isequal(realstates(j), i),
         		guess_counts(3,i) = guess_counts(3,i) + 1; %count true
         	else
         		guess_counts(4,i) = guess_counts(4,i) + 1; %count false
         	end
             end          		
      	   end
         end
       end

    for i = 1:datalength,
    	if isequal(realstates(i), stateguess(i)),
        	correct = correct+1;
    	end
    end
 end



%************************************************************************************************************************************************************
%If we are at the refreshlength 
%then we append our history matrics and write out
%************************************************************************************************************************************************************

if curind==refreshlength,
  
  start = samples-curind+1; 
  allmu(:,start:samples) = allmunew; 
  allsigma(:,start:samples) = allsigmanew;
  	alltransitions(:,:,start:samples) = alltransitionsnew;
  
  if opt == 0,
  	allstates(:,start:samples) = allstatesnew; 
  end
  
  sample_num = samples; 
  
  writeout; 
 
end


end



%************************************************************************************************************************************************************
%Final append at write out
%************************************************************************************************************************************************************


  start = samples-curind+1; 
  allmu(:,start:samples) = allmunew; 
  allsigma(:,start:samples) = allsigmanew;
  	alltransitions(:,:,start:samples) = alltransitionsnew; 
  
  if opt == 0,
  	allstates(:,start:samples) = allstatesnew; 
  end

  sample_num = samples; 
    writeout; 
