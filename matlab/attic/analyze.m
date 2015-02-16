
%***************************************************************************************************************************
%TAKES THE DATA AND MAKES PRETTY PICTURES AND STUFF
%PROVIDE AS AN ARGUMENT THE TRAJECTORY NUMBER, FOR EXAMPLE, to analyze trajectory 19 type 'analyze(10)' 
%***************************************************************************************************************************


function analyze(traj_index)

opt = 1

%***************************************************************************************************************************
%Load the trajectory to analyze
%***************************************************************************************************************************

	ind = int2str(traj_index); 
	s = strcat('load trajectory_', ind);
	eval(s); 

%***************************************************************************************************************************
%load the actual parameters and initialize (hopefully small) all state matrices if we are working with contrived data
%***************************************************************************************************************************
	
	if opt==0,
    		s = strcat('load actual_parameters_', int2str(origin_index)); 
    		eval(s);
                size(allstates)
		meanstates = mean(int16(allstates));
		stdstates = std(int16(allstates));
	end
%***************************************************************************************************************************
%Create means of mus, sigmas, standard deviations, transitions, and sigmas
%***************************************************************************************************************************


	meanmu = mean(allmu);
	stdmu = std(allmu);
	meansigma = mean(allsigma);
	stdsigma = std(allsigma);
	meantransition = squeeze(mean(alltransitions));
	stdtransition = squeeze(std(alltransitions));
	transgraph = zeros(size(alltransitions,3),size(alltransitions,2),size(alltransitions,1));

%************************************************************************************************************************************************************
	%format all transition matrices to be a readable plot   	
%************************************************************************************************************************************************************


for i = 1:numstates,
    for j = 1:numstates
        for k = 1:sample_num,
               transgraph(k,i,j) = alltransitions(j,i,k); 
        end
    end
end



%************************************************************************************************************************************************************
	%plot transition matrices    	
%************************************************************************************************************************************************************






for i = 1:numstates,
	  figure(i);  plot(transgraph(:,:,i)); 
          xlabel('sample count') ;
          ylabel('transition probability');
          titlename = 'transitprobabilities vs. sample for state';  
          title(titlename); 
          legend('state 1', 'state 2', 'state 3', 'state 4'); 
end



%************************************************************************************************************************************************************
	%plot mus	
%************************************************************************************************************************************************************


 figure(i+1); plot(allmu');
 xlabel('sample number');
 ylabel('mean value');
title('means vs. sample'); 


%************************************************************************************************************************************************************
	%plot observations  	
%************************************************************************************************************************************************************




figure(i+2); plot(observations);
axis([0,160000,0,5])
 legend('state 1', 'state 2', 'state 3', 'state 4'); 


%************************************************************************************************************************************************************
	%plot sigma values  	
%************************************************************************************************************************************************************



figure(i+3); plot(allsigma'); 
 xlabel('sample number'); 
ylabel('sigma value'); 
title('sigma vs. sample'); 
legend('state 1', 'state 2', 'state 3', 'state 4'); 

%************************************************************************************************************************************************************
	%remove dilimeters from observations so that we may display a histogram. 	
%************************************************************************************************************************************************************


histdata = observations; 
histdata(find(histdata == -1000)) = []; 

%************************************************************************************************************************************************************
	%plot histogram of FRET values in observations 	
%************************************************************************************************************************************************************

figure(i+4); hist(histdata, 1000); 
title('FRET spectrum'); 
xlabel('FRET value'); 
ylabel('Counts'); 
%************************************************************************************************************************************************************
	%now create probability density curve using markov model parameters for each state and plot.  	
%************************************************************************************************************************************************************



statesplot = zeros(4,1000);

%create statesplot data
for j = 1:numstates,
	temp = 0.001:.001:1;
	tmean = muguess(j);
	%observationnum;
	tsigma = sigmaguess(j);
	coeff = 1/(tsigma*sqrt(atan(1)*8));
	statesplot(j,:) =coeff*exp(-(power((temp-tmean),2)/(2*power(tsigma,2))));
	statesplot(j,:) = (statesplot(j,:)./sum(statesplot(j,:)))*datalength;

end

figure(i+5); plot(.001:.001:1,statesplot(:,:)'); 
title('Predicted spectrum'); 
legend('state 1','state 2', 'state 3', 'state 4'); 
xlabel('FRET value');
ylabel('Predicted counts'); 
%***********************************************************************************************************************************************************
	%if we are working with contrived data, find spcificity and senstivity of sampling 	
%************************************************************************************************************************************************************

figure(i+6); 
ax = plotyy((.001:.001:1),hist(histdata,1000),.001:.001:1,statesplot(:,:)',@bar,@plot);
axis(ax(1),[0,1,0,4500]);
axis(ax(2),[0,1,0,4500]); 

if opt==0,
    %calculate the total percent correct
    percentcorrect = 100 * (correct/(sample_num*datalength - length(find(observations == -1000))));  


    %calculate the sensitivity and specificity of each state prediction
    sensitivities = zeros(1, numstates); 
    specificities = zeros(1,numstates); 

    for i = 1:numstates,
        %sensitivity = true_positives / (true_positives + false_negatives)
        sensitivities(i) = guess_counts(1,i)/(guess_counts(1,i) + guess_counts(4,i)); 
        %specificity = true_negatives / (true_negatives + false_positives)	
        specificities(i) = guess_counts(3,i)/(guess_counts(3,i) + guess_counts(2,i));
    end

    sensitivities
    specificities 





    
end
%subplot(4,2,i+3); plotyy(1:50,meanstates(1:50),1:50,realstates(1:50));
%xlabel('observation number');
%ylabel('observed state'); 
%title('state vs. observation time'); 
%legend('mean states'); 


