function states = samplestateassigments(fretdata, mu, sigma, transitionmatrix)

  numstates = size(transitionmatrix,2);
  datalength = length(fretdata); 
  %del_ind = [0 find(fretdata == -1000) length(fretdata)+1]; 
 % numsets = length(del_ind); %THE NUMBER OF DATASETS WE WILL SAMPLESTATEASSIGNEMENTS FOR 
 

 %CUT UP FRETDATA INTO INDIVIDUAL TRAJECTORIES

% for i = 1:numsets,
%	s = strcat('subset_',int2str(i), '=[');
%	for j = del_ind(i)+1:del_ind(i+1)-1
%		s = strcat(s, ',',num2str(fretdata(j)),'%1.20f'); 
%	end
	

  % Calculate the stationary probability distribution of the transitionmatrix
  [fretevec, freteval] = eig(transitionmatrix);
  [sortedeigs, sortedorder] = sort(diag(freteval), 'descend');
  stationary = fretevec(:, sortedorder(1));
  stationary = stationary ./ sum(stationary);

  % Calculate the partial probabilities of observing state i at time t using observations 1...t

  partialprobs = zeros(numstates, datalength);

  for i = 1:numstates,
    %  normalprobability(fretdata(1),mu(i),sigma(i)) * stationary(i)
    partialprobs(i,1) = normalprobability(fretdata(1), mu(i), sigma(i)) * stationary(i);
  end
  % Normalize the partial probabilities at the first step
  partialprobs(:,1) = partialprobs(:,1) ./ sum(partialprobs(:,1));
 %partialprobs(:,1)
  for data = 2:datalength,
    for i = 1:numstates,
      partialprobs(i,data) = normalprobability(fretdata(data), mu(i), sigma(i));
      partsum = 0; 
      for j = 1:numstates,
       % i
        %j
      %  partialprobs(j,data-1)*transitionmatrix(i,j)
        %partialprobs(i,data) = partialprobs(i,data) + transitionmatrix(j,i) * partialprobs(j,data-1);
        partsum = partsum+(partialprobs(j,data-1)*transitionmatrix(j,i));
      end
     % partialprobs(:,data)
      partialprobs(i,data) = partialprobs(i,data)*partsum; 
    end
    % Normalize the partial probabilities at the 'data' step
   partialprobs(:,data) = partialprobs(:,data) ./ sum(partialprobs(:,data));
 %  partialprobs(:,data)
  end
   %partialprobs
  % Sample from state_t through state_1 using these partial probabilities
  cumprob = partialprobs(:,datalength);
  for i = 2:numstates,
    cumprob(i) = cumprob(i) + cumprob(i-1);
  end
  q = rand();
  [n, mybin] = histc(q, cumprob);
  states(datalength) = int8(mybin + 1);
  for data = datalength-1:-1:1,
    for i = 1:numstates,
      newprob(i) = transitionmatrix(i, states(data+1)) * partialprobs(i, data);
    end
    newprob = newprob ./ sum(newprob);
    cumprob = newprob;
    for i = 2:numstates,
      cumprob(i) = cumprob(i) + cumprob(i-1);
    end
    q = rand();
    [n, mybin] = histc(q, cumprob);
    states(data) = int8(mybin + 1);
  end
