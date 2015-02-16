% THIS FILE GENERATES SAMPLE OBSERVATIONS ASSUMING THERE ARE N STATES AND AN NxN TRANSITION MATRIX.  EACH STATE IS DEFINED BY A NORMAL DISTRIBUTION WITH MEAN MU AND STANDARD DEVIATION SIGMA

%accepts d_length = length of fret trajectory
%accepts data_index = index of data file to save to


function makefretdata(d_length,data_index)


outid = strcat(' data_set_', int2str(data_index)); 
outactid = strcat(' actual_parameters_', int2str(data_index)); 

% Definitions for mu, sigma, and the transition matrix
numstates = 4;
fretmu = [0.1 0.3 0.6 0.9];
fretsigma = [0.1 0.1 0.1 0.1];
% transitionmatrix = [0.6 0.1 0.2 0.1; 0.1 0.2 0.6 0.1; 0.5 0.1 0.3 0.1; 0.1 0.7 0.1 0.1];
transitionmatrix = [0.98 0.01 0.005 0.005; 0.01 0.975 0.01 0.05; 0.05 0.01 0.975 0.01; 0.05 0.05 0.01 0.98];
datalength = d_length;

% calculate the eigenvalues and eigenvectors
[fretevec, freteval] = eig(transitionmatrix');
[sortedeigs, sortedorder] = sort(diag(freteval), 'descend');

% Calculate the stationary eigenvector and the cumulative probabilities of the stationary eigenvector
stationary = fretevec(:, sortedorder(1));
stationary = stationary ./ sum(stationary);
cumstationary = stationary;
for i = 2:numstates,
  cumstationary(i) = cumstationary(i-1) + cumstationary(i);
end

% Calculate the cumulative transition probabilities
cumprobs = transitionmatrix;
for i = 1:numstates,
  for j = 2:numstates,
    cumprobs(i,j) = cumprobs(i,j-1) + cumprobs(i,j);
  end
end

%fid = fopen(strcat(filename, '.m', 'w');


fretdata = zeros(1, datalength); 


i=1;
% Define the initial position
q = rand();
[n, mybin] = histc(q, cumstationary);
mystate(i) = mybin+1;
% Pick the datapoint from the normal distribution corresponding to the state
mydata(i) = fretsigma(mystate(i)) * randn() + fretmu(mystate(i));

fretdata(i) = mydata(i); 
realstates(i) = mystate(i); 

for i = 2:datalength,
  % Pick the next point according to the transition probabilities
  q = rand();
  [n, mybin] = histc(q, cumprobs(mystate(i-1),:));
  mystate(i) = mybin+1;
  mydata(i) = fretsigma(mystate(i)) * randn() + fretmu(mystate(i));

  fretdata(i) = mydata(i); 
  realstates(i) = mystate(i); 
	
end

%NOW SAVE DATA AND ACTUAL MODEL PARAMETERS TO A DATA FILE

s = strcat('save',  outid, ' fretdata'); 
eval(s);
s = strcat('save', outactid, ' realstates fretmu fretsigma transitionmatrix numstates');
eval(s);
