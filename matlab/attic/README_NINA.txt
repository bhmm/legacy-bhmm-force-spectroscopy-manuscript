matlab files:
makefrets.m (THIS IS THE FILE TO CHANGE TO CHANGE FAKE DATA FRET PARAMETERS)
 - creates actual_parameters_"index".mat
 - creates data_set_"index".mat
 - you need to manually create data_list_"index".mat, the "index" numbers don't have to match, this somehow allows you to combine data input from multiple datasets, but not necessary right now 

initialize.m
% - requires filterfrets.m (filters trajectories and deletes values below 0, above 1, commented out for now)
 - requires writeout.m (writes out trajectory files into .mat file)
 - creates trajectory_"index".mat

sample.m
 - requires samplemusigma.m (samples mu and sigma)
 - requires sampletransitionmatrix.m (samples transition matrices)
   - requires dirichopt_rnd.m
 - requires sortBeta.m (sorts states by Aaron's algorithm)
   - requires state_merge_sort.m (???)
     - requires state_merge.m (???)
       - requires compare_states.m (???)
   - requires order_states.m (???)
 - requires samplestateassignments.m (this is the viterbi like algorithm)
   - requires normalprobability.m (calculates probabilities from a normal distribution)
 - updates trajectory_"index".mat

analyze.m
 I can't get this to quite run at the moment...


run:
> makefrets(data_length, data_index)
% data_length is how many observations, data_index is number for output files
% I ran "makefrets(500,2)"

> create a text file data_list_"data_index".m which contains the lines
 data_list = ["data_index"];
 origin_file = ' testfile';

> initialize(data_index, output_index, states, opt)
% data_index is same as above, output_index is output trajectory index, states is number of states in HMM, opt is 0 (don't really understand)
% I ran "initialize(2, 2, 4, 0)"

> sample(data_index, output_index, numsamples)
% data_index and output_index are same as above, num_samples is number of times to iterate sampling the stateassignments
% I ran "sample(2, 2, 1000)"


