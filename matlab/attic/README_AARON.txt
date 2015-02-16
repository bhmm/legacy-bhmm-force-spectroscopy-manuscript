Aaron Ewall-Wice
FRET ANALYSIS
Readme


1: list files, trajectory files, and data files. 

This program stores and retrieves data to and from three different files, data files, data lists, and trajectories. 

	The data file is where the actual trajectory is stored and read from. data files are  mat files named in the format data_set_'index'.mat so a valid example would be  data_set_0. So far, the only way to create these data files is with the contrived data generating function makefrets.m or by loading the actual data in matlab, storing it in the matrix "fretdata" and then typing "save data_set_'index' fretdata".

	 
	In order to sample from multiple trajectories, one must employ a list file. List files are named with the format data_list_'index'.m with the index setting each list apart from eachother. For example, the 0th list file would be named data_list_0.m. As i is read as a matlab script, it will be necessary for you to format the list of data files you wish to be sampled from accordingly. List files contain a single matrix called data_list, thus to have the sampler read from the data sets 1, 12, and 5 in a list file with index of 0 for example you will need to enter the following in a text editor and save as 'data_list_0.m'Lists are the files tat are actually read from regardless of how many data sets you wish to sample from. If you only wish to sample from one data set, then simply add only the on one data set you wish to analyze to the data_list. 

	data_list = [1 12 5]; 

	
	Trajectory files are themselves only written to and read by the various functions so you will probably very rarely if ever have to directly interact with them. In brief, they are named with the convention trajectory_'index'.mat and store all information regarding a sampling trajectory including the allstates, alltransitions, stateguess, muguess transitionguess arrays, and sample_num. Essentially, the trajectory files allow one too sample trajectories in multiple steps, perhaps choosing to initialize many unique trajectories from a single stopping point in an old trajectory. Because they are written in binary .mat format, the only way insofar to directly edit one of these trajectoryfiles is to load it into matlab, directly manipulate it's data members, and resave it. 



2: Functions

	There are four functions that divide up the sampling process into 3 to 4 steps. 

	makefrets:
	First there is makefrets which will generate a fake data trajectory and is called as makefrets(data_length, data_index). data_length specifies the number of observations one wishes to put in the trajectory and data_index is the number that the output data file will be indexed as. For example, makefrets(1000,0) will create a 1000 observation fret trajectory in data_set_0. For each data trajectory, makefrets will also publish a file that stores the actual parameters of the markov model used to generate the data stored in actual_parameters_'index'. Model parameters may be modified only by direct manipulation of the .m file contents. 

	initialize: 
	intilialize will read a data_list file and based on data_sets listed will initialize a sampling trajectory. Called as initialize(data_list_index, output_file_index, state count, opt) 

data_list_index indicates the data list you wish to retrieve data from, output_file_index is the index of the trajectory you wish to initialize, state_count is the number of state deffinitions you wish to sample with, opt indicates whether you are working with actual or contrived data(by makefrets), 0 = contrived, any other real indicates actual data. 

	sample:
	sample will samplethe trajectory for a number of times specified in it's parameters. Called as sample(sample_file_index, out_index, num_samples) sample_file_index specifies the index of the trajectory you wish to sample, out_index specifies index of the trajectory you wish to save to and num_samples indicates how many times you'd like to sample the given trajectory. samples add so when you sample a trajectory that has already been sampled 1000 times for say 200 times, sample counts will show up as beginning at 1001 and terminating at 1200. 

	analyze:
	analyze is called simply to read a trajectory and display any results one wishes to specify in the function. At present I've only made a few little plots of mu, sigma, and transition values over samples but we can add a lot more to it. it is called simply by analyze(traj_index) where traj_index is the index of the trajectory you wish to analyze.


These are simply the high level functions that I expect any users to be calling, there are many more and I will attempt to document them (as well as these in more detail) as time progresses.  




 

