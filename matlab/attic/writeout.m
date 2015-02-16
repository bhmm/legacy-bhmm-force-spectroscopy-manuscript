%***********************************************************************************************************************************************************
%Script will write out a trajectory 
%************************************************************************************************************************************************************

outid = strcat(' trajectory_',int2str(out_index));

s = strcat('save', outid, ' origin_file opt origin_index numstates stateguess muguess sigmaguess transitionguess allmu allsigma alltransitions allstates data_list_index correct observation_sets observations guess_counts sample_num datalength'); 


eval(s); 

