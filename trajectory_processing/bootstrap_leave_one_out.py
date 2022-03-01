from enspara.msm import MSM, builders
from enspara import ra
import ctypes
import itertools
import multiprocessing as mp
import numpy as np
import os
from functools import partial

#assignments:
#assignments is a [number of trajectories] x [number of trajectory frames] ragged array 
#of cluster centers to which each frame of each trajectory is assigned
#Each mutant has 7 trajectories and a total of 50,000 frames

#straps:
#straps is a [number of replicates] x [number of trajectories] x [number of trajectory frames] (ragged?) array

#draw n_trials random samples of size len(data) from data with replacement
def bootstrap(data):
	n_trials=len(data)
	straps = [ra.RaggedArray(assignments[np.arange(len(assignments))!=i]) for i in np.arange(n_trials)]      
	return straps

#build msm for one replicate
def eq_pop(assignments, MSM_lag_time, strapid):
	#checks input datatype and reads/converts data to a set type
	if isinstance(assignments, ra.RaggedArray):
		unique_states = np.unique(assignments._data)
	else:
		unique_states = np.unique(assignments.flatten())
		print("flattened")
		#I've never observed this case to occur

	#build msm
	print("building MSM %s" % strapid)
	b = partial(builders.normalize, prior_counts=1/nctrs)#nctrs is used in place of 1/unique_states.shape[0])
	msm_obj = MSM(lag_time=MSM_lag_time, method=b, max_n_states=nctrs) 
        #setting max_n_states makes sure that centers not visited by the trajectories used for a particular round of bootstrapping are still included in the MSM (connected only by prior counts and therefore lacking ergodic sampling) so that all TMPs/EPVs have the same dimensions 
 
	msm_obj.fit(assignments)

	#debugging
	if len(msm_obj.tprobs_)-len(msm_obj.tprobs_[0]) != 0:
		print(len(msm_obj.tprobs_)-len(msm_obj.tprobs_[0]))
	return msm_obj.eq_probs_, msm_obj.tprobs_

#define variables################
mutant_list = ["WT","Nep","Ser1","Ser2"]		#dataset to bootstrap
MSM_lag_time = 250		#lag time (in steps)
data_dir = "/project/bowmanlab/jjmiller/Simulations/PlasmepsinV/all_rmsd/"

#################################


for mutant in mutant_list:
	#load data and determine the total number of cluster centers by looking at the number of unique elements in the whole assignments file
	print("loading assignments")
	assignments = ra.load(data_dir+f"{mutant}_assignments.h5")
	print("assignments loaded")
	nctrs = len(np.unique(assignments.flatten()))
	print("%d cluster centers" % nctrs)

	#perform sampling for bootstrapping
	straps = bootstrap(assignments)

	#calculate equilibrium and transition probabilities for each bootstrapped replicate
	all_pops = [eq_pop(assig, MSM_lag_time, n) for n,assig in enumerate(straps)]
	#compile separate arrays of equilibrium and transition probabilities
	all_eq_pops = [x[0] for x in all_pops]
	all_t_probs = [x[1] for x in all_pops]

	#save equilibrium and transition probabilities
	print('eq_pops and t_probs calculated')
	np.save("/project/bowmanlab/jjmiller/Simulations/PlasmepsinV/all_rmsd/bootstrapping/"+f"strap_output_{mutant}_lagtime{MSM_lag_time}_eqpops.npy", all_eq_pops)
	np.save("/project/bowmanlab/jjmiller/Simulations/PlasmepsinV/all_rmsd/bootstrapping/"+f"strap_output_{mutant}_lagtime{MSM_lag_time}_tprobs.npy", all_t_probs)
