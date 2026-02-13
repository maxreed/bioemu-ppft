import numpy as np
import mdtraj

# i can replace these with the parameters for the state calculation.
residue_pairs = [[11,60],[12,31],[64,94]] # these are the three pairs of residues we care about
# recall that state 1 has larger distances (above thresholds), state 2 has smaller distances (below thresholds)
residue_pair_thresholds = np.array([6.5,11.5,13.5]) # giving distances in angstroms to be consistent with bioemu
residue_pair_k = np.array([4,4,4]) # i picked these arbitrarily

# this classifies ras sequences as being in states 1 or 2 (though state classification is continuous from 0 to 1)
def getState(traj: mdtraj.Trajectory, res_pairs=residue_pairs, thresholds=residue_pair_thresholds, k_vals=residue_pair_k):
	CA_indices = traj.topology.select('name CA')
	CA_coords = traj.xyz[:, CA_indices, :] # WT_CA_coords[0] gives all CA coordinates of the first frame
	distances = []
	for frame in CA_coords:
		these_distances =[]
		for residue_pair in residue_pairs:
			this_distance = np.linalg.norm(frame[residue_pair[0]] - frame[residue_pair[1]]) * 10 # * 10 to convert to angstroms
			these_distances.append(this_distance)
		distances.append(these_distances)
	distances = np.array(distances)
	return 1 / (1 + np.exp(np.sum((distances - residue_pair_thresholds) * residue_pair_k,axis=1)))

if __name__=="__main__":
	test_traj = mdtraj.load_xtc("test_RAS_Y40F/samples.xtc", top="test_RAS_Y40F/topology.pdb")
	state_classification = getState(test_traj)
	print("WT:\t" + str(np.sum(state_classification) / len(state_classification)))
