#Author: Hana Jaafari
#Date: 12/5/19
#Purpose: This code is developed by Nick and this author (https://github.com/xqding/DirectCouplingAnalysis_PottsModel_Tensorflow), and is for a "structure-based" DCA code using AWSEM predicted structures.

#HJ: I have made many packages that were being used in many functions into global libraries.
import os
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import Select


def structure_filtered_dca_find_sequence_identifiers_for_structures(pfam_id, structure_strings, pfam_msa_directory=os.path.join(".", "pfam_msas")):

	structure_sequence_dict = {}

	for line in open(os.path.join(pfam_msa_directory, "%s_full.txt" % pfam_id)):

		if not line[0] == "#":

			continue

		for structure_string in structure_strings:

			if structure_string in line:

				structure_sequence_dict[line.split()[1]] = structure_string

	return structure_sequence_dict



def structure_filtered_dca_find_most_aligned_structure(family_sequences, structure_sequence_dict, start_and_end_dict, pdb_directory=os.path.join(".", "pdbs")):

	aligned_residues = []

	max_aligned_residues = 0

	for sequence in family_sequences:

		if sequence.id in list(structure_sequence_dict.keys()):

			num_aligned_residues = np.sum([1 if x.isupper() else 0 for x in str(sequence.seq)])

			aligned_residues.append((sequence.id, num_aligned_residues, structure_sequence_dict[sequence.id]))



	sorted_aligned_residues = sorted(aligned_residues, key=lambda x: x[1], reverse=True) #HJ: This orders the list based on the number of aligned residues

#	print(sorted_aligned_residues)

	for sequence, num_aligned_residues, structure_id in sorted_aligned_residues:

		sequence_start = int(sequence.split('/')[1].split('-')[0]) #HJ: "sequence_start" and "sequence_end" refer to the range of the original sequence that is aligned to the MSA.

		sequence_end = int(sequence.split('/')[1].split('-')[1])

		aligned_structure = structure_sequence_dict[sequence]

		pdb_id = aligned_structure.split()[0]

		chain_id = aligned_structure.split()[1]

		start, end = start_and_end_dict[aligned_structure] #HJ: "start" and "end" refer to the resolved pdb structure of the aligned sequence. d

		structure_filtered_dca_download_pdb(pdb_id, pdb_directory=pdb_directory)

		structure_filtered_dca_truncate_pdb(pdb_id, chain_id, start, end)

		structure = structure_filtered_dca_parse_pdb(pdb_id+chain_id+'_'+str(start)+'-'+str(end)) #HJ: This generates a pdb of the aligned sequence by truncating the pdb file based on sequence's identified chain.

		num_residues = len([residue for residue in structure.get_residues() if not structure_filtered_dca_is_hetero(residue)]); print(sequence_end-sequence_start+1); print(num_residues)

		if num_residues == sequence_end-sequence_start+1: #HJ: if the number of residues in the truncated pdb structure is greater than the number of aligned sequenes.

			most_aligned_sequence = sequence

			most_aligned_structure = aligned_structure

			break



	pdb_id = most_aligned_structure.split()[0]

	chain_id = most_aligned_structure.split()[1]

	start, end = start_and_end_dict[most_aligned_structure]



	return most_aligned_sequence, most_aligned_structure, pdb_id, chain_id, start, end


#HJ: for this code, we provide the associated protein family (via "pfam_id"). I have modified "pfam_msa_directory" in this function's arguments.
def structure_filtered_dca_from_pdb_and_pfam_family(pfam_id, pfam_msa_directory=os.path.join(".", "pfam_msas"), pdb_directory=os.path.join(".", "pdbs"), pdb_fasta_directory=os.path.join(".", "pdb_fastas"), hmm_directory=os.path.join(".", "hmms"), hmm_binary_path="/Users/hanajaafari/Desktop/hmmer-3.2.1/src/hmmer-3.2.1/src", bash_command_string=r'C:\Windows\System32\bash.exe', pdb_alignments_directory=os.path.join(".", "pdb_alignments"), processed_alignments_directory=os.path.join(".", "processed_alignments"), model_directory=os.path.join(".", "models")):
	

	# create directories for storing files

	directory_list = [pfam_msa_directory, pdb_directory, pdb_fasta_directory, hmm_directory, pdb_alignments_directory, processed_alignments_directory, model_directory]

	for directory_name in directory_list:

		if not os.path.exists(directory_name):

			os.makedirs(directory_name)



	# parse pfam pdb file

	structures = structure_filtered_dca_parse_pfam_pdb_file() #HJ: "structures" is a dictionary with the protein families as the keys, with corresponding individual values including PDB ID, chain, PDB residue start, PDB residue end.
    


	# get structures for this family

	family_structures = structures[pfam_id]



	# construct structure strings to search the alignment file

	structure_strings = [pdb_id + ' ' + chain_id for (pdb_id, chain_id, start, end) in family_structures] #HJ: This produces a list of the pdb_ids and chain_id (seperated by space) for every value associated with the chosen protein family.



	# construct start and end dict

	start_and_end_dict = {}

	for (pdb_id, chain_id, start, end) in family_structures:

		start_and_end_dict[pdb_id+' '+chain_id] = (int(start), int(end))



	# download family msa

	structure_filtered_dca_download_pfam_msa(pfam_id, pfam_msa_directory=pfam_msa_directory)



	# read the full sequence alignment

	family_sequences = structure_filtered_dca_read_sto(os.path.join(pfam_msa_directory, "%s_full.txt" % pfam_id)); #HJ: This function reads the MSA of the protein family.



	# find the sequence identifiers that correspond to the structures

	structure_sequence_dict = structure_filtered_dca_find_sequence_identifiers_for_structures(pfam_id, structure_strings, pfam_msa_directory=pfam_msa_directory);print(structure_sequence_dict) #HJ: This function searches through the full alignment and finds the proteins with identified PDB structures. It produces a dictionary with the key being the protein name, and the PDB ID and chain being the corresponding values.



	# find the sequence with a structure that has the most aligned residues

	structure_filtered_dca_find_most_aligned_structure(family_sequences, structure_sequence_dict, start_and_end_dict)
#
#
#
#	# download pdb
#
#	structure_filtered_dca_download_pdb(pdb_id, pdb_directory=pdb_directory)
#
#
#
#	# truncate pdb
#
#	structure_filtered_dca_truncate_pdb(pdb_id, chain_id, start, end)
#
#
#
#	# # extract matrix of pairwise distances between all pairs of residues
#
#	# # extract sequence from pdb
#
#	pairwise_distances, sequence = structure_filtered_dca_get_pairwise_distances_and_sequence(pdb_id, chain_id, start, end, pdb_directory=pdb_directory)
#
#
#
#	# filter msa based on extracted and aligned sequence (only include columns where query sequence is not a gap)
#
#	structure_filtered_dca_process_alignment(pfam_id, most_aligned_sequence)
#
#
#
#	# calculate DCA parameters using the resulting filtered msa and while only allowing nonzero couplings where the pair of residues in the pdb are in contact
#
#	structure_filtered_dca_compute_dca_parameters_sparse(pfam_id, pairwise_distances)
#
#	#structure_filtered_dca_compute_dca_parameters(pfam_id, pairwise_distances)
#
#
#
#	# compare highest coupling norms to contact map
#
#	structure_filtered_dca_compare_norms_to_contact_map(pfam_id, pairwise_distances)



def structure_filtered_dca_compare_norms_to_contact_map(pfam_id, pairwise_distances, couplings_regularization_parameter=0.050, model_directory=os.path.join(".", "models"), contact_threshold=9.5, processed_alignments_directory=os.path.join(".", "processed_alignments"), num_predicted_contacts_to_show=200):


	mpl.rc('font', size = 16)

	mpl.rc('axes', titlesize = 'large', labelsize = 'large')

	mpl.rc('xtick', labelsize = 'large')

	mpl.rc('ytick', labelsize = 'large')

	import sys



	## load model

	with open(os.path.join(model_directory, "model_couplings_regularization_parameter_{:.3f}.pkl".format(couplings_regularization_parameter)), 'rb') as input_file_handle:

		model = pickle.load(input_file_handle)



	len_seq = model['len_seq']

	K = model['K']

	num_node = model['num_node']

	couplings_regularization_parameter = model['couplings_regularization_parameter']

	maxiter = model['max_iter']

	J = model['J']

	h = model['h']



	## calculate interaction scores

	J_prime_dict = {}

	score_FN = np.zeros([len_seq, len_seq])

	for i in range(len_seq):

		for j in range(i+1, len_seq):

			J_prime = J[(i*K):(i*K+K), (j*K):(j*K+K)]

			J_prime = J_prime - J_prime.mean(0).reshape([1,-1]) - J_prime.mean(1).reshape([-1,1]) + J_prime.mean()

			J_prime_dict[(i,j)] = J_prime

			score_FN[i,j] = np.sqrt(np.sum(J_prime * J_prime))

			score_FN[j,i] = score_FN[i,j]

	score_CN = score_FN - score_FN.mean(1).reshape([-1,1]).dot(score_FN.mean(0).reshape([1,-1])) / np.mean(score_FN)

			

	tmp = np.copy(score_CN).reshape([-1])

	tmp.sort()

	cutoff = tmp[-num_predicted_contacts_to_show*2]

	contact_plm = score_CN > cutoff

	for j in range(contact_plm.shape[0]):

		for i in range(j, contact_plm.shape[1]):

			contact_plm[i,j] = False



	contact_pdb = pairwise_distances <= contact_threshold



	for j in range(contact_pdb.shape[0]):

		for i in range(j, contact_pdb.shape[1]):

			contact_pdb[i,j] = False



	##### Plot contacts from both #####

	fig = plt.figure(figsize = (10,10))

	fig.clf()

	I,J = np.where(contact_pdb)

	plt.plot(I,J, 'bo', alpha = 0.2, markersize = 8, label = 'native contacts from PDB')

	plt.axes().set_aspect('equal')

	I,J = np.where(contact_plm)

	plt.plot(I,J, 'r^', markersize = 6, mew = 1.5, label = 'predicted contacts from Potts model')

	plt.title(pfam_id)

	plt.legend()

	output_dir_name = 'output'

	if not os.path.exists(output_dir_name):

		os.makedirs(output_dir_name)

	plt.savefig(os.path.join(output_dir_name, "contact_both.png"))

	plt.show()

	sys.exit()



def structure_filtered_dca_compute_dca_parameters(pfam_id, pairwise_distances, processed_alignments_directory=os.path.join(".", "processed_alignments"), structure_filtering=True, contact_threshold=9.5, couplings_regularization_parameter=0.050, model_directory=os.path.join(".", "models"), max_iter=200):

	import tensorflow as tf

	from sys import exit

	import timeit

	import argparse

	import subprocess




	## read msa

	msa_file_name = os.path.join(processed_alignments_directory, "%s_seq_msa_binary.pkl" % pfam_id)

	with open(msa_file_name, 'rb') as input_file_handle:

		seq_msa_binary = pickle.load(input_file_handle)



	msa_file_name = os.path.join(processed_alignments_directory, "%s_seq_msa.pkl" % pfam_id)

	with open(msa_file_name, 'rb') as input_file_handle:

		seq_msa = pickle.load(input_file_handle)



	weight_file_name = os.path.join(processed_alignments_directory, "%s_seq_weight.pkl" % pfam_id)

	with open(weight_file_name, 'rb') as input_file_handle:

		seq_weight = pickle.load(input_file_handle)



	## pseudolikelihood method for Potts model

	_, len_seq, K = seq_msa_binary.shape

	num_node = len_seq * K

	batch_size = tf.compat.v1.placeholder(tf.int32)

	data = tf.compat.v1.placeholder(tf.float32, shape = [None, num_node])

	data_weight = tf.compat.v1.placeholder(tf.float32, [None])

	half_J = tf.Variable(tf.zeros([num_node, num_node]))

	h = tf.Variable(tf.zeros([num_node]))

	J = half_J + tf.transpose(a=half_J)

	J_mask_value = np.ones((num_node, num_node), dtype = np.float32)

	for i in range(len_seq):

		J_mask_value[K*i:K*i+K, K*i:K*i+K] = 0



	# filter for contacts

	if structure_filtering:

		for i in range(len_seq):

			for j in range(len_seq):

				if not pairwise_distances[i][j] < contact_threshold:

					J_mask_value[K*i:K*i+K, K*j:K*j+K] = 0

			J_mask = tf.constant(J_mask_value)



	J = J * J_mask

	logits = tf.matmul(data, J) + h

	logits = tf.reshape(logits, [-1, K])

	cross_entropy = tf.nn.softmax_cross_entropy_with_logits( logits = logits, labels = tf.stop_gradient( tf.reshape(data, [-1,K])))

	cross_entropy = tf.reduce_sum(input_tensor=tf.reshape(cross_entropy, [-1, len_seq]), axis = 1)

	cross_entropy = tf.reduce_sum(input_tensor=cross_entropy * data_weight)

	couplings_regularization_parameter_holder = tf.compat.v1.placeholder(tf.float32)

	cross_entropy = cross_entropy + couplings_regularization_parameter_holder * tf.reduce_sum(input_tensor=tf.square(J))



	## create a session

	sess = tf.compat.v1.Session()

	init = tf.compat.v1.global_variables_initializer()

	sess.run(init)



	## trainging using L-BFGS-B algorithm

	feed_dict = {data: seq_msa_binary.reshape((-1,num_node)), data_weight: seq_weight, couplings_regularization_parameter_holder: couplings_regularization_parameter}

	print("Initial Cross Entropy: ", sess.run(cross_entropy, feed_dict = feed_dict))

	start_time = timeit.time.time()

	optimizer = tf.contrib.opt.ScipyOptimizerInterface(cross_entropy, var_list = [half_J,h], method = "L-BFGS-B", options={'maxiter': max_iter, 'disp': 1, 'iprint': 2})

	optimizer.minimize(sess, feed_dict = feed_dict)

	end_time = timeit.time.time()

	print("Time elapsed in seconds: ", end_time - start_time)

	print("Final Cross Entropy: ", sess.run(cross_entropy, feed_dict = feed_dict))



	## save J and h

	model = {}

	model['len_seq'] = len_seq

	model['K'] = K

	model['num_node'] = num_node

	model['couplings_regularization_parameter'] = couplings_regularization_parameter

	model['max_iter'] = max_iter

	model['J'] = sess.run(J)

	model['h'] = sess.run(h)



	output_file_name = os.path.join(model_directory, "model_couplings_regularization_parameter_{:.3f}.pkl".format(couplings_regularization_parameter))

	with open(output_file_name, 'wb') as output_file_handle:

		pickle.dump(model, output_file_handle)



def structure_filtered_dca_compute_dca_parameters_sparse(pfam_id, pairwise_distances, processed_alignments_directory=os.path.join(".", "processed_alignments"), structure_filtering=True, contact_threshold=9.5, couplings_regularization_parameter=0.050, model_directory=os.path.join(".", "models"), max_iter=200):

	import tensorflow as tf

	from sys import exit


	import timeit

	import argparse

	import subprocess



	## read msa

	msa_file_name = os.path.join(processed_alignments_directory, "%s_seq_msa_binary.pkl" % pfam_id)

	with open(msa_file_name, 'rb') as input_file_handle:

		seq_msa_binary = pickle.load(input_file_handle)



	msa_file_name = os.path.join(processed_alignments_directory, "%s_seq_msa.pkl" % pfam_id)

	with open(msa_file_name, 'rb') as input_file_handle:

		seq_msa = pickle.load(input_file_handle)



	weight_file_name = os.path.join(processed_alignments_directory, "%s_seq_weight.pkl" % pfam_id)

	with open(weight_file_name, 'rb') as input_file_handle:

		seq_weight = pickle.load(input_file_handle)



	## pseudolikelihood method for Potts model

	_, len_seq, K = seq_msa_binary.shape

	num_node = len_seq * K

	batch_size = tf.compat.v1.placeholder(tf.int32)

	data = tf.compat.v1.placeholder(tf.float32, shape = [None, num_node])

	data_weight = tf.compat.v1.placeholder(tf.float32, [None])

	#half_J = tf.Variable(tf.zeros([num_node, num_node]))

	#J_mask_value = np.ones((num_node, num_node), dtype = np.float32)

	# for i in range(len_seq):

	# 	J_mask_value[K*i:K*i+K, K*i:K*i+K] = 0



	# filter for contacts

	contact_indices = []

	mask_indices = []

	if structure_filtering:

		for i in range(len_seq):

			for j in range(len_seq):

				for k in range(K):

					for l in range(K):

						if pairwise_distances[i][j] < contact_threshold:

							if j > i:

								contact_indices.append([K*i+k,K*j+l])

						else:

							mask_indices.append([K*i+k,K*j+l])

	

	J_mask = tf.constant(tf.SparseTensor(indices=mask_indices, values=[0]*len(mask_indices), dense_shape=[num_node, num_node]))

	half_J = tf.Variable(tf.SparseTensor(indices=contact_indices, values=[0]*len(contact_indices), dense_shape=[num_node, num_node]))

	h = tf.Variable(tf.zeros([num_node]))

	J = half_J + tf.transpose(a=half_J)



	J = J * J_mask

	logits = tf.matmul(data, J) + h

	logits = tf.reshape(logits, [-1, K])

	cross_entropy = tf.nn.softmax_cross_entropy_with_logits( logits = logits, labels = tf.stop_gradient( tf.reshape(data, [-1,K])))

	cross_entropy = tf.reduce_sum(input_tensor=tf.reshape(cross_entropy, [-1, len_seq]), axis = 1)

	cross_entropy = tf.reduce_sum(input_tensor=cross_entropy * data_weight)

	couplings_regularization_parameter_holder = tf.compat.v1.placeholder(tf.float32)

	cross_entropy = cross_entropy + couplings_regularization_parameter_holder * tf.reduce_sum(input_tensor=tf.square(J))



	## create a session

	sess = tf.compat.v1.Session()

	init = tf.compat.v1.global_variables_initializer()

	sess.run(init)



	## trainging using L-BFGS-B algorithm

	feed_dict = {data: seq_msa_binary.reshape((-1,num_node)), data_weight: seq_weight, couplings_regularization_parameter_holder: couplings_regularization_parameter}

	print("Initial Cross Entropy: ", sess.run(cross_entropy, feed_dict = feed_dict))

	start_time = timeit.time.time()

	optimizer = tf.contrib.opt.ScipyOptimizerInterface(cross_entropy, var_list = [half_J,h], method = "L-BFGS-B", options={'maxiter': max_iter, 'disp': 1, 'iprint': 2})

	optimizer.minimize(sess, feed_dict = feed_dict)

	end_time = timeit.time.time()

	print("Time elapsed in seconds: ", end_time - start_time)

	print("Final Cross Entropy: ", sess.run(cross_entropy, feed_dict = feed_dict))



	## save J and h

	model = {}

	model['len_seq'] = len_seq

	model['K'] = K

	model['num_node'] = num_node

	model['couplings_regularization_parameter'] = couplings_regularization_parameter

	model['max_iter'] = max_iter

	model['J'] = sess.run(J)

	model['h'] = sess.run(h)



	output_file_name = os.path.join(model_directory, "model_couplings_regularization_parameter_{:.3f}.pkl".format(couplings_regularization_parameter))

	with open(output_file_name, 'wb') as output_file_handle:

		pickle.dump(model, output_file_handle)



def structure_filtered_dca_truncate_pdb(pdb_id, chain_id, start, end, pdb_directory=os.path.join(".", "pdbs")):

	from Bio.PDB import PDBIO

	io = PDBIO()

	structure = structure_filtered_dca_parse_pdb(pdb_id)

	structure_selector = structure_filtered_dca_pfam_is_in_structure_selection(chain_id, start, end)

	io.set_structure(structure)

	io.save(os.path.join(pdb_directory, "%s%s_%d-%d.pdb" % (pdb_id, chain_id, start, end)), select=structure_selector)



def structure_filtered_dca_download_pfam_msa(pfam_id, pfam_msa_directory=os.path.join(".", "pfam_msas")):

	import urllib3

	import gzip

	import os



	print("Downloading the full multiple sequence alignment for Pfam: {0} ......".format(pfam_id))

	if not os.path.exists(os.path.join(pfam_msa_directory, "{0}_full.txt".format(pfam_id))):

		http = urllib3.PoolManager()

		r = http.request('GET', 'http://pfam.xfam.org/family/{0}/alignment/full/gzipped'.format(pfam_id))

		data = gzip.decompress(r.data)

		data = data.decode()

		full_filename = os.path.join(pfam_msa_directory, "{0}_full.txt".format(pfam_id))

		with open(full_filename, 'w') as file_handle:

			print(data, file = file_handle)



	print("Downloading the seed multiple sequence alignment for Pfam: {0} ......".format(pfam_id))

	if not os.path.exists(os.path.join(pfam_msa_directory, "{0}_seed.txt".format(pfam_id))):

		http = urllib3.PoolManager()

		r = http.request('GET', 'http://pfam.xfam.org/family/{0}/alignment/seed/gzipped'.format(pfam_id))

		data = gzip.decompress(r.data)

		data = data.decode()

		full_filename = os.path.join(pfam_msa_directory, "{0}_seed.txt".format(pfam_id))

		with open(full_filename, 'w') as file_handle:

			print(data, file = file_handle)



def structure_filtered_dca_download_pdb(pdb_id, pdb_directory=os.path.join(".", "pdbs")):

	from Bio.PDB import PDBList

	if not os.path.isfile(os.path.join(pdb_directory, pdb_id+'.pdb')):

		pdbl = PDBList()

		pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_directory, obsolete=False, file_format="pdb")

		os.rename(os.path.join(pdb_directory, 'pdb'+pdb_id+'.ent'), os.path.join(pdb_directory, pdb_id+'.pdb'))



def structure_filtered_dca_get_interaction_atom(residue):

	try:

		if residue.resname == "GLY":

			return residue['CA']

		else:

			return residue['CB']

	except:

		raise



def structure_filtered_dca_get_interaction_distance(res1, res2):

	return structure_filtered_dca_get_interaction_atom(res1) - structure_filtered_dca_get_interaction_atom(res2)



def structure_filtered_dca_is_hetero(residue):

	if residue.id[0] != ' ':

		return True

	else:

		return False



def structure_filtered_dca_get_pairwise_distances_and_sequence(pdb_id, chain_id, start, end, pdb_directory='./pdbs/'):



	parser = PDBParser()

	structure = parser.get_structure(pdb_id, os.path.join(pdb_directory, "%s%s_%s-%s.pdb" % (pdb_id, chain_id, str(start), str(end))))

	structure = structure[0][chain_id]

	sequence = structure_filtered_dca_get_sequence_from_structure(structure)

	res_list = Selection.unfold_entities(structure, 'R')

	res_list = [residue for residue in res_list if not structure_filtered_dca_is_hetero(residue)]

	pairwise_distances = []

	for i, res1 in enumerate(res_list):

		pairwise_distances.append([])

		for j, res2 in enumerate(res_list):

			pairwise_distances[-1].append(structure_filtered_dca_get_interaction_distance(res1, res2))



	return np.array(pairwise_distances), sequence



def structure_filtered_dca_get_sequence_from_structure(structure):

	from Bio.PDB import PPBuilder

	sequence = ""

	ppb=PPBuilder(radius=10.0)

	for pp in ppb.build_peptides(structure, aa_only=False):

		sequence += '%s\n' % pp.get_sequence()

	return sequence.replace('\n', '')



def structure_filtered_dca_process_alignment(pfam_id, query_seq_id, pfam_msa_directory=os.path.join(".", "pfam_msas"), pdb_alignments_directory=os.path.join(".", "pdb_alignments"), maximum_percent_gaps=0.20, maximum_percent_sequences_with_gaps=0.20, processed_alignments_directory=os.path.join(".", "processed_alignments")):



	## read all the sequences into a dictionary

	seq_dict = {}

	alignment_file_name = os.path.join(pfam_msa_directory, pfam_id+"_full.txt")

	with open(alignment_file_name, 'r') as file_handle:

		for line in file_handle:

			line = line.strip()        

			if line == "" or line[0] == "#" or line[0] == "/" or line[0] == "":

				continue

			seq_id, seq = line.split()

			seq_dict[seq_id] = seq.upper()



	## remove gaps in the query sequence from the rest of the sequences

	query_seq = seq_dict[query_seq_id] # with gaps

	idx = [ s == "-" or s == "." for s in query_seq]

	for k in seq_dict.keys():

		seq_dict[k] = [seq_dict[k][i] for i in range(len(seq_dict[k])) if idx[i] == False]

	query_seq = seq_dict[query_seq_id] # without gaps



	## remove sequences that have too many gaps with respect to the query sequence

	len_query_seq = len(query_seq)

	seq_id = list(seq_dict.keys())

	for k in seq_id:

		if seq_dict[k].count("-") + seq_dict[k].count(".") >= len_query_seq * maximum_percent_gaps:

			seq_dict.pop(k)

	

	## convert aa type into num 0-20

	aa = ['R', 'H', 'K',

		'D', 'E',

		'S', 'T', 'N', 'Q',

		'C', 'G', 'P',

		'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']

	aa_index = {}

	aa_index['-'] = 0

	aa_index['.'] = 0



	i = 1

	for a in aa:

		aa_index[a] = i

		i += 1

	with open(os.path.join(processed_alignments_directory, "%s_aa_index.pkl" % pfam_id), 'wb') as file_handle:

		pickle.dump(aa_index, file_handle)



	seq_msa = []

	for k in seq_dict.keys():

		try:

			seq_msa.append([aa_index[s] for s in seq_dict[k]])

		except:

			continue

	seq_msa = np.array(seq_msa)



	with open(os.path.join(processed_alignments_directory, "%s_seq_msa.pkl" % pfam_id), 'wb') as file_handle:

		pickle.dump(seq_msa, file_handle)



	## reweighting sequences

	seq_weight = np.zeros(seq_msa.shape)

	for j in range(seq_msa.shape[1]):

		aa_type, aa_counts = np.unique(seq_msa[:,j], return_counts = True)

		num_type = len(aa_type)

		aa_dict = {}

		for a in aa_type:

			aa_dict[a] = aa_counts[list(aa_type).index(a)]

		for i in range(seq_msa.shape[0]):

			seq_weight[i,j] = (1.0/num_type) * (1.0/aa_dict[seq_msa[i,j]])

	tot_weight = np.sum(seq_weight)

	seq_weight = seq_weight.sum(1) / tot_weight 

	with open(os.path.join(processed_alignments_directory, "%s_seq_weight.pkl" % pfam_id), 'wb') as file_handle:

		pickle.dump(seq_weight, file_handle)



	## change aa numbering into binary

	K = 21 ## num of classes of aa

	D = np.identity(K)

	num_seq = seq_msa.shape[0]

	len_seq_msa = seq_msa.shape[1]

	seq_msa_binary = np.zeros((num_seq, len_seq_msa, K))

	for i in range(num_seq):

		seq_msa_binary[i,:,:] = D[seq_msa[i]]



	with open(os.path.join(processed_alignments_directory, "%s_seq_msa_binary.pkl" % pfam_id), 'wb') as file_handle:

		pickle.dump(seq_msa_binary, file_handle)



def structure_filtered_dca_read_sto(input_file_name):

	from Bio import AlignIO



	input_handle = open(input_file_name, "rU")

	alignments = AlignIO.read(input_handle, "stockholm")

	input_handle.close()



	return alignments



def structure_filtered_dca_parse_pfam_pdb_file(pdb_file="pdb_pfamA_reg.txt"):

	structures = {}

	for line in open(pdb_file, 'r'):

		auto_pdb_reg, auto_uniprot_reg_full, pdb_id, pfamA_acc, pfamseq_acc, chain, pdb_res_start, pdb_start_icode, pdb_res_end, pdb_end_icode, seq_start, seq_end, hex_colour = line.split('\t') #HJ: This loop is reading every line in the file, and seperating the values by tab.

		if not pfamA_acc in structures.keys(): 

			structures[pfamA_acc] = []

		structures[pfamA_acc].append((pdb_id, chain, pdb_res_start, pdb_res_end)) #HJ: A dictionary is created for each protein family. For every protein family, this loop progressively adds any corresponding structures (pdb_id) with their residue length.

	return structures



def structure_filtered_dca_parse_pdb(pdb_id, pdb_directory=os.path.join(".", "pdbs")):

	parser = PDBParser()

	return parser.get_structure(pdb_id, os.path.join(pdb_directory, "%s.pdb" % pdb_id))



class structure_filtered_dca_pfam_is_in_structure_selection(Select):

	def __init__(self, selected_chain, selected_start, selected_end):

		self.selected_chain = selected_chain

		self.selected_start = selected_start

		self.selected_end = selected_end

	def accept_atom(self, atom):

		accept = True

		pdb_id, model, chain, (het_flag, resid, _), (atomid, _) = atom.get_full_id()

		if model != 0 or chain != self.selected_chain or resid not in range(self.selected_start, self.selected_end+1) or het_flag != " ":

			accept = False

		return accept
    
structure_filtered_dca_from_pdb_and_pfam_family("PF05400")