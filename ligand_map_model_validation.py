from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager as MapModelManager
from libtbx import group_args
import pickle
import os
from collections import Counter
from mmtbx.maps.correlation import five_cc
from multiprocessing import Pool
import argparse
import logging
import requests

"""
This script contains functions for:
1. extracting ligands from models
2. extracting cryo-em map density from around the ligands
3. calculate CCmask[1], the cross correlation between the ligand and the map

The term "entry" generally refers to a libtbx.group_args object 
containing information related to a pair of PDB/EMDB depositions




[1] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6130467/
"""

def extract_ligand_models(model, desired_hetcodes=[], ignored_hetcodes=["UNK"]):
	"""

	Parameters
	----------
	model	:	mmtbx.model object, container for macromolecular structure
	desired_hetcodes	:	a list of strings, ligand codes to extract. If empty,
											extract all
	ignored_hetcodes	: a list of strings, ligand codes NOT to extract.

	The hetcode strings are compared to the unique_resnames property:
		model.get_hierarchy.residue_groups.unique_resnames


	Returns
	-------
	ligand_entries	: a group_args collection with the following properties

									resname: a string, the hetcode
									model: an mmtbx.model object containing the ligand only

	"""

	# The below is necessary for cryo-em models without crystal symmetry
	# Maybe move this to model.select() later?
	if model.crystal_symmetry() is None:
		from cctbx.maptbx.box import shift_and_box_model
		model = shift_and_box_model(model)

	sel_str_other = "not (water or nucleotide or protein)"  # "other_cnts" in model.composition()
	other_model = model.select(model.selection(sel_str_other))
	ligand_entries = []
	for rg in other_model.get_hierarchy().residue_groups():
		unique_resnames = list(rg.unique_resnames())
		assert(len(unique_resnames) == 1)
		resname = unique_resnames[0]
		if resname not in ignored_hetcodes:
			if (len(desired_hetcodes) == 0) or resname in desired_hetcodes:
				atoms = rg.atoms()
				atoms_sel = atoms.extract_i_seq()
				ligand_model = other_model.select(atoms_sel)
				ligand_id = resname+"_"+rg.parent().id
				group = group_args(resname=resname, model=ligand_model,
													 ligand_id=ligand_id)
				ligand_entries.append(group)

	return ligand_entries


def ligands_map_model_cc(entry):
	"""

	Parameters
	----------
	entry	: A libtbx.group_args object.
					Must contain the following properties:
					 	entry.map_file
						entry.model_file
						entry.resolution

					May contain folling properties:
						entry.output_directory
						entry.has_ligands

	Returns
	-------
	entry	:	The same object as the input, with new data attached.
					New data:
						entry.ligands	# a list of group_args objects for each ligand found
						entry.ligands[0].five_cc_obj	# the cross correlation information
	"""
	model_path = entry.model_file
	map_path = entry.map_file


	dm = DataManager()
	dm.process_model_file(model_path)
	dm.process_real_map_file(map_path)

	mm = dm.get_real_map()
	model = dm.get_model()
	composition = model.composition()
	entry.add(key="composition", value=composition)
	entry.add(key="has_ligands",value=(len(composition._result.other_cnts)>0))

	if entry.has_ligands:
		if hasattr(entry, "output_directory"):
			entry_output_path = os.path.join(entry.output_directory, entry.entry_id)
			if not os.path.exists(entry_output_path):
				os.mkdir(entry_output_path)
			entry.add(key="entry_output_directory",value=entry_output_path)

		ligand_model_entries = extract_ligand_models(model)

		for ligand_entry in ligand_model_entries:

			map_model_manager = MapModelManager(map_manager=mm.deep_copy(),
																					model=ligand_entry.model,
																					ignore_symmetry_conflicts=True)
			boxed_mmm = map_model_manager.extract_all_maps_around_model()
			ligand_mm = boxed_mmm.map_manager()
			ligand_model = boxed_mmm.model()

			five_cc_obj = five_cc(
				ligand_mm.map_data(),
				ligand_model.get_xray_structure(),
				entry.resolution,
				box=None,
				keep_map_calc=False,
				compute_cc_box=False,
				compute_cc_image=False,
				compute_cc_mask=True,
				compute_cc_volume=False,
				compute_cc_peaks=False)

			if hasattr(entry,"output_directory"):

				ligand_model_path = os.path.join(entry.entry_output_directory,
																				 "ligand_" + ligand_entry.ligand_id + ".pdb")
				ligand_map_path = os.path.join(entry.entry_output_directory,
																			 "ligand_" + ligand_entry.ligand_id + ".map")

				ligand_entry.add(key="ligand_model_path", value=ligand_model_path)
				ligand_entry.add(key="ligand_map_path", value=ligand_map_path)

				boxed_mmm.write_map(ligand_map_path)
				boxed_mmm.write_model(ligand_model_path)


			ligand_entry.add(key="five_cc", value=five_cc_obj)
			delattr(ligand_entry,"model") #for multiprocessing, model cannot be pickled

		entry.add(key="ligands", value=ligand_model_entries)
		if hasattr(entry, "output_directory"):
			entry_pkl_path = os.path.join(entry.entry_output_directory,
																		"entry_"+entry.entry_id+".pkl")
			with open(entry_pkl_path, "wb") as fh:
				pickle.dump(entry, fh)

	return entry




def read_data_directory(input_directory):
	"""

	Parameters
	----------
	input_directory	:	The directory of pdb/emdb entries.
										On the cci cluster:/net/cci/share/cryoem/maps_and_models/
										Expects one folder for each entry, and a pickle object in the directory containing a libtbx.group_args object.
										This in turn points to a map file, model file, and contains the resolution.

	nproc	:	The number of entries to process in parallel
	output_directory : The location to (optionally) write out ligand maps and
										 ligand models

	Returns
	-------
	results	: a list of libtbx.group_args objects, including CCmask
										information for each ligand model
	"""

	entry_ids = [entry for entry in os.listdir(input_directory) if
							 os.path.isdir(os.path.join(input_directory, entry))]

	entries = []  # this is a list of group_args objects

	# read the group_args objects in from the .pkl files in the maps_and_models directory
	for entry_id in entry_ids:
		entry_path = os.path.join(input_directory, entry_id)
		filenames = os.listdir(entry_path)
		extensions = [filename.split(".")[-1] for filename in filenames]
		extension_dict = Counter(extensions)
		if extension_dict["pkl"] == 1:
			pkl_path = entry_path + "/" + [filename for filename in filenames if
																		 filename.split(".")[-1] == "pkl"][0]
			with open(pkl_path, "rb") as fh:
				entry = pickle.load(fh)
			entry.add(key="entry_id", value=entry_id)
			entry.add(key="pdb_accession", value=entry_id[:4])
			entry.add(key="emdb_accession", value=entry_id[5:])
			entries.append(entry)

	# filter entries
	logging.info("Found "+str(len(entries))+" entries")
	entries_with_errors = [entry for entry in entries if entry.error is not None]
	logging.info("Dropped " + str(len(entries_with_errors)) + " due to "
																													 "pre-identified "
																													 "errors.")
	entries_without_errors = set(entries).difference(entries_with_errors)
	entries_with_files = [entry for entry in entries_without_errors if
												(entry.map_file is not None) and (
																	entry.model_file is not None)]
	entries_with_files = [entry for entry in entries_with_files if
												os.path.exists(entry.map_file) and os.path.exists(
													entry.model_file)]
	logging.info("Dropped " + str(len(entries_without_errors)-len(
		entries_with_files)) + " entries due to inaccessible files.")

	return entries_with_files



def process_data_directory(input_directory,
											output_directory=None,
											do_one_entry=None,
											overwrite=False,
											nproc=1):

	entries = read_data_directory(input_directory)
	if output_directory is not None:
		for entry in entries:
			entry.add(key="output_directory",value=output_directory)
		if not os.path.exists(output_directory):
			os.mkdir(output_directory)
		if (not overwrite):
			output_contents = os.listdir(output_directory)
			entry_ids = [entry.entry_id for entry in entries]
			entry_ids_to_process = set(entry_ids).difference(output_contents)
			entries = [entry for entry in entries if entry.entry_id in
								 entry_ids_to_process]
	if do_one_entry is not None:
		entries = [entry for entry in entries if entry.entry_id == do_one_entry]

	if nproc ==1:
		results = []
		for entry in entries:
			r = ligands_map_model_cc(entry)
			results.append(r)
	else:
		p = Pool(args.nproc)
		results = p.map(ligands_map_model_cc, entries)
	return results



if __name__ == '__main__':
	logfilepath = 'ligand_map_model_validation.log'
	if os.path.exists(logfilepath):
		pass
	logging.basicConfig(filename=logfilepath,
											level=logging.DEBUG)
	parser = argparse.ArgumentParser()
	parser.add_argument('--input_directory', help='Directory of PDB/EMDB '
																								'entries.')
	parser.add_argument('--nproc', help='The number of entries to process in '
																			'parallel.')
	parser.add_argument('--output_directory',
											help='Directory to write ligand maps and models.')
	parser.add_argument('--overwrite',
											help='If false (default), only process unfinished',
											default=False)
	parser.add_argument('--do_one_entry',
											help="Provide a single entry code.",
											default=None)

	args, extras = parser.parse_known_args()

	args.input_directory = os.path.abspath(args.input_directory)
	args.output_directory = os.path.abspath(args.output_directory)

	try:
		args.nproc = int(args.nproc)
	except:
		args.nproc = 1

	results = process_data_directory(**vars(args))



