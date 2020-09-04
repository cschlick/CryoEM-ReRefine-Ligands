#!/usr/bin/env python
# coding: utf-8

from iotbx.map_manager import map_manager as MapManager
from mmtbx.model import manager as ModelManager
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager as MapModelManager


import os
import pickle
from multiprocessing import Pool
from collections import Counter
import numpy as np
import gzip
import shutil


# read in nucleotide entries
nucleotide_entries_path = "data/nucleotide_entries.pkl"
working_directory = "../../maps_and_models_ligand_extracted"
error_log = "logs/generate_ligand_density.log"
if os.path.exists(error_log):
    os.remove(error_log)


with open(nucleotide_entries_path,"rb") as fh:
    nucleotide_entries = pickle.load(fh)
    
error_entries = [entry for entry in nucleotide_entries if entry.error != None]
nucleotide_entries = [entry for entry in nucleotide_entries if entry.error == None]

# because we can't deal with multiple ligands, filter for entries with only one type of ligand (and nucleotide)
nucleotide_entries_single_lig = [entry for entry in nucleotide_entries if len(entry.composition._result.other_cnts)==1]
print(str(len(nucleotide_entries_single_lig))+" entries")


def entry_to_ligand_density(entry):
    """
    1. read in entry
    2. find ligands
    3. extract density
    4. write density to disk
    
    """
    

    ligand_code = entry.composition._result.other_cnts.keys()[0]

    # check if we have a folder in the working directory
    entry_path = working_directory+"/"+entry.entry

    # uncomment to force making all entry directories, even if present
    if os.path.exists(entry_path):
        shutil.rmtree(entry_path)

    if not os.path.exists(entry_path): # only process if not already present
        os.mkdir(entry_path)

        # this needs to be fixed, we should have to copy the uncompressed file to open it
        with gzip.open(entry.map_file, 'rb') as fh_in:
            unzip_map_path = entry_path+"/"+entry.map_file.split("/")[-1].strip(".gz")
            with open(unzip_map_path, 'wb') as fh_out:
                shutil.copyfileobj(fh_in, fh_out)


        # All the below is ugly, but it is the only way I was able to get it to work. The map and model is ready in N+1 times where N is the number of ligands
        map_manager = MapManager(unzip_map_path)
        dm = DataManager()
        dm.process_model_file(entry.model_file)
        model_manager = dm.get_model()
        h0 = dm.get_model().get_hierarchy()
        sel0 = h0.atom_selection_cache().selection("not (water or nucleotide or protein)")
        h1 = h0.select(sel0)
        chain_ids = [chain.id for chain in h1.chains()]


        for chain_id in chain_ids:
            map_manager = MapManager(unzip_map_path)
            dm = DataManager()
            dm.process_model_file(entry.model_file)
            model_manager = dm.get_model()
            h0 = dm.get_model().get_hierarchy()
            sel0 = h0.atom_selection_cache().selection("not (water or nucleotide or protein)") # can we directly select the ligands?
            h1 = h0.select(sel0)
            sel1 = h1.atom_selection_cache().selection("chain "+chain_id)
            h2 = h1.select(sel1)
            model_manager_ligand = ModelManager(None,pdb_hierarchy=h2)


            map_model_manager=MapModelManager(map_manager=map_manager, model=model_manager_ligand)
            boxed_mmm = map_model_manager.extract_all_maps_around_model()
            small_map_manager=boxed_mmm.map_manager()
            small_map_manager.write_map(entry_path+"/"+"ligand_"+ligand_code+"_"+chain_id+".map") # save as gz?

            small_model=boxed_mmm.model()
            boxed_mmm.write_model(entry_path+"/"+"ligand_"+ligand_code+"_"+chain_id+".pdb") # write as cif?


from multiprocessing import Pool


# run the density extraction in parallel
p = Pool(21)
results = p.map(entry_to_ligand_density,nucleotide_entries_single_lig)




