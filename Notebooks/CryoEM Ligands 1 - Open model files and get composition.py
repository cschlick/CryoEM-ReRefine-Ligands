#!/usr/bin/env python
# coding: utf-8

from iotbx.map_manager import map_manager as MapManager
from mmtbx.model import manager as ModelManager
from iotbx.data_manager import DataManager
import numpy as np
from iotbx.map_model_manager import map_model_manager as MapModelManager


import os
import pickle
from multiprocessing import Pool
from collections import Counter


errors = []
error_log_path = "data/parsing_errors.txt"

maps_and_models_path = "/net/cci/share/cryoem/maps_and_models/" # this directory is curated EM structures


entry_ids = [entry for entry in os.listdir(maps_and_models_path) if os.path.isdir(maps_and_models_path+entry)]

entries = [] # this is a list of group_args objects

# read the group_args objects in from the .pkl files in the maps_and_models directory
for entry_id in entry_ids:
    entry_path = maps_and_models_path+entry_id
    filenames = os.listdir(entry_path)
    extensions = [filename.split(".")[-1] for filename in filenames]
    extension_dict = Counter(extensions)
    if extension_dict["pkl"] == 1:
        pkl_path = entry_path+"/"+[filename for filename in filenames if filename.split(".")[-1] == "pkl"][0]
        with open(pkl_path,"rb") as fh:
            group_args = pickle.load(fh)
        group_args.add(key="entry",value=entry)
        group_args.add(key="entry_pdb",value=entry[:4])
        group_args.add(key="entry_emdb",value=entry[5:])
        entries.append(group_args)

    else:
        errors.append("FAILED: "+entry_path+": .pkl occurrences != 1")

with open(error_log_path,"w") as fh:
    for error in errors:
        fh.write(error+"\n")
        


def extract_entry_composition(group_args):
    """
    1. takes an entry group_args
    2. reads it in, gets composition
    3. adds composition to group_args
    4. returns group_args
    """
    try:
        dm = DataManager()
        dm.process_model_file(group_args.model_file)
        group_args.add(key="composition",value=dm.get_model().composition())
    except:
        group_args.add(key="composition",value="None")
    return group_args





# Extract composition in parallel
p = Pool(21)
results = p.map(extract_entry_composition,entries)


# dump list of entries (group_arg objects) to disk
with open("data/entry_composition.pkl","wb") as fh:
    pickle.dump(results,fh)


# read it back in
entry_composition_path = "data/entry_composition.pkl" #list of group_args objects
with open(entry_composition_path,"rb") as fh:
    entry_composition_list = pickle.load(fh)


# prune for entries without a composition (todo: why do some not have composition?)
entries_no_composition = []
entries_yes_composition = []
ligand_list = []
for entry in entry_composition_list:
    if entry.composition == "None":
        entries_no_composition.append(entry)
    else:
        entries_yes_composition.append(entry)


# write list of entries to disk
with open("data/entry_composition.pkl","wb") as fh:
    pickle.dump(entries_yes_composition,fh)




