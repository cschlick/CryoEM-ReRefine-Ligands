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
from collections import OrderedDict
import operator


import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

plt.style.use("seaborn-notebook")
plot_scale = 2
font_scale = 3

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (16*plot_scale, 11*plot_scale),
         'axes.labelsize': 10*font_scale,
         'axes.titlesize':10*font_scale,
         'xtick.labelsize':10*font_scale,
         'ytick.labelsize':10*font_scale}
plt.rcParams.update(params)


# read in entries

entry_composition_path = "data/entry_composition.pkl" #list of group_args objects
with open(entry_composition_path,"rb") as fh:
    entry_composition_list = pickle.load(fh)


# get a list of all the ligands that appear in the entries
ligand_list = []
for entry in entry_composition_list:

    composition = entry.composition
    ligands = composition._result.other_cnts.keys()
    for ligand in ligands:
        ligand_list.append(ligand)





# Turn it into a Counter() (Kind of like a histogram)
ligand_counter = Counter(ligand_list)
ligand_counter_sorted = OrderedDict(sorted(ligand_counter.items(), key=operator.itemgetter(1),reverse=True))


#plot

labels, values = zip(*ligand_counter_sorted.items()[:20])

indexes = np.arange(len(labels))
width = 1
ret = plt.bar(indexes, values, width)
ret = plt.xticks(indexes + width * 0.1, labels)
plt.ylabel("Number")
plt.xlabel("Ligand Code")


# get entries with ADP, ATP, or GTP
nucleotide_entries = []
for entry in entry_composition_list:
    composition = entry.composition
    ligands = composition._result.other_cnts.keys()
    if ("ATP" in ligands) or ("ADP" in ligands) or ("GTP" in ligands):
        nucleotide_entries.append(entry)


# write nucleotide entries to disk
nucleotide_entries_path = "data/nucleotide_entries.pkl"
with open(nucleotide_entries_path,"wb") as fh:
    pickle.dump(nucleotide_entries,fh)

