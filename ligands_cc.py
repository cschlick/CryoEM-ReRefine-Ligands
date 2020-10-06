from __future__ import division, print_function
import os
import sys
import pickle

from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager as MapModelManager
from iotbx.pdb import common_residue_names_get_class
from libtbx import group_args
from iotbx.file_reader import splitext


"""
Summary
-------

This module contains functions to:
1. extract ligand models from an mmtbx.model object
2. extract density from around those models and calculate CCmask [1]



[1] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6130467/

"""


class LigandSelection:

  default_exclude_resnames = {"UNK"}
  default_exclude_residue_classes = {"common_amino_acid",
                                     "modified_amino_acid",
                                     "common_rna_dna",
                                     "modified_rna_dna",
                                     "ccp4_mon_lib_rna_dna",
                                     "common_water",
                                     "common_element"}

  def __init__(self,
               include_only_resnames=None,
               include_resnames=None,
               exclude_resnames=None,
               include_residue_classes=None,
               exclude_residue_classes=None
               ):
    self.include_resnames = include_resnames
    self.include_only_resnames = include_only_resnames
    self.exclude_resnames = exclude_resnames
    self.include_residue_classes = include_residue_classes
    self.exclude_residue_classes = exclude_residue_classes

  def query(self,resname):
    include_resnames = self.include_resnames
    include_only_resnames = self.include_only_resnames
    exclude_resnames = self.exclude_resnames
    include_residue_classes = self.include_residue_classes
    exclude_residue_classes = self.exclude_residue_classes

    if include_only_resnames is not None:
      if resname in include_only_resnames:
        return True
      else:
        return False
    else:
      residue_class_name = common_residue_names_get_class(name=resname)

      # first set up the disallowed list
      disallowed_residue_classes = self.default_exclude_residue_classes
      if exclude_residue_classes is not None:
        disallowed_residue_classes = disallowed_residue_classes.union(
          set(exclude_residue_classes))

      disallowed_resnames = self.default_exclude_resnames
      if exclude_resnames is not None:
        disallowed_resnames = disallowed_resnames.union(set(exclude_resnames))

      # if a resname or residue_class is explicitly included,
      # it is selected even if it is on the disallowed list
      retvalue = False
      if residue_class_name in disallowed_residue_classes:
        if include_residue_classes is not None:
          if residue_class_name in include_residue_classes:
            retvalue = True
        if include_resnames is not None:
          if resname in include_resnames:
            retvalue = True

      else:
        if resname in disallowed_resnames:
          if include_resnames is not None:
            if resname in include_resnames:
              retvalue = True
        else:
          retvalue = True

      return retvalue





def extract_ligand_models(model,ligand_selection=None):
  if ligand_selection is None:
    ligand_selection = LigandSelection()
  if model.crystal_symmetry() is None:
    from cctbx.maptbx.box import shift_and_box_model
    model = shift_and_box_model(model)

  ligand_models = []
  sel_str_other = "not (water or nucleotide or protein)"  # "other_cnts" in model.composition()
  other_model = model.select(model.selection(sel_str_other))

  for model in other_model.get_hierarchy().models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for resname in rg.unique_resnames():
          if ligand_selection.query(resname):
              iselection = rg.atoms().extract_i_seq()
              ligand_model = other_model.select(iselection)
              ligand_models.append(ligand_model)

  return ligand_models



class LigandsCC:

  @staticmethod
  def ligand_id_tuple(ligand_model):
    id_tuples = []
    for model in ligand_model.get_hierarchy().models():
      for chain in model.chains():
        for rg in chain.residue_groups():
          for resname in rg.unique_resnames():
            id_tuple = (model.id, chain.id, rg.resseq, resname)
            id_tuples.append(id_tuple)

    assert (len(id_tuples) == 1)
    return id_tuples[0]

  @staticmethod
  def ligand_name(ligand_model):
    ligand_id = LigandsCC.ligand_id_tuple(ligand_model)
    model_id, chain_id, resseq, resname = ligand_id
    name = ""
    name +=resname.strip()
    if len(model_id)>0:
      name += "_MODEL_"+model_id.strip()
    if len(chain_id)>0:
      name += "_CHAIN_"+chain_id.strip()
    if len(resseq)>0:
      name += "_RSEQ_"+resseq.strip()
    return name


  def __init__(self,map_model_manager,resolution):
    self.map_model_manager = map_model_manager
    self.resolution = resolution
    self.ligand_map_model_managers = None # populated on self.process_ligands()
    self.ccmask_ligands = None            # populated on self.process_ligands()
    self.ccmask_full_model = None         # populated upon self.validate()

  def validate(self):
    self.ccmask_full_model = self.map_model_manager.map_model_cc(
      resolution=self.resolution)

  def process_ligands(self,ligand_selection=None):
    if ligand_selection is None:
      ligand_selection = LigandSelection()
    ligand_models = extract_ligand_models(self.map_model_manager.model(),ligand_selection)
    self.ligand_map_model_managers= []
    for ligand_model in ligand_models:
      map_model_manager = MapModelManager(
        map_manager=self.map_model_manager.map_manager(),
                                          model=ligand_model,
                                          ignore_symmetry_conflicts=True)
      mmmlig = map_model_manager.extract_all_maps_around_model()
      self.ligand_map_model_managers.append(mmmlig)
    self.ccmask_ligands = []
    for mmmlig in self.ligand_map_model_managers:
      ccmask = mmmlig.map_model_cc(resolution=self.resolution)
      self.ccmask_ligands.append(ccmask)

