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


# def ligands_cc_from_mmm(map_model_manager, resolution=None,
#                         output_directory=None,
#                         log=sys.stdout): #TODO: move to program class
#   """
#   Parameters
#   ----------
#   map_model_manager	:	An iotbx.map_model_manager object
#   resolution        :   The nominal resolution. If not provided, will calculate
#   output_directory	:	string, optionally write out the ligand density
#
#   Returns
#   -------
#   (ligand_dict, ccmask_values) : tuple of  (dict,list)
#
#   ligand_dict keys    : A tuple describing the ligand's position in the
#             hierarchy (model.id, chain.id, rg.resseq, resname)
#   ligand_dict values  : aiotbx.map_model_manager objects for the ligand
#
#   ccmask_values       :map_model correlation values for each ligand
#   """
#   # TODO: Lower level functions should not do IO
#   dm = DataManager() #TODO: this does not belong here
#   mmm = map_model_manager
#   mm = mmm.map_manager()
#   model = mmm.model()
#   if resolution is not None:
#     mmm.map_manager().set_resolution(resolution)
#   if mmm.map_manager()._resolution is None:
#     print(
#       "The map_manager proved does not have a resolution property. Resolution " # TODO the user should provide resolution
#       "will be calculated using the d99 method.",
#       file=log)
#     resolution = mmm.map_manager().resolution()
#     print("Using resolution: %.3f A\n"%(resolution))
#
#   ligand_dict = extract_ligand_models(model)
#   print("Found " + str(len(ligand_dict)) + " ligands in the model\n", file=log)
#   return_ligand_dict = {}
#   ccmask_values = []
#   # TODO: check total CC before running ligands (no more validation for now)
#   # TODO: make a class for this, validation of proten fit goes in class,
#   #  offer users a cutoff for protein validation
#   for ligand_id, ligand_model in ligand_dict.items():
#
#     model_id, chain_id, rg_resseq, resname = ligand_id
#
#     map_model_manager = MapModelManager(map_manager=mm,
#                                         model=ligand_model,
#                                         ignore_symmetry_conflicts=True)
#     mmmlig = map_model_manager.extract_all_maps_around_model()
#     return_ligand_dict[ligand_id] = mmmlig
#
#     ccmask = mmmlig.map_model_cc(resolution=mm.resolution())
#     # This uses the five_cc object internally, and returns only CCmask
#
#     ccmask_values.append(ccmask)
#     if output_directory is not None:
#
#       if not os.path.exists(output_directory):
#         os.makedirs(output_directory)
#
#       ligand_name = resname.strip() + "_" + chain_id.strip()
#       ligand_map_file = os.path.abspath(
#         os.path.join(output_directory, ligand_name + ".map"))
#       ligand_model_file = os.path.abspath(
#         os.path.join(output_directory, ligand_name + ".cif"))
#
#       print("Writing map file: " + ligand_map_file, file=log)
#       print("Writing model file: " + ligand_model_file, file=log)
#       dm.write_model_file(mmmlig.model(),
#                           filename=ligand_model_file,
#                           extension="cif",
#                           overwrite=True)
#
#       dm.write_real_map_file(mmmlig.map_manager(),
#                              filename=ligand_map_file,
#                              overwrite=True)
#
#
#   return return_ligand_dict, ccmask_values
#
#
# def read_input_directory(input_directory, output_directory=None):
#   # TODO: this is really more something for a script
#   """
#   Parameters
#   ----------
#   input_directory  :  A directory to search for map/model pairs
#                    :  Each map/model pair should be together in
#                    :  a subdirectory.
#
#  output_directory  :  optional, the location to write ligand
#                    :  maps and models
#   """
#   input_directory = os.path.abspath(input_directory)
#   candidates = [os.path.join(input_directory, folder) for folder in
#                 os.listdir(input_directory)]
#   folders = [folder for folder in candidates if os.path.isdir(folder)]
#   entries = []
#   for folder in folders:
#     files = [os.path.join(folder, file) for file in os.listdir(folder) if
#              not os.path.isdir(file)]
#     model_files = []
#     model_extensions = []
#     map_files = []
#     pkl_files = []
#     skipping = False
#     resolution = None
#
#     for file in files:
#       basename, ext, ext_compress = splitext(file)
#       if ext in [".pdb", ".cif"]:
#         model_files.append(file)
#         model_extensions.append(ext)
#       elif ext in [".map", ".mrc", ".ccp4"]:
#         map_files.append(file)
#       elif ext == ".pkl":
#         pkl_files.append(file)
#
#
#
#       # pkl parsing takes advantage of possible existing output from Dorothee's maps_and_models folder
#       if len(pkl_files) == 1:
#         with open(pkl_files[0], "rb") as fh:
#           ga = pickle.load(fh)
#         if isinstance(ga, group_args):
#           if hasattr(ga, "error"):
#             if (ga.error is not None):
#               skipping = True
#               print("Skipping due to previous error")
#           if hasattr(ga, "resolution"):
#             if ga.resolution is not None:
#               resolution = ga.resolution
#
#     if not skipping:
#       if len(model_files) == 0 or len(map_files) == 0:
#         print("Skipping due to insufficient files: " + folder)
#         skipping = True
#       elif len(map_files) > 1:
#         print("Skipping due to multiple map files")
#         skipping = True
#       else:
#         # prefer cif file if possible
#         if ".cif" in model_files:
#           model_file = model_files[model_extensions.index(".cif")]
#         else:
#           model_file = model_files[0]
#         map_file = map_files[0]
#
#       if not os.path.exists(model_file) or not os.path.exists(map_file):
#         print("Skipping due to non-existent files.")
#         skipping = True
#
#       entry = group_args(
#         input_directory=input_directory,
#         entry_directory=os.path.join(input_directory, folder),
#         output_directory=output_directory,
#         model_file=model_file,
#         map_file=map_file,
#         resolution=resolution,
#         skipping=skipping)
#
#       entries.append(entry)
#
#   return entries
#
#
# def process_entry(entry):
#   """
#   Parameters
#   ----------
#   entry  :  a libtbx.group_args object
#
#   Expects the following fields:
#   entry_directory  :  the directory of input maps/models
#   input_directory  :  parent of entry directory resides
#   output_directory : optional, write ligand maps/models
#   map_file         : location of map file
#   model_file       : location nof model file
#   resolution       : resolution of map
#   skipping         : whether to skip due to errors
#
#   Returns
#   -------
#
#
#   """
#   if entry.output_directory is not None:
#     output_directory = os.path.join(entry.output_directory,os.path.split(
#       entry.entry_directory)[-1])
#   else:
#     output_directory = None
#   dm = DataManager()
#   mmm = dm.get_map_model_manager(model_file=entry.model_file,
#                                  map_files=entry.map_file)
#   ignore_ligands = None
#   if hasattr(entry,"ignore_ligands"):
#     ignore_ligands = entry.ignore_ligands
#     # TODO continue to propogate this
#
#   result = run_single_mmm(mmm,
#                resolution=entry.resolution,
#                output_directory=output_directory)
#   return result
#
#
# def run_full_directory(input_directory,
#                        output_directory=None,
#                        ignore_ligands = None,
#                        nproc=1):
#   entries = read_input_directory(input_directory,
#                                  output_directory=output_directory)
#   if ignore_ligands is not None:
#     for entry in entries:
#       entry.ignore_ligands = ignore_ligands
#   if nproc ==1:
#     results = []
#     for entry in entries:
#       result = process_entry(entry)
#       results.append(result)
#
#   return results
#
#
#
#
# def run_single_mmm(map_model_manager,
#         resolution=None,
#         output_directory=None,
#         ignore_ligands = None,
#         log=sys.stdout):
#   '''
#     A tool to extract ligand models and the surrounding density, and calculate map_model correlation.
#
#     Parameters
#     ----------
#       map_model_manager:   iotbx.map_model_manager holds map and model
#       resolution:          nominal resolution of map
#       output_directory:		 the optional path to write maps and models
#       ignore_ligands  :    list of ligand resnames to ignore
#       log:                 output stream
#
#     Returns
#     -------
#       group_args object containing value of map_model correlation
#     '''
#
#   print(
#     "Getting map-model correlation for all ligands in the map_model_manager:\n %s" % (
#       map_model_manager), file=log)
#   if resolution is None:
#     print("Resolution not provided.",file=log)
#   else:
#     print("Resolution will be: %.3f A " % (resolution), file=log)
#   map_model_manager.map_manager().set_resolution(resolution)
#   ligand_dict, ccmask_values = ligands_cc_from_mmm(map_model_manager,
#                                                    output_directory=output_directory,
#                                                    log=log)
#
#   print("Results:",file=log)
#   for i, ligand_id in enumerate(ligand_dict.keys()):
#     model_id, chain_id, rg_resseq, resname = ligand_id
#     ligand_name = resname.strip() + "_" + chain_id.strip()
#     print("\tLigand: " + ligand_name + "\tCCmask: %.3f" % (ccmask_values[i]),
#           file=log)
#   return [group_args(
#     ligand_dict=ligand_dict,
#     ccmask_values=ccmask_values
#   )]
#
#
