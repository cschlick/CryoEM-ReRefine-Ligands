from __future__ import division, print_function
import sys
import os
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
3. optionally write the ligand map/model data to disk

This was modeled on the tutorial found here [2]

[1] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6130467/
[2] http://cci2.lbl.gov:8080/docs/phenix/doc_program_template/

"""


def extract_ligand_models(model, include_resnames=[], exclude_resnames=[]):
  """
  Parameters
  ----------
  model:	mmtbx.model object, container for macromolecular structure
  include_resnames:	a list of strings, ligand codes to extract. If empty,
                      extract all
  ignored_resnames: a list of strings, ligand codes NOT to extract.

  The resnames are compared to the resnames in
  iotbx.pdb.hierarchy.residue_groups


  Returns
  -------
  ligand_dict: dict, tuple:mmtbx.model

  keys are a tuple describing the ligand's position in the
  hierarchy (model.id, chain.id, rg.resseq, resname)

  values are mmtbx.model objects for the ligand

  """
  # necessary for some cryo-em models
  if model.crystal_symmetry() is None:
    from cctbx.maptbx.box import shift_and_box_model
    model = shift_and_box_model(model)

  # In addition to using the resnames, we try to exclude ligands using
  # iotbx.pdb.common_residue_names_get_class
  exclude_classes = ["common_amino_acid",
                     "modified_amino_acid",
                     "common_rna_dna",
                     "modified_rna_dna",
                     "ccp4_mon_lib_rna_dna",
                     "common_water",
                     "common_element"]

  default_exclude_resnames = ["UNK"]

  exclude_resnames = default_exclude_resnames + exclude_resnames

  extract_all = True if len(include_resnames) == 0 else False

  ligand_dict = {}
  sel_str_other = "not (water or nucleotide or protein)"  # "other_cnts" in model.composition()
  other_model = model.select(model.selection(sel_str_other))
  # I profiled the above method compared to looping over the entire original model,
  # and this was >10X faster.

  for model in other_model.get_hierarchy().models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for resname in rg.unique_resnames():
          if (common_residue_names_get_class(
                  name=resname) not in exclude_classes
                  and resname not in exclude_resnames):
            if resname in include_resnames or extract_all:
              iselection = rg.atoms().extract_i_seq()
              id_tuple = (model.id, chain.id, rg.resseq, resname)
              assert (id_tuple not in ligand_dict.keys())
              ligand_model = other_model.select(iselection)
              ligand_dict[id_tuple] = ligand_model

  return ligand_dict


def ligands_cc_from_mmm(map_model_manager, resolution=None,
                        output_directory=None,
                        log=sys.stdout):
  """
  Parameters
  ----------
  map_model_manager	:	An iotbx.map_model_manager object
  resolution        :   The nominal resolution. If not provided, will calculate
  output_directory	:	string, optionally write out the ligand density

  Returns
  -------
  (ligand_dict, ccmask_values) : tuple of  (dict,list)

  ligand_dict keys    : A tuple describing the ligand's position in the
            hierarchy (model.id, chain.id, rg.resseq, resname)
  ligand_dict values  : aiotbx.map_model_manager objects for the ligand

  ccmask_values       :map_model correlation values for each ligand
  """
  dm = DataManager()
  mmm = map_model_manager
  mm = mmm.map_manager()
  model = mmm.model()
  if resolution is not None:
    mmm.map_manager().set_resolution(resolution)
  if mmm.map_manager()._resolution is None:
    print(
      "The map_manager proved does not have a resolution property. Resolution "
      "will be calculated using the d99 method.",
      file=log)
    resolution = mmm.map_manager().resolution(method="d99")
    print("Using resolution: %.3f A\n"%(resolution))

  ligand_dict = extract_ligand_models(model)
  print("Found " + str(len(ligand_dict)) + " ligands in the model\n", file=log)
  return_ligand_dict = {}
  ccmask_values = []

  for ligand_id, ligand_model in ligand_dict.items():

    model_id, chain_id, rg_resseq, resname = ligand_id

    map_model_manager = MapModelManager(map_manager=mm,
                                        model=ligand_model,
                                        ignore_symmetry_conflicts=True)
    mmmlig = map_model_manager.extract_all_maps_around_model()
    return_ligand_dict[ligand_id] = mmmlig

    ccmask = mmmlig.map_model_cc(resolution=mm.resolution())
    # This uses the five_cc object internally, and returns only CCmask

    ccmask_values.append(ccmask)
    if output_directory is not None:

      if not os.path.exists(output_directory):
        os.makedirs(output_directory)

      ligand_name = resname.strip() + "_" + chain_id.strip()
      ligand_map_file = os.path.abspath(
        os.path.join(output_directory, ligand_name + ".map"))
      ligand_model_file = os.path.abspath(
        os.path.join(output_directory, ligand_name + ".cif"))

      print("Writing model file: " + ligand_map_file, file=log)
      print("Writing map file: " + ligand_model_file, file=log)
      dm.write_model_file(mmmlig.model(),
                          filename=ligand_model_file,
                          extension="cif",
                          overwrite=True)

      dm.write_real_map_file(mmmlig.map_manager(),
                             filename=ligand_map_file,
                             overwrite=True)

    
  return return_ligand_dict, ccmask_values


def read_input_directory(input_directory, output_directory=None):
  """
  Parameters
  ----------
  input_directory  :  A directory to search for map/model pairs
                   :  Each map/model pair should be together in
                   :  a subdirectory.

 output_directory  :  optional, the location to write ligand
                   :  maps and models
  """
  input_directory = os.path.abspath(input_directory)
  candidates = [os.path.join(input_directory, folder) for folder in
                os.listdir(input_directory)]
  folders = [folder for folder in candidates if os.path.isdir(folder)]
  entries = []
  for folder in folders:
    files = [os.path.join(folder, file) for file in os.listdir(folder) if
             not os.path.isdir(file)]
    model_files = []
    model_extensions = []
    map_files = []
    pkl_files = []
    skipping = False
    resolution = None

    for file in files:
      basename, ext, ext_compress = splitext(file)
      if ext in [".pdb", ".cif"]:
        model_files.append(file)
        model_extensions.append(ext)
      elif ext in [".map", ".mrc", ".ccp4"]:
        map_files.append(file)
      elif ext == ".pkl":
        pkl_files.append(file)



      # pkl parsing takes advantage of possible existing output from Dorothee's maps_and_models folder
      if len(pkl_files) == 1:
        with open(pkl_files[0], "rb") as fh:
          ga = pickle.load(fh)
        if isinstance(ga, group_args):
          if hasattr(ga, "error"):
            if (ga.error is not None):
              skipping = True
              print("Skipping due to previous error")
          if hasattr(ga, "resolution"):
            if ga.resolution is not None:
              resolution = ga.resolution

    if not skipping:
      if len(model_files) == 0 or len(map_files) == 0:
        print("Skipping due to insufficient files: " + folder)
        skipping = True
      elif len(map_files) > 1:
        print("Skipping due to multiple map files")
        skipping = True
      else:
        # prefer cif file if possible
        if ".cif" in model_files:
          model_file = model_files[model_extensions.index(".cif")]
        else:
          model_file = model_files[0]
        map_file = map_files[0]

      if not os.path.exists(model_file) or not os.path.exists(map_file):
        print("Skipping due to non-existent files.")
        skipping = True

      entry = group_args(
        input_directory=input_directory,
        entry_directory=os.path.join(input_directory, folder),
        output_directory=output_directory,
        model_file=model_file,
        map_file=map_file,
        resolution=resolution,
        skipping=skipping)

      entries.append(entry)

  return entries


def process_entry(entry):
  """
  Parameters
  ----------
  entry  :  a libtbx.group_args object

  Expects the following fields:
  entry_directory  :  the directory of input maps/models
  input_directory  :  parent of entry directory resides
  output_directory : optional, write ligand maps/models
  map_file         : location of map file
  model_file       : location nof model file
  resolution       : resolution of map
  skipping         : whether to skip due to errors

  Returns
  -------


  """
  if entry.output_directory is not None:
    output_directory = os.path.join(entry.output_directory,os.path.split(
      entry.entry_directory)[-1])
  else:
    output_directory = None
  dm = DataManager()
  mmm = dm.get_map_model_manager(model_file=entry.model_file,
                                 map_files=entry.map_file)
  ignore_ligands = None
  if hasattr(entry,"ignore_ligands"):
    ignore_ligands = entry.ignore_ligands
    # TODO continue to propogate this

  result = run_single_mmm(mmm,
               resolution=entry.resolution,
               output_directory=output_directory)
  return result


def run_full_directory(input_directory,
                       output_directory=None,
                       ignore_ligands = None,
                       nproc=1):
  entries = read_input_directory(input_directory,
                                 output_directory=output_directory)
  if ignore_ligands is not None:
    for entry in entries:
      entry.ignore_ligands = ignore_ligands
  if nproc ==1:
    results = []
    for entry in entries:
      result = process_entry(entry)
      results.append(result)

  return results




def run_single_mmm(map_model_manager,
        resolution=None,
        output_directory=None,
        ignore_ligands = None,
        log=sys.stdout):
  '''
    A tool to extract ligand models and the surrounding density, and calculate map_model correlation.

    Parameters
    ----------
      map_model_manager:   iotbx.map_model_manager holds map and model
      resolution:          nominal resolution of map
      output_directory:		 the optional path to write maps and models
      ignore_ligands  :    list of ligand resnames to ignore
      log:                 output stream

    Returns
    -------
      group_args object containing value of map_model correlation
    '''

  print(
    "Getting map-model correlation for all ligands in the map_model_manager:\n %s" % (
      map_model_manager), file=log)
  if resolution is None:
    print("Resolution not provided.",file=log)
  else:
    print("Resolution will be: %.3f A " % (resolution), file=log)
  map_model_manager.map_manager().set_resolution(resolution)
  ligand_dict, ccmask_values = ligands_cc_from_mmm(map_model_manager,
                                                   output_directory=output_directory,
                                                   log=log)

  print("Results:",file=log)
  for i, ligand_id in enumerate(ligand_dict.keys()):
    model_id, chain_id, rg_resseq, resname = ligand_id
    ligand_name = resname.strip() + "_" + chain_id.strip()
    print("\tLigand: " + ligand_name + "\tCCmask: %.3f" % (ccmask_values[i]),
          file=log)
  return [group_args(
    ligand_dict=ligand_dict,
    ccmask_values=ccmask_values
  )]


