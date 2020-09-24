from __future__ import division, print_function
import sys
import os

from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager as MapModelManager
from iotbx.pdb import common_residue_names_get_class
from libtbx import group_args

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
      "The map_manager proved does not have a resolution property. This will be calculated now, but it is slow.",
      file=log)

  ligand_dict = extract_ligand_models(model)
  print("Found " + str(len(ligand_dict)) + " ligands in the model ", file=log)
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
        os.mkdir(output_directory)

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


def run(map_model_manager,
        resolution=None,
        output_directory=None,
        log=sys.stdout):
  '''
    A tool to extract ligand models and the surrounding density, and calculate map_model correlation.

    Parameters
    ----------
      map_model_manager:   iotbx.map_model_manager holds map and model
      resolution:          nominal resolution of map
      output_directory:		 the optional path to write maps and models
      log:                 output stream

    Returns
    -------
      group_args object containing value of map_model correlation
    '''

  print(
    "Getting map-model correlation for all ligands in the map_model_manager:\n %s" % (
      map_model_manager), file=log)
  print("Resolution will be: %.3f A " % (resolution), file=log)
  map_model_manager.map_manager().set_resolution(resolution)
  ligand_dict, ccmask_values = ligands_cc_from_mmm(map_model_manager,
                                                   output_directory=output_directory,
                                                   log=log)
  return group_args(
    ligand_dict=ligand_dict,
    ccmask_values=ccmask_values
  )


if __name__ == "__main__":
  running = False
  args = sys.argv[1:]
  if len(args) not in [3, 4]:
    print("Usage: 'ligands_cc.py <model_file> <map_file> <resolution> "
          "<optional: output_directory>'")
  elif len(args) == 3:
    model_file = args[0]
    map_file = args[1]
    resolution = float(args[2])
    output_directory = None
    running = True
  elif len(args) == 4:
    model_file = args[0]
    map_file = args[1]
    resolution = float(args[2])
    output_directory = args[3]
    running = True

  if running:
    log = sys.stdout

    dm = DataManager()
    mmm = dm.get_map_model_manager(model_file=model_file, map_files=map_file,
                                   log=log)
    result = run(mmm,
                 resolution=resolution,
                 output_directory=output_directory)

    ligand_dict = result.ligand_dict
    ccmask_values = result.ccmask_values
    print("Results:")
    for i, ligand_id in enumerate(ligand_dict.keys()):
      model_id, chain_id, rg_resseq, resname = ligand_id
      ligand_name = resname.strip() + "_" + chain_id.strip()
      print("\tLigand: " + ligand_name + "\tCC: %.3f" % (ccmask_values[i]))
