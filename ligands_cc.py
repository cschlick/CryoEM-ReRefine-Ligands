from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from libtbx.utils import Sorry, multi_out
import os


from iotbx.map_model_manager import map_model_manager as MapModelManager
from iotbx.pdb import common_residue_names_get_class



"""
Summary
-------

This module contains functionality to extract ligand models and calculate 
ligand CCmask, as described here.[1]


[1] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6130467/

"""


class LigandSelection:
  """
  A companion class of LigandsCC. The selection class stores preferences for
  how to extract ligands. LigandsCC extracts on the basis of:
    1) residue name :  iotbx_pdb_hierarchy_ext.residue_group.unique_resnames()
    2) residue class : iotbx.pdb.common_residue_names_get_class(resname)

  Instantiation is done with lists of strings. The query(resname) method
  returns True/False for whether the resname meets the selection criteria.
  """

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

    retvalue = False
    if include_only_resnames is not None:
      if resname in include_only_resnames:
        retvalue =  True
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
  """
  Extract ligands from a parent model.

  Params
  ------
  model  :  mmtbx.model.model_manager object
  ligand_selection : LigandSelection object

  Returns
  -------
  ligand_models : list of mmtbx.model.model_manager objects,
                  one for each ligand extracted
  """
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
  """
  A class to extract ligands from an iotbx.map_model_manager object.

  Vars
  ----
  map_model_manager : the parent map_model_manager
  resolution : resolution of the parent map
  ccmask_full_model : the ccmask of the full map/model
  ligand_map_model_managers : a list of ligand map_model_managers
  ccmask_ligands : a list of ccmask values for each ligand


  """

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
    self.ccmask_full_model = None         # populated upon self.validate()
    self.ligand_map_model_managers = None # populated on self.process_ligands()
    self.ccmask_ligands = None            # populated on self.process_ligands()

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

class Program(ProgramTemplate):

  description = """
  Program to calculate ligand map-model correlation
  Usage: 
    phenix.ligands_cc <model_file> <map_file> resolution=<resolution>

  """

  # Define the data types that will be used
  datatypes = ['model', 'phil', 'real_map']

  master_phil_str = """

  input_files {
    include scope iotbx.map_model_manager.map_model_phil_str
    }

  ligands_cc {
    resolution = None
      .type = float
      .help = Nominal resolution of map
    write_ligand_files = False
      .type = bool
      .help = Whether to write the map and model for each ligand  
    output_directory = None
      .type = str
      .help = Optional directory to write ligand maps and models
    include_only_resnames = None
      .type = str
      .multiple = True
      .help = Only extract ligands with this residue name
    include_resnames = None
      .type = str
      .multiple = True
      .help = Do extract ligands with this residue name
    exclude_resnames = None
      .type = str
      .multiple = True
      .help = Do NOT extract ligands with this residue name
    include_residue_classes = None
      .type = str
      .multiple = True
      .help = Do extract ligands of this residue class
    exclude_residue_classes = None
      .type = str
      .multiple = True
      .help = Do NOT extract ligands of this residue class
    }

  """
  # Define how to determine if inputs are ok
  def validate(self):
    print('Validating inputs', file=self.logger)

    # Model is provided
    self.data_manager.has_models(
      raise_sorry=True,
      expected_n=1,
      exact_count=True)
    # Map provided
    self.data_manager.has_real_maps(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

    if self.params.ligands_cc.resolution is None:
      raise Sorry("Resolution must be provided.")

    if self.params.ligands_cc.write_ligand_files:
      if self.params.ligands_cc.output_directory is None:
        self.params.ligands_cc.output_directory = os.getcwd()
      else:
        if not os.path.exists(self.params.ligands_cc.output_directory):
          os.mkdir(self.params.ligands_cc.output_directory)

    # set all ligand selection properties to lists, or if empty, None
    ligand_selection_params = [
      "include_only_resnames",
      "include_resnames",
      "exclude_resnames",
      "include_residue_classes",
      "exclude_residue_classes"
    ]
    for lsp in ligand_selection_params:
      setattr(self.params.ligands_cc,lsp,
              list(getattr(self.params.ligands_cc,lsp)))
      if len(getattr(self.params.ligands_cc,lsp)) == 0:
        setattr(self.params.ligands_cc,lsp,None)


  # Run the program
  def run(self):
    map_model_manager = self.data_manager.get_map_model_manager()
    ligands_cc = LigandsCC(map_model_manager,self.params.ligands_cc.resolution)

    ligand_selection = LigandSelection(
      include_resnames=self.params.ligands_cc.include_resnames,
      include_only_resnames=self.params.ligands_cc.include_only_resnames,
      exclude_resnames=self.params.ligands_cc.exclude_resnames,
      include_residue_classes=self.params.ligands_cc.include_residue_classes,
      exclude_residue_classes=self.params.ligands_cc.exclude_residue_classes)

    ligands_cc.process_ligands(ligand_selection)

    print("Number of ligands found: " + str(len(
      ligands_cc.ligand_map_model_managers)),file=self.logger)
    print("Results:",file=self.logger)

    for i,mmmlig in enumerate(ligands_cc.ligand_map_model_managers):
      ligand_name = LigandsCC.ligand_name(mmmlig.model())
      print("\tLigand: " + ligand_name + "\tCCmask: %.3f" % (
        ligands_cc.ccmask_ligands[i]),file=self.logger)

    if self.params.ligands_cc.write_ligand_files:
      print("\nWriting ligand maps and models:",file=self.logger)
      output_dir = self.params.ligands_cc.output_directory
      for mmmlig in ligands_cc.ligand_map_model_managers:
        ligand_name = LigandsCC.ligand_name(mmmlig.model())
        ligand_model_file = os.path.join(output_dir, ligand_name)
        print("\t" + ligand_name,file=self.logger)
        ligand_map_file = os.path.join(output_dir, ligand_name + ".map")
        print("\t\t"+"Model: "+ligand_model_file+".cif",file=self.logger)
        self.data_manager.write_model_file(mmmlig.model(),
                            filename=ligand_model_file,
                            extension="cif",
                            overwrite=True)
        print("\t\t" + "Map: " + ligand_map_file+"\n",file=self.logger)
        self.data_manager.write_real_map_file(mmmlig.map_manager(),
                               filename=ligand_map_file,
                               overwrite=True)


    self.result = group_args(
      ligand_selection = ligand_selection,
      ligands_cc = ligands_cc
    )
