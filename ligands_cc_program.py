from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from libtbx.utils import Sorry, multi_out
from phenix.ligands_cc.ligands_cc import  LigandsCC, LigandSelection
import os




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
    print("Number of ligands found: " + str(len(ligands_cc.ligand_map_model_managers)))
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
        print("\t" + ligand_name)
        ligand_map_file = os.path.join(output_dir, ligand_name + ".map")
        print("\t\t"+"Model: "+ligand_model_file+".cif")
        self.data_manager.write_model_file(mmmlig.model(),
                            filename=ligand_model_file,
                            extension="cif",
                            overwrite=True)
        print("\t\t" + "Map: " + ligand_map_file+"\n")
        self.data_manager.write_real_map_file(mmmlig.map_manager(),
                               filename=ligand_map_file,
                               overwrite=True)


    self.result = group_args(
      ligand_selection = ligand_selection,
      ligands_cc = ligands_cc
    )



  # Get the results
  def get_results(self):
    return self.result
