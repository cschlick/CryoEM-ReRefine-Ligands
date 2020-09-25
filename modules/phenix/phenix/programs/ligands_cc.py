from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from libtbx.utils import Sorry, multi_out

import os

# =============================================================================

class Program(ProgramTemplate):

  description = """
  Program to calculate ligand map-model correlation
  Usage: phenix.ligands_cc <model_file> <map_file> <optional: output_directory>

  """

  # Define the data types that will be used
  datatypes = ['model', 'phil', 'real_map']

  # Input parameters

  master_phil_str = """

  input_files {
    include scope iotbx.map_model_manager.map_model_phil_str
    }

  ligands_cc {
    resolution = None
      .type = float
      .help = Nominal resolution of map
      .short_caption = Resolution
    output_directory = None
      .type = str
      .help = Optional directory to write ligand maps and models
      .short_caption = Output directory
    input_directory = None
      .type = str
      .help = Optionally use a directory to search for map_model pairs 
      .short_caption = Input directory
    ignore_ligands = ["UNK"]
      .type = str
      .help = Do not extract these ligands
      .short_caption = Ignored ligands
    }

  """
  # Define how to determine if inputs are ok
  def validate(self):
    print('Validating inputs', file=self.logger)
    if self.params.ligands_cc.input_directory is None:
      # Model is provided
      if (not self.data_manager.has_models(
      raise_sorry=True,
      expected_n=1,
      exact_count=True)):
        raise Sorry("1 model file is required")
      # Map or set of map coefficients is provided
      if (not self.data_manager.has_real_maps(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)):
        raise Sorry("1 Real map is required.")
    else:
      if not os.path.exists(self.params.ligands_cc.input_directory):
        raise Sorry("Input directory not found.")

    ignore_string = self.params.ligands_cc.ignore_ligands
    ignore_string = ignore_string.strip()
    if ignore_string[0] in ["{","[","("]:
      ignore_string = ignore_string[1:]
    if ignore_string[-1] in ["}","]","]"]:
      ignore_string = ignore_string[:-1]
    ignore_list = ignore_string.split(",")
    self.params.ligands_cc.ignore_ligands = ignore_list



  # Run the program
  def run(self):
    if self.params.ligands_cc.input_directory is None:
      from phenix.ligands_cc import run_single_mmm
      self.result = run_single_mmm(
        map_model_manager = self.data_manager.get_map_model_manager(from_phil=True),
        resolution = self.params.ligands_cc.resolution,
        output_directory = self.params.ligands_cc.output_directory,
        ignore_ligands = self.params.ligands_cc.ignore_ligands,
        log        = self.logger)
    else:
      from phenix.ligands_cc import run_full_directory
      self.result = run_full_directory(self.params.ligands_cc.input_directory,
                                       output_directory=self.params.ligands_cc.output_directory,
                                       ignore_ligands=self.params.ligands_cc.ignore_ligands,
                                       nproc=1)

  # Get the results
  def get_results(self):
    return self.result
