from __future__ import absolute_import, division, print_function
import sys
from phenix.program_template import ProgramTemplate
from libtbx import group_args

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

  crystal_info {
    resolution = None
      .type = float
      .help = Nominal resolution of map
      .short_caption = Resolution
    }

  """
  # how to add output directory to Phil parameters?

  # Define how to determine if inputs are ok
  def validate(self):

    # Expect exactly one map and one model. Stop if not the case
    self.data_manager.has_real_maps(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

  # Set any defaults
  def set_defaults(self):
    pass


  # Run the program
  def run(self):
    sys.path.append("../")  # TODO: move import to phenix directory structure
    from ligands_cc import run
    self.result = run(
      map_model_manager = self.data_manager.get_map_model_manager(from_phil=True),
      resolution = self.params.crystal_info.resolution,
      log        = self.logger)

  # Get the results
  def get_results(self):
    return group_args(
      ligand_dict=self.result.ligand_dict,
      ccmask_values=self.result.ccmask_values
     )
