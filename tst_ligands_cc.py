# -*- coding: utf-8 -*-
from __future__ import division, print_function

from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager as MapModelManager
from libtbx.test_utils import approx_equal
from cctbx.development.create_models_or_maps import generate_model, \
  generate_map_coefficients
from cctbx.development.create_models_or_maps import generate_map \
  as generate_map_data

import time, sys, os

sys.path.append(".")  # TODO: move this import to phenix directories
from ligands_cc import extract_ligand_models, ligands_cc_from_mmm


def tst_01():
  # Excercise the ligands_cc.extract_ligand_models function
  # Also generate a map for tst_02 
  model_file = "data/6f1u_fragment.pdb"  # TODO: this should be in regression
  # dir

  dm = DataManager()
  dm.process_model_file(model_file)
  model = dm.get_model()

  # this is necessary for cryo-em models without crystal symmetry
  if model.crystal_symmetry() is None:
    from cctbx.maptbx.box import shift_and_box_model
    model = shift_and_box_model(model)

  ligand_dict = extract_ligand_models(model)
  assert (len(ligand_dict) == 1)
  assert (ligand_dict.keys()[0] == ('', 'O', '   1', 'ADP'))
  ligand_model = ligand_dict.values()[0]

  # try to generate a simulated map for the model
  # I have tried this multiple ways, and I still can't get it to work as expected.
  # TODO: ask about this
  d_min = 3.5
  map_coeffs = generate_map_coefficients(model=ligand_model,
                                         d_min=d_min,
                                         scattering_table="xray",
                                         # electron causes error
                                         log=None)
  new_mm = generate_map_data(
    map_coeffs=map_coeffs,
    d_min=d_min,
    gridding=(1, 1, 1),  # this is a placeholder
    wrapping=False,
    origin_shift_grid_units=None,
    high_resolution_real_space_noise_fraction=0,
    log=None)

  ligand_model.set_shift_cart(new_mm.shift_cart())
  ligand_model.set_unit_cell_crystal_symmetry(
    new_mm.unit_cell_crystal_symmetry())
  ############################################# end simulated map

  new_mmm = MapModelManager(map_manager=new_mm, model=ligand_model)
  new_mmm.write_model("data/6f1u_framgent_testing.pdb")
  new_mmm.write_map("data/6f1u_framgent_testing.map")


def tst_02():
  # excercise the ligands_cc.ligands_cc_from_mmm function
  model_file = "data/6f1u_framgent_testing.pdb"
  map_file = "data/6f1u_framgent_testing.map"
  assert(os.path.exists(model_file))
  assert(os.path.exists(map_file))
  dm = DataManager()
  dm.process_model_file(model_file)
  dm.process_real_map_file(map_file)
  model = dm.get_model()
  # this is necessary for cryo-em models without crystal symmetry
  if model.crystal_symmetry() is None:
    from cctbx.maptbx.box import shift_and_box_model
    model = shift_and_box_model(model)

  mmm = MapModelManager(map_manager=dm.get_real_map(), model=dm.get_model())
  ligand_dict, ccmask_values = ligands_cc_from_mmm(mmm, resolution=3.5)
  assert (len(ligand_dict) == 1)
  assert (len(ccmask_values) == 1)
  assert (ccmask_values[0] > 0.9)


if __name__ == "__main__":
  t0 = time.time()
  tst_01()
  tst_02()
  print("Time: %6.4f" % (time.time() - t0))
  print("OK")