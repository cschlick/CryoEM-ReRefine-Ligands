# -*- coding: utf-8 -*-
from __future__ import division, print_function

from iotbx.data_manager import DataManager
from iotbx.map_model_manager import map_model_manager as MapModelManager
from libtbx.test_utils import approx_equal


import time, sys, os


from phenix.programs.ligands_cc import LigandsCC, LigandSelection


def tst_01():
  # Excercise the ligand selection

  # test defaults only
  resname = "MG"
  ls = LigandSelection(include_resnames=None,
                       include_only_resnames=None,
                       exclude_resnames=None,
                       include_residue_classes=None,
                       exclude_residue_classes=None,
                       )
  assert (ls.query(resname) == False)

  resname = "ADP"
  ls = LigandSelection(include_resnames=None,
                       include_only_resnames=None,
                       exclude_resnames=None,
                       include_residue_classes=None,
                       exclude_residue_classes=None,
                       )
  assert (ls.query(resname) == True)

  # test include only

  resname = "MG"
  ls = LigandSelection(include_resnames=None,
                       include_only_resnames=["ATP"],
                       exclude_resnames=None,
                       include_residue_classes=None,
                       exclude_residue_classes=None,
                       )
  assert(ls.query(resname)== False)

  resname = "ATP"
  ls = LigandSelection(include_resnames=None,
                       include_only_resnames=["ATP"],
                       exclude_resnames=None,
                       include_residue_classes=None,
                       exclude_residue_classes=None,
                       )
  assert(ls.query(resname)==True)

  # test supplied include_resnames
  resname = "MG"
  ls = LigandSelection(include_resnames=["MG"],
                       include_only_resnames=None,
                       exclude_resnames=None,
                       include_residue_classes=None,
                       exclude_residue_classes=None,
                       )
  assert (ls.query(resname) == True)

  resname = "ATP"
  ls = LigandSelection(include_resnames=["MG"],
                       include_only_resnames=None,
                       exclude_resnames=None,
                       include_residue_classes=None,
                       exclude_residue_classes=None,
                       )
  assert (ls.query(resname) == True)

  # test supplied exclude_resnames
  resname = "ATP"
  ls = LigandSelection(include_resnames=None,
                       include_only_resnames=None,
                       exclude_resnames=["ATP"],
                       include_residue_classes=None,
                       exclude_residue_classes=None,
                       )
  assert (ls.query(resname) == False)

  # test supplied include_residue_classes
  resname = "MG"
  ls = LigandSelection(include_resnames=None,
                       include_only_resnames=None,
                       exclude_resnames=None,
                       include_residue_classes=["common_element"],
                       exclude_residue_classes=None)
  assert (ls.query(resname) == True)

  # test supplied exclude_residue_classes
  resname = "ATP"
  ls = LigandSelection(include_resnames=None,
                       include_only_resnames=None,
                       exclude_resnames=None,
                       include_residue_classes=None,
                       exclude_residue_classes=["other"])
  assert (ls.query(resname) == False)



def tst_02():
  # excercise the LigandsCC class
  #TODO: data paths should be existing files in regression directory
  model_file = "data/6qi8_4552/6qi8_4552.cif"
  map_file = "data/6qi8_4552/6qi8_4552.map"
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
  ligands_cc = LigandsCC(mmm, resolution=3.5)
  ligands_cc.validate()
  assert (approx_equal(ligands_cc.ccmask_full_model, 0.022970855944436454))
  ligands_cc.process_ligands()

  assert (len(ligands_cc.ligand_map_model_managers) == 6)

  expected_ligand_ccmask = [
    -0.080566400250404,
    -0.05039665095274361,
    0.02568153796222316,
    -0.058698045705798986,
    0.10242327846306232,
    -0.07061912170822986]

  for i, ccmask in enumerate(expected_ligand_ccmask):
    assert (approx_equal(ccmask, ligands_cc.ccmask_ligands[i]))


if __name__ == "__main__":
  t0 = time.time()
  tst_01()
  tst_02()
  print("Time: %6.4f" % (time.time() - t0))
  print("OK")