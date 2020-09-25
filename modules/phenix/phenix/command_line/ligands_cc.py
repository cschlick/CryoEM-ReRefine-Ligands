# LIBTBX_SET_DISPATCHER_NAME phenix.ligands_cc.py
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from phenix.programs import ligands_cc

if __name__ == '__main__':
  run_program(program_class=ligands_cc.Program)