"""
Here we put all system-dependent constants
"""
import numpy as np
import os

LIBRARY_PATH = os.path.dirname((os.path.abspath(__file__)))
THIRDPARTY_PATH = os.path.abspath(os.path.join(LIBRARY_PATH, '..', 'thirdparty'))

# Check for QuickBB
_quickbb_path = os.path.abspath(
    os.path.join(
        THIRDPARTY_PATH, 'quickbb', 'run_quickbb_64.sh')
    )
if os.path.isfile(_quickbb_path):
    QUICKBB_COMMAND = _quickbb_path
else:
    QUICKBB_COMMAND = None

# Check for Tamaki solver
_tamaki_solver_path = os.path.abspath(
    os.path.join(
        THIRDPARTY_PATH, 'pace2017_solvers', 'tamaki_treewidth')
    )
if os.path.isdir(_tamaki_solver_path):
    TAMAKI_SOLVER_PATH = _tamaki_solver_path
else:
    TAMAKI_SOLVER_PATH = None
    raise FileNotFoundError('Tamaki solver not found')
