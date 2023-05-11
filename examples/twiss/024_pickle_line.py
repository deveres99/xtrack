import time
import multiprocessing as mp

import numpy as np

import xtrack as xt
import xpart as xp

from xobjects.context_cpu import XobjectPointer

from cpymad.madx import Madx

# Load the line
line = xt.Line.from_json(
    '../../test_data/hllhc15_noerrors_nobb/line_w_knobs_and_particle.json')
# line._var_management = None

#line.particle_ref = xp.Particles(p0c=7e12, mass=xp.PROTON_MASS_EV)

import pickle
ss = pickle.dumps(line)

line2 = pickle.loads(ss)