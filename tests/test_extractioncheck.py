# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

import numpy as np

import xpart as xp
import xtrack as xt
from xobjects.test_helpers import for_all_test_contexts

@for_all_test_contexts
def test_extractioncheck(test_context):
    ext = xt.Extractioncheck(x_min=-0.001, 
                             x_max=0.005, 
                             px_min=1e-6, 
                             px_max=3e-5)
    line = xt.Line(elements=[ext], element_names=["ext"])
    line.build_tracker(_context=test_context)
    
    x_test  = [0.0, -0.001,    0.0, 1.3]
    px_test = [0.0,   1e-6, 2.5e-5, 2.9]
    state_test = [1, -1000, -1000, 1]

    p = xp.Particles(x=x_test, px=px_test, _context=test_context)

    line.track(p)

    x_out = test_context.nparray_from_context_array(p.x)[np.argsort(test_context.nparray_from_context_array(p.particle_id))]
    px_out = test_context.nparray_from_context_array(p.px)[np.argsort(test_context.nparray_from_context_array(p.particle_id))]
    state_out = test_context.nparray_from_context_array(p.state)[np.argsort(test_context.nparray_from_context_array(p.particle_id))]

    assert np.all(np.isclose(x_out, x_test, atol=0.0, rtol=1e-10))
    assert np.all(np.isclose(px_out, px_test, atol=0.0, rtol=1e-10))
    assert np.all(np.isclose(state_out, state_test, atol=0.0, rtol=1e-10))
