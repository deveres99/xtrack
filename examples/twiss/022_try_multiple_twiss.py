import time

import xtrack as xt
import xpart as xp

from cpymad.madx import Madx

# Load the line
line = xt.Line.from_json(
    '../../test_data/hllhc15_noerrors_nobb/line_w_knobs_and_particle.json')
line.particle_ref = xp.Particles(p0c=7e12, mass=xp.PROTON_MASS_EV)
line.build_tracker()

tw_ref = line.twiss()

ele_start_match = 's.ds.l1.b1'
ele_end_match = 'e.ds.r1.b1'
tw_init = tw_ref.get_twiss_init(ele_start_match)

ele_index_start = line.element_names.index(ele_start_match)
ele_index_end = line.element_names.index(ele_end_match)

ttt = line.twiss(
    #verbose=True,
    ele_start=ele_index_start,
    ele_stop=ele_index_end,
    twiss_init=tw_init,
    _keep_initial_particles=True,
    _keep_tracking_data=True,
    )

line._kill_twiss = False


tw = line.twiss(
    #verbose=True,
    ele_start=ele_index_start,
    ele_stop=ele_index_end,
    twiss_init=tw_init,
    _ebe_monitor=ttt.tracking_data,
    _initial_particles=ttt._initial_particles
        )

