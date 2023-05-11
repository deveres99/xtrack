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

tw = line.twiss(
    #verbose=True,
    ele_start=ele_index_start,
    ele_stop=ele_index_end,
    twiss_init=tw_init,
    _ebe_monitor=ttt.tracking_data,
    _initial_particles=ttt._initial_particles
        )

on_x1_values = [50, 100, 150]
buffers = []
for on_x1 in on_x1_values:
    line.vars['on_x1'] = on_x1
    buffers.append(line.tracker._buffer.buffer.copy())

twisses = []
for bb in buffers:
    line._buffer.buffer = bb
    twisses.append(line.twiss())

for on_x1, tw in zip(on_x1_values, twisses):
    print(f'on_x1 = {on_x1}, px["ip1"] = {tw["px", "ip1"]*1e6:f}e-6')