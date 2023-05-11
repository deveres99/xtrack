import time
import multiprocessing as mp

import numpy as np

import xtrack as xt
import xpart as xp

from xobjects.context_cpu import XobjectPointer

from cpymad.madx import Madx

# Load the line
if __name__ == '__main__':
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

    on_x1_values = np.linspace(50, 150, 18)

    init_part = ttt._initial_particles
    ebe_monitor = ttt.tracking_data

    n_proc = 3

    # slice on_x1_values in n_proc chunks
    inputs = []
    for i_proc in range(n_proc):
        inputs.append(on_x1_values[i_proc::n_proc])

    def f_for_pool(input):

        tw_chunk = []
        for on_x1 in input:
            line.vars['on_x1'] = on_x1

            tw = line.twiss(
                ele_start=ele_index_start,
                ele_stop=ele_index_end,
                twiss_init=tw_init,
                _initial_particles=init_part,
                _ebe_monitor=ebe_monitor
                )
            tw.particle_on_co = None
            tw_chunk.append(tw._data)

        return tw_chunk

    ip1_index = list(tw.name).index('ip1')

    n_repeat = 10

    t1 = time.perf_counter()
    print('Start serial')
    for i in range(n_repeat):
        tw_chunks_serial = list(map(f_for_pool, inputs))
    print('End serial')

    # Remerge chunks
    twisses_serial = []
    for i in range(len(tw_chunks_serial)):
        twisses_serial += tw_chunks_serial[i]

    t2 = time.perf_counter()
    print(f'Elapsed time serial: {t2-t1} s')


    pool = mp.Pool(processes=n_proc)
    t1 = time.perf_counter()
    print('Start parallel')
    for i in range(n_repeat):
        tw_chunks_parallel = pool.map(f_for_pool, inputs)

    # Remerge chunks
    twisses_parallel = []
    for i in range(len(tw_chunks_parallel)):
        twisses_parallel += tw_chunks_parallel[i]

    print('End parallel')

    t2 = time.perf_counter()
    print(f'Elapsed time parallel: {t2-t1} s')


    print('Result serial:')
    for on_x1, twdata in zip(on_x1_values, twisses_serial):
        tw._data = twdata
        print(f'on_x0 = {on_x1}, px["ip1"] = {tw["px", "ip1"]*1e6:f}e-6')

    print('Result parallel:')
    for on_x1, twdata in zip(on_x1_values, twisses_parallel):
        tw._data = twdata
        print(f'on_x0 = {on_x1}, px["ip1"] = {tw["px", "ip1"]*1e6:f}e-6')




    # td = line.tracker._tracker_data
    # td._element_ref_data = XobjectPointer(td._element_ref_data)
    # td._element_dict = None
    # td._elements = None
    # td.element_classes = None
    # td._ElementRefClass = None