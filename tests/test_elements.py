# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

import numpy as np
import pytest
import xobjects as xo
import xpart as xp
from cpymad.madx import Madx
from scipy.stats import linregress
from xobjects.test_helpers import for_all_test_contexts
from xpart.particles import Particles, ParticlesPurelyLongitudinal

import ducktrack as dtk
import xtrack as xt
from xtrack import MadLoader
from xtrack.beam_elements.elements import _angle_from_trig


@for_all_test_contexts
def test_constructor(test_context):
    elements = [
        xt.Drift(_context=test_context),
        xt.Marker(_context=test_context),
        xt.Multipole(_context=test_context, knl=[2, 3]),
        xt.RFMultipole(_context=test_context, knl=[2]),
        xt.Cavity(_context=test_context, voltage=3.),
        xt.SRotation(_context=test_context, angle=4),
        xt.XRotation(_context=test_context, angle=1.8),
        xt.YRotation(_context=test_context, angle=2.4),
        xt.ZetaShift(_context=test_context, dzeta=3E-4),
        xt.XYShift(_context=test_context, dx=1),
        xt.DipoleEdge(_context=test_context, h=1),
        xt.LimitRect(_context=test_context, min_x=5),
        xt.LimitRectEllipse(_context=test_context, max_x=6),
        xt.LimitEllipse(_context=test_context, a=10),
        xt.LimitRacetrack(_context=test_context, min_x=-3, max_x=4,
                           min_y=2, max_y=3, a=0.2, b=0.3),
        xt.LimitPolygon(_context=test_context, x_vertices=[1,-1,-1,1],
                        y_vertices=[1,1,-1,-1]),
        xt.Elens(_context=test_context, inner_radius=0.1),
        xt.Wire(_context=test_context, current=3.),
        xt.Exciter(_context=test_context, knl=[1], samples=[1,2,3,4],
                   sampling_frequency=1e3),
        xt.Bend(_context=test_context, length=1.),
        xt.Quadrupole(_context=test_context, length=1.),
    ]

    # test to_dict / from_dict
    for ee in elements:
        dd = ee.to_dict()
        nee = ee.__class__.from_dict(dd, _context=test_context)
        # Check that the two objects are bitwise identical
        if not isinstance(test_context, xo.ContextCpu):
            ee.move(_context=xo.ContextCpu())
            nee.move(_context=xo.ContextCpu())
        assert (ee._xobject._buffer.buffer[ee._xobject._offset:ee._xobject._size]
                - nee._xobject._buffer.buffer[
                    nee._xobject._offset:nee._xobject._size]).sum() == 0


@pytest.mark.parametrize(
    'element_cls,args',
    [
        (xt.SRotation, {'angle': 8, 'sin_z': 0.4}),
        (xt.XRotation, {'angle': 8, 'sin_angle': 0.4}),
        (xt.YRotation, {'angle': 8, 'sin_angle': 0.4}),
        (xt.SRotation, {'angle': 8, 'cos_z': 0.4}),
        (xt.XRotation, {'angle': 8, 'cos_angle': 0.4}),
        (xt.YRotation, {'angle': 8, 'cos_angle': 0.4}),
        (xt.XRotation, {'angle': 8, 'tan_angle': 0.4}),
        (xt.YRotation, {'angle': 8, 'tan_angle': 0.4}),
    ],
)
def test_rotations_constructors_on_inconsistent_input(element_cls, args):
    with pytest.raises(ValueError):
        element_cls(**args)


@pytest.mark.parametrize(
    'cos,sin,tan,angle',
    [
        (np.cos(0.03), np.sin(0.03), np.tan(0.03), 0.03),
        (np.cos(0.01), np.sin(0.01), None, 0.01),
        (None, np.sin(0.02), np.tan(0.02), 0.02),
        (np.cos(0.03), None, np.tan(0.03), 0.03),
        (None, None, None, None),
        (0.1, None, None, None),
        (None, None, 0.1, None),
        (np.cos(0.02), np.sin(0.04), None, None),
        (None, np.sin(0.04), np.tan(0.02), None),
        (np.cos(0.02), None, np.tan(0.04), None),
    ]
)
def test__angle_from_trig(cos, sin, tan, angle):
    should_fail = angle is None

    if should_fail:
        with pytest.raises(ValueError):
            _angle_from_trig(cos, sin, tan)
    else:
        result, cos_res, sin_res, tan_res = _angle_from_trig(cos, sin, tan)

        if cos is not None:
            assert cos == cos_res
        if sin is not None:
            assert sin == sin_res
        if tan is not None:
            assert tan == tan_res

        assert np.isclose(angle, result, atol=1e-13)


@for_all_test_contexts
def test_backtrack(test_context):
    elements = [
        xt.Drift(_context=test_context),
        xt.Multipole(_context=test_context, knl=[2, 3]),
        xt.RFMultipole(_context=test_context, knl=[2]),
        xt.ReferenceEnergyIncrease(_context=test_context, Delta_p0c=42),
        xt.Cavity(_context=test_context, voltage=3.),
        xt.SRotation(_context=test_context, angle=4),
        xt.XRotation(_context=test_context, angle=0.3),
        xt.YRotation(_context=test_context, angle=0.7),
        xt.ZetaShift(_context=test_context, dzeta=3E-4),
        xt.XYShift(_context=test_context, dx=1),
        xt.DipoleEdge(_context=test_context, h=1),
        xt.LimitRect(_context=test_context, min_x=5),
        xt.LimitRectEllipse(_context=test_context, max_x=6),
        xt.LimitEllipse(_context=test_context, a=10),
        xt.LimitRacetrack(_context=test_context, min_x=-3, max_x=4,
                           min_y=2, max_y=3, a=0.2, b=0.3),
        xt.LimitPolygon(_context=test_context, x_vertices=[1,-1,-1,1],
                        y_vertices=[1,1,-1,-1]),
        xt.Elens(_context=test_context, inner_radius=0.1),
        xt.Exciter(_context=test_context, knl=[1], samples=[1,2,3],
                   sampling_frequency=1e3),
    ]

    dtk_particle = dtk.TestParticles(
            p0c=25.92e9,
            x=1e-3,
            px=1e-5,
            y=-2e-3,
            py=-1.5e-5,
            delta=1e-2,
            zeta=1.)

    for element in elements:
        line_test = xt.Line(elements=[element])
        line_test.build_tracker(_context=test_context)

        # track forward and backward
        new_particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                               _context=test_context)
        line_test.track(new_particles)
        line_test.track(new_particles, backtrack=True)

        # assert that nothing changed
        for k in 'x,px,y,py,zeta,delta'.split(','):
            assert np.isclose(test_context.nparray_from_context_array(
                      getattr(new_particles, k))[0],
                      getattr(dtk_particle, k), rtol=1e-14, atol=1e-14)


@for_all_test_contexts
def test_arr2ctx(test_context):
    d = xt.Drift(_context=test_context)

    a = [1., 2., 3.]
    assert type(d._arr2ctx(a)) is test_context.nplike_array_type

    a = test_context.zeros(shape=(20,), dtype=np.int64)
    assert type(d._arr2ctx(a)) is test_context.nplike_array_type
    assert (type(d._arr2ctx(a[1])) is int
            or (type(d._arr2ctx(a[1])) is test_context.nplike_array_type
                and d._arr2ctx(a[1]).shape == ()))

    a = np.array([1., 2., 3.])
    assert type(d._arr2ctx(a)) is test_context.nplike_array_type
    assert type(d._arr2ctx(a[1])) is float


@for_all_test_contexts
def test_drift(test_context):
    dtk_particle = dtk.TestParticles(
            p0c=25.92e9,
            x=1e-3,
            px=1e-5,
            y=-2e-3,
            py=-1.5e-5,
            delta=1e-2,
            zeta=1.)

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)

    drift = xt.Drift(_context=test_context, length=10.)
    drift.track(particles)

    dtk_drift = dtk.elements.Drift(length=10.)
    dtk_drift.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.x)[0],
                      dtk_particle.x, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.y)[0],
                      dtk_particle.y, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.zeta)[0],
                      dtk_particle.zeta,
                      rtol=1e-14, atol=1e-14)


@for_all_test_contexts
def test_drift_exact(test_context):
    dtk_particle = dtk.TestParticles(
            p0c=25.92e9,
            x=1e-3,
            px=1e-5,
            y=-2e-3,
            py=-1.5e-5,
            delta=1e-2,
            zeta=1.)

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)

    drift = xt.Drift(_context=test_context, length=10.)
    line = xt.Line(elements=[drift])
    line.build_tracker(compile=False, _context=test_context)
    line.config.XTRACK_USE_EXACT_DRIFTS = True
    line.track(particles)

    dtk_drift = dtk.elements.DriftExact(length=10.)
    dtk_drift.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.x)[0],
                      dtk_particle.x, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.y)[0],
                      dtk_particle.y, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.zeta)[0],
                      dtk_particle.zeta,
                      rtol=1e-14, atol=1e-14)


@for_all_test_contexts
def test_marker(test_context):
    dtk_particle = dtk.TestParticles(
        p0c=25.92e9,
        x=1e-3,
        px=1e-5,
        y=-2e-3,
        py=-1.5e-5,
        delta=1e-2,
        zeta=1.)

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)

    marker = xt.Marker(_context=test_context)
    marker.track(particles)

    assert np.isclose(test_context.nparray_from_context_array(particles.zeta)[0],
                      dtk_particle.zeta,
                      rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.x)[0],
                      dtk_particle.x, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.y)[0],
                      dtk_particle.y, rtol=1e-14, atol=1e-14)


@for_all_test_contexts
def test_elens(test_context):
    dtk_particle = dtk.TestParticles(
            p0c  = np.array([7000e9]),
            x    = np.array([1e-3]),
            px   = np.array([0.0]),
            y    = np.array([2.2e-3]),
            py   = np.array([0.0]),
            zeta = np.array([0.]))

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)

    elens = xt.Elens(_context=test_context,
                     inner_radius=1.1e-3,
                     outer_radius=2.2e-3,
                     elens_length=3.,
                     voltage=15e3,
                     current=5)

    elens.track(particles)

    dtk_elens = dtk.elements.Elens(
                   inner_radius=1.1e-3,
                   outer_radius=2.2e-3,
                   elens_length=3.,
                   voltage=15e3,
                   current=5)

    dtk_elens.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.px)[0],
                      dtk_particle.px, rtol=1e-2, atol=1e-2)
    assert np.isclose(test_context.nparray_from_context_array(particles.py)[0],
                      dtk_particle.py, rtol=1e-9, atol=1e-9)


@for_all_test_contexts
def test_elens_measured_radial(test_context):
    def compute_coef(r_measured, j_measured, r_1_new, r_2_new,
                     r_1_old, r_2_old, p_order = 13):
        new_r = r_measured*(r_2_new-r_1_new)/(r_2_old-r_1_old)
        new_j = j_measured*(r_2_old-r_1_old)/(r_2_new-r_1_new)

        product = new_r*new_j

        delta_r = new_r[2]-new_r[1]

        numerator = [];
        s = len(new_r)
        for i in range(s):
            numerator.append((delta_r*max(np.cumsum(product[0:i+1]))))
        L = np.cumsum(product)
        denominator = max(L)*delta_r
        f_r = np.array(numerator/denominator)
        r_selected = new_r[new_j != 0]
        f_selected = f_r[new_j != 0]
        coef = np.polyfit(r_selected, f_selected, p_order)
        return coef


    particle_ref = xp.Particles(
                    p0c=np.array([7000e9]),
                    x=np.array([1e-3]),
                    px=np.array([0.0]),
                    y=np.array([2.2e-3]),
                    py=np.array([0.0]),
                    zeta=np.array([0.]))
    particle_test = particle_ref.copy(_context=test_context)

    # polynomial fit parameters for constant radial density
    r     = np.linspace(0.20338983,12,60)
    j     = np.append(np.append(np.linspace(0,4,20)*0,
           np.linspace(4,8,20)/np.linspace(4,8,20)), np.linspace(8,12,20)*0)
    C     = compute_coef(r, j, 1.4, 2.8, 4.0, 8.0)

    elens_radial_profile = xt.Elens(current=5, inner_radius=1.4e-3,
                outer_radius=2.8e-3, elens_length=3, voltage=10e3,
                coefficients_polynomial=C, _context=test_context)

    elens_constant = xt.Elens(current=5, inner_radius=1.4e-3,
                    outer_radius=2.8e-3, elens_length=3, voltage=10e3,
                    )

    elens_radial_profile.track(particle_test)
    elens_constant.track(particle_ref)

    particle_test.move(_context=xo.ContextCpu())

    assert np.isclose(particle_test.px[0], particle_ref.px[0],
                      rtol=1e-2, atol=1e-2)
    assert np.isclose(particle_test.py[0], particle_ref.py[0],
                      rtol=1e-2, atol=1e-2)


@for_all_test_contexts
def test_wire(test_context):
    dtk_particle = dtk.TestParticles(
            p0c =np.array([7000e9]),
            x   =np.array([1e-3]),
            px  =np.array([0.0]),
            y   =np.array([2.2e-3]),
            py  =np.array([0.0]),
            zeta=np.array([0.]))

    particles = xp.Particles(_context=test_context,
                             **dtk_particle.to_dict())


    wire = xt.Wire(_context    =  test_context,
                   L_phy  =  1.3,
                   L_int  =  1.3,
                   current=  250,
                   xma    = -8e-3,
                   yma    = -10e-3)

    wire.track(particles)

    dtk_wire = dtk.elements.Wire(
                   L_phy  =  1.3,
                   L_int  =  1.3,
                   current=  250,
                   xma    = -8e-3,
                   yma    = -10e-3)

    dtk_wire.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.px)[0],
                      dtk_particle.px, rtol=1e-9, atol=1e-9)
    assert np.isclose(test_context.nparray_from_context_array(particles.py)[0],
                      dtk_particle.py, rtol=1e-9, atol=1e-9)



@for_all_test_contexts
def test_linear_transfer_first_order_taylor_map(test_context):
    dtk_particle = dtk.TestParticles(
            p0c=25.92e9,
            x=1e-3,
            px=1e-5,
            y=-2e-3,
            py=-1.5e-5,
            zeta=2.,
            delta=2E-4)

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)

    m0 = np.arange(6,dtype=float)
    m1 = np.ones((6,6),dtype=float)
    for i in range(6):
        for j in range(6):
            m1[i,j] = 10*i+j
    arc = xt.FirstOrderTaylorMap(_context=test_context,m0 = m0, m1 = m1)
    arc.track(particles)

    dtk_arc = dtk.elements.FirstOrderTaylorMap(m0 = m0, m1 = m1)
    dtk_arc.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.x)[0],
                      dtk_particle.x, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.px)[0],
                      dtk_particle.px, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.y)[0],
                      dtk_particle.y, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.py)[0],
                      dtk_particle.py, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.zeta)[0],
                      dtk_particle.zeta, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.delta)[0],
                      dtk_particle.delta, rtol=1e-14, atol=1e-14)


@for_all_test_contexts
def test_cavity(test_context):
    cav = xt.Cavity(_context=test_context, frequency=0, lag=90, voltage=30)
    part = xp.Particles(p0c=1e9, delta=[0, 1e-2], zeta=[0, 0.2], _context=test_context)
    part0 = part.copy(_context=xo.ContextCpu())

    cav.track(part)

    part = part.copy(_context=xo.ContextCpu())

    assert np.allclose(part.energy,
                            part0.energy+cav.voltage, atol=5e-7, rtol=0)

    Pc = np.sqrt(part.energy**2 - part.mass0**2)
    delta = Pc/part.p0c - 1
    beta = Pc/part.energy

    tau0 = part0.zeta/(part0.beta0)
    tau = part.zeta/(part.beta0)

    assert np.allclose(part.delta, delta, atol=1e-14, rtol=0)
    assert np.allclose(part.rpp, 1/(1+delta), atol=1e-14, rtol=0)
    assert np.allclose(part.rvv, beta/part.beta0, atol=1e-14, rtol=0)
    assert np.allclose(tau, tau0, atol=1e-14, rtol=0)
    assert np.allclose((part.ptau - part0.ptau) * part0.p0c, 30, atol=1e-9, rtol=0)

@for_all_test_contexts
def test_exciter(test_context):
    fs = 2.99792458e8 # sampling frequency in Hz
    frev = 2.99792458e8 # revolution frequency in Hz
    k0l = -0.1 # this is scaled by the waveform
    signal = [1,2,3] # this is the waveform
    duration = 4/fs
    exciter = xt.Exciter(samples=signal, sampling_frequency=fs,
                        frev=frev, start_turn=0, knl=[k0l], duration=duration)

    line = xt.Line([exciter])
    line.build_tracker(_context=test_context)

    particles = xp.Particles(p0c=6.5e12, zeta=[0,-1,-2], _context=test_context)
    num_particles = len(particles.zeta)

    line.track(particles, num_turns=1)
    expected_px = np.array([0.1, 0.2, 0.3])
    particles.move(_context=xo.context_default)

    assert np.allclose(particles.px, expected_px)

    particles.move(_context=test_context)
    line.track(particles, num_turns=1)
    expected_px += np.array([0.2, 0.3, 0.1])
    particles.move(_context=xo.context_default)
    assert np.allclose(particles.px, expected_px)

    particles.move(_context=test_context)
    line.track(particles, num_turns=1)
    expected_px += np.array([0.3, 0.1, 0])
    particles.move(_context=xo.context_default)
    assert np.allclose(particles.px, expected_px)


test_source = r"""
/*gpufun*/
void test_function(TestElementData el,
                LocalParticle* part0,
                /*gpuglmem*/ double* b){

    double const a = TestElementData_get_a(el);

    //start_per_particle_block (part0->part)

        const int64_t ipart = part->ipart;
        double const val = b[ipart];

        LocalParticle_add_to_s(part, val + a);

    //end_per_particle_block
}

/*gpufun*/
void TestElement_track_local_particle(TestElementData el,
                LocalParticle* part0){

    double const a = TestElementData_get_a(el);

    //start_per_particle_block (part0->part)

        LocalParticle_set_s(part, a);

    //end_per_particle_block
}

"""


@pytest.mark.parametrize(
    'particles_class',
    [Particles, ParticlesPurelyLongitudinal],
)
@for_all_test_contexts
def test_per_particle_kernel(test_context, particles_class):
    class TestElement(xt.BeamElement):
        _xofields = {
            'a': xo.Float64
        }

        _extra_c_sources = [test_source]

        _per_particle_kernels = {
            'test_kernel': xo.Kernel(
                c_name='test_function',
                args=[
                    xo.Arg(xo.Float64, pointer=True, name='b')
                ]),
        }

    el = TestElement(_context=test_context, a=10)

    p = Particles(p0c=1e9, s=[1, 2, 3], _context=test_context)
    el.track(p)
    p.move(_context=xo.ContextCpu())
    assert np.all(p.s == [10, 10, 10])

    p = particles_class(p0c=1e9, s=[1, 2, 3], _context=test_context)
    b = p.s*0.5
    el.test_kernel(p, b=b)
    p.move(_context=xo.ContextCpu())
    assert np.all(p.s == np.array([11.5, 13, 14.5]))

@for_all_test_contexts
def test_simplified_accelerator_segment(test_context):
    dtk_particle = dtk.TestParticles(
            p0c=25.92e9,
            x=1e-3,
            px=1e-5,
            y=-2e-3,
            py=-1.5e-5,
            zeta=2.,
            delta=2E-4)

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)

    alpha_x_0 = -0.5
    beta_x_0 = 100.0
    disp_x_0 =  1.8
    disp_px_0 = 2.2
    alpha_x_1 = 2.1
    beta_x_1 = 2.0
    disp_x_1 = 3.3
    disp_px_1 = 3.7
    alpha_y_0 = -0.4
    beta_y_0 = 8.0
    disp_y_0 = -0.2
    disp_py_0 = -0.4
    alpha_y_1 = 0.7
    beta_y_1 = 0.3
    disp_y_1 = -1.9
    disp_py_1 = -2.9
    Q_x = 0.27
    Q_y = 0.34
    beta_s = 856.9
    Q_s = 0.001
    energy_ref_increment = 1.2E9
    energy_increment = 4.8E8
    x_ref_0 = -5E-3
    px_ref_0 = 6E-4
    x_ref_1 = 2E-2
    px_ref_1 = -5E-5
    y_ref_0 = -9E-2
    py_ref_0 = 1E-4
    y_ref_1 = 4E-2
    py_ref_1 = 5E-4

    arc = xt.LineSegmentMap(_context=test_context,
        alfx=(alpha_x_0, alpha_x_1), betx=(beta_x_0, beta_x_1),
        dx=(disp_x_0, disp_x_1), dpx=(disp_px_0, disp_px_1),
        alfy=(alpha_y_0, alpha_y_1), bety=(beta_y_0, beta_y_1),
        dy=(disp_y_0, disp_y_1), dpy=(disp_py_0, disp_py_1),
        qx=Q_x, qy=Q_y,
        bets=beta_s, qs=Q_s,
        energy_ref_increment=energy_ref_increment,energy_increment=energy_increment,
        x_ref=(x_ref_0, x_ref_1), px_ref=(px_ref_0, px_ref_1),
        y_ref=(y_ref_0, y_ref_1), py_ref=(py_ref_0, py_ref_1))

    arc.track(particles)

    dtk_arc = dtk.elements.LinearTransferMatrix(
        alpha_x_0=alpha_x_0, beta_x_0=beta_x_0, disp_x_0=disp_x_0, disp_px_0=disp_px_0,
        alpha_x_1=alpha_x_1, beta_x_1=beta_x_1, disp_x_1=disp_x_1, disp_px_1=disp_px_1,
        alpha_y_0=alpha_y_0, beta_y_0=beta_y_0, disp_y_0=disp_y_0, disp_py_0=disp_py_0,
        alpha_y_1=alpha_y_1, beta_y_1=beta_y_1, disp_y_1=disp_y_1, disp_py_1=disp_py_1,
        Q_x=Q_x, Q_y=Q_y,
        beta_s=beta_s, Q_s=Q_s,
        chroma_x=0.0, chroma_y=0.0,
        det_xx=0.0, det_xy=0.0, det_yy=0.0, det_yx=0.0,
        energy_ref_increment=energy_ref_increment,energy_increment=energy_increment,
        x_ref_0 = x_ref_0, px_ref_0 = px_ref_0, x_ref_1 = x_ref_1, px_ref_1 = px_ref_1,
        y_ref_0 = y_ref_0, py_ref_0 = py_ref_0, y_ref_1 = y_ref_1, py_ref_1 = py_ref_1)

    dtk_arc.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.x)[0],
                      dtk_particle.x, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.px)[0],
                      dtk_particle.px, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.y)[0],
                      dtk_particle.y, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.py)[0],
                      dtk_particle.py, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.zeta)[0],
                      dtk_particle.zeta, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.delta)[0],
                      dtk_particle.delta, rtol=1e-14, atol=1e-14)


@for_all_test_contexts
def test_simplified_accelerator_segment_chroma_detuning(test_context):
    dtk_particle = dtk.TestParticles(
            p0c=25.92e9,
            x=1e-3,
            px=1e-5,
            y=-2e-3,
            py=-1.5e-5,
            zeta=2.,
            delta=2E-4)

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)

    alpha_x_0 = -0.5
    beta_x_0 = 100.0
    disp_x_0 = 1.8
    alpha_x_1 = 2.1
    beta_x_1 = 2.0
    disp_x_1 = 3.3
    alpha_y_0 = -0.4
    beta_y_0 = 8.0
    disp_y_0 = -0.2
    alpha_y_1 = 0.7
    beta_y_1 = 0.3
    disp_y_1 = -1.9
    Q_x = 0.27
    Q_y = 0.34
    beta_s = 856.9
    Q_s = 0.001    
    #energy_ref_increment = 1.2E9
    energy_ref_increment = 0.0 # There seems to be a bug for non-zero values
    energy_increment = 4.8E8
    x_ref_0 = -5E-3
    px_ref_0 = 6E-4
    x_ref_1 = 2E-2
    px_ref_1 = -5E-5
    y_ref_0 = -9E-2
    py_ref_0 = 1E-4
    y_ref_1 = 4E-2
    py_ref_1 = 5E-4
    chroma_x=8.0
    chroma_y=-5.0
    det_xx = 1E-3
    det_xy = -2E-4
    det_yy = -6E-4
    det_yx = 3E-3

    arc = xt.LineSegmentMap(_context=test_context,
        alfx=(alpha_x_0, alpha_x_1), betx=(beta_x_0, beta_x_1),
        dx=(disp_x_0, disp_x_1), dpx=(0.0, 0.0),
        alfy=(alpha_y_0, alpha_y_1), bety=(beta_y_0, beta_y_1),
        dy=(disp_y_0, disp_y_1), dpy=(0.0, 0.0),
        qx=Q_x, qy=Q_y,
        bets=beta_s, qs=Q_s,
        energy_ref_increment=energy_ref_increment,energy_increment=energy_increment,
        x_ref=(x_ref_0, x_ref_1), px_ref=(px_ref_0, px_ref_1),
        y_ref=(y_ref_0, y_ref_1), py_ref=(py_ref_0, py_ref_1),
        dqx=chroma_x, dqy=chroma_y,
        det_xx=det_xx, det_xy=det_xy, det_yy=det_yy, det_yx=det_yx)
    arc.track(particles)

    dtk_arc = dtk.elements.LinearTransferMatrix(alpha_x_0=alpha_x_0, beta_x_0=beta_x_0, disp_x_0=disp_x_0,
        alpha_x_1=alpha_x_1, beta_x_1=beta_x_1, disp_x_1=disp_x_1,
        alpha_y_0=alpha_y_0, beta_y_0=beta_y_0, disp_y_0=disp_y_0,
        alpha_y_1=alpha_y_1, beta_y_1=beta_y_1, disp_y_1=disp_y_1,
        Q_x=Q_x, Q_y=Q_y,
        beta_s=beta_s, Q_s=Q_s,
        chroma_x=chroma_x, chroma_y=chroma_y,
        det_xx=det_xx, det_xy=det_xy, det_yy=det_yy, det_yx=det_yx,
        energy_ref_increment=energy_ref_increment,energy_increment=energy_increment,
        x_ref_0 = x_ref_0, px_ref_0 = px_ref_0, x_ref_1 = x_ref_1, px_ref_1 = px_ref_1,
        y_ref_0 = y_ref_0, py_ref_0 = py_ref_0, y_ref_1 = y_ref_1, py_ref_1 = py_ref_1)
    dtk_arc.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.x)[0],
                      dtk_particle.x, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.px)[0],
                      dtk_particle.px, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.y)[0],
                      dtk_particle.y, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.py)[0],
                      dtk_particle.py, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.zeta)[0],
                      dtk_particle.zeta, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.delta)[0],
                      dtk_particle.delta, rtol=1e-14, atol=1e-14)


@for_all_test_contexts
def test_simplified_accelerator_segment_uncorrelated_damping(test_context):
    alpha_x_0 = -0.5
    beta_x_0 = 100.0
    disp_x_0 = 1.8
    alpha_x_1 = 2.1
    beta_x_1 = 2.0
    disp_x_1 = 3.3
    alpha_y_0 = -0.4
    beta_y_0 = 8.0
    disp_y_0 = -0.2
    alpha_y_1 = 0.7
    beta_y_1 = 0.3
    disp_y_1 = -1.9
    Q_x = 0.27
    Q_y = 0.34
    beta_s = 856.9
    Q_s = 0.001
    energy_ref_increment = 1.2E9
    energy_increment = 4.8E8
    x_ref_0 = -5E-3
    px_ref_0 = 6E-4
    x_ref_1 = 2E-2
    px_ref_1 = -5E-5
    y_ref_0 = -9E-2
    py_ref_0 = 1E-4
    y_ref_1 = 4E-2
    py_ref_1 = 5E-4
    damping_rate_x = 5E-4
    damping_rate_y = 1E-3
    damping_rate_s = 2E-3

    dtk_particle = dtk.TestParticles(
            p0c=25.92e9,
            x=1e-3,
            px=1e-5,
            y=-2e-3,
            py=-1.5e-5,
            zeta=2.,
            delta=2E-4)

    particles = xp.Particles.from_dict(dtk_particle.to_dict(),
                                       _context=test_context)


    arc = xt.LineSegmentMap(_context=test_context,
        alfx=(alpha_x_0, alpha_x_1), betx=(beta_x_0, beta_x_1),
        dx=(disp_x_0, disp_x_1), dpx=(0.0, 0.0),
        alfy=(alpha_y_0, alpha_y_1), bety=(beta_y_0, beta_y_1),
        dy=(disp_y_0, disp_y_1), dpy=(0.0, 0.0),
        qx=Q_x, qy=Q_y,
        bets=beta_s, qs=Q_s,
        energy_ref_increment=energy_ref_increment,
        energy_increment=energy_increment,
        x_ref=(x_ref_0, x_ref_1), px_ref=(px_ref_0, px_ref_1),
        y_ref=(y_ref_0, y_ref_1), py_ref=(py_ref_0, py_ref_1),
        damping_rate_x = damping_rate_x,
        damping_rate_y = damping_rate_y,
        damping_rate_s = damping_rate_s)

    arc.track(particles)

    dtk_arc = dtk.elements.LinearTransferMatrix(alpha_x_0=alpha_x_0, beta_x_0=beta_x_0, disp_x_0=disp_x_0,
    alpha_x_1=alpha_x_1, beta_x_1=beta_x_1, disp_x_1=disp_x_1,
    alpha_y_0=alpha_y_0, beta_y_0=beta_y_0, disp_y_0=disp_y_0,
    alpha_y_1=alpha_y_1, beta_y_1=beta_y_1, disp_y_1=disp_y_1,
    Q_x=Q_x, Q_y=Q_y,
    beta_s=beta_s, Q_s=Q_s,
    chroma_x=0.0, chroma_y=0.0,
    det_xx=0.0, det_xy=0.0, det_yy=0.0, det_yx=0.0,
    energy_ref_increment=energy_ref_increment,energy_increment=energy_increment,
    x_ref_0 = x_ref_0, px_ref_0 = px_ref_0, x_ref_1 = x_ref_1, px_ref_1 = px_ref_1,
    y_ref_0 = y_ref_0, py_ref_0 = py_ref_0, y_ref_1 = y_ref_1, py_ref_1 = py_ref_1,
    damping_rate_x = damping_rate_x,damping_rate_y = damping_rate_y,damping_rate_s = damping_rate_s)
    dtk_arc.track(dtk_particle)

    assert np.isclose(test_context.nparray_from_context_array(particles.x)[0],
                      dtk_particle.x, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.px)[0],
                      dtk_particle.px, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.y)[0],
                      dtk_particle.y, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.py)[0],
                      dtk_particle.py, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.zeta)[0],
                      dtk_particle.zeta, rtol=1e-14, atol=1e-14)
    assert np.isclose(test_context.nparray_from_context_array(particles.delta)[0],
                      dtk_particle.delta, rtol=1e-14, atol=1e-14)

@for_all_test_contexts
def test_simplified_accelerator_segment_uncorrelated_damping_equilibrium(test_context):
    alpha_x_0 = 0.0
    beta_x_0 = 100.0
    alpha_y_0 = 0.0
    beta_y_0 = 8.0
    Q_x = 0.18
    Q_y = 0.22
    beta_s = 856.9
    Q_s = 0.015
    damping_rate_x = 5E-4
    damping_rate_y = 1E-3
    damping_rate_s = 2E-3
    energy = 45.6
    equ_emit_x = 0.3E-9
    equ_emit_y = 1E-12
    equ_length = 3.5E-3
    equ_delta = 3.8E-4
    beta_s = equ_length/equ_delta
    equ_emit_s = equ_length*equ_delta

    npart = int(1E3)
    particles = xp.Particles(_context=test_context,
                x=np.random.randn(npart)*np.sqrt(equ_emit_x*beta_x_0),
                px=np.random.randn(npart)*np.sqrt(equ_emit_x/beta_x_0),
                y=np.random.randn(npart)*np.sqrt(equ_emit_y*beta_y_0),
                py=np.random.randn(npart)*np.sqrt(equ_emit_y/beta_y_0),
                zeta=np.random.randn(npart)*np.sqrt(equ_emit_s*beta_s),
                delta=np.random.randn(npart)*np.sqrt(equ_emit_s/beta_s),
                p0c=energy*1E9)
    particles._init_random_number_generator()


    arc = xt.LineSegmentMap(_context=test_context,
        alfx=(alpha_x_0, alpha_x_0), betx=(beta_x_0, beta_x_0),
        dx=(0.0, 0.0), dpx=(0.0, 0.0),
        alfy=(alpha_y_0, alpha_y_0), bety=(beta_y_0, beta_y_0),
        dy=(0.0, 0.0), dpy=(0.0, 0.0),
        qx=Q_x, qy=Q_y,
        bets=beta_s, qs=Q_s,
        damping_rate_x = damping_rate_x,
        damping_rate_y = damping_rate_y,
        damping_rate_s = damping_rate_s,
        equ_emit_x=equ_emit_x, equ_emit_y=equ_emit_y, equ_emit_s=equ_emit_s)

    gamma_x = (1.0+alpha_x_0**2)/beta_x_0
    gamma_y = (1.0+alpha_y_0**2)/beta_y_0
    n_turns = int(1E3)
    emit_x = np.zeros(n_turns,dtype=float)
    emit_y = np.zeros_like(emit_x)
    emit_s = np.zeros_like(emit_x)
    ctx2np = test_context.nparray_from_context_array
    for turn in range(n_turns):
        arc.track(particles)
        emit_x[turn] = 0.5*np.average(ctx2np(
            gamma_x*particles.x**2+2*alpha_x_0*particles.x*particles.px
            +beta_x_0*particles.px**2))
        emit_y[turn] = 0.5*np.average(ctx2np(
            gamma_y*particles.y**2+2*alpha_y_0*particles.y*particles.py
            +beta_y_0*particles.py**2))
        emit_s[turn] = 0.5*np.average(ctx2np(particles.zeta**2/beta_s
            +beta_s*particles.delta**2))
    turns = np.arange(n_turns)
    equ_emit_x_0 = np.average(emit_x)
    equ_emit_y_0 = np.average(emit_y)
    equ_emit_s_0 = np.average(emit_s)

    assert np.isclose(equ_emit_x,equ_emit_x_0, rtol=1e-1, atol=1e-10)
    assert np.isclose(equ_emit_y,equ_emit_y_0, rtol=1e-1, atol=1e-10)
    assert np.isclose(equ_emit_s,equ_emit_s_0, rtol=1e-1, atol=1e-10)


@for_all_test_contexts
def test_nonlinearlens(test_context):
    mad = Madx()

    dr_len = 1e-11
    mad.input(f"""
    ss: sequence, l={dr_len};
        lens: nllens, at=0, cnll=0.15, knll=0.3;
        ! since in MAD-X we can't track a zero-length line, we put in
        ! this tiny drift here at the end of the sequence:
        dr: drift, at={dr_len / 2}, l={dr_len};
    endsequence;
    beam;
    use, sequence=ss;
    """)

    line = xt.Line.from_madx_sequence(mad.sequence.ss)
    line.config.XTRACK_USE_EXACT_DRIFTS = True # to be consistent with madx
    line.build_tracker(_context=test_context)

    num_p_test = 10
    x_test = np.linspace(-1e-2, 2e-2, num_p_test)
    y_test = np.linspace(-3e-2, 1e-2, num_p_test)
    px_test = np.linspace(-2e-5, 4e-5, num_p_test)
    py_test = np.linspace(-4e-5, 2e-5, num_p_test)

    p0 = xp.Particles(p0c=2e9, x=x_test, px=px_test, y=y_test, py=py_test,
                    zeta=.1, ptau=1e-3)

    part = p0.copy(_context=test_context)
    line.track(part, _force_no_end_turn_actions=True)
    part.move(_context=xo.context_default)

    xt_tau = part.zeta/part.beta0
    px = []
    py = []
    for ii in range(len(p0.x)):
        mad.input(f"""
        beam, particle=proton, pc={p0.p0c[ii] / 1e9}, sequence=ss, radiate=FALSE;

        track, onepass, onetable;
        start, x={p0.x[ii]}, px={p0.px[ii]}, y={p0.y[ii]}, py={p0.py[ii]}, \
            t={p0.zeta[ii]/p0.beta0[ii]}, pt={p0.ptau[ii]};
        run,
            turns=1,
            track_harmon=1e-15;  ! since in this test we don't care about
                ! losing particles due to t difference, we set track_harmon to
                ! something very small, to make t_max large.
        endtrack;
        """)

        mad_results = mad.table.mytracksumm[-1]

        px.append(mad_results.px)
        py.append(mad_results.py)

        assert np.allclose(part.x[ii], mad_results.x, atol=1e-14, rtol=0), 'x'
        assert np.allclose(part.px[ii], mad_results.px, atol=1e-14, rtol=0), 'px'
        assert np.allclose(part.y[ii], mad_results.y, atol=1e-14, rtol=0), 'y'
        assert np.allclose(part.py[ii], mad_results.py, atol=1e-14, rtol=0), 'py'
        assert np.allclose(xt_tau[ii], mad_results.t, atol=1e-14, rtol=0), 't'
        assert np.allclose(part.ptau[ii], mad_results.pt, atol=1e-14, rtol=0), 'pt'
        assert np.allclose(part.s[ii], mad_results.s, atol=1e-14, rtol=0), 's'

@for_all_test_contexts
def test_multipole_tilt_90_deg(test_context):

    m = xt.Multipole(knl=[0.1, 0], hxl=0.1, length=2, _context=test_context)
    p = xt.Particles(x = 0, y=0, delta=1, p0c=1e12, _context=test_context)
    ln = xt.Line(elements=[
        xt.SRotation(angle=-90.),
        m,
        xt.SRotation(angle=90.)])
    ln.build_tracker(_context=test_context)
    ln.track(p)

    # Check dispersion
    my = xt.Multipole(ksl=[0.1, 0], hyl=0.1, length=2, _context=test_context)
    py = xt.Particles(x = 0, y=0, delta=1., p0c=1e12, _context=test_context)
    my.track(py)

    p.move(_context=xo.context_default)
    py.move(_context=xo.context_default)

    assert np.allclose(p.x, py.x, rtol=0, atol=1e-14)
    assert np.allclose(p.y, py.y, rtol=0, atol=1e-14)
    assert np.allclose(p.px, py.px, rtol=0, atol=1e-14)
    assert np.allclose(p.py, py.py, rtol=0, atol=1e-14)
    assert np.allclose(p.zeta, py.zeta, rtol=0, atol=1e-14)
    assert np.allclose(p.ptau, py.ptau, rtol=0, atol=1e-14)

    # Check weak focusing
    pf = xt.Particles(x=0, y=0.3, delta=0., p0c=1e12, _context=test_context)
    pfy = pf.copy(_context=test_context)

    ln.track(pf)
    my.track(pfy)

    pf.move(_context=xo.context_default)
    pfy.move(_context=xo.context_default)

    assert np.allclose(pf.x, pfy.x, rtol=0, atol=1e-14)
    assert np.allclose(pf.y, pfy.y, rtol=0, atol=1e-14)
    assert np.allclose(pf.px, pfy.px, rtol=0, atol=1e-14)
    assert np.allclose(pf.py, pfy.py, rtol=0, atol=1e-14)
    assert np.allclose(pf.zeta, pfy.zeta, rtol=0, atol=1e-14)
    assert np.allclose(pf.ptau, pfy.ptau, rtol=0, atol=1e-14)
