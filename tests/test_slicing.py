# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2023.                 #
# ######################################### #
import numpy as np
import pytest
import math

import xtrack as xt
from xtrack.slicing import Strategy, Teapot, Uniform


def test_slicing_uniform():
    # Test for one slice
    slicing_3 = Uniform(1)
    assert slicing_3.element_weights() == [1.0]
    assert slicing_3.drift_weights() == [0.5] * 2
    assert [w for w in slicing_3] == [(0.5, True), (1.0, False), (0.5, True)]

    # Test for three slices
    slicing_3 = Uniform(3)
    assert slicing_3.element_weights() == [1/3] * 3
    assert slicing_3.drift_weights() == [0.25] * 4

    elem_info, drift_info = (1./3., False), (0.25, True)
    assert [w for w in slicing_3] == [
        drift_info, elem_info,
        drift_info, elem_info,
        drift_info, elem_info,
        drift_info,
    ]

    # Test error handling
    with pytest.raises(ValueError):
        Uniform(0)


def test_slicing_teapot():
    # Test for one slice
    slicing_3 = Teapot(1)
    assert slicing_3.element_weights() == [1.0]
    assert slicing_3.drift_weights() == [0.5] * 2
    assert [w for w in slicing_3] == [(0.5, True), (1.0, False), (0.5, True)]

    # Test for three slices
    slicing_3 = Teapot(3)
    assert slicing_3.element_weights() == [1/3] * 3
    assert slicing_3.drift_weights() == [0.125, 0.375, 0.375, 0.125]

    elem_info = (1./3., False)
    assert [w for w in slicing_3] == [
        (0.125, True), elem_info,
        (0.375, True), elem_info,
        (0.375, True), elem_info,
        (0.125, True),
    ]

    # Test error handling
    with pytest.raises(ValueError):
        Teapot(0)


def test_slicing_teapot_mode_thick():
    # Test for two slices
    slicing_3 = Teapot(2, mode='thick')
    assert slicing_3.drift_weights() == [0.5] * 2
    assert [w for w in slicing_3] == [(0.5, True), (0.5, True)]

    # Test for four slices
    slicing_3 = Teapot(4, mode='thick')
    assert slicing_3.drift_weights() == [0.125, 0.375, 0.375, 0.125]

    assert [w for w in slicing_3] == [
        (0.125, True),
        (0.375, True),
        (0.375, True),
        (0.125, True),
    ]

    # Test error handling
    with pytest.raises(ValueError):
        Teapot(0)


def test_slicing_strategy_matching():
    elements = [
        ('keep_this', xt.Quadrupole(length=1.0)),
        ('mb10', xt.Bend(length=1.0)),
        ('keep_drifts', xt.Drift(length=1.0)),
        ('mb11', xt.Quadrupole(length=1.0)),
        ('mq10', xt.Quadrupole(length=1.0)),
        ('something', xt.Bend(length=1.0)),
        ('mb20', xt.Bend(length=1.0)),
        ('keep_thin', xt.Multipole(length=1.0)),
        ('mb21', xt.Quadrupole(length=1.0)),
    ]

    slicing_strategies = [
        # Default: one slice
        Strategy(slicing=Uniform(1)),
        # All bends: two slices
        Strategy(slicing=Teapot(2), element_type=xt.Bend),
        # All CFDs: three slices
        Strategy(slicing=Uniform(3), element_type=xt.Quadrupole),
        # If the name starts with mb: five slices (the bend and the cfd 'mb11')
        Strategy(slicing=Teapot(5), name=r'mb1.*'),
        # If the name starts with some: four slices (the bend 'something')
        Strategy(slicing=Uniform(4), name=r'some.*'),
        # Keep the CFD 'keep'
        Strategy(slicing=None, name=r'keep'),
        # mq10: one slice
        Strategy(slicing=Teapot(1), name=r'mq10'),
    ]

    line = xt.Line(
        elements=dict(elements),
        element_names=[name for name, _ in elements],
    )

    line.slice_thick_elements(slicing_strategies)

    # Check that the slices are as expected:
    expected_names = [
        # Kept Quadrupole:
        'keep_this',
        # 5 slices for mb10:
        'mb10_entry',  # Marker
        'drift_mb10..0', 'mb10..0', 'drift_mb10..1', 'mb10..1', 'drift_mb10..2',
        'mb10..2', 'drift_mb10..3', 'mb10..3', 'drift_mb10..4', 'mb10..4',
        'drift_mb10..5',
        'mb10_exit',  # Marker
        # Kept drift:
        'keep_drifts',
        # Five slices for mb11:
        'mb11_entry',  # Marker
        'drift_mb11..0', 'mb11..0', 'drift_mb11..1', 'mb11..1', 'drift_mb11..2',
        'mb11..2', 'drift_mb11..3', 'mb11..3', 'drift_mb11..4', 'mb11..4',
        'drift_mb11..5',
        'mb11_exit',  # Marker
        # One slice for 'mq10':
        'mq10_entry',  # Marker
        'drift_mq10..0', 'mq10..0', 'drift_mq10..1',
        'mq10_exit',  # Marker
        # Four slices for 'something':
        'something_entry',  # Marker
        'drift_something..0', 'something..0', 'drift_something..1',
        'something..1', 'drift_something..2', 'something..2',
        'drift_something..3', 'something..3', 'drift_something..4',
        'something_exit',  # Marker
        # Two slices for 'mb20':
        'mb20_entry',  # Marker
        'drift_mb20..0', 'mb20..0', 'drift_mb20..1', 'mb20..1', 'drift_mb20..2',
        'mb20_exit',  # Marker
        # Keep thin:
        'keep_thin',
        # Three slices for 'mb21' (it's a Quadrupole!):
        'mb21_entry',  # Marker
        'drift_mb21..0', 'mb21..0', 'drift_mb21..1', 'mb21..1', 'drift_mb21..2',
        'mb21..2', 'drift_mb21..3',
        'mb21_exit',  # Marker
    ]
    assert line.element_names == expected_names

    # Check types:
    for name, element in line.element_dict.items():
        if name == 'keep_this':
            assert isinstance(element, xt.Quadrupole)
        elif name == 'keep_drifts' or name.startswith('drift_'):
            assert isinstance(element, xt.Drift)
        elif name == 'keep_thin':
            assert isinstance(element, xt.Multipole)
        elif name[-3:-1] == '..':
            assert isinstance(element, xt.Multipole)
        else:
            assert isinstance(element, xt.Marker)

    # Check the right scheme was used:
    def _lengths_of_drifts(name):
        return [
            line[nn].length for nn in line.element_names
            if nn.startswith(f'drift_{name}')
        ]

    # Check markers
    assert isinstance(line['mb10_entry'], xt.Marker)
    assert isinstance(line['mb11_entry'], xt.Marker)
    assert isinstance(line['mq10_entry'], xt.Marker)
    assert isinstance(line['something_entry'], xt.Marker)
    assert isinstance(line['mb20_entry'], xt.Marker)
    assert isinstance(line['mb21_entry'], xt.Marker)
    assert isinstance(line['mb10_exit'], xt.Marker)
    assert isinstance(line['mb11_exit'], xt.Marker)
    assert isinstance(line['mq10_exit'], xt.Marker)
    assert isinstance(line['something_exit'], xt.Marker)
    assert isinstance(line['mb20_exit'], xt.Marker)
    assert isinstance(line['mb21_exit'], xt.Marker)

    # Teapot
    expected_mb20_drift_lens = [1/6, 2/3, 1/6]
    assert _lengths_of_drifts('mb20') == expected_mb20_drift_lens

    expected_mb10_drift_lens = [1/12, 5/24, 5/24, 5/24, 5/24, 1/12]
    assert _lengths_of_drifts('mb10') == expected_mb10_drift_lens

    expected_mb11_drift_lens = [1/12, 5/24, 5/24, 5/24, 5/24, 1/12]
    assert _lengths_of_drifts('mb11') == expected_mb11_drift_lens

    # Uniform
    expected_mq10_drift_lens = [1/2] * 2
    assert _lengths_of_drifts('mq10') == expected_mq10_drift_lens

    expected_mb21_drift_lens = [1/4] * 4
    assert _lengths_of_drifts('mb21') == expected_mb21_drift_lens

    expected_something_drift_lens = [1/5] * 5
    assert _lengths_of_drifts('something') == expected_something_drift_lens

    # Test accessing compound elements
    assert line.get_compound_subsequence('mb10') == [
        'mb10_entry', 'drift_mb10..0', 'mb10..0', 'drift_mb10..1', 'mb10..1',
        'drift_mb10..2', 'mb10..2', 'drift_mb10..3', 'mb10..3', 'drift_mb10..4',
        'mb10..4', 'drift_mb10..5', 'mb10_exit',
    ]
    assert line.get_compound_subsequence('mb11') == [
        'mb11_entry', 'drift_mb11..0', 'mb11..0', 'drift_mb11..1', 'mb11..1',
        'drift_mb11..2', 'mb11..2', 'drift_mb11..3', 'mb11..3', 'drift_mb11..4',
        'mb11..4', 'drift_mb11..5', 'mb11_exit',
    ]
    assert line.get_compound_subsequence('mq10') == [
        'mq10_entry', 'drift_mq10..0', 'mq10..0', 'drift_mq10..1', 'mq10_exit',
    ]
    assert line.get_compound_subsequence('something') == [
        'something_entry', 'drift_something..0', 'something..0',
        'drift_something..1', 'something..1', 'drift_something..2',
        'something..2', 'drift_something..3', 'something..3',
        'drift_something..4', 'something_exit',
    ]
    assert line.get_compound_subsequence('mb20') == [
        'mb20_entry', 'drift_mb20..0', 'mb20..0', 'drift_mb20..1', 'mb20..1',
        'drift_mb20..2', 'mb20_exit',
    ]
    assert line.get_compound_subsequence('mb21') == [
        'mb21_entry', 'drift_mb21..0', 'mb21..0', 'drift_mb21..1', 'mb21..1',
        'drift_mb21..2', 'mb21..2', 'drift_mb21..3', 'mb21_exit',
    ]
    assert list(line.compound_container.compound_names) == [
        'mb10', 'mb11', 'mq10', 'something', 'mb20', 'mb21',
    ]


def test_slicing_thick_bend_simple():

    additional_kwargs = {}
    additional_kwargs['k1'] = 0.2

    bend = xt.Bend(
        length=3.0,
        k0=0.1,
        h=0.2,
        **additional_kwargs
    )
    line = xt.Line(elements=[bend], element_names=['bend'])
    line.slice_thick_elements([Strategy(slicing=Teapot(2))])

    assert len(line) == 7  # 2 markers + 2 slices + 3 drifts

    assert line['drift_bend..0'].length == 3.0 * 1/6
    assert line['drift_bend..1'].length == 3.0 * 2/3
    assert line['drift_bend..2'].length == 3.0 * 1/6

    bend0, bend1 = line['bend..0'], line['bend..1']
    assert bend0.length == bend1.length == 1.5

    expected_knl = [0.15, 0.3, 0, 0, 0]
    assert np.allclose(bend0.knl, expected_knl, atol=1e-16)
    assert np.allclose(bend1.knl, expected_knl, atol=1e-16)

    expected_hxl = 0.3
    assert np.allclose(bend0.hxl, expected_hxl, atol=1e-16)
    assert np.allclose(bend1.hxl, expected_hxl, atol=1e-16)

    # Make sure the order and the inverse factorial make sense:
    _fact = math.factorial
    assert np.isclose(_fact(bend0.order) * bend0.inv_factorial_order, 1, atol=1e-16)
    assert np.isclose(_fact(bend1.order) * bend0.inv_factorial_order, 1, atol=1e-16)

    # All else is zero:
    assert np.allclose(bend0.ksl, 0, atol=1e-16)
    assert np.allclose(bend0.hyl, 0, atol=1e-16)


def test_slicing_thick_bend_into_thick_bends_simple():

    additional_kwargs = {}
    additional_kwargs['k0'] = 0.1
    additional_kwargs['h'] = 0.2
    additional_kwargs['k1'] = 0.2

    bend = xt.Bend(
        length=3.0,
        knl=[0.1, 0.2, 0.3, 0.4, 0.5],
        ksl=[0.7, 0.6, 0.5, 0.4, 0.3],
        **additional_kwargs
    )
    line = xt.Line(elements=[bend], element_names=['bend'])
    line.slice_thick_elements([Strategy(slicing=Uniform(2, mode='thick'))])

    assert len(line) == 4  # 2 markers + 2 thick slices

    bend0, bend1 = line['bend..0'], line['bend..1']
    assert bend0.length == bend1.length == 1.5

    assert bend0.k0 == bend1.k0 == 0.1
    assert bend0.h == bend1.h == 0.2
    assert bend0.k1 == bend1.k1 == 0.2

    expected_knl = np.array([0.1, 0.2, 0.3, 0.4, 0.5]) / 2
    assert np.allclose(bend0.knl, expected_knl, atol=1e-16)
    assert np.allclose(bend1.knl, expected_knl, atol=1e-16)

    expected_ksl = np.array([0.7, 0.6, 0.5, 0.4, 0.3]) / 2
    assert np.allclose(bend0.ksl, expected_ksl, atol=1e-16)
    assert np.allclose(bend1.ksl, expected_ksl, atol=1e-16)

    # Make sure the order and the inverse factorial make sense:
    _fact = math.factorial
    assert np.isclose(_fact(bend0.order) * bend0.inv_factorial_order, 1, atol=1e-16)
    assert np.isclose(_fact(bend1.order) * bend0.inv_factorial_order, 1, atol=1e-16)


def test_slicing_xdeps_consistency():
    num_elements = 50000
    num_slices = 1

    line = xt.Line(
        elements=[xt.Bend(k0=1, length=100)] * num_elements,
        element_names=[f'bend{ii}' for ii in range(num_elements)],
    )
    line._init_var_management()

    for ii in range(num_elements):
        line.vars[f'k{ii}'] = 1
        line.element_refs[f'bend{ii}'].k0 = line.vars[f'k{ii}']

    sgy = xt.slicing.Strategy(
        element_type=xt.Bend,
        slicing=xt.slicing.Uniform(num_slices),
    )
    line.slice_thick_elements([sgy])
    assert len(line.to_dict()['_var_manager']) == num_elements * num_slices
