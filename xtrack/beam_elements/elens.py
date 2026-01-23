# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

from ..base_element import BeamElement
import xobjects as xo

import numpy as np

class Elens(BeamElement):
    '''Beam element modeling a hollow electron lens.

    Parameters
    ----------
    inner_radius : float
        Inner radius of the electron lens in meters. Default is ``0``.
    outer_radius : float
        Outer radius of the electron lens in meters. Default is ``0``.
    current : float
        Current of the electron lens in Ampere. Default is ``0``.
    elens_length : float
        Effective length of the electron lens in meters. Default is ``0``.
    voltage : float
        Voltage of the electron lens in Volts. Default is ``0``.
    residual_kick_x : float
        Residual kick in the horizontal plane in radians. Default is ``0``.
    residual_kick_y : float
        Residual kick in the vertical plane in radians. Default is ``0``.
    polynomial_coefficients : array
        Array of coefficients of the polynomial describing the radial profile
        of the electron beam if not uniform. Default is ``[0]``.
    polynomial_order : int
        Order of the polynomial describing the radial profile
        of the electron beam if not uniform. Default is ``0``.
    e_beam_dir : int
        Direction of the electron beam with respect to the main beam, 1 if
        same direction, -1 if opposite direction. Default is ``-1``.
    chebyshev_max_order : int
        Maximum order of the Chebyshev polynomial used to apply non-linear
        residual kicks. Default is ``0``.
    chebyshev_reference_radius : float
        Reference radius used to evaluate the Chebyshev fit of the residual
        kicks. Default is ``1``.
    chebyshev_coefficients : array
        Coefficients of the Chebyshev fit of the residual kicks. Default is
        ``[0]``.
    elens_on : int
        Variable to turn on or off the elens elements to be able to simulate
        pulsed operation, ``0`` for off and ``1`` for on. Default is ``1``.
    '''

    _xofields={
        'current': xo.Float64,
        'inner_radius': xo.Field(xo.Float64, default=0.),
        'outer_radius': xo.Field(xo.Float64, default=0.),
        'elens_length': xo.Field(xo.Float64, default=0.),
        'voltage': xo.Field(xo.Float64, default=0.),
        'residual_kick_x': xo.Field(xo.Float64, default=0.),
        'residual_kick_y': xo.Field(xo.Float64, default=0.),
        'polynomial_coefficients': xo.Field(xo.Float64[:], default=[0]),
        'polynomial_order': xo.Field(xo.Float64, default=0.),
        'e_beam_dir': xo.Field(xo.Float64, default=-1.),
        'chebyshev_max_order': xo.Field(xo.Float64, default=0.),
        'chebyshev_reference_radius': xo.Field(xo.Float64, default=1.),
        'chebyshev_coefficients': xo.Field(xo.Float64[:], default=[0]),
        'elens_on': xo.Field(xo.Float64, default=1.),
    }

    has_backtrack = True

    _extra_c_sources = [
        '#include "xtrack/beam_elements/elements_src/elens.h"',
    ]

    def __init__(self, **kwargs):

        if '_xobject' in kwargs and kwargs['_xobject'] is not None:
            self.xoinitialize(**kwargs)
            return

        super().__init__(**kwargs)

        # Check that radial electron density is correctly normalised
        polynomial_order = len(self.polynomial_coefficients) - 1
        self.polynomial_order = polynomial_order

        if polynomial_order != 0:
            EMASS = 510998.928
            C_LIGHT = 299792458.0
            etot_e = self.voltage + EMASS
            p_e = np.sqrt(etot_e*etot_e - EMASS*EMASS)
            beta_e = p_e / etot_e

            e_dens_sum = 0.0
            for n in range(len(self.polynomial_coefficients)):
                e_dens_sum += self.polynomial_coefficients[n] * (
                    self.out_radius**(n+2) - self.inner_radius**(n+2)
                ) / (n + 2)

            if not np.isclose(1.0, e_dens_sum):
                raise ValueError("The provided polynomial coefficient are not " \
                "correctly normalised!")

    def get_chebyshev_coeffs_from_potential(
        chebyshev_reference_radius,
        chebyshev_max_order,
        potential_file,
        sep=None,
        skiprows=0,
        col_names=['x', 'y', 'z', 'value'],
        si_conversion=[1., 1., 1., 1.],
        lim=None
    ):
        '''Function fitting a measured or simulated hollow electron lens
        electric potential with Chebyshev polynomials that can be used for
        tracking with non-linear residual kicks.

        Input
        -----
        chebyshev_reference_radius : float
            Reference radius within which all x and y values of interest
            can be found.
        chebyshev_max_order : int
            Maximum order of the Chebyshev fit.
        potential_file : str
            Path to the .txt file in which the measured/simulated electric
            potential is stored. All the data is assumed to be stored in a
            single .txt file in long format.
        sep : str | None
            Separator/delimiter in the data file. Default is ``None``, i.e.
            white spaces.
        skiprows : int
            Number of rows to skip when reading the file to account for any
            headers. Default is ``0``.
        col_names : list[str]
            Name of the columns in the data file, used to correctly identify
            coordinates in case they are randomly ordered. Default in
            ``['x', 'y', 'z', 'value']``
        si_conversion : list[float]
            Values to multiply data with in case they are not stored in SI
            units. E.g. if the x and y values are stored in mm, the z value is
            stored in m and the potential is stored in V, and
            col_names=["x", "z", "y", "value"], then
            si_conversion=[1e-3, 1., 1e-3, 1.].
        lim : float | None
            Absolute limit in x and y above which data points are discarded.
            Default is ``None``.

        Returns
        -------
        C : array[float]
            Chebyshev fit coefficients stored in such a way that the potential
            is approximated by
            $\\Phi\\approx\\sum_{n}^N\\sum_{m=0}^nC_{m,n-m}T_m(x/a)T_{n-m}(y/a)$
            and $C_{m,n-m}$ corresponds to C[k], where k=n*(n+1)/2+m.
        '''
        # Load potential data
        x_idx = col_names.index("x")
        y_idx = col_names.index("y")
        z_idx = col_names.index("z")
        value_idx = col_names.index("value")

        data = np.loadtxt(potential_file, dtype=float, delimiter=sep, skiprows=skiprows)
        data[:, x_idx] *= si_conversion[x_idx]
        data[:, y_idx] *= si_conversion[y_idx]
        data[:, z_idx] *= si_conversion[z_idx]
        data[:, value_idx] *= si_conversion[value_idx]

        if lim:
            valid_idx = np.where(
                np.logical_and(
                    np.abs(data[:, x_idx]) < lim,
                    np.abs(data[:, y_idx]) < lim,
                )
            )[0]
            data = data[valid_idx, :]

        x = np.unique(data[:, x_idx])
        y = np.unique(data[:, y_idx])
        z = np.unique(data[:, z_idx])

        nx, ny, nz = len(x), len(y), len(z)
        phi = data[:, value_idx].reshape(nz, ny, nx)
        phi = phi.transpose(2, 1, 0)

        # Integrate potential along z
        phi_int = np.trapezoid(phi, z, axis=2)

        # Get data back to long format
        y_, x_ = np.meshgrid(y, x)
        x = x_.flatten()
        y = y_.flatten()
        Phi = phi_int.flatten()

        # Rescale coordinates
        ux = x / chebyshev_reference_radius
        uy = y / chebyshev_reference_radius

        N = chebyshev_max_order

        # Construct Chebyshev Vandermode matrix
        V = np.polynomial.chebyshev.chebvander2d(ux, uy, [N, N])

        cols = []
        index_map = []
        k = 0
        for i in range(N+1):
            for j in range(N+1):
                if i + j <= N:
                    cols.append(V[:, k])
                    index_map.append((i, j))
                k += 1
        A = np.column_stack(cols)

        # Solve the least squares problem
        coeffs, *_ = np.linalg.lstsq(A, Phi, rcond=None)

        # Store coefficients in a single 1D array
        C_2d = np.zeros((N+1, N+1))
        for k, (i, j) in enumerate(index_map):
            C_2d[i, j] = coeffs[k]
        C = []
        for n in range(N+1):
            for m in range(n+1):
                C.append(C_2d[m, n-m])

        return np.asarray(C, dtype=float)
