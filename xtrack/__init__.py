# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

from .general import _pkg_root, _print

from .base_element import BeamElement
from .beam_elements import *
from .random import *
from .tracker_data import TrackerData
from .line import Line, Node, freeze_longitudinal, _temp_knobs, EnergyProgram
from .tracker import Tracker, Log
from .match import (Vary, Target, TargetList, VaryList, TargetInequality, Action,
                    TargetRelPhaseAdvance, TargetSet, START, END,
                    GreaterThan, LessThan)
from .targets import (TargetLuminosity, TargetSeparationOrthogonalToCrossing,
                      TargetSeparation)
from .twiss import TwissInit, TwissTable
from .loss_location_refinement import LossLocationRefinement
from .internal_record import (RecordIdentifier, RecordIndex, new_io_buffer,
                             start_internal_logging, stop_internal_logging)
from .pipeline import (PipelineStatus, PipelineMultiTracker, PipelineBranch,
                        PipelineManager)

from .monitors import *
from . import linear_normal_form
from .multiline import Multiline, MultiTwiss

from .mad_loader import MadLoader

from .multisetter import MultiSetter

from .footprint import Footprint, LinearRescale

# Flags and test functions
from .line import _is_drift, _behaves_like_drift, _is_aperture, _is_thick, _allow_backtrack
from .line import _lines_equal, _apertures_equal

from .slicing import Strategy, Uniform, Teapot
from .loss_location_refinement import _skip_in_loss_location_refinement

from xdeps import Table, FunctionPieceWiseLinear

from ._version import __version__

from xpart import(Particles, PROTON_MASS_EV, ELECTRON_MASS_EV,
                  enable_pyheadtail_interface, disable_pyheadtail_interface)

