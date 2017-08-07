# Licensed under a MIT style license - see LICENSE
"""minions - Science minions execute pipeline tasks."""

from .core import *

from .calibrate import Calibrate
from .background import Background
from .cometphot import CometPhot
from .calcomet import CalComet
from .plotcometphot import PlotCometPhot
from .plotcometloc import PlotCometLoc
from .megatable import MegaTable
from .source_mask import SourceMask

frame_minions = (Calibrate, SourceMask, Background, CometPhot, CalComet, PlotCometLoc)
table_minions = (PlotCometPhot, MegaTable)

frame_minion_names = [m.name for m in frame_minions]
table_minion_names = [m.name for m in table_minions]
