"""
Observation planning for the Magellanic Satellites Survey.
"""

try:
    from .version import __version__
except ImportError:
    from .utils.get_version import get_version
    __version__ = get_version()

from survey_simulator import Simulator
from survey_scheduler import Scheduler
