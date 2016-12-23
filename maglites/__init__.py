"""
Observation planning for the Magellanic Satellites Survey.
"""

#try:
#    from .version import __version__
#except ImportError:
#    from .utils.get_version import get_version
#    __version__ = get_version()

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from maglites.scheduler import Scheduler

