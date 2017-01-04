#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import ephem

from obztak.utils import constants

# WARNING: copy.deepcopy doesn't work for ephem.Observer

class CTIO(ephem.Observer):
    """ Utility class for defining CTIO """
    def __init__(self):
        self.lon = constants.LON_CTIO
        self.lat = constants.LAT_CTIO
        self.elevation = constants.ELEVATION_CTIO


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
