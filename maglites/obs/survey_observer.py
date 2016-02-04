#!/usr/bin/env python
"""
Module for running survey operations.
"""
from maglites.sim import Simulator
from datetime import datetime

class Observer(Simulator):

    def __init__(self,*args,**kwargs):

    def run(self, date=None, plot=False):

        if date is None: date = datetime.now().isoformat()
        

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
