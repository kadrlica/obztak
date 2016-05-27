#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import logging
import argparse
import dateutil
        
from maglites import __version__

class SpecialFormatter(logging.Formatter):
    """
    Class for overloading log formatting based on level.
    """
    FORMATS = {'DEFAULT'       : "%(message)s",
               logging.WARNING : "WARNING: %(message)s",
               logging.ERROR   : "ERROR: %(message)s",
               logging.DEBUG   : "DEBUG: %(message)s"}
 
    def format(self, record):
        self._fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        return logging.Formatter.format(self, record)

class VerboseAction(argparse._StoreTrueAction):
    """
    Class for setting logging level from verbosity.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        super(VerboseAction,self).__call__(parser, namespace, values, option_string)
        if self.const: logging.getLogger().setLevel(logging.DEBUG)

class DatetimeAction(argparse.Action):
    """
    Class for setting logging level from verbosity.
    """
    def __call__(self, parser, namespace, values, option_string=None):
        datetime = dateutil.parser.parse(values)
        setattr(namespace, self.dest, datetime)


class Parser(argparse.ArgumentParser):
    def __init__(self,*args,**kwargs):
        super(Parser,self).__init__(*args,**kwargs)
        self.add_argument('-v','--verbose',action=VerboseAction,
                          help='output verbosity')
        self.add_argument('--version', action='version',
                          version='maglites v'+__version__,
                          help="print version number and exit")

    def remove_argument(self,option_string):
        for i,action in enumerate(self._actions):
            if option_string in action.option_strings:
                self._handle_conflict_resolve(None, [(option_string,action)])
                #action.container._remove_action(action)

logger = logging.getLogger()
handler = logging.StreamHandler()
handler.setFormatter(SpecialFormatter())
logger.addHandler(handler)
logger.setLevel(logging.INFO)

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = Parser()
    args = parser.parse_args()
