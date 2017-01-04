#!/usr/bin/env python
"""
Code related to the Magellanic Satellites Survey (MagLiteS).
"""
import os,sys
import logging
import copy

import numpy as np

from obztak.field import FieldArray, SISPI_DICT, SEP
from obztak.survey import Survey
from obztak.scheduler import Scheduler

from obztak.utils import fileio

PROGRAM = 'bliss'
PROPID  = '2017A-0260'

class Bliss(Survey):
    """ Survey sublcass for BLISS. """

    # 2017A ACTUAL
    nights = [['2017/02/07', 'second'],
              ['2017/02/08', 'second'],
              ['2017/02/13', 'second'],
              ['2017/02/14', 'second'],
              ['2017/02/19', 'second'],
              ['2017/02/20', 'second'],
              ['2017/03/07', 'first'],
              ['2017/03/17', 'second'],
              ['2017/04/13', 'second'],
              ['2017/04/14', 'second'],
              ['2017/05/02', 'first'],
              ['2017/05/30', 'first'],
              ['2017/05/31', 'first'],              
              ['2017/06/01', 'first'],
              ['2017/06/02', 'full'],
              ['2017/06/03', 'full'],
              ['2017/06/04', 'full'],
              ['2017/07/01', 'second'],
              ['2017/07/02', 'second'],
              ['2017/07/15', 'second'],

    def prepare_fields(self):
        pass

class BlissFieldArray(FieldArray):
    """ Array of BLISS fields """
    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID

    OBJECT_FMT = 'BLISS field' + SEP + ' %s'
    SEQID_FMT = 'BLISS scheduled' + SEP + ' %(DATE)s'

class BlissScheduler(Scheduler):
    
    def write(self,filename):
        fields = BlissFieldArray(self.scheduled_fields)
        fields.write(filename)
