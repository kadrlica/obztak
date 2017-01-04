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

from obztak.utils import constants
from obztak.utils.date import datestring
from obztak.utils import fileio

PROGRAM = 'bliss'
PROPID  = '2017A-0260'

class Bliss(Survey):
    """ Survey sublcass for BLISS. """

    def prepare_fields(self):
        pass

class BlissFieldArray(FieldArray):
    """ Array of BLISS fields """
    SISPI_DICT = copy.deepcopy(SISPI_DICT)
    SISPI_DICT["program"] = PROGRAM
    SISPI_DICT["propid"] = PROPID

    OBJECT_FMT = 'BLISS field' + SEP + ' %s'
    SEQID_FMT = 'BLISS scheduled' + SEP + ' %(DATE)s'
