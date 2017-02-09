#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import os
import glob
import json

import numpy as np

from obztak import set_survey,get_survey
from obztak.maglites import MaglitesSurvey, MaglitesScheduler, MaglitesFieldArray
from obztak.utils.testing import call, make_options
from obztak.utils import fileio

set_survey('maglites')

def test_maglites_prepare_fields():
    survey = MaglitesSurvey()
    fields = survey.prepare_fields(plot=False)

    assert len(fields) == 4478

    idx = [17,583,3389]
    test=fields[idx]

    np.testing.assert_equal(test['HEX'],[2042,4947,2207])
    np.testing.assert_equal(test['TILING'],[1,1,4])
    np.testing.assert_allclose(test['RA'],[358.03, 104.729, 342.55892077])
    np.testing.assert_allclose(test['DEC'],[-70.2678, -65.2038, -68.50739044])


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
