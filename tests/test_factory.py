#!/usr/bin/env python
"""
Generic python script.
"""

import obztak.factory
import obztak.survey
import obztak.scheduler
import obztak.field

import obztak.maglites
import obztak.bliss

MSG = "\nOutput type: %s\nExpected type: %s"

def test_survey_factory():
    TEST = [
        ('Survey',obztak.survey.Survey),
        ('MaglitesSurvey',obztak.maglites.MaglitesSurvey),
        ('maglitessurvey',obztak.maglites.MaglitesSurvey),
        ('BlissSurvey',obztak.bliss.BlissSurvey),
    ]
    for string,cls in TEST:
        obj = obztak.factory.survey_factory(string)
        if not isinstance(obj,cls):
            raise TypeError(msg%(obj,cls))

def test_scheduler_factory():
    TEST = [
        (None,obztak.scheduler.Scheduler),
        ('Scheduler',obztak.scheduler.Scheduler),
        ('obztak',obztak.scheduler.Scheduler),
        ('MaglitesScheduler',obztak.maglites.MaglitesScheduler),
        ('maglitesscheduler',obztak.maglites.MaglitesScheduler),
        ('maglites',obztak.maglites.MaglitesScheduler),
        # UNCOMMENT ONCE BLISS FIELDS COMMITTED
        #('BlissScheduler',obztak.bliss.BlissScheduler),
        #('bliss',obztak.bliss.BlissScheduler),
    ]
    for string,cls in TEST:
        obj = obztak.factory.scheduler_factory(string)
        if not isinstance(obj,cls):
            raise TypeError(MSG%(obj,cls))

def test_field_factory():
    TEST = [
        (None,obztak.field.FieldArray),
        ('FieldArray',obztak.field.FieldArray),
        ('obztak',obztak.field.FieldArray),
        ('MaglitesFieldArray',obztak.maglites.MaglitesFieldArray),
        ('magliTEsFieldArray',obztak.maglites.MaglitesFieldArray),
        ('maglites',obztak.maglites.MaglitesFieldArray),
        ('BlissFieldArray',obztak.bliss.BlissFieldArray),
        ('bliss',obztak.bliss.BlissFieldArray),
    ]
    for string,cls in TEST:
        obj = obztak.factory.field_factory(string)
        if not isinstance(obj,cls):
            raise TypeError(MSG%(obj,cls))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
