#!/usr/bin/env python
"""
Factories for generating objects
"""
import os,sys
from collections import OrderedDict as odict
import inspect
import importlib
import logging

from obztak import get_survey

SURVEYS = odict([
    (None,'Survey'),
    ('obztak','Survey'),
    ('maglites','MaglitesSurvey'),
    ('bliss','BlissSurvey'),
])

SCHEDULERS = odict([
    (None,'Scheduler'),
    ('obztak','Scheduler'),
    ('maglites','MaglitesScheduler'),
    ('bliss','BlissScheduler'),
])

FIELDS = odict([
    (None,'FieldArray'),
    ('obztak','FieldArray'),
    ('maglites','MaglitesFieldArray'),
    ('bliss','BlissFieldArray'),
])

TACTICIANS = odict([
    (None,'Tactician'),
    ('tactician','tactician'),
    ('coverage','CoverageTactician'),
    ('good','ConditionTactician'),
    ('bad','ConditionTactician'),
    ('maglites','ConditionTactician'),
    ('bliss','BlissTactician'),
])

def factory(cls, modules=None, **kwargs):
    """
    Factory for creating objects. Arguments are passed directly to the
    constructor of the chosen class.
    """
    # Format modules into a list
    if modules is None: modules = [__name__]
    elif isinstance(modules,basestring): modules = [modules]

    # Import the requested modules
    for module in modules:
        importlib.import_module(module)

    # Define a preticate for selecting class members
    def fn(member):
        return inspect.isclass(member) and member.__module__ in modules

    # Fill a dictionary of classes
    classes = odict()
    for module in modules:
        classes.update(inspect.getmembers(sys.modules[module], fn))

    # Lowercase class names
    members = odict([(k.lower(),v) for k,v in classes.items()])

    # Select class (case-insensitive)
    lower = cls.lower()
    if lower not in members.keys():
        msg = "Unrecognized class: %s"%(cls)
        raise KeyError(msg)

    # Access the requested class and build the object
    return members[lower](**kwargs)

def scheduler_factory(cls=None, **kwargs):
    modules = ['obztak.scheduler','obztak.maglites','obztak.bliss']
    if cls is None: cls = get_survey()
    cls = SCHEDULERS.get(cls,cls)
    return factory(cls, modules=modules, **kwargs)

def survey_factory(cls=None, **kwargs):
    modules = ['obztak.survey','obztak.maglites','obztak.bliss']
    if cls is None: cls = get_survey()
    cls = SURVEYS.get(cls,cls)
    return factory(cls, modules=modules, **kwargs)

def field_factory(cls=None, **kwargs):
    modules = ['obztak.field','obztak.maglites','obztak.bliss']
    if cls is None: cls = get_survey()
    cls = FIELDS.get(cls,cls)
    return factory(cls, modules=modules, **kwargs)

def tactician_factory(cls=None, **kwargs):
    """Create a tactician for a survey.

    Parameters:
    -----------
    cls : The survey or tactician to create (default: $OBZTAK_SURVEY).

    Returns:
    --------
    tac : The Tactician
    """
    modules = ['obztak.tactician','obztak.maglites','obztak.bliss']
    if cls is None: cls = get_survey()
    cls = TACTICIANS.get(cls,cls)
    return factory(cls, modules=modules, **kwargs)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
