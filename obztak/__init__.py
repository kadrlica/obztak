"""
Bizarro observation tactician for DECam.
"""
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

def set_survey(survey):
    """Set the survey name in the environment variable $OBZTAK_SURVEY."""
    import os
    os.environ['OBZTAK_SURVEY'] = survey
    return get_survey()

def get_survey():
    """Get the survey name from the environment variable $OBZTAK_SURVEY."""
    import os
    survey = os.getenv('OBZTAK_SURVEY',None)
    if survey is None:
        msg = "$OBZTAK_SURVEY environment variable not set."
        raise Exception(msg)
    return survey.lower()

__survey__ = get_survey()
