"""
Adding some general utilities
"""

def setdefaults(kwargs,defaults):
    """ Set dictionary defaults """
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs

def isstring(obj):
    """ Python 2/3 compatible string check """
    import six
    return isinstance(obj, (six.string_types, bytes))
