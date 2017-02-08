#!/usr/bin/env python
"""
Testing code
"""
__author__ = "Alex Drlica-Wagner"

import subprocess
from subprocess import Popen, PIPE, check_call

def call(cmd,**kwargs):
    p = Popen(cmd,stdout=PIPE,stderr=PIPE,**kwargs)
    stdout,stderr = p.communicate()
    if p.returncode > 0:
        msg = str(cmd + '\n' + str(stdout) + '\n' + str(stderr))
        print(msg)
        raise subprocess.CalledProcessError(msg)
    print(stdout)
    return stdout

def check_dict(value,test):
    for k,v in value.items():
        val = test[k]
        if val != v:
            msg = '%s: %s (in), %s (out)'%(k,v,val)
            raise ValueError(msg)

    return test

def make_options(kwargs):
    return ' '.join(['--%s %s'%(k,v) for k,v in kwargs.items()])


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
