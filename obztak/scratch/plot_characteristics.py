#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import pylab as plt
from collections import OrderedDict as odict

import obztak.utils.ortho
import obztak.field
from obztak.utils.constants import COLORS

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fields')
    args = parser.parse_args()
    filename = args.fields
    fields = obztak.field.FieldArray.read(args.fields)

    g = fields[fields['FILTER']=='g']
    r = fields[fields['FILTER']=='r']
    i = fields[fields['FILTER']=='i']
    z = fields[fields['FILTER']=='z']

    bands = odict([('g',g),('r',r),('i',i),('z',z)])

    kwargs = dict(histtype='step',bins=35,lw=2)

    plt.figure()
    for b,f in bands.items():
        plt.hist(f['AIRMASS'],color=COLORS[b],label='%s-band'%b,**kwargs)
    plt.title('AIRMASS')
    plt.xlabel('AIRMASS')
    plt.legend()
    plt.savefig(filename.replace('.csv','_airmass.png'))

    plt.figure()
    for b,f in bands.items():
        plt.hist(f['SLEW'],color=COLORS[b],log=True,label='%s-band'%b,**kwargs)
    plt.title('SLEW')
    plt.xlabel('SLEW (deg)')
    plt.legend()
    plt.savefig(filename.replace('.csv','_slew.png'))

    plt.figure()
    for b,f in bands.items():
        plt.hist(f['MOONANGLE'],color=COLORS[b],label='%s-band'%b,**kwargs)
    plt.title('MOONANGLE')
    plt.xlabel('MOONANGLE (deg)')
    plt.legend()
    plt.savefig(filename.replace('.csv','_moonangle.png'))

    # plot bliss
    obztak.utils.ortho.plot_bliss_coverage(fields)
    plt.savefig(filename.replace('.csv','_skyplot.png'))
