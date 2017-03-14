#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import obztak.utils.ortho
import obztak.field

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fields')
    parser.add_argument('-i','--interact')
    args = parser.parse_args()
    fields = obztak.field.FieldArray.read(args.fields)
    obztak.utils.ortho.plot_bliss_coverage(fields)

    if args.interact:
        import pylab as plt
        plt.show()

