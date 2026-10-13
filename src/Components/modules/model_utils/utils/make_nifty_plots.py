#!/usr/bin/env python3
from optparse import OptionParser
from model_utils.utils import plot_utils
import os

if __name__ == "__main__":
    parser = OptionParser(usage="Usage: %prog [options] modelname date0",
                          version='xxx' )
    (options, args) = parser.parse_args()

#  GET OMI FILE FROM INPUT ARGUMENT LIST
    if len(args) == 1:
        fpdata     = args[0]
    else:
        parser.error("must have 1 argument: fpdata filename")


    fnameh = os.path.splitext(os.path.basename(fpdata))[0]
    fname  = f"{fnameh}.png"
    plot_utils.niftyplot(fpdata,fname=fname)
