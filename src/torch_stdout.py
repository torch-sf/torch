#!/usr/bin/env python

from __future__ import division, print_function

from datetime import datetime

def tprint(*args, **kwargs):
    #tstr = datetime.now().strftime("%H:%M:%S.%f")
    tstr = (datetime.now().strftime("%m-%d-%Y %H:%M:%S.%f"))[:-3]
    print("[torch {}] ".format(tstr), end='',flush=True)
    return print(*args, **kwargs)
