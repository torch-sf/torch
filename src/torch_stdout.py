#!/usr/bin/env python

from __future__ import division, print_function

from datetime import datetime

def tprint(*args, **kwargs):
    # option 1
    #print("[torch] ", end='')
    # option 2
    tstr = datetime.now().strftime("%H:%M:%S.%f")
    print("[torch {}] ".format(tstr), end='')
    return print(*args, **kwargs)
