#!/usr/bin/env python

from __future__ import division, print_function

from sys import version_info
PYTHON3 = version_info[0] == 3
PYTHON2 = version_info[0] == 2
assert PYTHON3 or PYTHON2

if PYTHON3:
    from collections.abc import MutableMapping, Sequence
elif PYTHON2:
    from collections import MutableMapping, Sequence

import numpy as np
import os
from os import path
import re
import warnings
from warnings import warn
warnings.simplefilter('always', UserWarning)


class FlashPar(MutableMapping):
    """Flash parameter file representation"""

    def __init__(self, fname):
        """Flash parameter file representation"""
        self.fname = fname
        if not path.isfile(fname):
            raise Exception("FLASH par file does not exist: {}".format(fname))
        self.read()

    # Five abstract methods required for MutableMapping
    # (superclasses: Mapping, Collection, Sized, Iterable, Container)

    def __getitem__ (self, key):
        return self._p[key]

    def __setitem__ (self, key, value):

        self._p[key] = value

        if key not in self.line_nums:
            line = "{} = {:s}\n".format(key, self.format_write(key, value))
            self.lines.append(line)
            self.line_nums[key] = len(self.lines) - 1

        else:
            line_num = self.line_nums[key]
            line = self.lines[line_num]

            # Preserve comment, try to preserve whitespace alignment
#            if '#' in line:
#                if PYTHON3:
#                    text, comment = line.split('#', maxsplit=1)
#                elif PYTHON2:
#                    text, comment = line.split('#', 1)
#            else:
#                text, comment = line, None
#
#            splitter = re.compile(r'^(.*?)(\s*=\s*)(.*)$')
#            match = splitter.match(text)
#            assert match, "Line# {} missing key/value pair".format(line_num)
#            assert match.group(1).strip() == key, "Got wrong key"
#
#            key_str = match.group(1)
#            assign_str = match.group(2)
#            n = len(match.group(3))  # old value_str
#
#            padfmt = "{{:<{}s}}".format(str(n))
#            value_str = padfmt.format(self.format_write(key, value))
#
#            # Perform update
#            line = key_str + assign_str + value_str
#            if comment:
#                line = line + '#' + comment  # comment carries newline
#            else:
#                line = line + '\n'  # re-attach newline

            # Don't try to preserve comment, but do try to keep whitespace
            splitter = re.compile(r'^(.*?)(\s*=\s*)(.*)$')
            match = splitter.match(line)
            assert match, "Line# {} missing key/value pair".format(line_num)
            assert match.group(1).strip() == key, "Got wrong key"

            key_str = match.group(1)
            assign_str = match.group(2)
            line = key_str + assign_str + self.format_write(key, value) + '\n'

            self.lines[line_num] = line

        # Finished updating; sanity check keys
        assert set(self._p.keys()) == set(self.line_nums.keys())

    def __delitem__ (self, key):
        line_num = self.line_nums[key]
        del self._p[key]
        del self.lines[line_num]
        del self.line_nums[key]
        for k in self._p:
            if self.line_nums[k] > line_num:
                self.line_nums[k] = self.line_nums[k] - 1

    def __iter__ (self):
        return iter(self._p)

    def __len__ (self):
        return len(self._p)

    # FLASH parameter handling

    def read(self):
        """Read in current parameters"""
        p = {}
        lines = []
        line_nums = {}

        with open(self.fname, 'r') as f:
            splitter = re.compile(r'^(.*?)=(.*)')
            for line_num, line in enumerate(f):

                if PYTHON3:
                    text = (line.split('#', maxsplit=1))[0]
                elif PYTHON2:
                    text = (line.split('#', 1))[0]

                match = splitter.match(text)
                if match:
                    key = match.group(1).strip()
                    value = match.group(2).strip()
                    fvalue = self.format_read(key, value)
                    if key in p:
                        warn("Duplicate param '{}' in {}; old {}, new {}".format(key, self.fname, p[key], fvalue),
                             stacklevel=3)
                    p[key] = fvalue
                    line_nums[key] = line_num

                # keep all lines, including empty or comments,
                # so that read() followed by write() does not modify the file
                # at all
                lines.append(line)

        assert set(p.keys()) == set(line_nums.keys())

        self._p = p
        self.lines = lines
        self.line_nums = line_nums

    def format_read(self, k, v):
        """Cast key's value to correct Python type"""

        # In FLASH source, see
        # source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_read.F90
        # to understand how parameters are parsed.

        # FLASH parses only double quotes '"' and not single quotes.
        dquoted = re.compile(r'^"(.*?)"$')
        integered = re.compile(r'^[\-+]?[0-9]+$')
        dfloated = re.compile(r'^[\-+0-9.]+[dD][\-+0-9]+$')

        if v.lower() in ['t', '.true.']:
            return True
        elif v.lower() in ['f', '.false.']:
            return False
        elif dquoted.match(v):
            return v[1:-1]  # remove the double quotes
        elif integered.match(v):
            return int(v)
        elif dfloated.match(v):
            return float(v.lower().replace('d', 'e'))

        return float(v)  # this will flag bad input for us

    def format_write(self, k, v):
        """Format key's value to string for file output; no whitespace padding"""
        # accepts key k so we can choose format for specific parameters, e.g.,
        # pt_maxPerProc represents int, but pt_maxPerProc=1e9 gets stored as float.
        # cvisc represents float, but cvisc=1 gets stored as int.
        # There are so many possible parameters, though, that I think it is
        # folly to try to format all of them nicely.
        if isinstance(v, bool):
            if v:
                return ".true."
            else:
                return ".false."
        elif isinstance(v, float):
            return "{:g}".format(v)
        elif isinstance(v, int):
            return "{:d}".format(v)
        elif isinstance(v, str):
            return '"' + v + '"'
        else:
            raise Exception("unexpected type for {}={}".format(k,v))

    def write(self, out=None, clobber=False):
        if out is None:
            out = self.fname
        if not clobber and path.exists(out):
            raise IOError("File exists, not clobbering: {}".format(out))
        else:
            with open(out, 'w') as f:
                f.writelines(self.lines)
                f.flush()
                os.fsync(f)

    def show(self):
        for line in self.lines:
            print(line, end='')
