#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Created on 18 Mar 2015                                                     #
#                                                                            #
# Copyright: Copyright (C) 2015 Liam Deacon                                  #
#                                                                            #
# License: MIT License                                                       #
#                                                                            #
# Permission is hereby granted, free of charge, to any person obtaining a    #
# copy of this software and associated documentation files (the "Software"), #
# to deal in the Software without restriction, including without limitation  #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,   #
# and/or sell copies of the Software, and to permit persons to whom the      #
# Software is furnished to do so, subject to the following conditions:       #
#                                                                            #
# The above copyright notice and this permission notice shall be included in #
# all copies or substantial portions of the Software.                        #
#                                                                            #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    #
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING    #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER        #
# DEALINGS IN THE SOFTWARE.                                                  #
#                                                                            #
##############################################################################
"""
eeasisss.py - calculate EEASiSSS phase shifts.

eeasisss calculates phase shifts using John Rundgren's EEASiSSS() subroutine.

Examples
--------
.. code:: bash

    eeasisss.py -i inputX -a ~/atlib/ -l ilogA


"""
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import sys
import os
import argparse

from ctypes import cdll, create_string_buffer
from ctypes.util import find_library

_ext = ".dll" if str(sys.platform).startswith("win") else ".so"
_lib = os.path.join(os.path.dirname(__file__), "lib")

os.environ["PATH"] = _lib + ";" + os.environ["PATH"]
_library_path = find_library("EEASiSSS") or os.path.join(_lib, "libEEASiSSS" + _ext)
libEEASiSSS = cdll.LoadLibrary(_library_path)


def eeasisss(input_file="inputX"):
    """Wrapper function to call EEASiSSS Fortran library using ctypes"""
    libEEASiSSS.hartfock_(create_string_buffer(str(input_file)), size=255)


if __name__ == "__main__":
    # use this file as a script
    main()
