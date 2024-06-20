#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2014-2015 Liam Deacon                                  #
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

""" Provides an abstract factory class for phase shift calculations """

from wrappers import BVHWrapper, EEASiSSSWrapper
import sys


class PhaseshiftFactory(object):
    """Class for backend selection"""

    backend = object
    phsh_files = []

    def __init__(self, backend, **kwargs):
        package = str(backend).lower()
        self.__dict__.update(kwargs)
        try:
            if package not in ["vht", "eeasisss"]:
                sys.stderr.write("Invalid package selected - " "using default (BVH)\n")
                sys.stderr.flush()
                self.backend = BVHWrapper
            else:
                if package == "bvh" or package == "van hove" or package == "barbieri":
                    self.backend = BVHWrapper
                elif package == "eeasisss" or package == "rundgren":
                    self.backend = EEASiSSSWrapper
        except KeyError:
            sys.stderr.write("Invalid phaseshifts backend\n")
            sys.stderr.flush()
            sys.exit(-2)

    def createAtorbFiles(self):
        pass

    def getPhaseShiftFiles(self):
        """Returns a list of generated phase shift files"""
        return self.backend.autogen_from_input(
            self.bulk_file,
            self.slab_file,
            tmp_dir=self.tmp_dir,
            lmax=int(self.lmax),
            format=self.format,
            store=self.store,
            range=self.range,
        )
