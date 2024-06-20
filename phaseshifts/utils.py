#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Created on 15 Apr 2015                                                     #
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

"""
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import sys
import os
from shutil import copy


def expand_filepath(path):
    """Expands the filepath for environment and user variables"""
    return os.path.normpath(
        os.path.expanduser(os.path.expandvars(os.path.expanduser(path)))
    )


class FileUtils(object):
    """Class for performing phase shift related file operations."""

    @staticmethod
    def copy_files(files, dst, verbose=False):
        """copy list of files into destination directory"""
        # check if using native Windows Python with cygwin
        env = ""
        if str(sys.platform).startswith("win") and dst.startswith("/cygdrive"):
            if os.environ["CLEED_PHASE"] == dst:
                env = "CLEED_PHASE="
            dst = '"%s"' % (
                dst.split("/")[2] + ":" + os.path.sep.join(dst.split("/")[3:])
            )

        # do check and create directory if needed
        if os.path.isfile(dst):
            dst = os.path.dirname(dst)
        os.makedirs(dst, exist_ok=True)

        # copy each phase shift file to directory
        if verbose:
            print("\nCopying files to %s'%s'" % (env, dst))
        for filename in files:
            try:
                copy(filename, dst)
                if verbose:
                    print(os.path.basename(filename))
            except IOError:
                sys.stderr.write("Cannot copy file '%s'\n" % filename)
                sys.stderr.flush()
