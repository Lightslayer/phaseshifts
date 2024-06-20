#!/usr/bin/env python
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Created on 26 Apr 2015                                                     #
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
hf.py - calculate atomic charge densities

hf calculates atomic charge densities using a modified version of
Eric Shirley's hartfock() subroutine.

Examples
--------
.. code:: bash

    hf.py -i inputA -a ~/atlib/ -l ilogA


"""

# cspell: ignore atlib ilog hartfock

from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import sys
import os
import argparse

from ..libhartfock import hartfock  # TODO: Handle import more gracefully

__date__ = "2015-04-26"
__updated__ = "2015-04-26"

DEBUG = 0
TEST_RUN = 0
PROFILE = 0


def main(argv=None):
    """Command line options to hf"""
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # display help if no arguments
    if len(argv) == 1:
        argv.append("--help")

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __date__
    program_build_date = str(__updated__)
    program_version_message = "%%(prog)s %s (%s)" % (
        program_version,
        program_build_date,
    )
    program_shortdesc = __import__("__main__").__doc__.split("\n")[1]
    program_license = """%s

        Please contact Eric Shirley <eric.shirley@nist.gov> for queries
        and comments.

    usage:-
    """ % (
        program_shortdesc
    )

    return_code = 0
    try:
        # Setup argument parser
        parser = argparse.ArgumentParser(
            description=program_license,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        parser.add_argument(
            "-i",
            "--input",
            dest="input",
            metavar="<inputA>",
            help="path to the atomic orbital input file which "
            "is an appended list of Barbieri/Van Hove atorb "
            "files (with some modifications).",
            required=True,
        )
        parser.add_argument(
            "-l",
            "--log",
            dest="log",
            metavar="<log_file>",
            default=None,
            help="Optionally redirects the "
            "hartfock() calculation log to the file specified "
            "by <log_file>. The default is to print to stdout "
            " if the argument is None. [default: %(default)s]",
        )
        parser.add_argument(
            "-a",
            "--atom_dir",
            dest="chgden_dir",
            metavar="<chgden_dir>",
            default=(
                os.path.expandvars("ATLIB")
                if os.path.expandvars("ATLIB") != ""
                else os.path.expanduser("~/atlib")
            ),
            help="Specifies the output directory to place "
            "the calculated atomic charge density files into. "
            "If no argument is given, the 'ATLIB' environment "
            "variable is used as the output directory, else "
            "'~/atlib' is used. Use '.' to place the output "
            "files into the current working directory. "
            "[default: %(default)s]",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            dest="verbose",
            action="count",
            help="Set verbosity level. [default: %(default)s]",
        )
        parser.add_argument(
            "-V", "--version", action="version", version=program_version_message
        )

        # Process arguments
        args, unknown = parser.parse_known_args()

        verbose = int(getattr(args, "verbose", False))

        if verbose > 0 and len(unknown) > 0:
            for arg in unknown:
                sys.stderr.write("hf - warning: Unknown option '%s'\n" % arg)
            sys.stderr.flush()

        args.log = "" if args.log is None else args.log

        if not os.path.isfile(args.input):
            sys.stderr.write("hf - error: input file '%s' not found\n" % args.input)
            sys.stderr.flush()

        args.chgden_dir = os.path.expanduser(os.path.expandvars(args.chgden_dir))
        if not os.path.isdir(args.chgden_dir):
            if verbose:
                sys.stderr.write(
                    "hf - warning: '%s' does not exist. "
                    "Creating directory..." % args.chgden_dir
                )
                sys.stderr.flush()
            os.makedirs(args.chgden_dir)

        # call underlying hartfock() subroutine to create charge density files
        hartfock(args.input, args.log, args.chgden_dir)

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        pass
    except Exception as err:  # pylint: disable=broad-except
        if DEBUG or TEST_RUN:
            raise (err)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(err) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return_code = 2

    return return_code


if __name__ == "__main__":
    main()
