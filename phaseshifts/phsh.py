#!/usr/local/bin/python2.7
# encoding: utf-8

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2013-2014 Liam Deacon                             #
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
phsh.py - quickly generate phase shifts

phsh provides convenience functions to create phase shifts files
suitable for input into LEED-IV programs such as SATLEED and CLEED.

Examples
--------
.. code:: bash

   phsh.py -i *.inp -b *.bul -f CLEED -S phase_dir


"""
from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import argparse
import datetime
import subprocess
import sys
import os
import platform
import tempfile
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from glob import glob
from shutil import copy

import phaseshifts
import phaseshifts.settings
from phaseshifts import model, atorb
from phaseshifts.conphas import Conphas
from phaseshifts.leed import Converter, CLEED_validator, CSearch

try:
    from phaseshifts.lib.libphsh import phsh_rel, phsh_wil, phsh_cav  # type: ignore [import-untyped]
except ImportError:
    # TODO: Create pure-python implementations of the fortran libphsh.f code
    def function_not_implemented(*args, **kwargs):
        """Simple placeholder function, which raises a NotImplementedError."""
        raise NotImplementedError

    phsh_rel = phsh_wil = phsh_cav = function_not_implemented  # type: ignore

__all__ = []


def required_length(nmin, nmax):
    """custom action to check range"""

    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                msg = 'argument "{arg}" requires between {nmin} and {nmax} arguments'.format(
                    arg=self.dest, nmin=nmin, nmax=nmax
                )
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)

    return RequiredLength


class Wrapper(object):
    """Wrapper class to easily generate phase shifts"""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @staticmethod
    def autogen_from_input(
        bulk_file, slab_file, tmp_dir=None, model_name=None, **kwargs
    ):
        """
        Description
        -----------
        Generate phase shifts from a slab/cluster input file.

        Parameters
        ----------
        slab_file : str
            Path to the cluster slab MTZ input file.
        bulk_file : str
            Path to the cluster bulk MTZ input file.
        tmp_dir : str
            Temporary directory for intermediate files.
        store : bool or int
            Specify whether to keep generated files.
        format : str
            Specify formatting of generated phase shift files
        range : tuple(float, float, float)
            Specify the energy of the start, stop and step in eV.
        model_name : str
            Name of model.

        Returns
        -------
        output_files : list(str)
           A list of phase shift output filenames
        """
        dummycell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        if model_name == None:
            model_name = "atomic"

        # check formatting
        if "format" in kwargs:
            out_format = kwargs["format"].lower()
        else:
            out_format = None

        # check lmax
        if "lmax" in kwargs:
            lmax = kwargs["lmax"]
        else:
            lmax = 10

        # check for intermediate storage directory, temp folder otherwise
        tmp_dir = str(tmp_dir)  # ensure string does not have escape chars
        if not os.path.isdir(tmp_dir):
            tmp_dir = tempfile.gettempdir()

        # Load bulk model and calculate MTZ
        bulk_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(bulk_file):
            bulk_mtz = Converter.import_CLEED(bulk_file, verbose=VERBOSE)
            full_fname = glob(os.path.expanduser(os.path.expandvars(bulk_file)))[0]
            bulk_file = os.path.join(
                tmp_dir, os.path.splitext(os.path.basename(full_fname))[0] + "_bulk.i"
            )
            bulk_mtz.gen_input(filename=bulk_file)
        else:
            bulk_mtz.load_from_file(bulk_file)

        # Load slab model and calculate MTZ
        slab_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(slab_file):
            slab_mtz = Converter.import_CLEED(slab_file)
            full_fname = glob(os.path.expanduser(os.path.expandvars(slab_file)))[0]
            slab_file = os.path.join(
                tmp_dir, os.path.splitext(os.path.basename(full_fname))[0] + "_slab.i"
            )
            slab_mtz.gen_input(filename=slab_file)
        else:
            slab_mtz.load_from_file(slab_file)

        # generate atomic charge densities for each element in bulk model
        if not isinstance(bulk_mtz, model.MTZ_model):
            raise AttributeError("bulk_mtz is not an MTZ_model")

        # get unique elements in bulk and slab
        atomic_dict = {}
        bulk_elements = [atom.element.symbol for atom in bulk_mtz.atoms]
        slab_elements = [atom.element.symbol for atom in slab_mtz.atoms]
        for elem in set(bulk_elements + slab_elements):
            at_file = os.path.join(tmp_dir, "at_%s.i" % elem)
            if not os.path.isfile(at_file):
                print("\nCalculating atomic charge density for %s..." % elem)
                atomic_dict[elem] = atorb.Atorb.calculate_Q_density(
                    element=elem, output_dir=tmp_dir
                )
            else:
                atomic_dict[elem] = at_file

        # prepare at files for appending into atomic file
        bulk_at_files = [
            atomic_dict[atom.element.symbol] for atom in set(bulk_mtz.atoms)
        ]

        # create atomic.i input file from mtz model
        bulk_model_name = os.path.basename(os.path.splitext(bulk_file)[0])
        bulk_atomic_file = bulk_mtz.gen_atomic(
            input_files=bulk_at_files,
            output_file=os.path.join(tmp_dir, bulk_model_name + "_bulk.i"),
        )

        if VERBOSE:
            print("\nModel")
            print("bulk atoms: %s" % [s for s in bulk_mtz.atoms])
            print("slab atoms: %s" % [s for s in slab_mtz.atoms])

        # calculate muffin-tin potential for bulk model
        print("\nCalculating bulk muffin-tin potential...")
        if VERBOSE:
            print("\tcluster file: '%s'" % bulk_file)
            print("\tatomic file: '%s'" % bulk_atomic_file)
            print("\tslab calculation: '%s'" % str(False))
            print(
                "\toutput file: '%s'" % os.path.join(tmp_dir, bulk_model_name + ".bmtz")
            )
            print(
                "\tmufftin file: '%s'"
                % os.path.join(tmp_dir, bulk_model_name + "_mufftin.d")
            )

        bulk_mtz_file = bulk_mtz.calculate_MTZ(
            cluster_file=bulk_file,
            atomic_file=bulk_atomic_file,
            slab=False,
            output_file=os.path.join(tmp_dir, bulk_model_name + ".bmtz"),
            mufftin_file=os.path.join(tmp_dir, bulk_model_name + "_mufftin.d"),
        )
        print("Bulk MTZ = %f" % bulk_mtz.mtz)

        # prepare at files for appending into atomic file
        slab_at_files = [
            atomic_dict[atom.element.symbol] for atom in set(slab_mtz.atoms)
        ]

        # create atomic.i input file from mtz model
        slab_model_name = os.path.basename(os.path.splitext(slab_file)[0])
        slab_atomic_file = slab_mtz.gen_atomic(
            input_files=slab_at_files,
            output_file=os.path.join(tmp_dir, slab_model_name + "_slab.i"),
        )

        # calculate muffin-tin potential for slab model
        mufftin_filepath = os.path.join(tmp_dir, slab_model_name + "_mufftin.d")
        print("\nCalculating slab muff-tin potential...")
        if VERBOSE:
            print("\tcluster file: '%s'" % slab_file)
            print("\tatomic file: '%s'" % slab_atomic_file)
            print("\tslab calculation: %s" % str(True))
            print(
                "\toutput file: '%s'" % os.path.join(tmp_dir, slab_model_name + ".bmtz")
            )
            print("\tmufftin file: '%s'" % os.path.join(tmp_dir, mufftin_filepath))
            print("\tmtz value: %s" % str(bulk_mtz.mtz))

        slab_mtz_file = slab_mtz.calculate_MTZ(
            cluster_file=slab_file,
            output_file=os.path.join(tmp_dir, slab_model_name + ".mtz"),
            atomic_file=slab_atomic_file,
            mufftin_file=mufftin_filepath,
            mtz_string=str(bulk_mtz.mtz),
            slab=True,
        )

        # create raw phase shift files
        print(
            "\nGenerating phase shifts from '%s'..."
            % os.path.basename(mufftin_filepath)
        )
        filepath = os.path.join(tmp_dir, slab_model_name)
        phasout_filepath = filepath + "_phasout.i"
        dataph_filepath = filepath + "_dataph.d"

        phaseshifts = [atom.tag for atom in set(slab_mtz.atoms)]
        phasout_files = [
            os.path.join(tmp_dir, atom.tag + ".ph") for atom in set(slab_mtz.atoms)
        ]
        phsh_files = []

        # assign phase shift specific lmax values
        lmax_dict = {}
        for atom in set(slab_mtz.atoms + bulk_mtz.atoms):
            try:
                lmax_dict[atom.tag] = atom.lmax
            except AttributeError:
                print(atom.tag, "default lmax used:", lmax)
                lmax_dict[atom.tag] = lmax

        try:
            if slab_mtz.nform == 0 or str(slab_mtz.nform).lower().startswith("cav"):
                # calculate phase shifts
                phsh_cav(
                    mufftin_filepath,
                    phasout_filepath,
                    dataph_filepath,
                    filepath + "_zph.o",
                )

                # split phasout
                phasout_files = Conphas.split_phasout(
                    filename=phasout_filepath, output_filenames=phasout_files
                )

                # eliminate pi-jumps
                for i, phaseshift in enumerate(phaseshifts):
                    filename = os.path.splitext(phasout_files[i])[0]
                    if out_format == "curve":
                        filename += ".cur"
                    else:
                        filename += ".phs"
                    phsh_files.append(filename)
                    print(
                        "\nRemoving pi/2 jumps in '%s':\n" % os.path.basename(filename)
                    )
                    phsh = Conphas(
                        input_files=[phasout_files[i]],
                        output_file=filename,
                        formatting=out_format,
                        lmax=lmax_dict[phaseshift],
                    )
                    phsh.calculate()

            if slab_mtz.nform == 1 or str(slab_mtz.nform).lower().startswith("wil"):
                # calculate phase shifts
                phsh_wil(
                    mufftin_filepath,
                    phasout_filepath,
                    dataph_filepath,
                    filepath + "_zph.o",
                )

                # split phasout
                phasout_files = Conphas.split_phasout(
                    filename=phasout_filepath, output_filenames=phasout_files
                )

            if slab_mtz.nform == 2 or str(slab_mtz.nform).lower().startswith("rel"):
                # check energy range
                if "range" in kwargs:
                    try:
                        # get values
                        with open(mufftin_filepath, "r") as f:
                            lines = [line for line in f]

                        ei, de, ef, lsm, vc = [
                            t(s)
                            for t, s in zip(
                                (float, float, float, int, float),
                                lines[1].replace("D", "E").split()[:5],
                            )
                        ]

                        # assign new values
                        (ei, ef, de) = [
                            t(s) for t, s in zip((float, float, float), kwargs["range"])
                        ]

                        # edit energy range
                        lines[1] = str(
                            "%12.4f%12.4f%12.4f    %3i    %12.4f\n"
                            % (ei, de, ef, lsm, vc)
                        ).replace("e", "D")
                        #                         lines[1] = ff.FortranRecordWriter(
                        #                                         '(3D12.4,4X,I3,4X,D12.4)'
                        #                                         ).write([ei, de, ef, lsm, vc]) + '\n'

                        with open(mufftin_filepath, "w") as f:
                            f.write("".join([str(line) for line in lines]))

                    except any as e:
                        sys.stderr.write(
                            "Unable to change phase shift energy "
                            "range - using Barbieri/Van Hove "
                            "default of 20-300eV in 5eV steps\n"
                        )
                        sys.stderr.flush()

                # calculate phase shifts
                # print("Current time " + time.strftime("%X"))
                phsh_rel(
                    mufftin_filepath,
                    phasout_filepath,
                    dataph_filepath,
                    filepath + "_inpdat.txt",
                )
                # print("Current time " + time.strftime("%X"))

                # split phasout
                phasout_files = Conphas.split_phasout(
                    filename=phasout_filepath, output_filenames=phasout_files
                )

                # eliminate pi-jumps
                for i, phaseshift in enumerate(phaseshifts):
                    filename = os.path.splitext(phasout_files[i])[0]
                    if out_format == "curve":
                        filename += ".cur"
                    else:
                        filename += ".phs"
                    phsh_files.append(filename)
                    print(
                        "\nRemoving pi/2 jumps in '%s':\n" % os.path.basename(filename)
                    )
                    phsh = Conphas(
                        input_files=[phasout_files[i]],
                        output_file=filename,
                        formatting=out_format,
                        lmax=lmax_dict[phaseshift],
                    )
                    phsh.calculate()

        except AttributeError:
            raise AttributeError("MTZ_model has no NFORM (0-2) specified!")

        # copy files to storage location
        if "store" in kwargs and out_format != "cleed":
            if kwargs["store"] != ".":
                dst = os.path.abspath(
                    os.path.expanduser(os.path.expandvars(kwargs["store"]))
                )
            else:
                dst = os.path.abspath(".")
            Wrapper._copy_files(phsh_files, dst, verbose=True)

        elif "CLEED_PHASE" in os.environ and out_format == "cleed":
            dst = os.path.abspath(
                os.path.expanduser(os.path.expandvars("$CLEED_PHASE"))
            )
            Wrapper._copy_files(phsh_files, dst, verbose=True)

        else:
            Wrapper._copy_files(phsh_files, os.path.abspath("."), verbose=True)

        return phsh_files

    @staticmethod
    def _copy_files(files, dst, verbose=False):
        """copy list of files into destination directory"""
        # check if using native Windows Python with cygwin
        env = ""
        if platform.system() == "Windows" and dst.startswith("/cygdrive"):
            if os.environ["CLEED_PHASE"] == dst:
                env = "CLEED_PHASE="
            dst = '"%s"' % (
                dst.split("/")[2] + ":" + os.path.sep.join(dst.split("/")[3:])
            )

        # do check and create directory if needed
        if os.path.isfile(dst):
            dst = os.path.dirname(dst)
        if not os.path.exists(dst):
            try:
                os.makedirs(dst)
            except (PermissionError, WindowsError):
                pass

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


class CLIError(Exception):
    """Generic exception to raise and log different fatal errors."""

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def main(argv=None):
    """Command line options."""

    global VERBOSE

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    # display help if no arguments
    if len(argv) == 1:
        argv.append("--help")

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % phaseshifts.__version__
    program_version_message = "%%(prog)s %s" % program_version
    program_shortdesc = __import__("__main__").__doc__.split("\n")[1]
    program_license = """%s

      Created by Liam Deacon using LEEDPACK code (kindly permitted by Micheal Van Hove).
      Copyright 2013-%s Liam Deacon. All rights reserved.

      Licensed under the MIT license (see LICENSE file for details)

      Please send your feedback, including bug notifications
      and fixes, to: %s

    usage:-
    """ % (
        program_shortdesc,
        datetime.date.today().year,
        phaseshifts.__contact__,
    )

    try:
        # Setup argument parser
        parser = ArgumentParser(
            description=program_license, formatter_class=RawDescriptionHelpFormatter
        )
        parser.add_argument(
            "-b",
            "--bulk",
            dest="bulk",
            metavar="<bulk_file>",
            help="path to MTZ bulk or CLEED *.bul input file",
        )
        parser.add_argument(
            "-i",
            "--slab",
            dest="slab",
            metavar="<slab_file>",
            help="path to MTZ slab or CLEED *.inp input file",
            required=True,
        )
        parser.add_argument(
            "-t",
            "--tmpdir",
            dest="tmpdir",
            metavar="<temp_dir>",
            help="temporary directory for intermediate " "file generation",
        )
        parser.add_argument(
            "-l",
            "--lmax",
            dest="lmax",
            metavar="<lmax>",
            default=10,
            type=int,
            help="Maximum angular momentum " "quantum number [default: %(default)s]",
        )
        parser.add_argument(
            "-f",
            "--format",
            dest="format",
            metavar="<format>",
            default="CLEED",
            help="Use specific phase shift format "
            "i.e. 'cleed' or 'curve' "
            "[default: %(default)s]",
        )
        parser.add_argument(
            "-r",
            "--range",
            dest="range",
            nargs="+",
            action=required_length(2, 3),
            type=float,
            metavar="<energy>",
            default=(20.0, 600.0, 5.0),
            help="Energy range in eV with the format: "
            "'<start> <stop> [<step>]'. The <step> "
            "value is optional. Valid for relativistic "
            "calculations only. [default: %(default)s]",
        )
        parser.add_argument(
            "-g",
            "--generate-only",
            dest="generate",
            action="store_true",
            default=False,
            help="Exit after generating "
            "phaseshifts; do not launch subprocess using "
            "PHASESHIFTS_LEED environment variable. "
            "[default: %(default)s]",
        )
        parser.add_argument(
            "-a",
            "--atorbs-only",
            dest="atorbs_only",
            action="store_true",
            default=False,
            help="Only generate atomic orbitals of elements "
            "found in the input files using Eric Shirley's "
            "hartfock routine, then exit. "
            "[default: %(default)s]",
        )
        parser.add_argument(
            "-p",
            "--package",
            dest="package",
            metavar="<package>",
            default="VHT",
            help="Selects package to use for phase shift "
            "calculations. Choices are 'VHT' (van Hove-Tong) "
            "or 'Rundgren' (EEASiSSS). [default: %(default)s]",
        )
        parser.add_argument(
            "-S",
            "--store",
            dest="store",
            metavar="<subdir>",
            default=False,
            help="Keep intermediate files in subdir when done",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            dest="verbose",
            action="count",
            help="Set verbosity level. Note this will also "
            "produce postscript graphs when using the EEASiSSS"
            "backend. [default: %(default)s]",
        )
        parser.add_argument(
            "-V", "--version", action="version", version=program_version_message
        )

        # Process arguments
        args, unknown = parser.parse_known_args()

        verbose = False
        try:
            verbose = args.verbose
            VERBOSE = verbose
        except AttributeError:
            pass

        if verbose > 0 and len(unknown) > 0:
            for arg in unknown:
                sys.stderr.write("phsh - warning: Unknown option '%s'\n" % arg)
            sys.stderr.flush()

        if args.bulk is None:
            args.bulk = str(os.path.splitext(args.slab)[0] + ".bul")

        if args.store is False:
            args.store = "."

        if args.lmax < 1 or args.lmax > 18:
            raise argparse.ArgumentError("lmax is not between 1 and 18")

        if len(args.range) < 3:  # add default step to list
            args.range = list(args.range).append(5)

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as err:
        if DEBUG or TESTRUN:
            raise
        indent = len(program_name) * " "
        sys.stderr.write("{}: '{}'\n".format(program_name, err))
        sys.stderr.write("{} for help use --help".format(indent))
        return 2

    # create phase shifts (warning: black magic within - needs testing)
    if verbose:
        print("Phase shift auto-generation parameters")
        print("\tbulk input file: %s" % args.bulk)
        print("\tslab input file: %s" % args.slab)
        print("\tformat: %s" % args.format)
        print("\tlmax: %s" % args.lmax)
        print("\trange: %s eV" % [s for s in args.range])

    phsh_files = Wrapper.autogen_from_input(
        args.bulk,
        args.slab,
        tmp_dir=args.tmpdir,
        lmax=int(args.lmax),
        format=args.format,
        store=args.store,
        range=args.range,
    )

    # chain loop commands to next program
    if "PHASESHIFTS_LEED" in os.environ and not args.generate:
        # copy files into subdirectory
        csearch = CSearch(os.path.splitext(args.slab)[0])
        last_iteration = csearch.getIteration(-1)
        if last_iteration is not None:
            it = (
                str(last_iteration)
                .split("par:")[0]
                .replace(" ", "")
                .replace("#", "")
                .rjust(3, "0")
            )
            model = os.path.splitext(os.path.basename(args.slab))[0]
            parent = os.path.dirname(args.slab)
            name, ext = os.path.splitext(os.path.basename(args.slab))
            dest = os.path.join(parent, "phsh_" + model, name + "_" + it + ext)
            Wrapper._copy_files(phsh_files, dest, verbose)

        # create subprocess
        leed_cmd = [os.environ["PHASESHIFTS_LEED"]]
        # check if using native Windows Python with cygwin
        if platform.system() == "Windows" and leed_cmd.startwith("/cygdrive"):
            leed_cmd = '"%s"' % (
                leed_cmd.split("/")[2] + ":" + os.path.sep.join(leed_cmd.split("/")[3:])
            )

        for arg in argv:
            leed_cmd.append(arg)

        if verbose:
            print("phsh - starting subprocess: '{}'...".format(" ".join(leed_cmd)))

        # execute subprocess
        subprocess.check_call(leed_cmd)


if __name__ == "__main__":
    if phaseshifts.settings.DEBUG:
        sys.argv.append("-v")

    if phaseshifts.settings.TESTRUN:
        import doctest

        doctest.testmod()

    if phaseshifts.settings.PROFILE:
        import cProfile
        import pstats

        profile_filename = "wrapper_profile.txt"
        cProfile.run("main()", profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats("cumulative")
        stats.print_stats()
        statsfile.close()
        sys.exit(0)

    sys.exit(main())
