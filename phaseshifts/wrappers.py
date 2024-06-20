#!/usr/bin/env python
# encoding: utf-8
from phaseshifts.utils import FileUtils

##############################################################################
# Author: Liam Deacon                                                        #
#                                                                            #
# Contact: liam.deacon@diamond.ac.uk                                         #
#                                                                            #
# Copyright: Copyright (C) 2014-2015 Liam Deacon                             #
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
Provides wrapper classes for phaseshift calculation backends
"""

from __future__ import print_function, unicode_literals
from __future__ import absolute_import, division, with_statement

import sys
import os
import tempfile

from glob import glob
from abc import ABCMeta, abstractmethod, abstractproperty
from getpass import getuser
from time import gmtime, strftime
from re import compile, findall
from copy import deepcopy

from . import model, atorb
from .leed import Converter, CLEED_validator
from .lib.libphsh import phsh_rel, phsh_wil, phsh_cav
from .conphas import Conphas
from .utils import FileUtils


class Wrapper(object):
    """Abstract base wrapper class for generating phase shifts"""

    __metaclass__ = ABCMeta

    @abstractmethod
    def autogen_from_input(self):
        pass

    @abstractmethod
    def autogen_atorbs(self, elements=(), output_dir="."):
        """
        Abstract base method for generating atomic orbital input for
        Eric Shirley's hartfock program.

        Parameters
        ----------
        elements : set of elements.Element objects.
        output_dir : destination directory for generated files.
        """
        pass

    @staticmethod
    def _copy_phsh_files(phsh_files, store=".", out_format=None):
        """
        Description
        -----------
        Copies a list of phase shift files
        """
        out_format = None if out_format is None else out_format.lower()

        # change file extension for CLEED phase shifts
        if out_format == "cleed":
            phsh_files = [os.path.splitext(f)[0] + ".phs" for f in phsh_files]

        if out_format != "cleed":
            dst = FileUtils.expand_filepath(str(store)) if store != "." else "."
            dst = dst if os.path.isdir(dst) else "."
            dst = os.path.abspath(dst)
            FileUtils.copy_files(phsh_files, dst, verbose=True)

        elif "CLEED_PHASE" in os.environ and out_format == "cleed":
            dst = os.path.abspath(FileUtils.expand_filepath("$CLEED_PHASE"))
            FileUtils.copy_files(phsh_files, dst, verbose=True)

        else:
            FileUtils.copy_files(phsh_files, os.path.abspath("."), verbose=True)

    def _add_header(self, phsh_file=None, fmt=None):
        """
        Description
        -----------
        Prepends a header line to the beginning of the phase shift file.
        It assumes that the number of lines in the file are the number of phase
        shifts (excluding any trailing blank lines).

        Parameters
        ----------
        phsh_file : str
            Path to raw phase shift file to prepend header to. If no file is
            given then the function simply returns a string containing the
            header line.
        fmt : str
            Format of the added header line. Currently only 'cleed' and None
            are supported.

        Returns
        -------
        Header line string.

        Raises
        ------
        IOError
            If phsh_file cannot be read or is empty.

        Notes
        -----
        If neng and lmax cannot be found in the object dictionary then a best
        estimate will be made from the input file given.

        """
        if str(self.format).lower() == "cleed":
            # add formatted header for Held CLEED package
            header = ""
            neng = self.neng if "neng" in self.__dict__ else 0
            lmax = self.lmax if "lmax" in self.__dict__ else 0
            if phsh_file is not None:
                lines = []

                # get old content from file
                with open(phsh_file, mode="r", encoding="ascii") as file_ptr:
                    lines = file_ptr.readlines()

                if lines.count("\n") == 0:
                    raise IOError("phase shift file '%s' is empty" % phsh_file)

                # remove trailing lines
                while lines[-1] == "\n":
                    lines.pop(-1)  # removes last element

                # remove comments and empty lines for input
                stripped_lines = [
                    s
                    for s in lines
                    if s != "\n" and s not in compile("[a-zA-Z]").findall(lines)
                ]

                if lmax == 0:
                    # very crude count of l elements in input line
                    # uses multiple input lines to get best estimate
                    n = 10 if len(lines) > 10 else len(lines)
                    for i in range(n):
                        num = self._get_phaseshift_input(stripped_lines[i])
                        num = len(num.split())
                        if num > lmax:
                            lmax = num

                if neng == 0:
                    # guess that neng is half number of non-empty input lines
                    neng = "".join(stripped_lines).count("\n") / 2

                header = "%i %i neng lmax (calculated by %s on %s)\n" % (
                    neng,
                    lmax,
                    getuser(),
                    strftime("%Y-%m-%d at %H:%M:%S", gmtime()),
                )
                lines = lines.insert(0, header)

                # write new contents to file
                with open(phsh_file, mode="w", encoding="ascii") as output_file_ptr:
                    output_file_ptr.writelines(lines)

            return header

        return ""  # zero length string

    def _get_phaseshift_input(self, line):
        """
        Description
        -----------
        Parses a Fortran formatted string containing (multiple) decimal numbers
        and returns a Python-friendly line string to be processed further.
        """
        regex = "[-+0-9]{1,5}\.\d{1,6}"  # basic decimal number
        decimal = compile("({})".format(regex))
        decimal_and_exponent = compile("({}[eED][+-\d]\d{0,6})".format(regex))

        output = decimal.findall(line)
        output = output if output != [] else decimal_and_exponent.findall(line)
        output = [s.replace("D", "E") for s in output]
        return " ".join(output)


class EEASiSSSWrapper(Wrapper):
    """Wrapper class to easily generate phase shifts using EEASiSSS backend"""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def autogen_atorbs(self, elements=[], output_dir="."):
        """
        Generates chgden input files for a set of elements and calculates their
        atomic charge densities using the EEASiSS variant of Eric Shirley's
        hartfock program.

        Parameters
        ----------
        elements : set of elements.Element objects.
        output_dir : destination directory for generated files.

        Returns
        -------
        dictionary of atorb filepaths for all elements
        """
        EEASiSSSWrapper.calculate_Q_density(elements, output_dir=output_dir)
        atomic_dict = {}
        for symbol in [element.symbol for element in elements]:
            atomic_dict[symbol] = os.path.join(output_dir, "chgden%s" % symbol)
        return atomic_dict

    @staticmethod
    def autogen_from_input(
        bulk_file,
        slab_file,
        tmp_dir=None,
        model_name=None,
        lmax=10,
        verbose=False,
        **kwargs
    ):
        """
        Description
        -----------
        Generate phase shifts from a slab/cluster input file.

        Parameters
        ----------
        slab_file : str
            Path to the cluster slab input file.
        bulk_file : str
            Path to the cluster bulk input file.
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
        verbose : bool
            Determines whether to print output to screen.

        Returns
        -------
        output_files : list(str)
            A list of phase shift output filenames
        """
        dummycell = model.Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

        # check formatting of arguments
        out_format = kwargs["format"].lower() if "format" in kwargs else None
        model_name = model_name or "atomic"
        lmax = lmax if isinstance(lmax, int) else 10

        # check for intermediate storage directory, temp folder otherwise
        tmp_dir = tmp_dir or tempfile.gettempdir()

        # Load bulk model and calculate MTZ
        bulk_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(bulk_file):
            bulk_mtz = Converter.import_CLEED(bulk_file, verbose=verbose)
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
            raise AttributeError("bulk_mtz is not an MTZ_model() instance")

        # get unique elements in bulk and slab
        bulk_elements = [atom.element.symbol for atom in bulk_mtz.atoms]
        slab_elements = [atom.element.symbol for atom in slab_mtz.atoms]
        atomic_dict = EEASiSSSWrapper.autogen_atorbs(
            elements=set(bulk_elements + slab_elements), output_dir=tmp_dir
        )

        if verbose:
            print("\nModel")
            print("bulk atoms: %s" % [s for s in bulk_mtz.atoms])
            print("slab atoms: %s" % [s for s in slab_mtz.atoms])

        phsh_files = EEASiSSSWrapper.autogen_from_input(
            inputX, tmp_dir=tmp_dir, **kwargs
        )

        # copy files to storage location
        if "store" in kwargs and out_format != "cleed":
            if kwargs["store"] != ".":
                dst = os.path.abspath(
                    os.path.expanduser(os.path.expandvars(kwargs["store"]))
                )
            else:
                dst = os.path.abspath(".")
            EEASiSSSWrapper._copy_files(phsh_files, dst, verbose=True)

        elif "CLEED_PHASE" in os.environ and out_format == "cleed":
            dst = os.path.abspath(
                os.path.expanduser(os.path.expandvars("$CLEED_PHASE"))
            )
            EEASiSSSWrapper._copy_files(phsh_files, dst, verbose=True)

        else:
            EEASiSSSWrapper._copy_files(phsh_files, os.path.abspath("."), verbose=True)

        return phsh_files


class BVHWrapper(object):
    """
    Wrapper class to easily generate phase shifts
    using the Barbieri - Van Hove backend.
    """

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def autogen_atorbs(self, elements=[], output_dir="."):
        """
        Generates atomic orbital input files for a set of elements and
        calculates their atomic charge densities according to the Barbieri /
        Van Hove variant of Eric Shirley's hartfock program.

        Parameters
        ----------
        elements : set of elements.Element objects.
        output_dir : destination directory for generated files.

        Returns
        -------
        dictionary of atorb filepaths for all elements
        """
        atomic_dict = {}
        for elem in [element.symbol for element in elements]:
            at_file = os.path.join(output_dir, "at_%s.i" % elem)
            if not os.path.isfile(at_file):
                print("\nCalculating atomic charge density for %s..." % elem)
                atomic_dict[elem] = atorb.Atorb.calculate_Q_density(
                    element=elem, output_dir=output_dir
                )
            else:
                atomic_dict[elem] = at_file
        return atomic_dict

    @staticmethod
    def autogen_from_input(
        bulk_file, slab_file, tmp_dir=None, model_name=None, verbose=False, **kwargs
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

        model_name = model_name or "atomic"

        # check formatting
        store = kwargs["store"] if "store" in kwargs else "."
        out_format = kwargs["format"].lower() if "format" in kwargs else None
        lmax = kwargs["lmax"] if "lmax" in kwargs else 10

        # check for intermediate storage directory, temp folder otherwise
        tmp_dir = str(tmp_dir)  # ensure string does not have escape chars
        if os.path.isdir(str(tmp_dir)):
            tmp_dir = tempfile.gettempdir()

        # Load bulk model and calculate MTZ
        bulk_mtz = model.MTZ_model(dummycell, atoms=[])
        if CLEED_validator.is_CLEED_file(bulk_file):
            bulk_mtz = Converter.import_CLEED(bulk_file, verbose=verbose)
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
        atomic_dict = BVHWrapper.autogen_atorbs(
            elements=set(bulk_elements + slab_elements), output_dir=tmp_dir
        )

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

        if verbose:
            print("\nModel")
            print("bulk atoms: %s" % [s for s in bulk_mtz.atoms])
            print("slab atoms: %s" % [s for s in slab_mtz.atoms])

        # calculate muffin-tin potential for bulk model
        print("\nCalculating bulk muffin-tin potential...")
        if verbose:
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
        if verbose:
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

        phaseshifts = [atom.tag for atom in set(slab_mtz.atoms + bulk_mtz.atoms)]
        phasout_files = [
            os.path.join(tmp_dir, atom.tag + ".ph")
            for atom in set(slab_mtz.atoms + bulk_mtz.atoms)
        ]
        phsh_files = []

        # assign phase shift specific lmax values
        lmax_dict = {}  # FIXME: artifact from cherry-pick?
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

                        with open(mufftin_filepath, "w") as f:
                            f.write("".join([str(line) for line in lines]))

                    except Exception as err:
                        sys.stderr.write(
                            "Unable to change phase shift energy "
                            "range (due to '{}') - using Barbieri/Van Hove "
                            "default of 20-300eV in 5eV steps\n".format(err)
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
        Wrapper._copy_phsh_files(phsh_files, store, out_format)

        return phsh_files
