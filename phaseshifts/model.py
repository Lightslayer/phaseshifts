#!/usr/bin/env python
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
**model.py**

Provides convenience functions for generating input and calculating
atomic charge densities for use with the Barbieri/Van Hove phase
shift calculation package.

"""

from __future__ import print_function
from __future__ import division

import os
from copy import deepcopy
from shutil import move
from glob import glob

import phaseshifts.atorb as atorb
import phaseshifts.elements as elements
from phaseshifts.lib import libphsh


class Atom(object):
    """
    Atom class for input into cluster model for muffin-tin potential
    calculations.
    """

    def __init__(self, element, coordinates=[0.0, 0.0, 0.0], **kwargs):
        """
        Constructor for Atom class.

        Parameters
        ----------

        element : str or int
            This is either the elemental symbol, name or atomic number.
        coordinates : list[float, float, float] or ndarray
            The fractional coordinates within the unitcell in terms of the
            basis vector a.
        tag : str, optional
            Add a name tag to this element (useful if dealing with multiple
            atoms of the same element in a given model). Default is the
            symbol for that element - numeric ids may be appended in the model
            class.
        radius : float, optional
            The muffin-tin radius of the atom in Angstroms (default is to
            lookup 'atmrad' in the element dictionary).
        valence : int, optional
            The valency of the atom (default is to assume neutral atom).
        occupancy : float, optional
            The fractional occupancy of the atom in the given site.

        """
        self.element = elements.ELEMENTS[element]
        self.coordinates = None  # dummy
        self.set_coordinates(coordinates)
        self.name = self.element.name.title()
        self.tag = self.element.symbol.title()
        self.Z = self.element.protons
        self.radius = self.element.atmrad
        self.valence = 0
        if "valence" in kwargs:
            if "radius" not in kwargs:
                # assume covrad for non-zero valency
                self.radius = self.element.covrad
        self.__dict__.update(kwargs)

    # checks whether two atoms are equal w.r.t. name, radius and valence
    def __eq__(self, other):
        if isinstance(other, Atom):
            return (
                self.name == other.name
                and self.radius == other.radius
                and self.valence == other.valence
            )
        else:
            return False

    # checks whether two atoms are not equal w.r.t. name, radius and valence
    def __neq__(self, other):
        return not self.__eq__(other)

    # reprinting of Atom object
    def __repr__(self):
        return str("Atom(%s, tag='%s', radius=%s, " "valence=%s)") % (
            self.name,
            self.tag,
            self.radius,
            self.valence,
        )

    # redefine hash method for checking uniqueness of class instance
    def __hash__(self):
        return hash(self.__repr__())

    # set coordinates of atom within unitcell in terms of a
    def set_coordinates(self, coordinates):
        try:
            self.coordinates = coordinates
            self._coordinates = [r / 0.529 for r in coordinates]
        except any as e:
            raise e

    # set valence of atom
    def set_valence(self, valency):
        """Sets the valency of the atom"""
        self.valence = float(valency)

    # set muffin-tin radius of atom
    def set_mufftin_radius(self, radius):
        """
        Sets the muffin-tin radius of the atom in Angstroms.
        """
        try:
            self.radius = float(radius)
            self._radius = self.radius / 0.529  # in Bohr radii
        except:
            pass


class Unitcell(object):
    """
    Unitcell class
    """

    def __init__(self, a, c, matrix_3x3, **kwargs):
        """
        Constructor for the Unitcell class

        Parameters
        ----------
        a : float
            The in-plane lattice vector in Angstroms
        c : float
            The out-of-plane lattice vector in Angstroms. For cubic systems this
            will be equal to a.
        matrix_3x3: ndarray
            A 3x3 matrix describing the x,y,z construction of a,b,c basis vectors
            of the unitcell. Units for x, y & z should be in terms of fractional
            coordinates.
        alpha : float, optional
            Angle alpha in degrees.
        beta : float, optional
            Angle beta in degrees.
        gamma : float, optional
            Angle gamma in degrees.

        """
        # Convert Angstrom input to Bohr radii
        self.set_a(a)
        self.set_c(c)

        # Set basis vectors
        self.set_vectors(matrix_3x3)

        # Set xstal system
        self.alpha = 90.0
        self.beta = 90.0
        self.gamma = 90.0

        # Update additional information
        self.__dict__.update(kwargs)

    # checks if two class instances are equal w.r.t. name, radius & valence
    def __eq__(self, other):
        if isinstance(other, Atom):
            return (
                self.a == other.a
                and self.c == other.c
                and self.alpha == other.alpha
                and self.beta == other.beta
                and self.gamma == other.gamma
            )
        else:
            return False

    # checks if two class instances are not equal w.r.t. name, radius & valence
    def __neq__(self, other):
        return not self.__eq__(other)

    # reprinting of class object
    def __repr__(self):
        return str("Unitcell(a=%s, c=%s, alpha=%s, beta=%s, gamma=%s, basis=%s)") % (
            self.a,
            self.c,
            self.alpha,
            self.beta,
            self.gamma,
            self.basis,
        )

    # redefine hash method for checking uniqueness of class instance
    def __hash__(self):
        return hash(self.__repr__())

    # set basis vectors from (3x3) matrix in fractional coordinates
    def set_vectors(self, m3x3):
        self.basis = m3x3

    # set a lattice parameter
    def set_a(self, a):
        """
        Description
        -----------
        Set the magnitude of the in-plane lattice vector a in Angstroms.

        Parameters
        ----------
        a: float
            The magnitude of the in-plane lattice vector in Angstroms

        Notes
        -----
        To retrieve a in terms of Angstroms use 'unitcell.a', whereas the
        internal parameter 'unitcell._a' converts a into Bohr radii
        (1 Bohr = 0.529Å), which is used for the muffin-tin potential
        calculations in libphsh (CAVPOT subroutine).

        """
        self.a = float(a)
        self._a = self.a / 0.529  # (1 Bohr = 0.529Å)

    # set c lattice parameter
    def set_c(self, c):
        """
        Description
        -----------
        Set the magnitude of the out-of-plane lattice vector a.

        Parameters
        ----------
        c : float
            The magnitude of the in-plane lattice vector in Angstroms

        Notes
        -----
        To retrieve c in terms of Angstroms use 'unitcell.c', whereas the
        internal parameter 'unitcell._c' converts c into Bohr radii
        (1 Bohr = 0.529Å), which is used for the muffin-tin potential
        calculations in libphsh (CAVPOT subroutine).

        """
        self.c = float(c)
        self._c = self.c / 0.529  # (1 Bohr = 0.529Å)

    # set angle alpha in degrees
    def set_alpha(self, alpha):
        self.alpha = float(alpha) % 360.0

    # set angle beta in degrees
    def set_beta(self, beta):
        self.beta = float(beta) % 360.0

    # set angle gamma in degrees
    def set_gamma(self, gamma):
        self.gamma = float(gamma) % 360.0


class CoordinatesError(Exception):
    """Coordinate exception to raise and log duplicate coordinates."""

    def __init__(self, msg):
        super(CoordinatesError).__init__(type(self))
        self.msg = "CoordinatesError: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


class Model(object):
    """
    Generic model class.
    """

    def __init__(self, unitcell, atoms, **kwargs):
        """
        Constructor for Model class.

        Parameters
        ----------
        unitcell : Unitcell
            An instance of the Unitcell class.
        atoms : list
            Array of Atom class instances which constitute the model.

        """
        self.atoms = []
        self.set_atoms(atoms)
        self.set_unitcell(unitcell)
        self.__dict__.update(kwargs)

    # checks if two models are equal
    def __eq__(self, other):
        if isinstance(other, Atom):
            return self.atoms == other.atoms and self.unitcell == other.unitcell
        else:
            return False

    # checks if two models are not equal
    def __neq__(self, other):
        return not self.__eq__(other)

    # reprinting of Atom object
    def __repr__(self):
        return str("Model(atoms=%s, unitcell=%s)") % (self.atoms, self.unitcell)

    # redefine hash method for checking uniqueness of class instance
    def __hash__(self):
        return hash(self.__repr__())

    # estimate number of inequivalent atoms
    def _nineq_atoms(self):
        """
        Description
        -----------
        Internal method for estimating the number of inequivalent atoms

        Returns
        -------
        nineq_atoms, element_dict : tuple
            nineq_atoms : The estimated number of inequivalent atoms based on
                the valence and radius of each atom.
            element_dict : a dictionary of each element in the atom list where
                each element contains an atom dictionary of 'nineq_atoms',
                'n_atoms' and a complete 'atom_list'

        Example
        -------
        >>> C1 = Atom('C', [0, 0, 0])
        >>> Re1 = Atom('Re', [0, 0, 0], valence=2.0)
        >>> Re2 = Atom('Re', [0, 0, 0], radius=1)
        >>> uc = Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> mtz = MTZ_model(uc, atoms=[C1, Re1, Re2])
        >>> print(mtz._nineq_atoms())
         (3, {'Carbon': {'n_atoms': 1, 'atom_list': [Atom(Carbon, tag='C',
         coordinates=[0, 0, 0], Z=6, radius=0.91, valence=0)], 'nineq_atoms':
         1,'nineq_atoms_list': set([Atom(Carbon, tag='C',
         coordinates=[0, 0, 0], Z=6, radius=0.91, valence=0)])}, 'Rhenium': {
         'n_atoms': 2, 'atom_list': [Atom(Rhenium, tag='Re',
         coordinates=[0, 0, 0], Z=75, radius=1.97, valence=2.0), Atom(Rhenium,
         tag='Re', coordinates=[0, 0, 0], Z=75, radius=1, valence=0)],
         'nineq_atoms': 2, 'nineq_atoms_list': set([Atom(Rhenium, tag='Re',
         coordinates=[0, 0, 0], Z=75, radius=1.97, valence=2.0),
         Atom(Rhenium, tag='Re', coordinates=[0, 0, 0], Z=75, radius=1,
         valence=0)])}})

        """
        nineq_atoms = 0
        element_dict = {}
        atom_dict = {}
        # loop through atom list, testing each element for duplicates
        # get list of elements
        elements = set([atom.name for atom in self.atoms])
        for element in elements:
            atoms = [atom for atom in self.atoms if atom.name == element]
            n_atoms = len(set(atoms))
            nineq_atoms += n_atoms
            atom_dict = {
                "nineq_atoms": n_atoms,
                "n_atoms": len(atoms),
                "atom_list": atoms,
            }
            element_dict[element] = atom_dict
        return nineq_atoms, element_dict

    def add_atom(self, element, position, **kwargs):
        """
        Append an Atom instance to the model

        Parameters
        ----------
        element : str or int
            Either an element name, symbol or atomic number.
        position : list(float) or ndarray
            (1x3) array of the fractional coordinates of the atom
            within the unit cell in terms of the lattice vector a.

        """
        self.atoms.append(Atom(element, position, kwargs))

    def check_coordinates(self):
        """
        Check for duplicate coordinates of different atoms in model.

        Raises
        ------
        CoordinateError : exception
          If duplicate positions found.

        """
        positions = [str(atom.coordinates) for atom in self.atoms]
        info = ""
        for position in set(
            [position for position in positions if positions.count(position) > 1]
        ):
            for i, atom in enumerate(
                [atom for atom in self.atoms if str(atom.coordinates) == position]
            ):
                info += "%s, coordinates=%s, index=%i\n" % (
                    str(atom),
                    atom.coordinates,
                    i,
                )
        if len(set(positions)) < len(self.atoms):
            raise CoordinatesError(
                "Not every atom position in model is unique!\n%s\n" % info
            )

    @property
    def name(self):
        """Returns the model name"""
        return self.name or ""

    @name.setter
    def name(self, name):
        """Sets the name of the model"""
        self.name = str(name) if not isinstance(name, str) else name

    @property
    def elements(self):
        """Returns a list of unique elements within the model"""
        return {atom.element for atom in self.atoms}

    @property
    def atoms(self):
        """Returns a list of atoms within this model"""
        return self.atoms or []

    @atoms.setter
    def atoms(self, atoms):
        """
        Set the atoms for the model.

        Parameters
        ----------
        atoms : list(Atoms)
            Array of Atom instances. Entries in the list which are
            not of type Atom will be ignored.

        Raises
        ------
        TypeError : exception
          If atoms parameter is not a list.

        """
        if isinstance(atoms, list):
            self.atoms = [atom for atom in atoms if isinstance(atom, Atom)]
        else:
            raise TypeError

    def set_unitcell(self, unitcell):
        """
        Set the unitcell for the model

        Parameters
        ----------
        unitcell : Unitcell
            Instance of Unitcell class to set to model.

        Raises
        ------
        TypeError : exception
          If unitcell parameter is not a Unitcell.

        """
        if isinstance(unitcell, Unitcell):
            self.unitcell = unitcell
        else:
            raise TypeError


class MTZ_model(Model):
    """
    Muffin-tin potential Model subclass for producing input file
    for muffin-tin calculations in the Barbieri/Van Hove phase
    shift calculation package.
    """

    def __init__(self, unitcell, atoms, **kwargs):
        """
        Constructor for Model class.

        Parameters
        ----------
        unitcell : Unitcell
            An instance of the Unitcell class.
        atoms : list
            Array of Atom class instances which constitute the model.
        nh : int, optional
            Parameter for estimating the muffin-tin zero (default: 10).
        exchange : float, optional
            Hartree type exchange term alpha (default: 0.7200).
        c : float, optional
            The height of the slab in Angstroms - if this is much larger
            than the bulk c distance then there will be a large vacuum
            and therefore should be used when calculating a thin slab
            rather than a bulk muffin-tin potential. Default is to lookup
            the out-of-plane basis vector bulk value.
        nform : int or str, optional
            The phase shift calculation type, which can be 0 or 'cav' for
            using the cavpot subroutine, 1 or 'wil' for using the William's
            method, and 2 or 'rel' for using relativistic calculations suitable
            for the CLEED package.

        """
        self.atoms = []
        self.set_atoms(atoms)
        self.set_unitcell(unitcell)
        self.set_exchange(0.72)
        self.set_nh(10)
        self.mtz = None
        self.__dict__.update(kwargs)

    def set_nh(self, nh):
        """Sets the nh muffin-tin zero estimation parameter"""
        self.nh = int(nh)  # check this is not float

    def set_exchange(self, alpha):
        """Sets the alpha exchange term for muffin-tin calculation"""
        try:
            self.exchange = float(alpha)
        except:
            pass

    # set form of muffin-tin calculation: 0=cav, 1=wil, 2=rel
    def set_nform(self, nform):
        """
        Sets form of muffin-tin calculation

        Parameters
        ----------
        nform : int or str
          This governs the type of calculation, where nform can be:

          **1.** "cav" or 0 - use Cavendish method

          **2.** "wil" or 1 - use William's method

          **3.** "rel" or 2 - Relativistic calculations

        """
        try:
            if isinstance(nform, int):
                if nform >= 0 and nform <= 2:
                    self.nform = nform
            elif isinstance(nform, str):
                if nform == "cav":
                    self.nform = 0
                elif nform == "wil":
                    self.nform = 1
                elif nform == "rel":
                    self.nform = 2
                elif nform in ["0", "1", "2"]:
                    self.nform = int(nform)
        except any as e:
            raise TypeError(e.msg)

    def set_slab_c(self, c):
        """
        Description
        -----------
        Set the maximum height of the slab in Angstroms - if this is
        much larger than the bulk c distance then there will be a large
        vacuum and therefore should be used when calculating a thin slab
        rather than a bulk muffin-tin potential.

        Examples
        --------
        For Re the bulk c distance is 2.76Å, whereas a possible slab c distance
        could be ~10Å.

        """
        try:
            self.c = float(c)
            self._c = self.c / 0.529
        except:
            pass

    def _load_input_file(self, filename):
        """
        Description
        -----------
        Load an cluster/slab input file and update the class instance variables

        Parameters
        ----------
        filename : str
            The path of the input file (e.g. cluster*.i or *slab*.i)

        Raises
        ------
        IOError : exception
          If the file cannot be read.
        TypeError : exception
          If a input line cannot be parsed correctly.

        """

        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]

        try:
            with open(filename, "r") as f:
                self.header = f.readline()
                a = float(f.readline().split("#")[0].split()[0]) * 0.529
                a1 = [
                    t(s)
                    for (t, s) in zip(
                        (float, float, float), f.readline().split("#")[0].split()[:3]
                    )
                ]
                a2 = [
                    t(s)
                    for (t, s) in zip(
                        (float, float, float), f.readline().split("#")[0].split()[:3]
                    )
                ]
                a3 = [
                    t(s)
                    for (t, s) in zip(
                        (float, float, float), f.readline().split("#")[0].split()[:3]
                    )
                ]
                basis = [a1, a2, a3]
                c = float(a3[-1]) * 0.529  # change to Angstroms from Bohr
                self.set_unitcell(Unitcell(a, c, basis))
                self.nineq_atoms = int(f.readline().split("#")[0].split()[0])

                self.atoms = []
                line = f.readline().split("#")[0]
                while line.split()[0].isalpha():
                    n, Z, val, rad = [
                        t(s)
                        for (t, s) in zip(
                            (int, float, float, float),
                            f.readline().split("#")[0].split()[:4],
                        )
                    ]

                    for i in range(0, n):
                        position = [
                            t(s)
                            for (t, s) in zip(
                                (float, float, float),
                                f.readline().split("#")[0].split()[:3],
                            )
                        ]

                        atom = Atom(
                            line.split()[0].capitalize(),
                            coordinates=position,
                            tag=line.split()[1],
                            valence=val,
                            radius=rad,
                            Z=int(Z),
                        )

                        self.atoms.append(atom)
                    line = f.readline().split("#")[0]

                # print(line)
                self.set_nform(line.split()[0])
                self.set_exchange(f.readline().split("#")[0].split()[0])
                self.set_nh(f.readline().split("#")[0].split()[0])

        except IOError:
            raise IOError("cannot read '%s'" % filename)

        except ValueError:
            raise ValueError("malformatted input in '%s'" % filename)

    def load_from_file(self, filename):
        """
        Description
        -----------
        Load an input file and update the class instance variables

        Parameters
        ----------
        filename : str
            The path of the input file (e.g. cluster*.i or *slab*.i)

        Raises
        ------
        IOError : exception
            If the file cannot be read.
        TypeError : exception
            If a input line cannot be parsed correctly.

        """

        filename = glob(os.path.expanduser(os.path.expandvars(filename)))[0]
        if CLEED_validator.is_CLEED_file(filename):
            self = Converter.import_CLEED(filename)
        else:
            self._load_input_file(filename)

    def create_atorbs(self, **kwargs):
        """
        Description
        -----------
        Create Atorb input files for each element in model.

        Arguments
        ---------
        output_dir : str
            Path to output directory for files
        library_dir : str
            Path to library of input files. Use this if you've previously
            created input files for a given element as it doesn't need to
            recalculate the radial charge density every time - note this
            is also a workaround to allow precalculated files if your machine's
            output is drastically different from what is expected.
        config : dict
            See help(atorb.gen_input) for list of keywords and values.

        Returns
        -------
        output_files : dict
            Dictionary list of atorb*.i input files for the Atorb class to
            calculate the charge density from.

        """

        # check output path
        if "output_dir" in kwargs:
            output_dir = kwargs["output_dir"]
        else:
            if "library_dir" in kwargs:
                output_dir = kwargs["library_dir"]
            else:
                output_dir = os.path.abspath(".")

        atorb_files = []
        at_files = []

        # generate input files for atomic orbitals and charge densities
        for element in set([str(atom.element.symbol) for atom in self.atoms]):
            if not os.path.isdir(os.path.join(output_dir, "Atorb")):
                try:
                    os.makedirs(os.path.join(output_dir, "Atorb"))
                except WindowsError:
                    pass  # directory already exists on Windows

            atorb_file = os.path.join(output_dir, "Atorb", "atorb_%s.i" % element)
            if not os.path.isfile(atorb_file):  # create new atorb input file
                atorb_file = atorb.Atorb.gen_input(element, filename=atorb_file)

            at_file = os.path.join(output_dir, "Atorb", "at_" + element + ".i")
            if not os.path.isfile(at_file):
                at_file = atorb.Atorb.calculate_Q_density(
                    input=atorb_file, output_dir=os.path.join(output_dir, "Atorb")
                )

            # update lists
            atorb_files.append(atorb_file)
            at_files.append(at_file)

        return {"atorb_files": atorb_files, "at_files": at_files}

    def gen_atomic(self, **kwargs):
        """
        Description
        -----------

        Create atomic*.i input file for MTZ input based on model or
        list of files.

        Parameters
        ----------
        input_dir : str
            Input directory where at*.i files are kept.
        input_files : tuple
            List of input files to generate atomic input file from.
        output_file : str
            The filename of the resulting atomic*.i output file, which is
            simply a superimposed set of the radial charge densities from
            the individual input files.

        Returns
        -------
        output_file : str
            Returns the output file path string.

        Raises
        ------
        IOError : exception
            If either input directory or files do not exist.

        Notes
        -----
        If 'input_files' is not given then the default list of input files are
        inferred from the list of atoms in the model.

        """

        # input
        if "input_dir" in kwargs:
            input_dir = os.path.abspath(
                glob(os.path.expanduser(os.path.expandvars(kwargs["input_dir"])))[0]
            )
        else:
            input_dir = os.path.abspath(".")

        if os.path.isfile(input_dir):
            input_dir = os.path.dirname(input_dir)

        if not os.path.isdir(input_dir):
            raise IOError(
                "'%s' is not a valid input directory - " "does not exist!" % input_dir
            )

        # output filename
        if "output_file" in kwargs:
            output_file = os.path.abspath(kwargs["output_file"])
        else:
            output_file = os.path.abspath("atomic.i")

        # get list of input
        if "input_files" in kwargs:
            files = kwargs["input_files"]
        else:  # assume using atoms from model
            files = [
                os.path.join(input_dir, "at_" + atom.element.symbol + ".i")
                for atom in self.atoms
            ]

        # generate atomic.i input file by appending multiple at.i files
        with open(output_file, "w") as f:
            # loop through each atomic charge density file in list
            for input_file in files:
                if not os.path.isfile(str(input_file)) or input_file == None:
                    raise IOError(
                        "Radial charge density file "
                        "'%s' does not exist!" % input_file
                    )

                # append next input file to output
                with open(input_file) as infile:
                    f.write(infile.read())

        return output_file

    def get_MTZ(self, filename):
        """Retrieves muffin-tin potential from file"""
        try:
            with open(filename, "r") as f:
                self.mtz = float([line for line in f][0])  # read first line
        except IOError:
            raise IOError

    def calculate_MTZ(self, mtz_string="", **kwargs):
        """
        Description
        -----------
        Calculate the muffin-tin potential (MTZ) for a given cluster file

        Parameters
        ----------
        atomic_file : str
            The path to the atomic input file. If this is omitted the default
            is generate one using the MTZ_model.gen_atomic() method.
        cluster_file : str
            The path to the cluster input file which may be a bulk or slab
            model.
        slab : int or bool
            Determines whether the MTZ calculation is for a slab model (True).
            The default is a bulk calculation.
        output : dict
            Dictionary output of 'mtz' - muffin-tin potential & 'output_file'
            - the path to the MTZ output file.


        Returns
        -------
        output_files : list(str)
            Paths to the MTZ output file.

        """

        # check to see if cluster input exists
        if "cluster_file" in kwargs:
            cluster_file = os.path.abspath(
                glob(os.path.expanduser(os.path.expandvars(kwargs["cluster_file"])))[0]
            )
        else:
            cluster_file = os.path.abspath("cluster.i")

        if not os.path.isfile(cluster_file):
            raise IOError("MTZ cluster file '%s' does not exist!" % cluster_file)

        if (
            not os.access(os.path.dirname(cluster_file), os.W_OK)
            and "atomic_file" not in kwargs
        ):
            raise IOError(
                "Do not have write access to '%s'" % os.path.dirname(cluster_file)
            )

        # determine type of calculation - bulk or slab
        if "slab" in kwargs:
            slab = bool(kwargs["slab"])
            fid = "mtz"
        else:
            slab = False
            fid = "bmtz"
            mtz_string = ""

        if "atomic_file" in kwargs:
            atomic_file = os.path.abspath(kwargs["atomic_file"])
            if not os.path.isfile(atomic_file):
                raise IOError(
                    "Appended radial charge densities file "
                    "'%s' does not exist!" % atomic_file
                )
        else:  # generate on the fly
            input_dir = os.path.abspath(os.path.dirname(cluster_file))
            self.create_atorbs(output_dir=input_dir)
            self.gen_atomic(input_dir=input_dir)

        if "output_file" in kwargs:
            output_file = os.path.abspath(kwargs["output_file"])
        else:
            output_file = "%s.i" % fid

        try:
            os.makedirs(os.path.dirname(output_file))
        except WindowsError:
            pass

        if os.path.isfile(output_file):
            move(output_file, output_file + ".bak")

        # create mufftin debug file
        mufftin_file = os.path.splitext(output_file)[0] + "_mufftin.d"
        info_file = os.path.splitext(output_file)[0] + "_info.txt"

        # call cavpot routine
        self.mtz = libphsh.cavpot(
            mtz_string,
            int(slab),
            atomic_file,
            cluster_file,
            mufftin_file,
            output_file,
            info_file,
        )

        # check to see if new file has been written
        if not os.path.isfile(output_file):
            raise IOError(
                "Failed to write muffin-tin potential file '%s'" % output_file
            )

        return output_file

    def gen_input(self, **kwargs):
        """
        Description
        -----------
        Generate input file for use in the Barbieri/Van Hove phase shift
        calculation package (phsh1 or libphsh)

        Arguments
        ---------
        bulk : bool
            If True the c value is set to the bulk value...
        filename : str
            The path of the generated file (e.g. 'cluster.i')
        header : str
            A file header string.
        bulk : bool
            Specify whether generated file is a bulk or slab model
        pos_check : bool
            If true, the positions of each atom will be checked and
            any duplicates will cause a CoordinateError to be raised.
        nform : int or str
            Choose calculation type from either:
              - Cavendish (0)
              - Williams (1)
              - Relativistic (2)
        exchange : float
            Hartree-Fock exchange term alpha.
        nh : int
            Term used to estimate the muffin-tin zero.

        Returns
        -------
        filename on success

        Raises
        ------
        CoordinatesError : exception
          if the model atoms have duplicate coordinates and the 'pos_check'
          kwarg is set to True.

        """
        # create backup of atoms
        mtz_atoms = deepcopy(self.atoms)

        # get file header - should be single line (force if not)
        if "header" in kwargs:
            header = kwargs["header"].replace("\n", " ").replace("\r", "")
        else:
            header = "MTZ cluster input file"

        if "bulk" in kwargs:
            if kwargs["bulk"]:
                fid = "_bulk"
            else:
                fid = "_slab"
        else:
            fid = ""

        if "nform" in kwargs:
            self.set_nform(kwargs["nform"])
        else:
            self.nform = 2  # rel

        if "exchange" in kwargs:
            self.set_exchange(kwargs["exchange"])
        elif not isinstance(self.exchange, float):
            self.set_exchange(0.72)

        if "nh" in kwargs:
            self.set_nh(kwargs["nh"])
        elif not isinstance(self.nh, int):
            self.set_nh(10)

        # get filename or make educated guess
        if "filename" in kwargs:
            filename = kwargs["filename"]
        else:
            filename = "cluster%s.i" % fid

        # determine if coordinate duplicates
        if "pos_check" in kwargs:
            if bool(kwargs["pos_check"]):
                self.check_coordinates()

        # auto-generate file
        a = float(self.unitcell._a) * 0.529
        with open(filename, "w") as f:
            f.write(header + "\n")
            f.write(
                str(" %7.4f" % self.unitcell._a).ljust(33)
                + "# a lattice parameter distance in Bohr radii\n"
            )
            f.write(
                str(
                    " %7.4f %7.4f %7.4f"
                    % (
                        self.unitcell.basis[0][0] / a,
                        self.unitcell.basis[0][1] / a,
                        self.unitcell.basis[0][2] / a,
                    )
                ).ljust(33)
                + "# SPA coordinates of unit cell\n"
            )
            f.write(
                str(
                    " %7.4f %7.4f %7.4f"
                    % (
                        self.unitcell.basis[1][0] / a,
                        self.unitcell.basis[1][1] / a,
                        self.unitcell.basis[1][2] / a,
                    )
                ).ljust(33)
                + "\n"
            )
            f.write(
                str(
                    " %7.4f %7.4f %7.4f"
                    % (
                        self.unitcell.basis[2][0] / a,
                        self.unitcell.basis[2][1] / a,
                        self.unitcell.basis[2][2] / a,
                    )
                ).ljust(33)
                + "# Notice the value %.2f (%s calculation)\n"
                % (self.unitcell.basis[2][2] / a, fid.replace("_", ""))
            )

            # TODO: better nineq_atoms prediction
            try:
                nineq_atoms = int(kwargs["nineq_atoms"])
            except TypeError:
                nineq_atoms = self.nineq_atoms()[0]
            except KeyError:
                nineq_atoms = len(set(self.atoms))

            tags = []

            # check to see if nineq_atoms is estimated in code
            if isinstance(nineq_atoms, tuple):
                f.write(
                    str("%4i" % nineq_atoms[0]).ljust(33)
                    + "# number of ineq. atoms in this file (NINEQ)\n"
                )

                # now loop through each inequivalent atom and add to file
                elements_dict = nineq_atoms[1]
                for element in elements_dict:
                    atoms = elements_dict.get(element)["atom_list"]

                    # Loop through each inequivalent atom type
                    for ineq_atom in set(atoms):
                        # get list of atoms of this type
                        # i.e. same element, radius & valence
                        ineq_atoms = [atom for atom in atoms if atom == ineq_atom]
                        ineq_tags = set(
                            [atom.tag for atom in ineq_atoms if atom.tag not in tags]
                        )

                        # select first unused tag from list
                        for tag in ineq_tags:
                            if tag not in tags:
                                ineq_atom.tag = tag
                                break

                        # avoid duplicate tags for different atoms
                        while ineq_atom.tag in tags:
                            number = "".join([ch for ch in atom.tag if ch.isdigit()])
                            try:
                                number = int(number)
                                number += 1
                            except ValueError:
                                number = 1
                            ineq_atom.tag = (
                                "".join(
                                    [
                                        ch
                                        for ch in ineq_atom.tag
                                        if ch.isalpha() or ch in ["_", "-", "+"]
                                    ]
                                )
                                + "_"
                                + str(number)
                            )

                        f.write(
                            "{0} {1}".format(
                                ineq_atom.element.name.capitalize(), ineq_atom.tag
                            ).ljust(33)
                            + "# element, name tag\n"
                        )
                        f.write(
                            str(
                                "%4i %7.4f %7.4f %7.4f"
                                % (
                                    len(ineq_atoms),
                                    ineq_atom.Z,
                                    ineq_atom.valence,
                                    ineq_atom.radius / 0.529,
                                )
                            ).ljust(33)
                            + "# atoms in unit cell, "
                            "Z, valence, Muffin-tin radius (Bohr radii)\n"
                        )

                        # write each inequivalent atom list to file
                        for atom in ineq_atoms:
                            f.write(
                                str(
                                    " %7.4f %7.4f %7.4f"
                                    % (
                                        atom.coordinates[0] / a,
                                        atom.coordinates[1] / a,
                                        atom.coordinates[2] / a,
                                    )
                                ).ljust(33)
                                + "# coordinates in SPA units\n"
                            )

            else:  # assume each element is an inequivalent atom
                f.write(
                    str("%4i" % nineq_atoms).ljust(33)
                    + "# number of ineq. atoms in this file (NINEQ)\n"
                )

                for atom in set(self.atoms):
                    while atom.tag in tags:
                        number = "".join([ch for ch in atom.tag if ch.isdigit()])
                        try:
                            number = int(number)
                            number += 1
                        except ValueError:
                            number = 1
                        atom.tag = (
                            "".join(
                                [
                                    ch
                                    for ch in atom.tag
                                    if ch.isalpha() or ch in ["_", "-", "+"]
                                ]
                            )
                            + "_"
                            + str(number)
                        )

                    tags.append(atom.tag)
                    f.write(
                        "{0} {1}".format(
                            atom.element.name.capitalize(), atom.tag
                        ).ljust(33)
                        + "# element, name tag\n"
                    )
                    f.write(
                        str(
                            "%4i %7.4f %7.4f %7.4f"
                            % (
                                tags.count(atom.tag),
                                atom.Z,
                                atom.valence,
                                atom.radius / 0.529,
                            )
                        ).ljust(33)
                        + "# atoms in unit cell, Z, valence, "
                        "Muffin-tin radius (Bohr radii)\n"
                    )
                    f.write(
                        str(
                            " %7.4f %7.4f %7.4f"
                            % (
                                atom.coordinates[0] / a,
                                atom.coordinates[1] / a,
                                atom.coordinates[2] / a,
                            )
                        ).ljust(33)
                        + "# coordinates in SPA units\n"
                    )

            f.write(
                str("%4i" % self.nform).ljust(33)
                + "# nform=2|1|0 (for rel, will or cav)\n"
            )
            f.write(
                str(" %7.4f" % self.exchange).ljust(33)
                + "# alpha (for Hartree type exchange term)\n"
            )
            f.write(
                str("%4i" % self.nh).ljust(33)
                + "# nh (for estimating muffin-tin zero)\n"
            )
        # restore backup
        self.atoms = mtz_atoms

    def get_elements(self):
        """Return the unique elements in model"""
        return set([atom.name for atom in self.atoms])


# #==============================================================================
# # Testing
# #==============================================================================
# at = Atom('C', [0, 0, 0])
# ab = Atom('Re', [0, 0, 0], tag='Re2')
# ac = Atom('Re', [0, 0, 0], tag='Re1')
# uc = Unitcell(1, 2, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
# print(uc)
# print(set([at, ab, ac]))
#
# mtz = MTZ_model(uc, atoms=[at, ab, ac])
# mtz.load_from_file('C:\\Users\\Liam\\Dropbox\\Programming\\Python\\LEED-PyV\\phaseshifts\\test\\Re0001\\cluster_Re_bulk.i')
# print(mtz.get_elements())
# mtz.load_from_file('C:\\Users\\kss07698\\Desktop\\test_cluster.bak.i')
# mtz.gen_input(filename='C:\\Users\\kss07698\\Desktop\\test_cluster.bak.i')
