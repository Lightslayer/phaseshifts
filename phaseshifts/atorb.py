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
**atorb.py**

Provides convenience functions for generating input and calculating
atomic charge densities for use with the Barbieri/Van Hove phase
shift calculation package.

:See: http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/

:Requires: f2py (for libphsh fortran wrapper generation)

.. note::
   To generate libphsh fortran wrappers (libphsh.pyd) for your platform
   then use 'python setup.py' in the lib directory of this package to
   install into your python distribution. Alternatively, use::

     f2py -c -m libphsh libphsh.f

   Windows users may have to add appropriate compiler switches, e.g. ::

    # 32-bit
    f2py -c -m libphsh --fcompiler=gfortran --compiler=mingw-32 libphsh.f

    # 64-bit
    f2py -c -m libphsh --fcompiler=gfortran --compiler=mingw-64 libphsh.f

"""

from configparser import ConfigParser
import os
import re
from collections import OrderedDict
import sys

from phaseshifts.utils import expand_filepath

# TODO: Clean this up when the project officially drops support for Python 2.7
try:
    from typing import Optional  # noqa: F401
    from typing_extensions import Literal

    ElementBackendType = Literal["mendeleev", "elementy", "periodictable"]
except ModuleNotFoundError:
    # we are likely on an old version of python
    ElementBackendType = str  # type: ignore

try:
    import periodictable  # type: ignore [import-untyped]
except ImportError:
    periodictable = None

try:
    import elementy  # type: ignore [import-not-found]
except ImportError:
    elementy = None

try:
    import mendeleev  # type: ignore [import-not-found]
except ImportError:
    mendeleev = None  # Need to install mendeleev

from phaseshifts import elements

elements_dict = OrderedDict(
    [
        ("H", "Hydrogen"),
        ("He", "Helium"),
        ("Li", "Lithium"),
        ("Be", "Beryllium"),
        ("B", "Boron"),
        ("C", "Carbon"),
        ("N", "Nitrogen"),
        ("O", "Oxygen"),
        ("F", "Fluorine"),
        ("Ne", "Neon"),
        ("Na", "Sodium"),
        ("Mg", "Magnesium"),
        ("Al", "Aluminium"),
        ("Si", "Silicon"),
        ("P", "Phosphorus"),
        ("S", "Sulfur"),
        ("Cl", "Chlorine"),
        ("Ar", "Argon"),
        ("K", "Potassium"),
        ("Ca", "Calcium"),
        ("Sc", "Scandium"),
        ("Ti", "Titanium"),
        ("V", "Vanadium"),
        ("Cr", "Chromium"),
        ("Mn", "Manganese"),
        ("Fe", "Iron"),
        ("Co", "Cobalt"),
        ("Ni", "Nickel"),
        ("Cu", "Copper"),
        ("Zn", "Zinc"),
        ("Ga", "Gallium"),
        ("Ge", "Germanium"),
        ("As", "Arsenic"),
        ("Se", "Selenium"),
        ("Br", "Bromine"),
        ("Kr", "Krypton"),
        ("Rb", "Rubidium"),
        ("Sr", "Strontium"),
        ("Y", "Yttrium"),
        ("Zr", "Zirconium"),
        ("Nb", "Niobium"),
        ("Mo", "Molybdenum"),
        ("Tc", "Technetium"),
        ("Ru", "Ruthenium"),
        ("Rh", "Rhodium"),
        ("Pd", "Palladium"),
        ("Ag", "Silver"),
        ("Cd", "Cadmium"),
        ("In", "Indium"),
        ("Sn", "Tin"),
        ("Sb", "Antimony"),
        ("Te", "Tellurium"),
        ("I", "Iodine"),
        ("Xe", "Xenon"),
        ("Cs", "Cesium"),
        ("Ba", "Barium"),
        ("La", "Lanthanum"),
        ("Ce", "Cerium"),
        ("Pr", "Praseodymium"),
        ("Nd", "Neodymium"),
        ("Pm", "Promethium"),
        ("Sm", "Samarium"),
        ("Eu", "Europium"),
        ("Gd", "Gadolinium"),
        ("Tb", "Terbium"),
        ("Dy", "Dysprosium"),
        ("Ho", "Holmium"),
        ("Er", "Erbium"),
        ("Tm", "Thulium"),
        ("Yb", "Ytterbium"),
        ("Lu", "Lutetium"),
        ("Hf", "Hafnium"),
        ("Ta", "Tantalum"),
        ("W", "Tungsten"),
        ("Re", "Rhenium"),
        ("Os", "Osmium"),
        ("Ir", "Iridium"),
        ("Pt", "Platinum"),
        ("Au", "Gold"),
        ("Hg", "Mercury"),
        ("Tl", "Thallium"),
        ("Pb", "Lead"),
        ("Bi", "Bismuth"),
        ("Po", "Polonium"),
        ("At", "Astatine"),
        ("Rn", "Radon"),
        ("Fr", "Francium"),
        ("Ra", "Radium"),
        ("Ac", "Actinium"),
        ("Th", "Thorium"),
        ("Pa", "Protactinium"),
        ("U", "Uranium"),
        ("Np", "Neptunium"),
        ("Pu", "Plutonium"),
        ("Am", "Americium"),
        ("Cm", "Curium"),
        ("Bk", "Berkelium"),
        ("Cf", "Californium"),
        ("Es", "Einsteinium"),
        ("Fm", "Fermium"),
        ("Md", "Mendelevium"),
        ("No", "Nobelium"),
        ("Lr", "Lawrencium"),
        ("Rf", "Rutherfordium"),
        ("Db", "Dubnium"),
        ("Sg", "Seaborgium"),
        ("Bh", "Bohrium"),
        ("Hs", "Hassium"),
        ("Mt", "Meitnerium"),
        ("Ds", "Darmstadtium"),
        ("Rg", "Roentgenium"),
        ("Cn", "Copernicium"),
        ("Uut", "Ununtrium"),
        ("Fl", "Flerovium"),
        ("Uup", "Ununpentium"),
        ("Lv", "Livermorium"),
        ("Uus", "Ununseptium"),
        ("Uuo", "Ununoctium"),
    ]
)


def get_element(
    element, backend=None
):  # type: (str, Optional[ElementBackendType]) -> object
    """Obtain an element object for querying information using `backend`."""
    ele_obj = elements.ELEMENTS.get(element)
    if mendeleev and not ele_obj and backend in ("mandeleev", None):
        elements_data = mendeleev.get_all_elements()
        elements_data = {e.protons: e for e in elements_data}
        elements_data.update({e.symbol: e for e in elements_data})
        elements_data.update({e.name: e for e in elements_data})
        ele_obj = elements_data.get(element)
    if elementy and not ele_obj and backend in ("elementy", None):
        periodic_table = elementy.PeriodicTable()
        elements_data = {e.protons: e for e in periodic_table.elements}
        elements_data.update({e.symbol: e for e in periodic_table.elements})
        elements_data.update({e.name: e for e in periodic_table.elements})
        ele_obj = elements_data.get(element)
    if periodictable and not ele_obj and backend in ("periodtable", None):
        ele_obj = getattr(periodictable, element, None) or getattr(
            periodictable, str(element).lower()
        )
    if ele_obj is None:
        raise LookupError("Unable to match element {}".format(element))
    return ele_obj


def get_electron_config(element_obj):
    """Obtain the electronic orbital configuration for the given `element_obj`."""
    electron_config = None
    if hasattr(element_obj, "orbitals"):
        electron_config = " ".join(
            ["{orbital}{electrons}".format(**x) for x in element_obj.orbitals]
        )
    else:
        electron_config = getattr(element_obj, "eleconfig", None) or getattr(
            element_obj, "econf", None
        )
    return electron_config


class Atorb(object):
    r"""
    Description
    -----------
    A python wrapper for the `atorb` program by Eric Shirley for use in
    calculating atomic scattering for different elements

    Notes
    -----
    Original author: Eric Shirley

    There are nr grid points, and distances are in Bohr radii
    :math:`a_0 \simeq 0.539 \mathrm{\AA}`

    :math:`r(i) = r_{min} \cdot (r_{max} / r_{min})^{(i/n_r)}`,
    :math:`i=1,2,3,...n_r-1,n_r`

    The orbitals are stored in phe(), first index goes :math:`1...n_r`, the
    second index is the orbital index (:math:`i...n_{el}`)

    Look at the atomic files after printing this out to see everything...
    Suffice it to say, that the charge density at radius :math:`r(i)`
    in units of electrons per cubic Bohr radius is given by:

    :math:`\sum_{j-1}^{n_el}{occ(j) \cdot phe(i,j)^2 / (4.0\,\pi\,{r(i)^2)}}`

    Think of the phe functions as plotting the radial wave-functions
    as a function of radius on a logarithmic mesh...

    The Dirac equation is solved for the orbitals, whereas their density
    is treated by setting :math:`phe(i,j)` to Dirac's
    :math:`\sqrt{F(i,j)^2 + G(i,j)^2}` times the sign of :math:`G(i,j)`...

    So we are doing Dirac-Fock, except that we are not treating exchange
    exactly, in terms of working with major and minor components of the
    orbitals, and the phe's give the CORRECT CHARGE DENSITY...

    The above approximation ought to be very small for valence states,
    so you need not worry about it...

    The Breit interaction has been neglected altogether...it should not
    have a huge effect on the charge density you are concerned with...

    """

    atlib = "$ATLIB" if "ATLIB" in os.environ else "~/atlib/"
    userhome = "~/hf.conf"
    datalib = (
        os.path.join(os.environ["APPDATA"], "phaseshifts")
        if sys.platform.startswith("win")
        else "~/.phaseshifts"
    )

    def __init__(
        self,
        ngrid=1000,
        rel=True,
        exchange=0.0,
        relic=0,
        mixing_SCF=0.05,
        tolerance=0.0005,
        xnum=100,
        ifil=0,
        **kwargs,
    ):
        """
        Constructor
        """
        # set private data members
        self.ngrid = ngrid if isinstance(ngrid, int) else 1000
        self.rel = rel if isinstance(rel, bool) else True
        self.exchange = (
            exchange
            if isinstance(exchange, float) or isinstance(exchange, int)
            else 0.0
        )
        self.relic = relic if isinstance(relic, int) else 0
        self.mixing_SCF = mixing_SCF if isinstance(mixing_SCF, float) else 0.05
        self.tolerance = tolerance if isinstance(tolerance, float) else 0.0005
        self.xnum = xnum if isinstance(xnum, int) else 100

        # set other (compatibility) kwargs
        self.ifil = ifil if isinstance(ifil, int) else 0
        self.fmt = kwargs.get("fmt") or "bhv"
        self.__dict__.update(kwargs)

    @property
    def ngrid(self):
        """Returns the number of points in the radial charge grid"""
        return self._ngrid

    @ngrid.setter
    def ngrid(self, ngrid):
        """
        Description
        -----------
        Sets the number of points in the radial charge grid

        Parameters
        ----------
        ngrid : int
            Number of points in the radial grid.
        """
        try:
            self._ngrid = abs(int(ngrid))
        except ValueError:
            pass

    @property
    def rel(self):
        """Returns boolean value of whether to consider relativistic effects"""
        return self._rel

    @rel.setter
    def rel(self, rel):
        """Sets flag to consider relativistic effects"""
        self._rel = rel == "rel" or rel is True or rel == 1

    @property
    def exchange(self):
        """
        Returns the exchange correlation value,
        where 0.0=Hartree-Fock, 1.0=LDA or <float>=-alpha"""
        return self._exchange

    @exchange.setter
    def exchange(self, exchange):
        """
        Description
        -----------
        Sets the exchange correlation value.

        Parameters
        ----------
        exchange : float
            Exchange value for calculation. The value determines the
            calculation method, where :code:`0.0` = Hartee-Fock,
            :code:`1.0` = LDA or :code:`<float>` = -alpha .
        """
        try:
            self._exchange = float(exchange)
        except ValueError:
            pass

    @property
    def tolerance(self):
        """Returns the eigenvalue tolerance"""
        return self._tolerance

    @tolerance.setter
    def tolerance(self, tolerance):
        """
        Description
        -----------
        Sets the eigenvalue tolerance

        Parameters
        ----------
        tolerance : float
            The eigenvalue tolerance value to set to.
        """
        try:
            self._tolerance = float(tolerance)
        except ValueError:
            pass

    @property
    def relic(self):
        """Returns the relic value for calculation"""
        return self._relic

    @relic.setter
    def relic(self, relic):
        """
        Description
        -----------
        Sets the relic value for calculation.

        Parameters
        ----------
        relic : int
            Relic value for calculation.
        """
        try:
            self._relic = int(relic)
        except ValueError:
            pass

    @property
    def mixing_SCF(self):
        """Returns the self-consisting field value"""
        return self._mixing_SCF

    @mixing_SCF.setter
    def mixing_SCF(self, mixing):
        """Sets the self-consisting field value"""
        try:
            self._mixing_SCF = mixing
        except ValueError:
            pass

    @property
    def xnum(self):
        """Returns xnum value"""
        return self._xnum

    @xnum.setter
    def xnum(self, xnum):
        """
        Description
        -----------
        Sets the xnum value.

        Parameters
        ----------
        xnum : float
            ???
        """
        try:
            self._xnum = float(xnum)
        except ValueError:
            pass

    def gen_conf_file(self, conf_file="hf.conf"):
        """
        Description
        -----------
        Generates conf file from Atorb() object instance.

        Parameters
        ----------
        conf_file : str
            Filepath for conf output file (default: 'hf.conf').

        """
        conf_file = expand_filepath(conf_file)
        if not os.path.isdir(os.path.dirname(conf_file)):
            os.makedirs(os.path.dirname(conf_file))

        config = ConfigParser(allow_no_value=True)
        config.set("DEFAULT", "# parameters common to all backends")
        config.set("DEFAULT", "ngrid", str(self.ngrid))
        config.set("DEFAULT", "rel", str(self.ngrid))
        config.set("DEFAULT", "exchange", str(self.exchange))
        config.set("DEFAULT", "exchange", str(self.exchange))
        config.set("DEFAULT", "relic", str(self.relic))
        config.set("DEFAULT", "mixing_SCF", str(self.mixing_SCF))
        config.set("DEFAULT", "tolerance", str(self.tolerance))
        config.set("DEFAULT", "xnum", str(self.xnum))

        # write to file
        header = "hartfock config file"
        with open(conf_file, "w") as f:
            f.write(str("#").ljust(len(header) + 3, "#") + "\n")
            f.write("# {} #\n".format(header))
            f.write(str("#").ljust(len(header) + 3, "#") + "\n")
            config.write(f)

    def update_config(self, conf):
        """
        Description
        -----------
        Updates Atorb() instance with arguments found from ``conf``.

        Parameters
        ----------
        conf : str or dict
            Either filepath to the user-specified ``*.conf`` file containing
            the atomic charge density calculation parameters or else
            a dictionary of the keyword arguments to update.

        Raises
        ------
        ValueError if ``conf`` is neither a str or dict instance.
        """
        if isinstance(conf, str):
            if os.path.isfile(conf):
                self.__dict__.update(self._get_conf_parameters(conf))
        elif isinstance(conf, dict):
            self.__dict__.update(conf)
        else:
            raise ValueError(
                "conf argument '{}' is neither a " "str or dict instance".format(conf)
            )

    def _get_conf_lookup_dirs(self):
        """
        Description
        -----------
        Returns a list of lookup locations for configuration files.

        Locations include (in order):

                1. :envvar:`ATLIB` or ``~/atlib/``
                2. ``~/`` or ``%USERPROFILE%/``
                3. ``~/.phaseshifts/`` or ``%APPDATA%/phaseshifts``
                4. './'

        """
        filenames = [
            os.path.join(directory, "hf.conf")
            for directory in list(
                self.atlib, self.userhome, self.datalib, os.path.curdir
            )
        ]
        return [os.path.abspath(expand_filepath(f)) for f in filenames]

    def _get_conf_parameters(self, conf_file="hf.conf"):
        """
        Description
        -----------
        Reads ``*.conf`` file for Atorb.gen_input() user-specified defaults and
        returns a dictionary of the relevant keyword arguments.

        Parameters
        ----------
        conf_file : str
            Path to ``*.conf`` file to read from. If the file does not exist
            then the function will attempt to read ``hf.conf`` from normal
            lookup storage locations, including (in order):

                1. :envvar:`ATLIB` or ``~/atlib/``
                2. ``~/`` or ``%USERPROFILE%/``
                3. ``~/.phaseshifts/`` or ``%APPDATA%/phaseshifts``
                4. './'

        Returns
        -------
        Dictionary of Atorb.gen_input() keyword arguments.

        """
        conf_file = expand_filepath(conf_file)

        config = ConfigParser()
        config.read(list(conf_file) + self._get_conf_lookup_dirs())

        return config.items("DEFAULT")

    @staticmethod
    def get_quantum_info(shell):  # (str) -> Tuple[int|float|List[int|float], ...]
        r"""
        Description
        -----------
        Get a tuple of quantum information for a given orbital 's', 'p', 'd'
        or 'f' from a given subshell string.

        Returns
        -------
        tuple : (int, int, list[float, float], list[float, float])
            (n, l, j=[l-s, l+s], occ=[:math:`n^-_r`, :math:`n^+_r`])

        Notes
        -----
        - *n* is the principle quantum number (:math:`n > 0`).
        - *l* is the azimuthal quantum number (:math:`0 \leq l \leq n-1`).
        - *s* is the spin quantum number (:math:`s \pm \frac{1}{2}`).
        - *j* is the total angular momentum quantum numbers for both
          :math:`l-s` or :math:`l+s`, respectively.
        - :math:`n_r` is the occupancy of the spin-split :math:`l-s`
          and :math:`l+s` levels, respectively.

        Example
        -------
        >>> Atorb.get_quantum_info('3d6')
        (3, 2, [1.5, 2.5], [2.4, 3.6])

        """

        subshell = "".join([s for s in shell if s.isalpha()])
        try:
            (n, nelectrons) = [
                t(s) for t, s in zip((int, int), shell.replace(subshell, " ").split())
            ]
        except ValueError:  # assume 1 electron in shell
            n = int(shell.replace(subshell, " ").split()[0])
            nelectrons = 1

        s = 0.5
        shell_info = None
        if subshell == "s":
            l = 0
            occ = [nelectrons / 1.0]
            j = [l + s]
            shell_info = (n, l, j, occ)
        elif subshell == "p":
            # 3 subshells
            l = 1
            max_occ = 6
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            shell_info = (n, l, [l - s, l + s], occ)
        elif subshell == "d":
            # 5 subshells
            l = 2
            max_occ = 10
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            shell_info = (n, l, [l - s, l + s], occ)
        elif subshell == "f":
            # 7 subshells!
            l = 3
            max_occ = 14
            occ = []
            for j in [l - s, l + s]:
                occ.append(((2.0 * j) + 1) * nelectrons / max_occ)
            shell_info = (n, l, [l - s, l + s], occ)
        else:
            raise NotImplementedError(
                "Exotic sub-shells beyond f-block have not been implemented"
            )
        return shell_info

    @staticmethod
    def replace_core_config(electron_config):
        """
        Description
        -----------
        Replace nobel gas core with equivalent electronic shell configuration

        Parameters
        ----------
        electron_config : str
            String containing the electronic configuration of the given
            element.

        Returns
        -------
        str :
            A substituted string where the nobel gas core has been replaced.

        Examples
        --------
        >>> Atorb.replace_core_config('[Ar] 4s2')
         '1s2 2s2 2p6 3s2 3p6 4s2'

        >>> Atorb.replace_core_config('[Xe] 6s2 5d1')
         '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6 6s2 5d1'

        """
        cores = {
            "[He]": "1s2",
            "[Ne]": "1s2 2s2 2p6",
            "[Ar]": "1s2 2s2 2p6 3s2 3p6",
            "[Kr]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6",
            "[Xe]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6",
            "[Rn]": "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2 4d10 5p6" "4f14 5d10 6s2 6p6",
        }
        core = electron_config.split()[0]

        if core in cores:
            config = electron_config.replace(core, cores.get(core))
        else:
            config = electron_config
        return config

    def _gen_input(self, element, conf_file=None):
        if conf_file is not None:
            self.update_config(self.get_conf_parameters(conf_file))

        Atorb.gen_input(
            element,
            ngrid=self.ngrid(),
            rel=self.rel(),
            exchange_method=self.exchange(),
            relic=self.relic(),
            mixing_SCF=self.mixing_SCF(),
            tolerance=self.tolerance(),
            xnum=self.xnum(),
            atorb_file=(
                self.__dict__["atorb_file"] if "atorb_file" in self.__dict__ else None
            ),
            output=(self.__dict__["output"] if "output" in self.__dict__ else None),
            header=(self.__dict__["header"] if "header" in self.__dict__ else None),
            ifil=(0 if "ifil" not in self.__dict__ else int(self.__dict__["ifil"])),
            fmt=("vht" if "fmt" not in self.__dict__ else self.__dict__["fmt"]),
        )

    @staticmethod
    def gen_input(element, **kwargs):
        """
        Description
        -----------
        Generate hartfock atorb input file from <element> and optional **kwargs
        arguments. The generated input can then be inputted into
        :py:meth:`Atorb.calculate_Q_density()`.

        Parameters
        ----------
        element : int or str
            Either the atomic number, symbol or name for a given element
        output : str, optional
            File string for atomic orbital output (default: 'at_<symbol>.i')
        ngrid : int, optional
            Number of points in radial grid (default: 1000)
        rel : bool, optional
            Specify whether to consider relativistic effects
        filename : str, optional
            Name for generated input file (default: 'atorb')
        header : str, optional
            Comment at beginning of input file
        method : str, optional
            Exchange correlation method using either 0.0=Hartree-Fock,
            1.0=LDA, -alpha = xalpha (default: 0.0)
        relic : float, optional
            Relic value for calculation (default: 0)
        mixing_SCF : float, optional
            Self consisting field value (default: 0.5)
        tolerance : float, optional
            Eigenvalue tolerance (default: 0.0005)
        xnum : float, optional
            ??? (default: 100)
        ifil : int, optional
            flag to read :code:`vpert` array from :file:`vvalence` - possibly
            redundant. Only used when fmt='rundgren' or 'eeasisss' (default: 0)
        fmt : str, optional
            Format of generated atorb input file; can be either 'vht' for the
            van Hove-Tong package or 'rundgren' for the EEASiSSS package
            (default: 'vht')

        Returns
        -------
        Filename of input file once generated or else instance of StringIO
        object containing written input text.

        Notes
        -----
        output can also be a StringIO() object to avoid saving to file.

        Example
        -------
        >>> Atorb.gen_input('H',rel=False,filename="atorb.txt",output='at_H.i')
        >>> with open('atorb_H.txt', 'r') as f: print("".join(f.readlines())
        C*********************************************************************
        C  atorb input file: atorb_H.txt.
        C*********************************************************************
        i
        1 1000                   ! Z NR (number of points in radial grid)
        d
        0                        ! 1=rel, 0=n.r.
        x
        0.d0                     ! 0.d0=HF, 1.d0=LDA, -alfa = xalfa...
        a
        0 1 0.5 0.0005 100       ! relic,levels,mixing SCF, eigen. tol,for ech
        1 0 0 -0.5 1 1.0         ! n, l, l, -j, <1>, occupation
        w
        at_H.i
        q


        """
        ele = get_element(element, backend=None)
        Z = ele.protons

        # get full electronic configuration
        electron_config = get_electron_config(ele)
        config = Atorb.replace_core_config(electron_config)

        # get quantum numbers & occupancy for each electronic orbital in atom
        electrons = []
        nlevels = 0
        for shell in config.split():
            (n, l, J, occ) = Atorb.get_quantum_info(shell)
            for i, j in enumerate(J):
                electrons.append((n, l, l, -j, 1, occ[i]))
                nlevels += 1

        # test kwargs and generate output arguments
        if "output" in kwargs:
            output = kwargs.get("output")
        else:
            output = "at_{0}.i".format(ele.symbol)

        if "ngrid" in kwargs:
            NR = kwargs.get("ngrid")
        else:
            NR = 1000  # default grid resolution

        if "rel" in kwargs:
            rel = int(kwargs.get("rel"))
        else:
            rel = 1  # default is relativistic

        if "filename" in kwargs:
            filename = kwargs.get("filename")
        else:
            filename = "atorb_{0}.txt".format(ele.symbol)

        if "header" in kwargs:
            header = kwargs.value("header")
        else:
            header = "  atorb input file: {0}.".format(os.path.basename(filename))

        if "method" in kwargs:
            method = str(kwargs.value("method"))
            methods_dict = {"HF": "0.d0", "LDA": "1.d0", "xalpha": "-alpha"}
            if method in methods_dict:
                method = methods_dict.get("method")
            else:
                method = "0.d0"
        else:
            method = "0.d0"

        if "relic" in kwargs:
            relic = float(kwargs.get("relic"))
        else:
            relic = 0

        if "mixing_SCF" in kwargs:
            mixing_SCF = float(kwargs.get("mixing_SCF"))
        else:
            mixing_SCF = 0.5

        if "tolerance" in kwargs:
            eigen_tol = float(kwargs.get("tolerance"))
        else:
            eigen_tol = 0.0005

        if "ech" in kwargs:
            ech = int(kwargs.get("ech"))
        else:
            ech = 100

        fmt = kwargs.get("fmt") or "bvh"

        # produce output file
        atorb_file = kwargs["atorb_file"]
        if isinstance(atorb_file, str):
            f = open(atorb_file, "w")
        else:
            f = atorb_file

        ifil_str = ", ifil"
        if fmt.lower() not in ["rundgren", "eeasisss"]:
            ifil = ""
            ifil_str = ""

        if header is not None:
            f.write("!".ljust(70, "*") + "\n")
            f.write("! " + str(header) + "\n")
            f.write("!".ljust(70, "*") + "\n")
        elif f.tell() == 0:
            f.write("!".ljust(70, "*") + "\n")
            f.write(
                "! %s hartfock input auto-generated by phaseshifts\n"
                % (fmt.upper() if fmt.lower() in ["vht", "eeasisss"] else fmt.title())
            )
            f.write("!".ljust(70, "*") + "\n")

        f.write("i\n")
        if fmt.lower() in ["rundgren", "eeasisss"]:
            # add line for element symbol
            f.write("{}\n".format(ele.symbol))
        f.write(
            "{0} {1}".format(Z, int(kwargs["ngrid"])).ljust(30, " ")
            + "! Z NR (number of points in radial grid)\n"
        )
        f.write("d\n")
        f.write("{0}".format(int(rel)).ljust(30) + "! 1=rel, 0=n.r.\n")
        f.write("x\n")
        f.write(
            "{0}".format(method).ljust(30)
            + "! 0.d0=HF, 1.d0=LDA, -(float) = xalfa...\n"
        )
        f.write("a\n")
        f.write(
            "{} {} {} {} {} {}".format(
                relic, nlevels, mixing_SCF, kwargs["tolerance"], kwargs["xnum"], ifil
            ).ljust(30)
            + "! relic, levels, mixing SCF, eigen. tolerance, xnum%s\n" % ifil_str
        )
        for i in range(0, nlevels):
            f.write(
                "{0} {1} {2} {3} {4} {5}".format(*electrons[i]).ljust(30)
                + "! n, l, l, -j, <1>, occupation\n"
            )
        f.write("w\n")
        f.write("{0}\n".format(output))
        f.write("q\n")

        return filename  # return output filename for further use

    @staticmethod
    def calculate_Q_density(**kwargs):
        """
        Description
        -----------
        Calculate the radial charge density of a given element or atorb input
        file.

        Usage
        -----
        Atorb.calculate_Q_density(**kwargs)

        Parameters
        ----------
        kwargs may be any of the following.

        element : int or str, optional
            Generate element atorb input file on the fly. Additional
            kwargs may be used to govern the structure of the input
            file - please use ``help(phaseshifts.Atorb.gen_input)``
            for more information.
        input : str, optional
            Specify atorb input file otherwise will use the class
            instance value.
        output_dir : str, optional
            Specify the output directory for the `at_*.i` file
            generated, otherwise the default current working directory
            is used.

        Returns
        -------
        str : filename
            Path to calculated charge density file or :code:`None` if failed.

        Examples
        --------
        >>> Atorb.calculate_Q_density(input='atorb_C.txt')
                18.008635    -33.678535
                4.451786    -36.654271
                1.569616    -37.283660
                0.424129    -37.355634
                0.116221    -37.359816
                0.047172    -37.360317
                0.021939    -37.360435
                0.010555    -37.360464
                0.005112    -37.360471
                0.002486    -37.360473
                0.001213    -37.360473
                0.000593    -37.360473
                0.000290    -37.360474
            N L M J S OCC.
            1   0 0  -1/2   1    2.0000        -11.493862
            2   0 0  -1/2   1    2.0000         -0.788618
            2   1 1  -1/2   1    0.6667         -0.133536
            2   1 1  -3/2   1    1.3333         -0.133311
        TOTAL ENERGY =      -37.360474  -1016.638262

        >>> Atorb.calculate_Q_density(element='H')
                0.500007     -0.343752
                0.152392     -0.354939
                0.065889     -0.357254
                0.028751     -0.357644
                0.012732     -0.357703
                0.005743     -0.357711
                0.002641     -0.357712
                0.001236     -0.357713
                0.000587     -0.357713
                0.000282     -0.357713
        N L M J S OCC.
            1   0 0  -1/2   1    1.0000         -0.229756
        TOTAL ENERGY =       -0.357713     -9.733932


        """
        inp = None
        if "input" in kwargs:
            inp = os.path.abspath(kwargs.pop("input"))

        if "element" in kwargs:
            inp = os.path.abspath(
                Atorb.gen_input(element=kwargs.pop("element"), **kwargs)
            )

        current_dir = os.path.curdir
        if "output_dir" in kwargs:
            output_dir = kwargs.pop("output_dir")
            if os.path.isdir(output_dir):
                os.chdir(output_dir)
            else:
                try:
                    os.makedirs(output_dir, exist_ok=True)
                    os.chdir(output_dir)
                except (OSError, IOError) as err:
                    raise IOError(
                        "Unable to create output directory due to {!r}".format(err)
                    )  # noqa

        if not inp:
            raise ValueError("Input file not specified")

        # do lazy loading due to documentation not needing compiled code
        import phaseshifts.lib.libphsh  # noqa

        phaseshifts.lib.libphsh.hartfock(
            inp
        )  # calculates atomic orbital charge densities for atom

        # get output filename
        lines = []
        output_filename = "atorb"

        try:
            with open(inp, "r") as f:
                lines = [line for line in f]
        except IOError:
            raise IOError

        lines = [
            line.replace("\n", "")
            .replace("\r", "")
            .lstrip(" ")
            .split("!")[0]
            .split("#")[0]
            .rstrip(" ")
            for line in lines
        ]

        for i in range(len(lines)):
            if lines[i].startswith("w"):
                output_filename = lines[i + 1]
                break

        os.chdir(current_dir)  # return to original directory

        return (
            os.path.join(output_dir, output_filename)
            if output_dir is not None
            else output_filename
        )


class EEASiSSSAtorb(Atorb):
    def __init__(self, ifil=0, **kwargs):
        # set higher default number of grid points for EEASiSSS
        kwargs["ngrid"] = 2000 if "ngrid" not in kwargs else kwargs
        kwargs["fmt"] = "eeasisss"
        Atorb.__init__(self, kwargs)
        self.ifil = int(ifil) if isinstance(ifil, bool) or isinstance(ifil, int) else 0

    @property
    def ifil(self):
        """
        Returns flag for reading :code:`vpert` array from file :file:`vvalence`
        """
        return self._ifil

    @ifil.setter
    def ifil(self, ifil):
        """
        Sets whether to read :code:`vpert` array from :file:`vvalence`
        """
        try:
            self._ifil = int(ifil)
        except ValueError:
            pass

    def gen_conf_file(
        self,
        conf_file=("$ATLIB/hf.conf" if "ATLIB" in os.environ else "~/atlib/hf.conf"),
    ):
        """
        Description
        -----------
        Generates hartfock conf file from EEASiSSSAtorb() object

        Parameters
        ----------
        conf_file : str
            Filepath for conf output file (default: '~/atlib/hf.conf').

        Examples
        --------
        >>> from phaseshifts.atorb import EEASiSSSAtorb
        >>> atorb = EEASiSSSAtorb()  # create an object instance
        >>> # create a config file in the default location
        >>> atorb.gen_conf_file()
        """
        conf_file = expand_filepath(conf_file)

        # call parent method
        Atorb.gen_conf_file(self, conf_file)

        # add new section specific to EEASiSSS
        config = ConfigParser(allow_no_value=True)
        config.add_section("EEASiSSS")
        config.set("EEASiSSS", "# parameters specific to EEASiSSS backend")
        config.set("EEASiSSS", "ngrid", str(self.ngrid))  # override of base
        config.set("EEASiSSS", "ifil", str(self.ifil))

        # append new configuration data
        with open(conf_file, "a") as f:
            config.write(f)

    def _get_conf_parameters(self, conf_file="~/atlib/hf.conf"):
        """
        Description
        -----------
        Reads ``*.conf`` file for Atorb.gen_input() user-specified defaults and
        returns a dictionary of the relevant keyword arguments.

        Parameters
        ----------
        conf_file : str
            Path to ``*.conf`` file to read from. If the file does not exist
            then the function will attempt to read ``hf.conf`` from normal
            storage locations, including (in order):

                1. :envvar:`ATLIB` or ``~/atlib/``
                2. ``~/`` or ``%USERPROFILE%/``
                3. ``~/.phaseshifts/`` or ``%APPDATA%/phaseshifts``
                4. './'

        Returns
        -------
        Dictionary of keyword arguments for :py:meth:`Atorb.gen_input()`.

        """
        conf_file = expand_filepath(conf_file)

        config = ConfigParser()
        config.read(list(conf_file) + self._get_conf_lookup_dirs())

        conf_dict = {}
        conf_dict.update(config.items("DEFAULT"))
        conf_dict.update(config.items("EEASiSSS"))
        return conf_dict

    def _gen_input(self, element, conf_file=None):
        """Internal bound version of :py:meth:`EEEASiSSSAtorb.gen_input`"""
        if conf_file is not None:
            self.update_config(self.get_conf_parameters(conf_file))

        EEASiSSSAtorb.gen_input(
            element,
            ngrid=self.ngrid(),
            rel=self.rel(),
            exchange_method=self.exchange(),
            relic=self.relic(),
            mixing_SCF=self.mixing_SCF(),
            tolerance=self.tolerance(),
            xnum=self.xnum(),
            ifil=self.ifil(),
            atorb_file=(
                "inputA"
                if "atorb_file" in self.__dict__
                else self.__dict__["atorb_file"]
            ),
            output=(self.__dict__["output"] if "output" in self.__dict__ else None),
            header=(self.__dict__["header"] if "header" in self.__dict__ else None),
            fmt=("eeasisss" if "fmt" not in self.__dict__ else self.__dict__["fmt"]),
        )

    @staticmethod
    def gen_input(elements=[], atorb_file="inputA", **kwargs):
        """
        Description
        -----------
        :py:class:`EEASiSSSAtorb` override of :py:class:`Atorb` base class
        method which produces an input file for a set of elements rather
        than just individual ones.

        Parameters
        ----------
        elements : list
            List of elements to include in the generated input file. Each
            element can be either the atomic number, symbol, name or an
            :py:class:`phaseshifts.Element` instance.
        output : str, optional
            File string for atomic orbital output (default: 'at_<symbol>.i')
        ngrid : int, optional
            Number of points in radial grid (default: 1000)
        rel : bool, optional
            Specify whether to consider relativistic effects (default: True)
        atorb_file : str, optional
            Name for generated input file (default: 'inputA')
        header : str, optional
            Comment at beginning of input file (default: None)
        method : str or float, optional
            Exchange correlation method using either 0.0=Hartree-Fock,
            1.0=LDA, -alpha = float (default: 0.0)
        relic : float, optional
            Relic value for calculation (default: 0)
        mixing_SCF : float, optional
            Self consisting field value (default: 0.5)
        tolerance : float, optional
            Eigenvalue tolerance (default: 0.0005)
        xnum : float, optional
            ??? (default: 100)
        ifil : int, optional
            flag to read :code:`vpert` array from :file:`vvalence` - possibly
            redundant. Only used when fmt='rundgren' or 'eeasisss' (default: 0)

        Returns
        -------
        Filename of input file once generated or else instance of StringIO
        object containing written input text. :code:`None` will be returned if
        the method failed to generate the input.

        Notes
        -----
        output can also be a StringIO() object to avoid saving to file.

        Examples
        --------
        >>> from phaseshifts.elements import ELEMENTS, SERIES
        >>>
        >>> # generate hartfock input file for InGaAs
        >>> InGaAs = ['In', 31, 'Arsenic']
        >>> EEASiSSS.gen_input(elements=InGaAs, atorb_file='~/atlib/InGaAs.hf')
        >>>
        >>> # generate hartfock input file for all halogens
        >>> halogens = [e for e in ELEMENTS if SERIES[e.series] == 'Halogens']
        >>> EEASiSSS.gen_input(halogens, atorb_file='$ATLIB/halogens.hf')
        >>>
        >>> # do likewise for all non-metals, but using hf.conf file parameters
        >>> non_metals = [e for e in ELEMENTS if SERIES[e.series] == 'Nonmetals']
        >>> EEASiSSS.gen_input(non_metals, atorb_file='./nonmetals.hf')
        """
        io = StringIO()
        successful = False
        try:
            kwargs = kwargs.pop("fmt") if "fmt" in kwargs else kwargs
            # generate buffer string of input for each element
            for element in set(elements):
                Atorb.gen_input(element, atorb_file=io, fmt="eeasisss", **kwargs)
            # write buffered string to disk
            with open(atorb_file, "w") as f:
                f.write(io.getvalue())
            successful = True
        except Exception as e:
            raise e
        finally:
            # clean up
            io.close()
        return atorb_file if successful else None

    @staticmethod
    def calculate_Q_density(
        elements=[],
        atorb_input="inputA",
        output_dir=(
            expand_filepath("ATLIB")
            if "ATLIB" in os.environ
            else expand_filepath("~/atlib/")
        ),
        **kwargs,
    ):
        """
        Description
        -----------
        :py:class:`EEASiSSS` override of
        :py:class:`Atorb.calculate_Q_density()`
        base method to produce hartfock input files for calculating the
        atomic charge densities of different elements.

        Parameters
        ----------
        elements : list of Element, int or str
            Generate atorb input file on the fly. If the list is empty then
            the function will return an empty list.
        atorb_input : str, optional
            Specifies the path to the atorb input file. If :code:`None` then
            a temporary file will be created.
        output_dir : str, optional
            Specifies the output directory for the `at_*.i` file
            generated. (default: :envar:`$ATLIB` or :file:`~/atlib`)
        subroutine : function, optional
            Specifies the hartfock function to use (default: eeasisss_hartfock)

            .. warning:: Do not modify the :option:`subroutine` value without
                         good cause - here be dragons!

        Returns
        -------
        List of filepaths to the calculated atomic charge density files.

        Notes
        -----
        This method implicitly calls :py:meth:`EEASiSSSAtorb.gen_input` for
        generating a suitable input file for the charge density calculations.

        For the patient, it may be worth having a initial one-time atomic
        charge density calculation for every element.
        This can be done like so::

            from phaseshifts.elements import ELEMENTS, SERIES
            from phaseshifts.atorb import EEASiSSSAtorb
            # calculate chgden* files for EVERY element
            # and place them in the default location
            EEASiSSSAtorb.calculate_Q_density(elements=ELEMENTS)

        .. note::
            The above calculation is very useful and can be customised using a
            user generated ``hf.conf`` file.
            For more details see
            :py:meth:`EEASiSSSAtorb.gen_conf_file()`

        Examples
        --------
        >>> from phaseshifts.elements import ELEMENTS, SERIES
        >>> from phaseshifts.atorb import EEASiSSSAtorb
        >>>
        >>> # generate hartfock input file for InGaAs
        >>> InGaAs = ['In', 31, 'Arsenic']
        >>> input = '~/atlib/InGaAs.hf'
        >>> EEASiSSSAtorb.calculate_Q_density(elements=InGaAs,
        >>>                                   atorb_input=input)
        >>>
        >>> # generate hartfock input file for all halogens
        >>> halogens = [e for e in ELEMENTS if SERIES[e.series] == 'Halogens']
        >>> input = '$ATLIB/halogens.hf'
        >>> EEASiSSSAtorb.calculate_Q_density(halogens, atorb_file=input)
        >>>
        >>> # do likewise for all non-metals, but using hf.conf file parameters
        >>> series = 'Nonmetals'
        >>> non_metals = [e for e in ELEMENTS if SERIES[e.series] == series]
        >>> input = './nonmetals.hf'
        >>> EEASiSSSAtorb.calculate_Q_density(non_metals, atorb_file=input)
        """
        # do not do anything if no elements given, otherwise get the set
        if elements == []:
            return []
        else:
            elements = set(
                [
                    ELEMENTS[element] if not isinstance(element, Element) else element
                    for element in elements
                ]
            )

        output_dir = expand_filepath(output_dir) or os.curdir
        atorb_input = expand_filepath(atorb_input) or os.path.join(
            gettempdir(), "".join([e.symbol for e in elements]) + ".hf"
        )

        EEASiSSSAtorb.gen_input(elements=elements, atorb_file=atorb_input, **kwargs)

        Atorb.calculate_Q_density(
            atorb_input=atorb_input,
            output_dir=output_dir,
            subroutine=eeasisss_hartfock,
            **kwargs,
        )

        elements = [
            ELEMENTS[element] if not isinstance(element, Element) else element
            for element in set(elements)
        ]

        return [
            os.path.join(output_dir, "chgden" + element.symbol)
            if output_dir != os.path.curdir
            else "chgden" + element.symbol
            for element in set(elements)
        ]


if __name__ == "__main__":
    atorb = EEASiSSSAtorb()
    atorb.gen_conf_file("~/atlib/hf.conf")
    for i in range(1, 112):
        atorb.calculate_Q_density(elements=[ELEMENTS[i].symbol])


def get_substr_positions(string, substring="\n"):
    return [m.start() for m in re.finditer(substring, string)]
