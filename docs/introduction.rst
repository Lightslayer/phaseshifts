.. _introduction:

************
Introduction
************

This package is a Python package which produces atomic phase shifts for
various LEED packages (including CLEED), as well as for certain XPD packages.

Currently, it uses the Barbieri/Van Hove phase shift calculation package and
preliminary support for John Rundgren's EEASiSSS package as backends
that are wrapped using f2py and a few Python modules to provide
a user-friendly object-orientated interface for the end user.

The aim of this package is to both automate and simplify the generation of
phase shift files in a manner that is easy for the computational hitch-hiker,
but powerful for those that wish to extend the package for particular needs.

The :ref:`phsh` script unifies many of the steps needed for the phase
shift calculations into a single command intended for the end-user. For more
information please read the documentation at
`<http://pythonhosted.org//phaseshifts/>`_

-------------------------
Barbieri/Van Hove backend
-------------------------

The original phase shift package developed by A. Barbieri & M. A. Van Hove
during the 1970's & 1980's. To quote the authors' site:

"The phase shift calculation is performed in several steps:

1. Calculation of the radial charge density for a free atom.

2. Calculation of the radial muffin-tin potential for atoms embedded in a
   surface defined by the user (the surface is represented by a slab that
   is periodically repeated in 3 dimensions, within vacuum between the
   repeated slabs); various approximations to the exchange potential
   are available; relativistic effects are taken into account.

3. Calculation of phase shifts from the muffin-tin potential.

4. Elimination of pi-jumps in the energy dependence of the phase shifts."

.. note:: You can get the original Fortran source (& learn more about the
   *phshift* programs) from `here
   <http://www.icts.hkbu.edu.hk/surfstructinfo/SurfStrucInfo_files/leed/leedpack.html>`_

Please contact Michel Van Hove <vanhove@hkbu.edu.hk> regarding this package.

A local copy of the source files can be found under ``phaseshifts/lib/.phsh.orig/phsh[0-2].f``.

----------------
EEASiSSS backend
----------------

package [#]_ [#]_ [#]_ can also calculate phase shifts and is used
in a number of recent works on oxides [#]_ [#]_. In the words of the package's
author, John Rundgren, the main qualifications of the program are:

+ The program accepts three sources of atomic potentials:

    1. E. L. Shirley�s atomic program [#]_ applied together with Mattheiss�s
    superposition method.

    2. The DFT-SCF program of V. Eyert using the full-potential
    Augmented Spherical Wave method [#]_ .

    3. The DFT-SCF program `WIEN2k <http://www.wien2k.at/>`_ using the
    full-potential Augmented Plane Wave method.

+ The exchange-correlation interaction between the scattered electron and the
  crystal�s electron gas generates an energy-dependent inner potential.
  The phase shifts are referred to the in-crystal kinetic energy, and it is
  supposed that an associated LEED code uses the same standard.
+ The crystal potential is everywhere continuous so as to exclude fortuitous
  standing-wave electron resonances in the muffin-tin spheres and pertaining
  fortuitous wobblings in the phase shift versus energy curves.
+ The optimization of the muffin-tin radii is made using the method of
  `Differential Evolution <http://www.physik.uni-augsburg.de/~eyert/ASWhome.shtml/>`_
  , an extremely efficient minimizer.

.. note:: Please refer to the short :ref:`EEASiSSS_Phase_Shift_Package_Guide` for more information
          relating to using EEASiSSS.

.. [#] J Rundgren, Phys. Rev. B 68, 125405 (2003).
.. [#] J Rundgren, Phys. Rev. B 76, 195441 (2007).
.. [#] E. A. Soares, C. M. C. De Castillho, and V. E. Carvalho, J. Phys.: Condens. Matter
   23,303001 (2011).
.. [#] R. Pentcheva, W. Moritz, J. Rundgren, S. Frank, D. Schrupp, and M. Scheffler, Surf. Sci
   602, 1299 (2008).
.. [#] V.B. Nascimento, R.G. Moore, J. Rundgren, J. Zhang, L. Cai, R. Jin, D.G. Mandrus,
   and E.W. Plummer, Phys. Rev. B 75, 035408 (2007).
.. [#] S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles, and C. W. Clark, Phys. Rev. B
   55, 191 (1997).
.. [#] R Storn and K. Price, J. Global Optimization 11, 341 (1997).

Please contact John Rundgren <jru@kth.se> for queries, comments or suggestions
related to this package.

.. include:: included/about_the_code.rst

.. include:: included/alternatives.rst
