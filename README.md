# PHASESHIFTS PACKAGE

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Liam-Deacon/phaseshifts/HEAD?labpath=tutorial.ipynb)
![PyPI - Version](https://img.shields.io/pypi/v/phaseshifts?logo=pypi&logoColor=white)
![Python](https://img.shields.io/badge/python-2.7%20%7C%203.5%20--%203.11-blue?logo=python&logoColor=white)
![PyPI - Status](https://img.shields.io/pypi/status/phaseshifts?logo=pypi&logoColor=white)
![GitHub License](https://img.shields.io/github/license/Liam-Deacon/phaseshifts?logo=github)
[![Read the Docs](https://img.shields.io/readthedocs/phaseshifts?logo=readthedocs&logoColor=white)](https://phaseshifts.readthedocs.io/en/latest/)
![GitHub Release](https://img.shields.io/github/v/release/Liam-Deacon/phaseshifts?logo=github)
![GitHub (Pre-)Release Date](https://img.shields.io/github/release-date-pre/Liam-Deacon/phaseshifts?logo=github)
[![GitHub issues](https://img.shields.io/github/issues/Liam-Deacon/phaseshifts?logo=github)](https://github.com/Liam-Deacon/phaseshifts/issues)
![Codacy grade](https://img.shields.io/codacy/grade/ed204d9cc22f4b37b0ad6612b42e8b1e?logo=codacy)
![Code Climate technical debt](https://img.shields.io/codeclimate/tech-debt/Liam-Deacon/phaseshifts?logo=codeclimate)
![Code Climate maintainability](https://img.shields.io/codeclimate/maintainability/Liam-Deacon/phaseshifts?logo=codeclimate)
![Code Climate issues](https://img.shields.io/codeclimate/issues/Liam-Deacon/phaseshifts?logo=codeclimate)
![Codecov](https://img.shields.io/codecov/c/github/Liam-Deacon/phaseshifts?logo=codecov&logoColor=white)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/Liam-Deacon/phaseshifts/total?logo=github)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/phaseshifts?logo=pypi&logoColor=white)](https://pypi.org/project/phaseshifts/)
[![GitHub closed issues](https://img.shields.io/github/issues-closed/Liam-Deacon/phaseshifts?logo=github)](https://github.com/Liam-Deacon/phaseshifts/issues?q=is%3Aissue+is%3Aclosed+)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FLiam-Deacon%2Fphaseshifts.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FLiam-Deacon%2Fphaseshifts?ref=badge_shield)
<!-- ![Code Climate coverage](https://img.shields.io/codeclimate/coverage/Liam-Deacon/phaseshifts) -->
<!-- ![Codacy coverage](https://img.shields.io/codacy/coverage/phaseshifts) -->
<!-- ![Coveralls branch](https://img.shields.io/coverallsCoverage/github/Liam-Deacon/phaseshifts) -->

[![Test Package](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/test-package.yml/badge.svg)](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/test-package.yml)
[![GitHub Pages](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/publish-docs-to-github-pages.yaml/badge.svg)](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/publish-docs-to-github-pages.yaml)
[![Publish Package](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/publish-to-pypi.yaml/badge.svg)](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/publish-to-pypi.yaml)
[![Publish Docker Image(s)](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/publish-to-dockerhub.yaml/badge.svg)](https://github.com/Liam-Deacon/phaseshifts/actions/workflows/publish-to-dockerhub.yaml)

[![Star on GitHub](https://img.shields.io/github/stars/Liam-Deacon/phaseshifts.svg?style=social)](https://github.com/Liam-Deacon/phaseshifts/stargazers)
[![Watch on GitHub](https://img.shields.io/github/watchers/Liam-Deacon/phaseshifts.svg?style=social)](https://github.com/Liam-Deacon/phaseshifts/watchers)

This package is a Python package which produces elastic electron atom scattering
(EEAS) phase shifts for modelling within various LEED packages (including CLEED),
as well as for certain XPD packages.

The aim of this package is to both automate and simplify the generation of
phase shift files in a manner that is easy for the computational hitch-hiker,
but powerful for those that wish to extend the package for particular needs.
The eventual goal is to integrate with LEED-IV/XPD computational packages to
iteratively generate phase shifts during model optimisation.

## Backends

Currently, it primarily uses the Barbieri/Van Hove phase shift calculation package
as the calculation backend wrapped using f2py, however, there is also preliminary
support for John Rundgren's EEASiSSS package as an alternative backend.
Both backends benefit from a few Python modules to provide a more user-friendly
interface to automate a lot of the work for the end user.

### Barbieri/Van Hove backend

The original phase shift package developed by A. Barbieri & M. A. Van Hove
during the 1970's & 1980's. To quote the authors' site:

"The phase shift calculation is performed in several steps:

  1. Calculation of the radial charge density for a free atom.
  2. Calculation of the radial muffin-tin potential for atoms embedded in
     a surface defined by the user (the surface is represented by a slab
     that is periodically repeated in 3 dimensions, within vacuum between
     the repeated slabs); various approximations to the exchange
     potential are available; relativistic effects are taken into
     account.
  3. Calculation of phase shifts from the muffin-tin potential.
  4. Elimination of pi-jumps in the energy dependence of the phase
     shifts."

<!--lint disable no-unused-definitions-->
> [!NOTE]
> You can get the original Fortran source (& learn more about the
> *phshift* programs) from Michel Van Hove's LEED Calculation Home Page:
>
> <https://www.icts.hkbu.edu.hk/VanHove_files/leed/leedpack.html>
>
> A local copy of the source files can be found under
> `phaseshifts/lib/.phsh.orig/phsh[0-2].f`.
<!--lint enable no-unused-definitions-->

Please contact [Michel Van Hove <vanhove@hkbu.edu.hk>](mailto://vanhove@hkbu.edu.hk) regarding this package.

### EEASiSSS backend

The Elastic Electron-Atom Scattering in Solids and Surface Slabs (EEASiSSS)
package [^1] [^2] [^3] can also calculate phase shifts and is used
in a number of recent works on oxides [^4] [^5]. In the words of the package's
author, John Rundgren, the main qualifications of the program are:

+ The program accepts three sources of atomic potentials:
    1. E. L. Shirley's atomic program [^6] applied together with Mattheiss's
    superposition method.
    2. The DFT-SCF program of V. Eyert using the full-potential
    Augmented Spherical Wave method [^7].
    3. The DFT-SCF program [WIEN2k](http://www.wien2k.at/) using the
    full-potential Augmented Plane Wave method.

+ The exchange-correlation interaction between the scattered electron and the
  crystal's electron gas generates an energy-dependent inner potential.
  The phase shifts are referred to the in-crystal kinetic energy, and it is
  supposed that an associated LEED code uses the same standard.
+ The crystal potential is everywhere continuous so as to exclude fortuitous
  standing-wave electron resonances in the muffin-tin spheres and pertaining
  fortuitous wobblings in the phase shift versus energy curves.
+ The optimization of the muffin-tin radii is made using the method of
  [Differential Evolution](https://serp.ai/differential-evolution/),
  an extremely efficient minimizer.

>[!NOTE]
> A short EEASiSSS users guide is appended to the input template files
>`inputA` and `inputX` distributed with the program package.

[^1]: J Rundgren, Phys. Rev. B 68, 125405 (2003).
[^2]: J Rundgren, Phys. Rev. B 76, 195441 (2007).
[^3]: E. A. Soares, C. M. C. De Castillho, and V. E. Carvalho, J. Phys.: Condens. Matter
   23,303001 (2011).
[^4]: R. Pentcheva, W. Moritz, J. Rundgren, S. Frank, D. Schrupp, and M. Scheffler, Surf. Sci
   602, 1299 (2008).
[^5]: V.B. Nascimento, R.G. Moore, J. Rundgren, J. Zhang, L. Cai, R. Jin, D.G. Mandrus,
   and E.W. Plummer, Phys. Rev. B 75, 035408 (2007).
[^6]: S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles, and C. W. Clark, Phys. Rev. B
   55, 191 (1997).
[^7]: R Storn and K. Price, J. Global Optimization 11, 341 (1997).

Please contact John Rundgren <jru@kth.se> for queries, comments or suggestions
related to this package.


## Running

The <span class="title-ref">phsh.py</span> script (available after
installing the package) aims to simplify these steps with a single
command. For more information please read the documentation at
<http://pythonhosted.org/phaseshifts/> (latest PyPI release) or [GitHub
Pages](https://liam-deacon.github.io/phaseshifts/) (latest master)

The simplest and most reliable cross-platform way to run
<span class="title-ref">phsh.py</span> is through docker:

```bash
# obtain the image
docker pull ghcr.io/Liam-Deacon/phaseshifts:latest  # should only need to do this once

# run phsh.py via the docker image
docker run ghcr.io/Liam-Deacon/phaseshifts:latest  # will display usage

# or more generally (adjust as needed)
docker run ghcr.io/Liam-Deacon//phaseshifts:latest -v /path/to/host/input/data:/data [<phsh-args> ...]
```

<!--lint disable no-unused-definitions-->
> [!TIP]
> Development docker images can be built locally, e.g.
> `DOCKER_TAG=dev make docker`
<!--lint enable no-unused-definitions-->

<!--lint disable no-unused-definitions-->
> [!WARNING]
> There is a [known possible
> bug](https://github.com/Liam-Deacon/phaseshifts/issues/6) where the
> compiled `libphsh.f` is not thread-safe (as ascertained by the fortran
> compiler), as such if you anticipate using this library in concurrent
> environments then it is advised to run `phsh.py` via
> `docker run ghcr.io/Liam-Deacon/phaseshifts:latest` as this works
> around this limitation due to the emphereal nature of container
> instances created using `docker run`.
<!--lint enable no-unused-definitions-->

## Install

### TDLR;

For python 3.11 or older:

```bash
#  install latest release
pip install phaseshifts

# development install
pip install wheel numpy setuptools  # needed for older python/pip versions
pip install -e .
phsh --help
```

### Details

The [phaseshifts](http://https://pypi.python.org/pypi/phaseshifts/)
package requires CPython 2.7 or later and also uses the
[numpy](http://www.scipy.org/scipylib/download.html),
[scipy](http://www.scipy.org/scipylib/download.html) and
[periodictable](http://https://pypi.python.org/pypi/periodictable)
packages. Currently, it has only been tested most extensively with Python 2.7
on Windows, so there are no guarantees with other platforms. To perform
a setup follow the steps below.

> 1.  Install the numpy, scipy and periodictable packages.
>
>     On systems compatible with PyPI this can be done using the
>     command:
>
>     ```bash
>     pip install numpy scipy periodictable
>     ```
>
>     Or if you have the easy_install package:
>
>     ```bash
>     easy_install install numpy scipy periodictable
>     ```
>
>     Older versions of numpy & scipy did not allow simultaneous
>     installation -if you experience problems then try first installing
>     numpy before attempting to install scipy. The periodictable
>     package allows lookup of the most common crystal structure for a
>     given element and is instrumental in many of the convenience
>     functions contained in the model module.
>
>     Alternatively download and install these packages manually
>     following the instructions provided for the respective packages.
>
> 2.  To install the phaseshifts package:
>
>      ```bash
>      python setup.py install
>      ```
>
>     With any luck the package has been installed successfully. A set
>     of test scripts are provided, however a simple check may suffice
>     using an interactive session of the python interpreter:
>
>     ```python
>     >>> import phaseshifts
>     >>> from phaseshifts.lib import libphsh # compiled FORTRAN .pyd or .so using f2py
>     ```
>
>     If these execute without errors then it is likely that all is
>     well, but in case of problems or bugs please use the contact
>     provided below and I will do my best to address the problem
>     quickly.

<!--lint disable no-unused-definitions-->
> [!TIP]
> On Windows systems it may be easier to install a scientific python
> distibution rather than install the dependencies from source -
> [Python(x,y)](http://code.google.com/p/pythonxy) or
> [Anaconda](https://www.anaconda.com/download) with mingw (gcc &
> gfortran) installed is highly recommended. Mac OS X users can simply
> do `brew install gfortran` and Debian/Ubuntu users can do
> `sudo apt-get install -y gfortran`.
<!--lint enable no-unused-definitions-->

<!--lint disable no-unused-definitions-->
> [!NOTE]
> On unix systems, setup the virtualenv on Python 3.10 or lower,
> activate it and run <span class="title-ref">make</span>.
<!--lint enable no-unused-definitions-->

<!--lint disable no-unused-definitions-->
> [!WARNING]
> Python 3.12 compatibility is a work in progress due to the removal of
> `numpy.distuils` build backend for `f2py` preventing simple
> installation via `pip install`, [this github
> issue](https://github.com/Liam-Deacon/phaseshifts/issues/8) tracks
> progress on fixing this known issue.
<!--lint enable no-unused-definitions-->

## About the code

The example source codes provided in this package are intended to be
instructional in calculating phase shifts. While it is not recommended
to use the example code in production, the code should be sufficient to
explain the general use of the library.

If you aren't familiar with the phase shift calculation process, you can
read further information in `doc/` folder:

  - [`phshift2007`](https://phaseshifts.readthedocs.io/en/latest/phshift2007.html) - a brief user guide/documentation concerning the
    input files (& details of the original fortran
    <span class="title-ref">phshift</span> package).
  - [phaseshifts API](https://phaseshifts.readthedocs.io/en/latest/modules.html) - a more detailed overview of the library functions
    and how to calculate phase shifts using the convenience functions in
    this package. This is not yet finished and so the reader is referred
    to the above document for the time being.

For those wanting a crash course of the Van Hove / Tong programs, I
advise reading the [phsh2007](https://phaseshifts.readthedocs.io/en/latest/phshift2007.html) document. See the `examples/` directory
to get an idea of the structure of the input files (for a random
selection of models & elements). In particular see the `cluster_Ni.i`
file for helpful comments regarding each line of input.

Those of you who are eager to generate phase shifts - first look at the
example cluster files for a bulk and slab calculation, noting that the
atoms in the model are in fractional units of the *a* basis vector for
the unit cell (SPA units). Next, after creating a bulk and slab model in
the `cluster.i` format, simply use the following python code:

```python
>>> from phaseshifts.phsh import Wrapper as phsh
>>> phsh.autogen_from_inputs(bulk_file, slab_file)
```

This will hopefully produce the desired phase shift output files (at
least for simple models) and works by assessing the two models to
determine what output to produce. For more detailed documentation and
function use refer to the pdf manual.

<!--lint disable no-unused-definitions-->
> [!TIP]
> A standalone command line utility **phsh.py** is provided as a way of
> automating the generation of phase shifts as part of a script. For
> more information use:
>
> ``` bash
> phsh.py --help
> ```
<!--lint enable no-unused-definitions-->

<!--lint disable no-unused-definitions-->
> [!NOTE]
> The <span class="title-ref">phaseshifts.leed</span> module provides a
> conversion class for CLEED `.inp` and `.bul` files. This is included
> as part of the <span class="title-ref">phsh.py</span> module, however
> the file extension is important (needs `.inp`, `.pmin`, `.bul`, or
> `.bmin`) and error checking is limited. There are also plans to
> include a validator to check the files for malformatted input at some
> point in the future.
<!--lint enable no-unused-definitions-->

### Alternatives

A number of alternatives are available, notably the following:

  1. [AQuaLEED](https://physics.mff.cuni.cz/kfpp/povrchy/files/) (with
     a useful [poster overview of phaseshifts
     calculations](https://physics.mff.cuni.cz/kfpp/povrchy/files/1179-Poster.pdf)).
     This is an officially mentioned piece of software on Michel Van
     Hove's [LEED Calculation Homepage](https://www.icts.hkbu.edu.hk/VanHove_files/leed/leedpack.html). Although the
     poster mentions that the software is written in python, this
     software is not (currently) distributed on <https://PyPI.org> (or via alternative means such as a docker image on [DockerHub](https://www.docker.com/products/docker-hub/))  and
     therefore harder to integrate with other python LEED-related
     projects such as [CLEED](https://github.com/Liam-Deacon/CLEED) and
     [cleedpy](https://github.com/empa-scientific-it/cleedpy).
  2. Elastic Electron-Atom Scattering in Solids and Solid Surfaces
     [(EEASiSSS)](https://www.researchgate.net/profile/John-Rundgren-2/publication/235583683_Optimized_surface-slab_excited-state_muffin-tin_potential_and_surface_core_level_shifts/links/5a266f89a6fdcc8e866bd7e5/Optimized-surface-slab-excited-state-muffin-tin-potential-and-surface-core-level-shifts.pdf)
     is authored by John Rundgren and first described in the paper: *"J. Rundgren Phys. Rev. B 68 125405 (2003)"*.
     This program takes a different approach to calculating phase shifts by using optimised muffin-tin potentials
     for surface slabs with preassigned surface core-level shifts.
     Whilst the source code is not publicly available online (to this author's best knowledge), John Rundgren
     has been more than happy to assist when approached in the past.
  3. A fortran program is described in "[McGreevy, E., & Stewart, A.L. (-Apr
     1978).](https://inis.iaea.org/search/search.aspx?orig_q=RN:9399501)
     A program for calculating elastic scattering phase shifts for an
     electron colliding with a one-electron target using perturbation
     theory. Computer Physics Communications, 14(1-2), 99-107.", however
     this code is not publicly available online (pay-walled by journal).

<!--lint disable no-unused-definitions-->
> [!NOTE]
> It would be fantastic to include this software (and document it's use) as part of the phaseshifts python package
> allowing the user to choose the backend they wish to use for calculating phase shifts (e.g. ``EEASiSSS`` or ``phshift2007``).
> As such [John Rundgren](https://www.researchgate.net/profile/John-Rundgren-2) should be contacted to see if
> he would be happy to collaborate on making this possible. This is being tracked by
> [this item](https://github.com/Liam-Deacon/phaseshifts/issues/92).
<!--lint enable no-unused-definitions-->

<!--lint disable no-unused-definitions-->
> [!IMPORTANT]
> Should you know of alternatives, please either [open an
> issue](https://Liam-Deacon/phaseshifts/issues) or (better yet) create
> a PR with changes to this documentation to keep this list up to date.
<!--lint enable no-unused-definitions-->

## Acknowledgements

As with all scientific progress, we stand on the shoulders of giants. If
this package is of use to you in publishing papers then please
acknowledge the following people who have made this package a reality:

  - **A. Barbieri** and **M.A. Van Hove** - who developed most of the
    original fortran code. Use *A. Barbieri and M.A. Van Hove, private
    communication.* (see `doc/phsh2007.txt` for further details).
  - **E.L. Shirley** - who developed part of the fortran code during
    work towards his PhD thesis (refer to the thesis: *E.L. Shirley,
    "Quasiparticle calculations in atoms and many-body core-valence
    partitioning", University of Illinois, Urbana, 1991*).
  - **Christoph Gohlke** - who developed the elements.py module used
    extensively throughout for the modelling convenience functions (see
    'elements.py' for license details).

I would also be grateful if you acknowledge this python package
(*phaseshifts*) as: *L.M. Deacon, private communication.*


[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FLiam-Deacon%2Fphaseshifts.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FLiam-Deacon%2Fphaseshifts?ref=badge_large)

### Thanks

I wish to personally add a heart-felt thanks to both Eric Shirley and
Michel Van Hove who have kindly allowed the use of their code in the
[libphsh.f](https://github.com/Liam-Deacon/phaseshifts/blob/master/phaseshifts/lib/libphsh.f) file needed for the underlying low-level functions in this
package.

## Contact

This package is developed/maintained in my spare time so any bug
reports, patches, or other feedback are very welcome.

The project is (still) in the early developmental stages and so anyone who
wishes to get involved are most welcome.  Please either
[create an issue](https://github.com/Liam-Deacon/phaseshifts/issues/new) or (better yet) submit a [pull request](https://github.com/Liam-Deacon/phaseshifts/pulls).

<!--lint disable no-unused-definitions-->
> [!TIP]
> Please [star](https://github.com/Liam-Deacon/phaseshifts) it on GitHub as this will help
> to easily indicate that others find the package useful.
<!--lint enable no-unused-definitions-->

## To Do

> 1.  Documentation - the manual has been started, but is not complete
>     and thus is a high priority. The current aim is to use sphinx to
>     generate html and latex documents for semi-automated generation of
>     both the tutorial and supporting website. If you have the
>     phaseshifts source and the
>     [sphinx](https://pypi.python.org/pypi/Sphinx) and the
>     [numpydoc](https://pypi.python.org/pypi/numpydoc) PyPi packages
>     then you can try making html or latex manuals using `make html` or
>     `make latexpdf` commands from the `doc/` directory.
> 2.  Test suit to verify the package is working as expected.
> 3.  GUI frontend (Qt ui files are provided in the `gui/` directory for
>     anyone wishing to undertake this challenge). Other front ends are
>     welcome (I use Qt due to familiarity/experience). For those
>     wishing a sneak preview, try executing `main.pyw`

See either [todo issues](https://github.com/Liam-Deacon/phaseshifts/issues?q=is%3Aopen+is%3Aissue+label%3A%22todo+%3Aspiral_notepad%3A%22) or [TODO.rst](https://github.com/Liam-Deacon/phaseshifts/blob/master/TODO.rst) for more information.

## Contacts

  - [Liam Deacon](mailto://liam.m.deacon@gmail.com) - *current maintainer*
  - [Michel Van Hove](mailto://vanhove@cityu.edu.hk) - Contact for
  original LEEDPACK `phsh[0-3].f` programs
