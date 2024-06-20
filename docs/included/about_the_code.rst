About the code
==============

The example source code provided in this package are intended to be
instructional in calculating phase shifts. While it is not recommended to
use the example code in production, the code
should be sufficient to explain the general use of the library.

If you aren't familiar with the phase shift calculation process, you can
read further information in ``doc/`` folder:

+ :ref:`phshift2007.rst <Van_Hove_Phase_Shift_Package_Guide_Overview>` -
  a brief user guide/documentation concerning the input files
  (& details of the original fortran `phshift` package).
+ :ref:`phaseshifts API reference <api>` - a more detailed overview of the library
  functions and how to calculate phase shifts using the convenience functions in this package.
  This is not yet finished and so the reader is referred to the above document for the time being.

For those wanting a crash course of the Van Hove / Tong programs, I advise reading the
``phsh2007.txt`` document.
See the ``examples/`` directory to get an idea of the structure of the input files
(for a random selection of models & elements). In particular see the ``cluster_Ni.i``
file for helpful comments regarding each line of input.

.. tip:: There is also a nice diagram overview (ignoring the LEED parts) contained within the
         `AQuaLEED poster <https://physics.mff.cuni.cz/kfpp/povrchy/files/1179-Poster.pdf>`_
         available online.

Those of you who are eager to generate phase shifts - first look at the example
cluster files for a bulk and slab calculation, noting that the atoms in the model
are in fractional units of the *a* basis vector for the unit cell (SPA units). Next,
after creating a bulk and slab model in the ``cluster.i`` format, simply use
the following python code:

   >>> from phaseshifts.phsh import Wrapper as phsh
   >>> phsh.autogen_from_inputs(bulk_file, slab_file)

This will hopefully produce the desired phase shift output files (at least for
simple models) and works by assessing the two models to determine what output to
produce. For more detailed documentation and function use refer to the pdf manual.

.. tip:: A standalone command line utility **phsh.py** is provided as a way of
         automating the generation of phase shifts as part of a script. For more
         information use:

         .. code:: bash

            phsh.py --help

.. note:: The `phaseshifts.leed` module provides a conversion class for CLEED ``.inp`` and
          ``.bul`` files. This is included as part of the `phsh.py` module,
          however the file extension is important (needs ``.inp``, ``.pmin``, ``.bul``,
          or ``.bmin``) and error checking is limited. There are also plans to include a
          validator to check the files for mal-formatted input at some point in the
          future.
