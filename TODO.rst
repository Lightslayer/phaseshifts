====
TODO
====

This is a highly experimental pre-alpha release and as such changes are being
made frequently to the source code. There are a number of features that will
be included as this package matures, which include:

 1. Full implementation of python convenience functions and classes to easily
    generate phase shift output files. This really boils down to the following:

    - Provide function to automatically generate phase shifts from a single
      function by providing the slab and bulk cluster input files. Python
      will then handle all the tedius atorb, muffin-tin potential (MTZ) and
      phase shift calculations (including removing pi jumps). This will be in
      ``phsh.py`` for ease of access.

    - Classes and modules to support the above.

    - Additional functionality not originally envisaged can be added to the
      ``plugins/`` directory.

 2. Documentation - the manual has not yet been started and so is a high priority
    after #1. Let the docstrings be a guide, but this should also include greater
    detail concerning models and the calculations themselves (i.e. background info
    on radial charge distributions, muffin-tin potentials and phase shifts). The aim
    is to also include a GUI section to help guide the casual scientific reader.

    .. important:: *STARTED* - Liam Deacon (will accept others to take this role).

 3. Test suite to determine whether the functions (particularly the low-level fortran
    functions in the ``libphsh.pyd`` or ``libphsh.so`` modules work as expected. There may be
    scenarios where users copy dynamic libraries that are incompatible with their
    system (for instance if they do not compile ``libphsh.f`` themselves) and therefore
    ``libphsh`` will need verification. Other potential problems may arise from the use of
    Jython, which (currently) does not support `scipy` or `numpy`. The software has been
    written in Python 2.7.5 and so checks should also be made to ensure that the CPython
    version is >=2.6 (Python 3 may be supported later using *e.g.* 2to3).

    .. important:: See https://github.com/Liam-Deacon/phaseshifts/issues/16 for current status.

 4. GUI frontend (Qt ui files are provided in the ``gui/`` directory for anyone
    wishing to undertake this challenge). Other frontends are welcome (I use Qt due
    to familiarity). [NOT YET ASSIGNED - any takers?].

    .. important:: See https://github.com/Liam-Deacon/phaseshifts/issues/29 for current status

 5. REST API to allow some sort of web stack usage, see https://github.com/Liam-Deacon/phaseshifts/issues/9

 6. Miscellanea:

    Address limitations of libphsh.f. In particular:
      + Allow dynamic allocation of key arrays for input files with any number of atoms
      + Address truncation of phase shift labels over 9 characters in length
      + Improve readability (it would be nice to eliminate *goto* statements)
      + Possibly form into a FORTRAN 90 module to allow optional arguments to
        functions and subroutines.
      + Pure-Python implementation of functionality, which would be particularly useful for users
        with restricted environments and/or to avoid specific compiler quirks
        (needs rigorous regression testing to ensure faithfulness of port)

    Ensure mutual compatibility between hartfock and model input files from
    EEASiSSS and Barbieri/Van Hove backends.
