.. _scripts:

*******
Scripts
*******

.. _phsh:

phsh.py
=======

Command line usage
------------------

The *phsh.py* script is placed into the system ``PATH`` during installation of the
phaseshifts package. It can then be used from the command line, e.g. ``phsh.py --help``
will produce a list of command line options::

  usage: phsh.py [-h] -b <bulk_file> -i <slab_file> [-t <temp_dir>] [-l <lmax>]
               [-r <start_energy> <final_energy> <step>] [-f <format>]
               [-S <subdir>] [-v] [-V]


  phsh - quickly generate phase shifts

        Created by Liam Deacon on 2013-11-15.
        Copyright 2013-2015 Liam Deacon. All rights reserved.

        Licensed under the MIT license (see LICENSE file for details)

        Please send your feedback, including bugs notifications
        and fixes, to: liam.deacon@diamond.ac.uk

      usage:-

    optional arguments:
    -h, --help            show this help message and exit
    -b <bulk_file>, --bulk <bulk_file>
                          path to MTZ bulk or CLEED *.bul input file
    -i <slab_file>, --slab <slab_file>
                          path to MTZ slab or CLEED *.inp input file
    -t <temp_dir>, --tmpdir <temp_dir>
                          temporary directory for intermediate file generation
    -l <lmax>, --lmax <lmax>
                          Maximum angular momentum quantum number. [default: 10]
    -f <format>, --format <format>
                          Use specific phase shift format i.e. 'cleed', 'curve'
                          or 'none'. Choose 'curve' if you wish to produce
                          XYY... data for easy plotting. <format> is case
                          in-sensitive. [default: 'cleed']
    -r <energy> [<energy> ...], --range <energy> [<energy> ...]
                          Energy range in eV with the format:
                          '<start> <stop> [<step>]', where the <step> value is
                          optional.  Valid for relativistic calculations
                          only. [default: (20, 600, 5)]
    -g, --generate-only
                          Exit after generating phaseshifts; do not launch
                          subprocess using PHASESHIFTS_LEED environment
                          variable. [default: False]
    -a, --atorbs-only
                          Only generate atomic orbitals of elements found in
                          the input files using Eric Shirley's hartfock
                          routine, then exit. [default: False]
    -p, --package
                          Selects package to use for phase shift calculations.
                          Choices are 'VHT' (van Hove-Tong) or 'Rundgren'
                          (EEASiSSS). [default: 'VHT')

    -S <subdir>, --store <subdir>
                          Keep intermediate files in subdir when done
    -v, --verbose         Set verbosity level [default: None].
    -V, --version         Show program's version number and exit

CLEED compatibility
-------------------
It is possible to use this script to generate phase shift files iteratively
during a geometry search for the CLEED package. In this manner phase shifts
will be generated at the beginning of each cycle of the search.

For this to work, the environment variable :envvar:`CSEARCH_LEED` must point to the
:code:`phsh.py` script, which will invoke the LEED program in :envvar:`PHASESHIFT_LEED`
after execution. When operating in this mode, the following assumptions are made:

 1. `-b <bulk_file>` option is not needed and the filename is assumed by
    changing the file extension of `<slab_file>` to '.bul'
 2. `-f CLEED` format is implied.
 3. The generated phase shifts are stored in the directory set by the
    :envvar:`CLEED_PHASE` environment variable, however a named copy with the
    iteration number (read from the matching '.log' file) will be placed in the
    same directory as the <slab_file>.
 4. `<lmax>` is equal to 10, unless additional parameter syntax is given in the CLEED
    `\*.inp` file. To use phase shift specific lmax values, then add a new line with::

        lmax:  <phase_shift> <lmax>

    for each phase shift you wish to have a different lmax to that of the default.
 5. The element and oxidation of each atom in a model is guessed by reading the phase
    shift tag from the CLEED input file. For example::

        po:  O_-2_COOH ...

    will be interpreted as a Oxygen with a -2 oxidation state and with a unique name
    tag of "O_-2_COOH" to show it is in a carboxylic group. Note the '-' must
    be at the beginning the oxidation sub-string. If no oxidation state is
    given then the atom is assumed to have zero charge.
 6. The muffin-tin radius of the phase shift species is guessed from lines with::

        rm:  <phase_shift> <radius>

    However, if no value is found the radius is guessed from the
    ::code::`ELEMENTS` dictionary within :py:mod:`phaseshifts.elements`
    depending on the valency of the given phase shift element.

A full list of additional syntax to customise the generation of the phase shifts
when using CLEED input files can be found in
:py:meth:`phaseshifts.leed.Converter.import_CLEED`.

.. note::
  If the :envvar:`PHASESHIFT_LEED` environment variable is not found, but
  :envvar:`CLEED_PHASE` is, however, found then the program will place the generated
  files in this directory unless a specific :code:`-S <subdir>` is provided.
