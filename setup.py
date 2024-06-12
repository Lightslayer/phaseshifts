#!/usr/bin/env python
from __future__ import print_function  # noqa

import os
import sys
import sysconfig

sys.path.append(os.path.dirname(__file__))

try:
    import phaseshifts
    import phaseshifts.phshift2007
except ModuleNotFoundError:
    phaseshifts = None  # type: ignore [assignment]

try:
    import setuptools  # type: ignore [import-untyped]
    from setuptools import find_packages, setup, Extension  # type: ignore [import-untyped]
except ImportError:
    from distutils.core import find_packages  # type: ignore [attr-defined,no-redef]

INCLUDE_DIRS = []

# WARNING: numpy.distutils is completely removed in python 3.12 and deprecated for removal in python 3.11 by Oct 2025
# The project will therefore need to be migrated to use a different build backend, see
# https://numpy.org/doc/stable/reference/distutils_status_migration.html#distutils-status-migration
try:
    import skbuild
except ImportError:
    skbuild = None  # type: ignore [assignment]

BUILD_BACKEND = None
try:
    from numpy.distutils.core import Extension, setup  # type: ignore  # noqa

    BUILD_BACKEND = "numpy.distutils"
except ModuleNotFoundError as npy_err:
    if skbuild and tuple(sys.version_info[:2]) >= (3, 11) and "bdist_wheel" in sys.argv:
        # TODO: Need to use migrate to a new f2py build backend, but have issues with pyproject.toml PEP-517 install
        # FIXME: Currently skbuild with CMakeLists.txt does not work, check with `make libphsh.cmake`
        try:
            from skbuild import setup

            BUILD_BACKEND = "skbuild"
        except ImportError:
            raise NotImplementedError(
                "TODO: Generate binary wheels correctly using pyproject.toml, scikit-build and cmake"
            )
    if tuple(sys.version_info[:2]) >= (3, 11):
        # We can invoke f2py and compile manually
        # HACK: Workaround missing distutils by invoking f2py directly
        # FIXME: Generated wheels unaware of native extension & deemed universal; *.so must be included in MANIFEST.in
        import subprocess

        setup = setuptools.setup  # noqa: F811

        BUILD_BACKEND = "numpy.f2py"
        if not any(x in sys.argv for x in ("sdist", "wheel")):
            args = [
                sys.executable,
                "-m",
                "numpy.f2py",
                "libphsh.f",
                "-m",
                "libphsh",
                "-c",
                "-f77flags='-frecursive'",
            ]
            try:
                args += phaseshifts.phshift2007.COMPILER_FLAGS["gfortran"]
            except NameError:
                pass
            subprocess.check_call(args, cwd="./phaseshifts/lib")  # nosec
        try:
            import numpy.f2py

            INCLUDE_DIRS += [
                "-I{}".format(x)
                for x in (numpy.get_include(), numpy.f2py.get_include())
            ]
        except ImportError:
            print(
                "WARNING: Unable to import numpy.f2py module for build", file=sys.stderr
            )
    else:
        raise npy_err

if len(sys.argv) == 1:
    sys.argv.append("install")

CMAKE_ARGS = {}

if BUILD_BACKEND == "skbuild":
    CMAKE_ARGS = {
        "cmake_args": [
            '-DPYTHON_INCLUDE_DIR="{}"'.format(sysconfig.get_path("include")),
            '-DPYTHON_LIBRARY="{}"'.format(sysconfig.get_config_var("LIBDIR")),
        ]
    }

BUILD_EXT_INPLACE_ARGS = ["build_ext", "--inplace"]
if len(sys.argv) < 2 and BUILD_BACKEND == "numpy.distutils":
    # add build_ext and --inplace flags when performing: `python setup.py`
    sys.argv.extend(BUILD_EXT_INPLACE_ARGS)
elif BUILD_BACKEND == "numpy.f2py":
    # FIXME: The following is a crude workaround to build_ext to avoid compiled libphsh*.so with undefined symbols
    for arg in filter(lambda x: x in sys.argv, BUILD_EXT_INPLACE_ARGS):
        sys.argv.remove(arg)
    if "build" not in sys.argv:
        sys.argv.append("build")

# build f2py extensions
f2py_exts_sources = {
    "libphsh": [
        os.path.join(
            "phaseshifts",
            "lib",
            "libphsh" + (".f" if BUILD_BACKEND == "numpy.distutils" else "module.c"),
        ),
    ]
}

try:
    gfortran_compiler_args = phaseshifts.phshift2007.COMPILER_FLAGS["gfortran"]
except NameError:
    gfortran_compiler_args = []

f2py_platform_extra_args = {
    "darwin": {"extra_link_args": [], "extra_compile_args": []},
    "win32": {
        "extra_link_args": gfortran_compiler_args,
        "extra_compile_args": gfortran_compiler_args,
    },
    "linux": {
        "extra_link_args": gfortran_compiler_args + ["-lgomp"],
        "extra_compile_args": gfortran_compiler_args + ["-fopenmp"],
    },
    "linux2": {
        "extra_link_args": gfortran_compiler_args + ["-lgomp"],
        "extra_compile_args": gfortran_compiler_args + ["-fopenmp"],
    },
}[sys.platform]

f2py_exts = (
    [
        # NOTE: When hacking the build process for Python 3.12, we still want to force wheel to be platform specific
        #       See https://stackoverflow.com/a/53463910/22567758
        Extension(
            name="phaseshifts.lib._native_build",
            sources=[os.path.join("phaseshifts", "lib", "_native_build.c")],
        )
    ]
    if BUILD_BACKEND != "numpy.distutils"
    else [
        Extension(
            name="phaseshifts.lib.libphsh",
            extra_compile_args=f2py_platform_extra_args["extra_compile_args"],
            extra_link_args=f2py_platform_extra_args["extra_link_args"],
            sources=f2py_exts_sources["libphsh"],
        )
    ]
)

print("BUILD_BACKEND: {}".format(BUILD_BACKEND))

README = "README.md"

dist = setup(
    name="phaseshifts",
    packages=find_packages(),
    version=getattr(phaseshifts, "__version__", "0.1.8-dev"),
    author="Liam Deacon",
    author_email="liam.m.deacon@gmail.com",
    license="MIT License",
    url="https://github.com/Liam-Deacon/phaseshifts",
    description=(
        "Python-based version of the Barbieri/Van Hove phase shift calculation package for LEED/XPD modelling"
    ),
    long_description=(None if not os.path.exists(README) else open(README).read()),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: X11 Applications :: Qt",  # The end goal is is to have Qt or other GUI frontend
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    extras_require={
        "atorb": ["mendeleev", "elementy"],
        "gui": [
            "six",
            "qtpy",
            "pyqt5; python_version < '3.7'",
            "pyside6; python_version >= '3.7'",
            "pyqtgraph",
        ],
        "dev": [
            "black",
            "isort",
            "numpy",
            "pre-commit",
            "ruff",
            "scikit-build; python_version > '3.11'",
            "wheel",
        ],
        "test": ["pytest", "pytest-cov", "pytest-qt"],
        "doc": ["sphinx>=7,<8", "sphinx_rtd_theme", "numpydoc", "ipykernel"],
    },
    keywords="phaseshifts atomic scattering muffin-tin diffraction",
    include_package_data=True,
    package_data=(
        {}
        if BUILD_BACKEND == "skbuild"
        else {
            # If any package contains *.f or *.pyd files, include them:
            "": ["*.f", "*.pyd", "*.so", "*.dll"],
            # If any package contains *.txt or *.rst files, include them:
            ".": ["*.txt", "*.rst", "*.pyw", "ChangeLog"],
            "lib": ["lib/*.f", "lib/*.c", "lib/*.h"] + ["*.so", "*.pyd"],
            "gui": ["gui/*.ui", "gui/*.bat"],
            "gui/res": ["gui/res/*.*"],
        }
    ),
    scripts=["phaseshifts/phsh.py"],
    install_requires=[
        "scipy >= 0.7",
        "numpy >= 1.3",
        "periodictable",
        "typing_extensions",
    ],
    ext_modules=f2py_exts,
    console=[os.path.join("phaseshifts", "phsh.py")],
)
