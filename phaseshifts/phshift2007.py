"""A collection of utilities to assist with extracting phshift2007 in a cross-platform manner."""

import argparse
import os.path
import sys
import subprocess  # nosec
import urllib.request
import zipfile

PHSHIFT2007_DOWNLOAD_URL = (
    "https://www.icts.hkbu.edu.hk/VanHove_files/leed/phshift2007.zip"
)
DEFAULT_DIRPATH = "."

COMPILER_FLAGS = {
    "gfortran": [
        "-pie",
        "-Wall",
        "-static-libgcc",
        "-static-libgfortran",
        "-frecursive",
        "-fcheck=bounds",
        "-std=legacy",
    ]
}


def download_phshift2007(dirpath=DEFAULT_DIRPATH):
    # type: (str) -> str
    """Download phshift2007 from the web."""
    os.makedirs(dirpath, exist_ok=True)
    zipfile_path = os.path.join(dirpath, "phshift2007.zip")
    urllib.request.urlretrieve(PHSHIFT2007_DOWNLOAD_URL, zipfile_path)  # nosec: B310
    return zipfile_path


def unzip_phshift2007(zipfile_path, dirpath=DEFAULT_DIRPATH):
    # type: (str, str) -> None
    """Unzip phshift2007 contents from the zipfile."""
    with zipfile.ZipFile(zipfile_path, mode="r") as zip_ref:
        zip_ref.extractall(dirpath)


def _extract_fortran_programs_from_ab_file(ab_filepath, output_dirpath=None):
    """Extracts the Fortran programs from the .ab files in the phshift2007.zip archive."""
    with open(ab_filepath, encoding="ascii") as fp_in:
        lines = fp_in.readlines()[3:]  # skip first three comment lines
        prog_marker = "C  program "
        program_start_lines = [
            i for i, line in enumerate(lines) if line.startswith(prog_marker)
        ]
        programs = {}
        for i, marker in enumerate(program_start_lines):
            prog_name = (
                lines[marker].lstrip(prog_marker).lower().strip("\n").strip("\r")
            )
            filename = os.path.join(
                output_dirpath or os.path.dirname(ab_filepath), prog_name
            )
            with open(filename, mode="w", encoding="ascii") as fp_out:
                if i == len(program_start_lines) - 1:
                    fp_out.writelines(lines[marker - 1 :])
                else:
                    fp_out.writelines(
                        lines[marker - 1 : program_start_lines[i + 1] - 1]
                    )
            programs[prog_name] = filename
        return programs


def extract_phshift2007(dirpath=DEFAULT_DIRPATH):
    """Extracts the phshift2007 Fortran programs from the phshift2007.zip archive."""
    for filename in ("psprog.ab3", "psprog.ab4"):
        filepath = os.path.join(dirpath, filename)
        print(_extract_fortran_programs_from_ab_file(filepath, dirpath))


def _compile_fortran_program(source_filepath, output_filepath=None):
    # type: (str, str) -> str
    """Compiles the Fortran program from the .for file in the phshift2007.zip archive."""
    if not output_filepath:
        output_filepath = os.path.join(
            os.path.dirname(source_filepath),
            os.path.basename(source_filepath).replace(".for", ""),
        )
        if sys.platform == "win32":
            output_filepath += ".exe"
    args = (
        ["gfortran"]
        + COMPILER_FLAGS["gfortran"]
        + [
            "-o",
            output_filepath,
            source_filepath,
        ]
    )
    print(" ".join(args))
    with subprocess.Popen(args) as proc:  # nosec
        status = proc.wait()
        if status != 0:
            raise subprocess.CalledProcessError(status, " ".join(proc.args))
    return output_filepath


def compile_phshift2007(dirpath=DEFAULT_DIRPATH):
    """Compiles the phshift2007 Fortran programs."""
    for program in ("phsh0", "phsh1", "phsh2cav", "phsh2wil", "phsh2rel", "phsh3"):
        source_filepath = os.path.join(dirpath, program + ".for")
        print("Compiled %s" % _compile_fortran_program(source_filepath))


def do_action(action, output_dirpath=DEFAULT_DIRPATH):
    """Perform the specified action."""
    if action in ("download", "all"):
        zipfile_path = download_phshift2007(output_dirpath)
        print("Downloaded phshift2007 to {}".format(zipfile_path))

    if action in ("unzip", "all"):
        unzip_phshift2007(zipfile_path, output_dirpath)
        print("Unzipped phshift2007 to {}".format(output_dirpath))

    if action in ("extract", "all"):
        extract_phshift2007(output_dirpath)
        print("Extracted phshift2007 Fortran programs to {}".format(output_dirpath))

    if action in ("compile", "all"):
        compile_phshift2007(output_dirpath)
        print("Compiled phshift2007 Fortran programs to {}".format(output_dirpath))


def get_parser():
    # type: () -> argparse.ArgumentParser
    """Get the parser object for the script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dirpath",
        type=str,
        default=DEFAULT_DIRPATH,
        help="The directory path to download and extract phshift2007 to.",
    )
    parser.add_argument(
        "--action",
        type=str,
        default="all",
        choices=["download", "unzip", "extract", "compile", "all"],
        help="The action to perform [default %(default)s].",
    )
    return parser


def main(argv=None):
    """Main entry point for the script."""
    parser = get_parser()
    args, _ = parser.parse_known_args(argv or sys.argv)
    do_action(action=args.action, output_dirpath=args.output_dirpath)


if __name__ == "__main__":
    main()
