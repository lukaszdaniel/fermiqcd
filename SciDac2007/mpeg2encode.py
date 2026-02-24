#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil

def setup_environment():
    """Sets up the environment variables."""
    os.environ["PATH"] = ":".join([
        os.environ.get("PATH", ""),
        "/bin", "/usr/bin", "/usr/sbin",
        "/usr/local/bin", "/usr/bsd", "/usr/ucb"
    ])

def get_program_paths():
    """Returns program directory and current working directory."""
    progname = os.path.basename(sys.argv[0])
    progdir = os.path.dirname(sys.argv[0])
    cwd = os.getcwd()

    # Normalize the path for program directory
    tmpdir = progdir if progdir.startswith("/") else os.path.join(cwd, progdir)
    tmpdir = os.path.normpath(tmpdir)
    
    return progname, progdir, tmpdir, cwd

def parse_arguments():
    """Parse command line arguments and return relevant data."""
    want_version = False
    ver = ""
    beta = False
    using_dev = False
    visitargs = []
    legacyvisitargs = sys.argv[1:]

    for arg in sys.argv[1:]:
        if arg == "-v":
            ver = sys.argv[sys.argv.index(arg) + 1]
            if not ver:
                sys.exit("The '-v' version argument requires a value.")
        elif arg == "-beta":
            beta = True
        elif arg == "-dv":
            using_dev = True
        elif arg == "-version":
            want_version = True
        elif arg in [
            "-xml2atts", "-xml2avt", "-xml2info", "-xml2makefile",
            "-xml2projectfile", "-xml2plugin", "-xml2python", "-xml2window",
            "-xml2java", "-xmltest", "-xmledit", "-makemili", "-convert",
            "-visitconvert", "-silex", "-surfcomp", "-text2polys",
            "-time_annotation"
        ]:
            print(f"NOTE: Specifying tools as an argument to VisIt is no longer necessary.\n"
                  f"In the future, you should just run '{arg[1:]}' instead.\n")
        elif arg in ["-mpeg2encode", "-mpeg_encode"]:
            progname = arg[1:]
        elif arg == "-composite":
            progname = "visit_composite"
        elif arg == "-transition":
            progname = "visit_transition"
        else:
            visitargs.append(arg)

    return want_version, ver, beta, using_dev, visitargs, legacyvisitargs

def determine_version(visitdir, ver, want_version, using_dev, progname):
    """Determine the correct version of VisIt."""
    if os.path.isdir(f"{visitdir}/exe"):
        if want_version:
            print(f"The version of VisIt in the directory {visitdir}/exe/ will be launched.")
            sys.exit(0)

        ver = ""
        if using_dev or progname == "visit":
            visitargs.append("-dv")
    else:
        exe = [os.path.join(visitdir, d, "bin", "internallauncher") for d in os.listdir(visitdir) if os.path.isdir(os.path.join(visitdir, d))]
        exeversions = [os.path.basename(os.path.dirname(os.path.dirname(e))) for e in exe]
        current_version = os.readlink(f"{visitdir}/current") if os.path.islink(f"{visitdir}/current") else None
        beta_version = os.readlink(f"{visitdir}/beta") if os.path.islink(f"{visitdir}/beta") else None

        if want_version:
            if not current_version:
                sys.exit("There is no current version of VisIt.")
            print(f"The current version of VisIt is {current_version}.")
            sys.exit(0)

        if not ver:
            if beta:
                if not beta_version:
                    sys.exit("There is no beta version of VisIt.")
                ver = beta_version
            else:
                if not current_version:
                    sys.exit("There is no current version of VisIt.")
                ver = current_version

        if ver not in exeversions:
            if not os.path.isdir(f"{visitdir}/{ver}"):
                sys.exit(f"There is no version '{ver}' of VisIt.")
            legacyvisitargs.append(f"-{progname}")
            legacycmd = [f"{visitdir}/bin/legacylauncher"] + legacyvisitargs
            os.execvp(legacycmd[0], legacycmd)

    return ver

def set_visit_environment(ver, progname, visitdir):
    """Set the environment variables for the internal visit launcher."""
    os.environ["VISITVERSION"] = ver
    os.environ["VISITPROGRAM"] = progname
    os.environ["VISITDIR"] = visitdir

def launch_visit(visitdir, visitargs):
    """Launch the internal VisIt launcher."""
    visitcmd = [f"{visitdir}/bin/internallauncher"] + visitargs
    os.execvp(visitcmd[0], visitcmd)

###############################################################################
#
# Purpose:
#   This is the primary front end launcher for programs in the VisIt toolchain.
#   It is separate from the pieces than change on a per-version basis.
#
#   More specifically:
#     Parse out version arguments and determine the version to run.
#     Fall back to legacy launcher if needed.
#     Use argv[0] (i.e. $0) to determine if which program is being run.
#     Determine if someone is trying to use the old "visit -prog" method
#         of launching tools, and if so both warn them and fix things.
#     Keep track of mixing of public and private versions of VisIt.
#     Finally, if all goes well, launch the "internallauncher" which is
#         allowed to contain version-specific pieces.
#
# Note: Place NO version specific code here!
#       Place NOTHING that requires backwards compatibility code here.
#       ...
#       The only exception in this file is that which is just enough to
#       fall back to the pre-version-specific visit launcher script, and
#       by the time you read this note, that should not require changes.
#
# Programmer:  Jeremy Meredith
# Date      :  December  8, 2004
#
# Modifications:
#    Jeremy Meredith, Mon Dec 13 13:38:05 PST 2004
#    It was erroneously grepping for partial version name matches.
#    E.g. it would think 1.4 had an internallauncher if it found
#    one for version 1.4.1.  This caused old versions to fail to run.
#
#    Jeremy Meredith, Fri Jan  7 13:48:45 PST 2005
#    My last fix wasn't completely portable across perl versions.  I fixed it.
#
#    Jeremy Meredith, Tue Apr 19 08:36:36 PDT 2005
#    Made sure both "-convert" and "-visitconvert" work.
#
#    Brad Whitlock, Tue Sep 19 17:48:12 PST 2006
#    Added -mpeg2encode, -composite, -transition.
#
###############################################################################
def main():
    # Setup environment
    setup_environment()

    # Get program paths and current working directory
    progname, progdir, tmpdir, cwd = get_program_paths()

    # Add program directory to PATH
    os.environ["PATH"] = ":".join([progdir, os.environ["PATH"]])

    # Parse command line arguments
    want_version, ver, beta, using_dev, visitargs, legacyvisitargs = parse_arguments()

    # Determine the correct VisIt version
    visitdir = os.path.dirname(tmpdir)
    ver = determine_version(visitdir, ver, want_version, using_dev, progname)

    # Set environment variables for internal VisIt launcher
    set_visit_environment(ver, progname, visitdir)

    # Launch VisIt
    launch_visit(visitdir, visitargs)

if __name__ == "__main__":
    main()
