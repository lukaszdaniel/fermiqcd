#!/bin/sh
# Wrapper to run the embedded Perl script from any path.

# Get the number of lines in this script
lines=$(wc -l < "$0")
# Subtract lines above the Perl script; adjust if Perl part moves
lines=$((lines - 22))

# Set TMPDIR if not already set
TMPDIR="${TMPDIR:-$HOME}"

# Extract the embedded Perl script
outfile="$TMPDIR/visitfelperl$$"
tail -n "$lines" "$0" > "$outfile" 2>/dev/null || tail --lines="$lines" "$0" >> "$outfile"

# Mark end of script and pass original command
echo "__END__" >> "$outfile"
echo "$0 $*" >> "$outfile"

# Execute the Perl script
exec perl "$outfile" "$0" "$@"

#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd 'cwd';

# Remove this script from disk
unlink $0;

# Update $0 to the actual program being run
$0 = shift @ARGV;

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

# -----------------------------------------------------------------------------
#                                 Setup
# -----------------------------------------------------------------------------

# Set a secure and reasonable path
$ENV{PATH} = join(':', $ENV{PATH}, qw(/bin /usr/bin /usr/sbin /usr/local/bin /usr/bsd /usr/ucb));

# Get script name and path
my $progname = basename($0);
my $progdir  = dirname($0);
my $cwd      = cwd();

# Add the program directory to the PATH
$ENV{PATH} = join(':', $progdir, $ENV{PATH});

# Normalize progdir into absolute path
my $tmpdir;
if ($progdir =~ m|^/|) {
    $tmpdir = $progdir;
}
elsif ($progdir eq '.') {
    $tmpdir = $cwd;
}
else {
    $tmpdir = "$cwd/$progdir";
}

# Normalize path by removing /./ and compressing /dir/..
$tmpdir =~ s|/\./|/|g;
$tmpdir =~ s|(/[^/]+/\.\.)||g;

my $visitdir = dirname($tmpdir);

# -----------------------------------------------------------------------------
#                            Parse the arguments
# -----------------------------------------------------------------------------

my $want_version = 0;
my $ver          = "";
my $beta         = 0;
my $using_dev    = 0;
my $ver_set      = 0;

my @legacyvisitargs = @ARGV;
my @visitargs;

while (@ARGV) {
    my $arg = shift @ARGV;

    if ($arg eq "-v") {
        $ver = shift @ARGV // die "The '-v' version argument requires a value.\n";
        $ver_set = 1;
    }
    elsif ($arg eq "-beta") {
        $beta = 1;
    }
    elsif ($arg eq "-dv") {
        $using_dev = 1;
    }
    elsif ($arg eq "-version") {
        $want_version = 1;
    }
    elsif ($arg eq "-xml2atts"        or
           $arg eq "-xml2avt"         or
           $arg eq "-xml2info"        or
           $arg eq "-xml2makefile"    or
           $arg eq "-xml2projectfile" or
           $arg eq "-xml2plugin"      or
           $arg eq "-xml2python"      or
           $arg eq "-xml2window"      or
           $arg eq "-xml2java"        or
           $arg eq "-xmltest"         or
           $arg eq "-xmledit"         or
           $arg eq "-makemili"        or
           $arg eq "-convert"         or
           $arg eq "-visitconvert"    or
           $arg eq "-silex"           or
           $arg eq "-surfcomp"        or
           $arg eq "-text2polys"      or
           $arg eq "-time_annotation") {
        $progname = substr($arg, 1);
        warn "NOTE: Use of 'visit $arg' is deprecated. Run '$progname' directly.\n";
    }
    elsif ($arg =~ /^-(mpeg2encode|mpeg_encode)$/) {
        $progname = substr($arg, 1);
    }
    elsif ($arg eq "-composite") {
        $progname = "visit_composite";
    }
    elsif ($arg eq "-transition") {
        $progname = "visit_transition";
    }
    else {
        push @visitargs, $arg;
    }
}

# -----------------------------------------------------------------------------
#                          Find the right version
# -----------------------------------------------------------------------------
# If we have a top-level "exe" directory, then don't bother looking
# for versions to use; this is a development executable.
if (-d "$visitdir/exe") {
    if ($want_version) {
        print "The version of VisIt in the directory $visitdir/exe/ will be launched.\n";
        exit 0;
    }

    # They may have specified the version, but we need to ignore it
    # because development trees don't have version directories.
    $ver = "";

    # We want to make sure we know if we are trying to launch a public
    # version from under a development version.  Keep track of this.
    # (Note -- don't add it on our own if we're launching a tool.)
    push @visitargs, "-dv" if $using_dev || $progname eq "visit";
}
else {
    # look for the version-specific visit script to determine viable versions
    my @exe = glob("$visitdir/*/bin/internallauncher");
    my @exeversions = map { m|^$visitdir/(.*?)/|; $1 } @exe;

    my $current_version = readlink("$visitdir/current");
    my $beta_version    = readlink("$visitdir/beta");

    if ($want_version) {
        if (!$current_version) {
            die "There is no current version of VisIt.\n";
        }
        print "The current version of VisIt is $current_version.\n";
        exit 0;
    }

    if (!$ver_set) {
        if ($beta) {
            die "There is no beta version of VisIt.\n" unless $beta_version;
            $ver = $beta_version;
        } else {
            die "There is no current version of VisIt.\n" unless $current_version;
            $ver = $current_version;
        }
    }

    # If there was no internal laucher for that version, then either
    # that version wasn't installed, or it is an old version that
    # didn't have a version-specific launcher script.
    if (!grep { $_ eq $ver } @exeversions) {
        unless (-d "$visitdir/$ver") {
            die "There is no version '$ver' of VisIt.\n";
        }

        # Fall back to legacy here; we set a version and didn't have a
        # version-specific launcher script.

        # Legacy tools were launched using the visit script with
        # an argumnet specifying which tool to run.  Fake it here.
        if ($progname ne "visit" && !grep { /^-\Q$progname\E$/ } @legacyvisitargs) {
            push @legacyvisitargs, "-$progname";
        }

        my @legacycmd = ("${visitdir}/bin/legacylauncher", @legacyvisitargs);
        exec @legacycmd or die "Can't execute legacy launcher: $!\n";

        # If we ever stop supporting versions before 1.4.1, we can
        # change to the "right" behavior here.
        #print STDERR "Version $ver of does not exist.\n";
        #exit 1;
    }

    # Warn if we mixed public and private development versions.
    if ($using_dev) {
        warn "\nWARNING: You are launching a public version of VisIt from within a development version!\n\n";
        push @visitargs, "-dv";
    }

    # The actual visit directory is now version-specific
    $visitdir = "$visitdir/$ver";
}

# -----------------------------------------------------------------------------
#     Set the environment variables needed for the internal visit launcher
# -----------------------------------------------------------------------------

$ENV{VISITVERSION} = $ver;
$ENV{VISITPROGRAM} = $progname;
$ENV{VISITDIR}     = $visitdir;

# -----------------------------------------------------------------------------
#                       Run the internal launcher!
# -----------------------------------------------------------------------------

my @visitcmd = ("${visitdir}/bin/internallauncher", @visitargs);
exec @visitcmd or die "Can't execute internal launcher: $!\n";
