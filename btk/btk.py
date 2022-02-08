#!/usr/bin/env python3

"""
BlobToolKit.

Usage: btk [<command>] [<subcommand>] [<args>...] [-h|--help] [--version]

Commands:
    blobtools       Blobtools commands
    pipeline        Pipeline commands
    viewer          Viewer commands
    -h, --help      Show this
    -v, --version   Show version number
See 'btk <command> --help' for more information on a specific command.

"""

import sys

from docopt import DocoptExit
from docopt import docopt
from pkg_resources import working_set

from .lib.version import __version__


def main():
    """Entry point."""
    if len(sys.argv) > 1:
        try:
            args = docopt(__doc__, help=False, version=__version__)
        except DocoptExit:
            args = {"<command>": sys.argv[1]}
        if args["<command>"]:
            # load <command> from entry_points
            for entry_point in working_set.iter_entry_points("btk.subcmd"):
                if entry_point.name == args["<command>"]:
                    subcommand = entry_point.load()
                    sys.exit(subcommand())
    print(__doc__)
    raise DocoptExit
