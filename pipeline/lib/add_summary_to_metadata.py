#!/usr/bin/env python3

"""
Add summary data to metadata.

Usage: add-summary-to-metadata --config FILE --out FILE --gitdir PATH

Options:
    --config FILE   YAML format config file
    --out FILE      Output file
    --gitdir PATH   Path to pipeline .git directory
"""

import logging
import os
import re
import subprocess
import sys
import traceback
from collections import OrderedDict

import git
import yaml
from docopt import DocoptExit
from docopt import docopt
from functions import read_similarity_settings
from functions import reads_by_prefix
from tolkein import tofile

logger_config = {
    "level": logging.INFO,
    "format": "%(asctime)s [%(levelname)s] line %(lineno)d %(message)s",
    "filemode": "w",
}
try:
    logger_config.update({"filename": snakemake.log[0]})
except NameError as err:
    pass
logging.basicConfig(**logger_config)
logger = logging.getLogger()

logger.info("Starting script: " + __file__)


def add_git_meta(meta, args):
    """Add git info to metadata."""
    p = subprocess.Popen(
        ["git", "--git-dir", args["--gitdir"], "rev-parse", "--short", "HEAD"],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )
    (output, err) = p.communicate()
    p_status = p.wait()
    meta["settings"]["commit"] = output.rstrip("\n")

    p = subprocess.Popen(
        ["git", "--git-dir", args["--gitdir"], "remote", "-v"],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )
    (output, err) = p.communicate()
    p_status = p.wait()
    try:
        meta["settings"]["pipeline"] = output.split()[1]
    except IndexError:
        meta["settings"]["pipeline"] = "UNKNOWN"

    p = subprocess.Popen(
        ["git", "--git-dir", args["--gitdir"], "describe", "--tags", "--abbrev=0"],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )
    (output, err) = p.communicate()
    p_status = p.wait()
    try:
        meta["settings"]["release"] = output.split()[0]
    except IndexError:
        pass


def add_software_versions(meta, args):
    """Add software versions to meta."""
    programs = {
        "blastn": {"flag": "-version", "regex": r":\s*(\S+)$"},
        "blobtools": {"regex": r"\s+v(\S+)$"},
        "busco": {"regex": r"\s+(\S+)$"},
        "diamond": {"regex": r"version\s*(\S+)$"},
        "minimap2": {"regex": r"(\S+)"},
        "mosdepth": {"regex": r"\s+(\S+)$"},
        "python": {"regex": r"\s+(\S+)$"},
        "samtools": {"regex": r"\s+(\S+)$"},
        "seqtk": {"flag": "", "line": 2, "status": 1, "regex": r":\s*(\S+)$"},
        "snakemake": {"regex": r"(\S+)"},
    }
    versions = {}
    logger.info("Adding software version details")
    for key, opts in programs.items():
        cmd = [key]
        flag = opts.get("flag", "--version")
        if flag:
            cmd.append(flag)
        try:
            p = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8",
            )
            (output, err) = p.communicate()
            status = p.wait()
            expected_status = opts.get("status", 0)
            if expected_status == 1:
                output = err
                status = 0
            if status == 0:
                for idx, line in enumerate(output.split("\n")):
                    if idx == opts.get("line", 0):
                        match = re.search(opts["regex"], line)
                        version = match.groups()[0]
            else:
                version = "unknown"
        except FileNotFoundError:
            version = "not found"
        logger.info("%s: %s", key, version)
        versions.update({key: version})
    meta["settings"]["software_versions"] = versions


def parse_args():
    """Parse snakemake args if available."""
    try:
        sys.argv["--config"] = snakemake.config
        sys.argv["--out"] = str(snakemake.output)
        sys.argv["--gitdir"] = (
            git.Repo(
                os.path.dirname(os.path.abspath(__file__)),
                search_parent_directories=True,
            ).working_tree_dir
            + "/.git"
        )
    except NameError as err:
        pass


def main():
    """Entry point."""
    try:
        parse_args()
        args = docopt(__doc__)
    except DocoptExit:
        raise DocoptExit
    try:
        meta = {}
        config = args["--config"]
        if not isinstance(config, dict):
            config = tofile.load_yaml(config)
        meta["assembly"] = config["assembly"]
        meta["taxon"] = config["taxon"]
        meta["settings"] = {
            "blast_chunk": 100000,
            "blast_max_chunks": 10,
            **config["settings"],
        }
        meta["similarity"] = {
            "diamond_blastx": read_similarity_settings(config, "diamond_blastx"),
            "diamond_blastp": read_similarity_settings(config, "diamond_blastp"),
        }
        meta["reads"] = reads_by_prefix(config)
        for key, value in meta["reads"].items():
            if "url" in value and not isinstance(value["url"], list):
                value["url"] = value["url"].split(";")
            elif "url" not in value:
                value["url"] = []
        add_git_meta(meta, args)
        add_software_versions(meta, args)

        yaml.add_representer(
            OrderedDict,
            lambda dumper, data: dumper.represent_mapping(
                "tag:yaml.org,2002:map", data.items()
            ),
        )
        yaml.add_representer(
            tuple,
            lambda dumper, data: dumper.represent_sequence(
                "tag:yaml.org,2002:seq", data
            ),
        )
        yaml.Dumper.ignore_aliases = lambda *args: True
        with open(args["--out"], "w") as fh:
            fh.write(yaml.dump(meta, indent=1))

    except Exception as err:
        logger.error(traceback.format_exc())
        exit(13)


if __name__ == "__main__":
    main()
