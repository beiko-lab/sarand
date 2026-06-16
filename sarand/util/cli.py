"""argparse type helpers and startup dependency checks."""
import argparse
import sys
from pathlib import Path

from sarand.config import PROGRAM_VERSION_NA
from sarand.external.bandage import Bandage
from sarand.external.blastn import Blastn
from sarand.util.logger import LOG


def check_file(path: str) -> Path:
    """argparse type: check an input file exists and is readable."""
    path = Path(path)
    if path.exists() and path.is_file():
        return path.resolve()
    raise argparse.ArgumentTypeError(
        f"{path} can't be found, please double check the path"
    )


def validate_range(value_type, minimum, maximum):
    """argparse type factory: ensure the value parses and is within [min, max]."""

    def range_checker(arg):
        if value_type is float:
            try:
                val = float(arg)
            except ValueError:
                raise argparse.ArgumentTypeError("Must be a float")
        elif value_type is int:
            try:
                val = int(arg)
            except ValueError:
                raise argparse.ArgumentTypeError("Must be an int")

        if val < minimum or val > maximum:
            raise argparse.ArgumentTypeError(f"must be in range [{minimum}-{maximum}]")
        return val

    return range_checker


def assert_dependencies_exist(blastn=True, bandage=True):
    """Check the external tools sarand shells out to exist and report versions."""
    versions = list()
    missing = list()
    if blastn:
        blastn_v = Blastn.version()
        versions.append(f'Blastn v{blastn_v}')
        if blastn_v is PROGRAM_VERSION_NA:
            missing.append('Blastn')
    if bandage:
        ba_v = Bandage.version()
        versions.append(f'Bandage v{ba_v}')
        if ba_v is PROGRAM_VERSION_NA:
            missing.append('Bandage')

    if len(versions) > 0:
        LOG.info('All dependencies found: ' + ', '.join(versions))
    if len(missing) > 0:
        LOG.error(f'The following tools are missing: {", ".join(missing)}')
        sys.exit(1)
