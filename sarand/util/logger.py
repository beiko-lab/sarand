import logging
import ntpath
import sys
from typing import Optional

from sarand import __version__, __program_name__
from sarand.config import SARAND_LOGGER_NAME


def get_logger() -> logging.Logger:
    """Returns the default logger for sarand."""
    return logging.getLogger(SARAND_LOGGER_NAME)


def create_logger(output: Optional[str] = None, verbose: bool = False):
    """Initialize the logger."""

    level = logging.DEBUG if verbose else logging.INFO

    # create logger
    logger = get_logger()
    if logger.hasHandlers():
        return logger
    logger.setLevel(level)

    # set the log format
    formatter = logging.Formatter('%(asctime)s [%(levelname)-5.5s]	%(message)s')

    # create console handler
    stream_h = logging.StreamHandler()
    stream_h.setLevel(level)
    stream_h.setFormatter(formatter)
    logger.addHandler(stream_h)

    # create the file handler
    if output is not None:
        file_h = logging.FileHandler(output)
        file_h.setLevel(level)
        file_h.setFormatter(formatter)
        logger.addHandler(file_h)

    logger.info(f'{__program_name__} v{__version__}')
    base_name = ntpath.basename(sys.argv[0])
    if base_name == '__main__.py':
        prog_name = __name__.split('.')[0]
    else:
        prog_name = base_name
    logger.info(f'{prog_name} {" ".join(sys.argv[1:])}')

    return logger


# Default class to export the logger for easy use
LOG = get_logger()
