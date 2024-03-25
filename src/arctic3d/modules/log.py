"""Manage logging."""
import logging
import sys
from functools import partial
from logging import FileHandler, StreamHandler


log_file_name = "log"


info_formatter = "[%(asctime)s %(module)s %(levelname)s] %(message)s"
debug_formatter = (
    "[%(asctime)s] "
    "%(filename)s:%(name)s:%(funcName)s:%(lineno)d: "
    "%(message)s"
)

log_formatters = {
    "DEBUG": debug_formatter,
    "INFO": info_formatter,
    "WARNING": info_formatter,
    "ERROR": info_formatter,
    "CRITICAL": info_formatter,
}

log_levels = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


def add_handler(
    log_obj,
    handler,
    stream,
    log_level="INFO",
    formatter=info_formatter,
):
    """Add a logging Handler to the log object."""
    ch = handler(stream)
    ch.setLevel(log_levels[log_level.upper()])
    ch.setFormatter(logging.Formatter(formatter))
    log_obj.addHandler(ch)
    return ch


def add_log_for_CLI(log, log_level, logfile):
    """Configure log for command-line clients."""
    llu = log_level.upper()

    params = {
        "log_level": llu,
        "formatter": log_formatters[llu],
    }

    log.handlers.clear()
    add_sysout_handler(log, **params)
    add_logfile_handler(log, stream=logfile, **params)
    return


add_sysout_handler = partial(
    add_handler, handler=StreamHandler, stream=sys.stdout
)  # noqa: E501
add_logfile_handler = partial(
    add_handler, handler=FileHandler, stream=log_file_name
)  # noqa: E501
