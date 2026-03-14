"""
Centralized logging configuration for the scheduler package.

Configures the 'scheduler' root logger and provides get_logger() for
module-level loggers. Log level can be set via SCHEDULER_LOG_LEVEL
(DEBUG, INFO, WARNING, ERROR, CRITICAL); default is INFO.
"""

import logging
import os

# Package root logger name
PACKAGE_LOGGER_NAME = "scheduler"

# Default format: timestamp, level, logger name, message
_DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

_configured = False


def get_logger(name: str) -> logging.Logger:
    """
    Return a logger for the given name (typically __name__).

    The logger is a child of the 'scheduler' package logger. On first
    call, the package logger is configured with a stream handler and
    level from SCHEDULER_LOG_LEVEL (default INFO).

    Parameters
    ----------
    name : str
        Logger name; use __name__ in modules.

    Returns
    -------
    logging.Logger
        Logger instance.
    """
    global _configured
    if not _configured:
        _configure()
        _configured = True
    if name.startswith(PACKAGE_LOGGER_NAME + ".") or name == PACKAGE_LOGGER_NAME:
        return logging.getLogger(name)
    return logging.getLogger(PACKAGE_LOGGER_NAME + "." + name)


def _configure() -> None:
    """Configure the scheduler package root logger."""
    root = logging.getLogger(PACKAGE_LOGGER_NAME)
    level_name = os.environ.get("SCHEDULER_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    root.setLevel(level)
    if not root.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(level)
        handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT, datefmt=_DATE_FORMAT))
        root.addHandler(handler)
