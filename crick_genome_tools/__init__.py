#!/usr/bin/env python
"""
Set app version
"""

from importlib.metadata import version


# Set version from package information
__version__ = version("crick-genome-tools")
