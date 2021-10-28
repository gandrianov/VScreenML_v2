"""Modification of original VScreenML module (https://www.pnas.org/content/117/31/18477)"""

# Add imports here
from .vscreenml_v2 import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
