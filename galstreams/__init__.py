from ._version import version as __version__  # noqa


def _check_dependencies():
    try:
        import gala
    except ImportError as err:
        raise ImportError(
            "\n"
            "galstreams requires the 'gala' package.\n\n"
            "Installation instructions:\n\n"
            "  pip install gala\n"
            "or (recommended):\n"
            "  conda install -c conda-forge gala\n"
        ) from err


_check_dependencies()

from .core import *  # noqa
