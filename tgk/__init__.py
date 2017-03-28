# Licensed under a MIT style license - see LICENSE
"""tgk --- 41P/Tuttle-Giacobini-Kresak outburst detection with LCO."""

from .sync import Sync
from .science import Science

def show_config(config_file):
    """Show the TGKSync/Science configuration file.
    
    Parameters
    ----------
    config_file : string
      The configuration file to show or `None` to show the defaults.

    """
    from .core import TGKMaster
    TGKMaster.show_config(config_file)

