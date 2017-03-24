import logging

class TGKSync:
    """Sync program data with LCO.

    Parameters
    ----------
    refresh : int
      Check for new data every `refresh` hours after last download.
      Must be >=2.

    """

    def __init__(self, refresh=2):
        from datetime import datetime
        from astropy.time import Time
        import astropy.units as u

        now = Time(datetime.now())

        self.refresh = max(refresh, 2)
        self.last_download = None

        
    def sync(self, **kwargs):
        import requests

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Master control program for monitoring LCO images of 41P.')
    args = parser.parse_args()
    sync = TGKsync()
    sync.run()

