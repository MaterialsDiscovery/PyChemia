__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2016"
__version__ = "0.1.2"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "May 13, 2016"


class Version:
    @staticmethod
    def full_version():
        return 'PyChemia Version=' + __version__ + ' from=' + __date__

    def __init__(self):
        pass
