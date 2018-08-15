# Cube Class
# Spaxel Class

##################################################################################
class Spaxel():


    __slots__ = ['flux', 'error', 'flux_weight', 'iflux']

    def __init__(self):
        self.flux = 0.0
        self.flux_weight = 0.0
        self.iflux = 0
        self.error = 0

class SpaxelAB():

    __slots__ = ['flux', 'error', 'flux_weight', 'iflux']

    def __init__(self):

        self.flux = 0
        self.error = 0
        self.flux_weight = 0.0
        self.iflux = 0.0
