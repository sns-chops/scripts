import numpy as np

class FermiChopper:

    radius = None # mm
    slit_thickness = None #mm
    frequency = None # Hz
    
    def delta_t(self):
        """compute time window size
        """
        period = 1./self.frequency
        return self.slit_thickness/2/np.pi/self.radius * period * 1e6
