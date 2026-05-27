"""
stalk.star
----------
Stellar model with quadratic limb darkening for transit simulations.
"""


import numpy as np
from dataclasses import dataclass

# Constants
#-------------
R_SUN_M = 6.957e+8




# Quadratic limb darkening coefficients by stectral type of the Star
# Source Claret & Bloemen (2011), Kepler bandpass

LDCs = {
    "O": (0.27, 0.22),
    "B": (0.34, 0.28),
    "A": (0.42, 0.33),
    "F": (0.48, 0.30),
    "G": (0.44, 0.23),
    "K": (0.55, 0.20),
    "M": (0.62, 0.14),
}


@dataclass 
class Star: 
    radius: float = 1.0  # Radius in Solar radii
    teff: float = 5778.0 # Effective temperature in Kelvin
    u1: float = 0.44     # First limb darkening coefficient
    u2: float = 0.23     # Second limb darkening coefficient
    name: str = "Star"   # Star name



    @property
    def radius_m(self):
        """
        Stellar radius (meters)
        """
        
        return self.radius * R_SUN_M



    def limb_darkening(self,mu):
        """
        Limb darkening model

        Params
        ------
        mu : float or array
            cos(theta), where theta is the angle between the line of sight
            and the normal to a given point of stellar surface.
            mu = 1 at centre, mu = 0 at limb.


        Returns 
        -------
        Normalized intesity profile of the star I(mu)/I(1)
        """


        mu = np.asarray(mu)
        return 1.0 - self.u1*(1.0 - mu) - self.u2*(1.0 - mu)**2
    



    def intensity_grid(self,n):
        """
        2D Intensity map of the stellar disk

        Params
        ------
        n : int 
            Pixel resolution (n x n)


        Returns
        -------
        xx, yy : np.ndarray
                Coordinates normalised to the star radius (-1 to +1)
        intensity : np.ndarray
                Intensity map

        """
        # Grid from -1 to +1 in units of stellar radii
        x = np.linspace(-1,1,n)
        y = np.linspace(-1,1,n)

        # 2 2D meshgrid of coordinates
        xx, yy = np.meshgrid(x,y)

        # Square distance form the center
        r2 = (xx**2 + yy**2)

        # Calculate mu for every pixel ( z = sqrt(1-r^2) )
        mu = np.sqrt(np.maximum(0.0, 1.0 - r2))

        # Inside disk we use limb_darkening, outside 0 
        intensity = np.where (r2 <= 1.0, self.limb_darkening(mu), 0.0)
        return xx, yy, intensity
    


    @classmethod
    def ldcs_from_spectral_type(cls, spectral_type, radius=1.0, teff=5778.0, name="Star"):
        """
        Creates a Star from spectral type and Limb Darkening Coefficients

        Params
        ------

        cls: class
            existing class

        spectral_type : str
                        Choose from "O", "B", "A", "F", "G", "K", "M"
        """
        sp = spectral_type.upper()[0]

        if sp not in LDCs:
            raise ValueError(f"Wrong spectral type: {spectral_type}")
        
        u1, u2 = LDCs[sp]
        return cls(radius = radius, teff = teff, u1 = u1, u2 = u2, name = name)