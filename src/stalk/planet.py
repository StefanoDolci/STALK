"""
stalk.planet
------------
Planetary model for transiting exoplanets
"""

from dataclasses import dataclass
import numpy as np 

# Constants
#----------------
R_EARTH_M = 6.371e6
AU_M = 1.496e11


#------------------


@dataclass
class Planet: 
    radius : float = 1.0            # Radius of the planet in Earth radii
    period : float = 365.25         # Period od the orbit in days
    semi_maj_axis : float = 1.0     # Distance from the Star in AU
    inclination : float = 90.0      # Orbital inclination in degrees ° 
    name : str = "Planet"



    @property
    def radius_m(self):
        """
        Planet radius (meters)
        """
        
        return self.radius * R_EARTH_M
    

    @property
    def distance_m(self):
        """
        Planet distance from the Star (meters)

        """
        
        return self.semi_maj_axis * AU_M
    
    


    
    def transit_depth(self,star):
        """
        Geometric depth of the transit in ppm
        

        Params
        ------
        Star : class
            Star class (radius, teff, u1, u2, name)
        
            
        Returns
        -------
        Flux depth given by the ratio of the Planet_radius with the Star_radius
        ppm
        """

        return (self.radius_m /star.radius_m)**2 * 1e+6
    

