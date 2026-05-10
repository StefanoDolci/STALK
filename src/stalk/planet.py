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
    period : float = 365.25         # Period od the orbit [days]
    semi_maj_axis : float = 1.0     # Distance from the Star [AU]
    inclination : float = 90.0      # Orbital inclination [degrees] 
    eccentricity : float = 0.0      # Eccentricity 0 = circle
    arg_periastron: float = 90.0    # Omega [degrees]
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
    

    def impact_parameter(self,star):
        """"
        Formula for impact parameter

        Params
        -------
        star : class
        
        """
        inclination = np.radians(self.inclination)

        return ((self.distance_m / star.radius_m) * np.cos(inclination))
    

    def transit_duration(self,star):
        """
        Total duration of the transit (from firts to fourth contact) in hours

        Params
        ------
        star : class
        
        """
        # Ratio of the radii
        k = self.radius_m / star.radius_m

        # Impact parameter
        b = self.impact_parameter(star)

        arg = (star.radius_m / self.distance_m ) * np.sqrt((1+k)**2 - b**2)
        T_days = (self.period / np.pi) * np.arcsin(arg)

        ecc_factor = np.sqrt(1 - self.eccentricity**2) / (1 + self.eccentricity * np.sin(np.radians(self.arg_periastron)))

        T_days = T_days * ecc_factor
        return T_days * 24.0
    