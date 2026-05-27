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
    radius : float = 1.0                # Radius of the planet in Earth radii
    period : float = 365.25             # Period od the orbit [days]
    semi_maj_axis : float = 1.0         # Distance from the Star [AU]
    inclination : float = 90.0          # Orbital inclination [degrees] 
    eccentricity : float = 0.0          # Eccentricity 0 = circle
    arg_periastron: float = 90.0        # omega [degrees]
    long_ascending_node : float = 0.0   # Omega [degrees]
    t_periastron : float = 0.0          # Time of the periastron passage [days]
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
    



    def solve_kepler(self, M, tol=1e-10, max_iter=100):
        """
        Solve Kepler's equation for the eccentric anomaly E,

        with Netwon-Raphson iteration.

        M = E - e*sin(E)

        Params
        ------
        M : float or array
            Mean anomaly [radians]  
            M(t) = 2pi/P * (t - tp)
            where P is the period and tp is the time of passage from periastron
        
        tol : float
            Convergence tollerance

        max_iter : int
            Number of maximum iteration 

        Returns
        -------
        E : shape as M 
            Eccentric anomaly [radians]    
        """
        M_was_scalar=np.isscalar(M)
        M = np.atleast_1d(np.asarray(M, dtype=float))


        e = self.eccentricity

        # M to [-pi,pi]
        M  = np.mod(M + np.pi, 2.0* np.pi) -np.pi

        E = M + 0.85 * e * np.sign(np.sin(M))

        for _ in range(max_iter):
            f = E - e * np.sin(E) - M
            fp = 1.0 - e * np.cos(E)
            delta = f/fp
            E = E - delta 

            if np.all(np.abs(delta)<tol):
                break

        return E.item() if M_was_scalar else E
    


    def sky_position(self, times, star):

        """
        Projected sky-plane position of the planet at given time 


        Params:
        -------

        times : float / array
            Times [days]. Reference is t_periastron

        star: Star
            Host star

        
        Returns: 
        --------
        x_sky, y_sky, z_sky: np.ndarray
            Sky-plane coord. in stellar radii 

            z_sky > 0 means planet in front star
            z_sky < 0 means planet behind star

        """

        times = np.atleast_1d(np.asarray(times, dtype=float))

        e = self.eccentricity

        # 1. Mean anomaly
        M = 2.0 * np.pi * (times - self.t_periastron) / self.period

        # 2. Eccentric anomaly 
        E = np.atleast_1d(self.solve_kepler(M))

        # 3. True anomaly nu (
        nu = 2.0 * np.arctan2(
            np.sqrt(1.0 + e) * np.sin(0.5 * E),
            np.sqrt(1.0 - e) * np.cos(0.5 * E),
        )

        # 4. Heliocentric distance in units of stellar radii
        a_over_R = self.distance_m / star.radius_m
        r = a_over_R * (1.0 - e * np.cos(E))

        # 5. Position in the orbital plane (periastron along + x_orb axis)
        x_orb = r * np.cos(nu)
        y_orb = r * np.sin(nu)

        # 6. Rotation angles
        omega = np.radians(self.arg_periastron)
        inc   = np.radians(self.inclination)
        Omega = np.radians(self.long_ascending_node)

        cos_w, sin_w = np.cos(omega), np.sin(omega)
        cos_i, sin_i = np.cos(inc),   np.sin(inc)
        cos_O, sin_O = np.cos(Omega), np.sin(Omega)

        # 7. Rotation orbital plane -> observer frame

        x_sky = ( cos_O * cos_w - sin_O * sin_w * cos_i) * x_orb \
              + (-cos_O * sin_w - sin_O * cos_w * cos_i) * y_orb
        y_sky = ( sin_O * cos_w + cos_O * sin_w * cos_i) * x_orb \
              + (-sin_O * sin_w + cos_O * cos_w * cos_i) * y_orb
        z_sky = (sin_w * sin_i) * x_orb + (cos_w * sin_i) * y_orb

        return x_sky, y_sky, z_sky

