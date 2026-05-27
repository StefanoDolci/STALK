"""
stalk.transit

----------
Transit model 

"""




import numpy as np
import matplotlib.pyplot as plt 
from dataclasses import dataclass



@dataclass
class LightCurve:

    time : np.ndarray
    flux : np.ndarray




    def plot(self):
        """
        Plot of the light curve

        """
        fig, ax = plt.subplots(figsize=(10,5))
        ax.plot(self.time, self.flux, color = 'steelblue', lw=1.5)
        ax.set_xlabel("Time from mid transit [days]", fontsize = 12)
        ax.set_ylabel("Normalized flux", fontsize = 12)
        ax.set_title("STALK Light Curve", fontsize = 13)
        ax.grid(True, alpha = 0.3)
        return fig, ax
    



def simulate(star, planet, n_points=500, grid_size=300, duration_factor=2.0):
    """
    Simulate the transit light curve

    Params
    ------
    star : Star
    
    planet : Planet

    n_points : int 
        Number of time samples

    grid_size : int 
        Pixel resolution of the stellar grid

    duration_dactor : float
        Total time window = duration_factor * transit_duration


    Returns
    -------
    LightCurve    
    """

    # Stellar grid
    xx, yy, intensity = star.intensity_grid(grid_size)
    total_flux = intensity.sum()

    # Time array, centered mid-transit [hours]
    transit_duration = planet.transit_duration(star)
    half_window = 0.5 * duration_factor * transit_duration
    times = np.linspace(-half_window, half_window, n_points)

    # Planet position (in stellar radii)

    









