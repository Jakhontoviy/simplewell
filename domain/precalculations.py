import numpy as np
import client.services.log_processing_service as lps
import client.domain.charts_processing as chp

class precalculations:
    '''
    Petrophysics precalculations 
    '''

    def temp_gradient(bottom_hole_temperature, surface_temperature, bottom_hole_depth):
        """
        Temperature gradient calculation.
        Parameters
        ----------
        bottom_hole_temperature : float
            Bottom hole temperature (deg F or deg C)
        surface_temperature : float
            Surface temperature (deg F or deg C)
        bottom_hole_depth : float
            Bottom hole depth (ft or m)
        Returns
        -------
        float
            Returns temperature gradient in deg per depth unit (degF/ft or deg C/m)
        """
        gradient = (bottom_hole_temperature - surface_temperature) / bottom_hole_depth
        return gradient

    def formation_temperature_gradient(surface_temperature, gradient, depth):
        """
        Calculates formation temperature based on a gradient.
        Parameters
        ----------
        surface_temperature : float
            Surface temperature (deg F or deg C)
        gradient : float
            Temperature gradient (degF/ft or degC/m)
        depth : float
            Depth at which temperature is required (ft or m)
        Returns
        -------
        float
            Returns formation temperature at a entered depth
        """
        formation_temperature = surface_temperature + gradient * depth
        return formation_temperature
    
    def formation_temperature_2points(temp_surface=20, temp_bottom = 80, temp_bottom_depth = 2500, depth = 2300):
        """
        Calculates the formation temperature based on two points.

        Parameters:
            temp_surface (float): The surface temperature.
            temp_bottom (float): The temperature at the bottom.
            temp_bottom_depth (float): The depth at which the temperature at the bottom is measured.
            depth (float): The depth at which the formation temperature is required.

        Returns:
            float: The formation temperature at the specified depth.
        """
        temp_formation = temp_surface + (temp_bottom - temp_surface) / temp_bottom_depth * depth
        return temp_formation
    
    def resistivity_water_from_salinity(salinity = 30, temp_formation = 80):
        """
        Converts salinity to resistivity
        Parameters
        ----------
        salinity : float
            Salinity (parts-pers-thousand, ppt), example: 30
        temperature : float
            Temperature (degC)
        Returns
        -------
        float
            Returns resistivity in ohm-m
        """
        water_resisitivity = (222/(temp_formation+17.7)/salinity)**0.88
        return water_resisitivity
    
    

    def resisitivity_from_induction(ind_val, zond_type):
        """
        Calculate the resistivity value from the given induction value and zone type.

        Parameters:
            ind_val (float): The induction value.
            zond_type (str): The zond type, options availabe are:
                "4F0.75", "6F1", "3I1", "4I1", "7I1.6A", "7I1.6R", "8I1.4", "6E1", "IC4A".

        Returns:
            float or None: The calculated resistivity value, or None if the zone type is not recognized.
        """
        if zond_type == "4F0.75":
            rt_val = 10**(3.01074-0.99234*np.log10(ind_val)-0.01659*(np.log10(ind_val))**2)
        elif zond_type == "6F1":
            rt_val = 10**(2.9587-0.91352*np.log10(ind_val)-0.04752*(np.log10(ind_val))**2)
        elif zond_type == "3I1":
            rt_val = 1383.1*(ind_val)**(-1.1054)
        elif zond_type == "4I1":
            rt_val = 10**(3.01074-0.99234*np.log10(ind_val)-0.01659*(np.log10(ind_val))**2)
        elif zond_type == "7I1.6A":
            rt_val = 10**(3.755-2.872*np.log10(ind_val)+1.37307*(np.log10(ind_val))**2-0.33413*(np.log10(ind_val))**3)
        elif zond_type == "7I1.6R":
            rt_val = 10**(2.1522-1.17306*np.log10(ind_val)+0.32776*(np.log10(ind_val))**2-0.07849*(np.log10(ind_val))**3)
        elif zond_type == "8I1.4":
            rt_val = 0.663925-0.0002356*ind_val-29.4654/(ind_val**0.5)+1100.9287/ind_val-72.129/(ind_val**1.5)
        elif zond_type == "6E1":
            rt_val = -0.23231+869.755/ind_val+416.1173/(ind_val**1.5)-285.6378/(ind_val**2)
        elif zond_type == "IC4A":
            rt_val = 2507.8*ind_val**(-1.2952)
        else:
            rt_val = None
        return rt_val
    
    def resistivity_water_computation_archi(rt, phit, a = 1, m = 2, n = 2, sw = 1):
        """
        Calculate the water resistivity based on Dakhnov and Archi equation.

        Parameters:
            rt (float): Formation Resisitivity
            phit (float): Formation Porosity
            a (float, optional): The constant term in the equation. Defaults to 1.
            m (float, optional): The exponent in the equation. Defaults to 2.
            n (float, optional): The exponent in the equation. Defaults to 2.
            sw (float, optional): The initial value of the water resistance. Defaults to 1.

        Returns:
            float: The calculated water resistance.

        References:
            Dakhnov, A. and Archimedes, R. C. (1994) Shale fluid properties. 3rd ed. London: John Wiley and Sons.

        """
        return (sw**n*rt*phit**m)/a
    
    
    def aps_from_sp(sp, sp_sand = 20, sp_shale=180):
        """
        Calculate the Spontaneous Potential Alpha (APS) log from the Spontaneous Potential (SP) log.

        Parameters:
            sp (float): The Spontaneous Potential log value.
            sp_sand (float, optional): The Spontaneous Potential log value for sand. Defaults to 20.
            sp_shale (float, optional): The Spontaneous Potential log value for shale. Defaults to 180.

        Returns:
            float: The calculated Spontaneous Potential Alpha log.
        """
        return (sp_shale - sp)/(sp_shale - sp_sand)
    