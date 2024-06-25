import numpy as np

class petrophysics:
    '''
    Petrophysics basic equations
    '''

    def shale_volume_from_gammaray(log, min_value, max_value, method="linear", limit_result=True, low_limit=0, high_limit=1):
        """
        Calculates the volume of clay or shale from a gamma ray log.

        Maxvalue can be set to represent the 100% shale or clay value when working with either VShale or VClay.
        Non-linear equations can be selected using the method argument:
            larionov-young
            larionov-old
            steiber
            clavier

        Parameters
        ----------
        minvalue : float
            Value representing a 100% clean interval.
        maxvalue : float
            Value representing either 100% clay or 100% shale.
        inputvalue : float
            Gamma ray value from log measurements.
        method : string
            Select method for calculating VClay or VShale:
            - linear
            - larionov-young
            - larionov-old
            - steiber
            - clavier

                By default: linear
        limit_result : bool, optional
            Apply limits to the result value.
            By default False
        low_limit : int, optional
            Low limit. If value falls below this limit it will be set to this value.
            By default 0
        high_limit : float, optional
            High limit. If value falls above this limit it will be set to this value.
            By default: 1
        Returns
        -------
        float
            Returns a VShale or VClay in decimal unit.

        References
        -------
        Bhuyan, K. and Passey, Q. R. (1994) Clay estimation from GR and neutron-density porosity logs, SPWLA 35th Annual Logging Symposium, pp. 1–15.
        Larionov VV (1969) Borehole radiometry: Moscow, U.S.S.R., Nedra
        Steiber RG (1973) Optimization of shale volumes in open hole logs. J Petrol Technol 1973(31):147 – 162
        Asquith and Krygowski (2004) Basic Well Log Analysis. Second Edition.
        Original paper reference cannot be found for Clavier (1971), however, the following discussion provides insight to its origins.
        https://www.researchgate.net/post/Volume_fraction_of_shale_full_reference_for_Clavier_1971
        """

        vsh = (log - min_value)/(max_value - min_value)

        if method == "linear":
            result = vsh
        elif method == "larionov-young":
            result = 0.083 * ((2**(3.71 * vsh)) - 1)
        elif method == "larionov-old":
            result = 0.33 * ((2**(2 * vsh)) - 1)
        elif method == "steiber":
            result = igr / (3 - 2 * vsh)
        elif method == "clavier":
            result = 1.7 - ((3.38-(vsh + 0.7)**2)**0.5)
        else:
            raise Exception(
                "Enter a valid method value: linear, larionov-young, larionov-old, steiber, clavier")

        if limit_result is True:
            return np.clip(result, low_limit, high_limit)
        else:
            return result

    def shale_vol_from_sp(log, min_value, max_value, limit_result=True, low_limit=0, high_limit=1):
        """
        Calculates a clay or shale volume from SP log.

        Maxvalue can be set to represent 100% shale value when working with VShale, or it can
        be set to 100% clay when working with VClay.
        Parameters
        ----------
        minvalue : float
            Value representing a 100% clean interval.
        maxvalue : float
            Value representing either 100% clay or 100% shale.
        inputvalue : float
            Gamma ray value from log measurements.
        limit_result : bool, optional
            Apply limits to the result value.
            By default False
        low_limit : int, optional
            Low limit. If value falls below this limit it will be set to this value.
            By default 0
        high_limit : float, optional
            High limit. If value falls above this limit it will be set to this value.
            By default: 1
        Returns
        -------
        float
            Returns a VShale or VClay in decimal unit.

        References
        ----------
        Dewan, J. T., 1983, Essentials of modern open- hole log interpretation: PennWell Books, Tulsa, Oklahoma.
        """

        result = (log - min_value)/(max_value - min_value)

        if limit_result is True:
            return np.clip(result, low_limit, high_limit)
        else:
            return result

    def shale_vol_from_den_neut(
            neut_porosity, dens_porosity, neut_shale_porosity, dens_shale_porosity, limit_result=False, low_limit=0,
            high_limit=1):
        """
        Calculates a shale volume from density and neutron porosity logs.
        Parameters
        ----------
        neut_porosity : float
            Neutron porosity
        dens_porosity : float
            Density porosity
        neut_shale_porosity : [type]
            [description]
        dens_shale_porosity : [type]
            [description]
    limit_result : bool, optional
            Apply limits to the result value.
            By default False
        low_limit : int, optional
            Low limit. If value falls below this limit it will be set to this value.
            By default 0
        high_limit : float, optional
            High limit. If value falls above this limit it will be set to this value.
            By default: 1
        Returns
        -------
        float
            Returns a shale volume in decimal unit.
        References
        ----------
        Bhuyan, K. and Passey, Q. R. (1994) Clay estimation from GR and neutron-density porosity logs, SPWLA 35th Annual Logging Symposium, pp. 1–15.
        Dewan, J. T., 1983, Essentials of modern open- hole log interpretation: PennWell Books, Tulsa, Oklahoma.
        """
        result = (neut_porosity - dens_porosity) / (neut_shale_porosity - dens_shale_porosity)

        if limit_result is True:
            return np.clip(result, low_limit, high_limit)
        else:
            return result

    def shale_vol_from_rt(rt_log, rt_clean, rt_shale):
        """
        Calculate the shale volume from the resistivity logs.

        Parameters:
            rt_log (float): The resistivity log value.
            rt_clean (float): The clean resistivity log value.
            rt_shale (float): The shale resistivity log value.

        Returns:
            float: The calculated shale volume.

        Calculates the shale volume based on the resistivity logs. The formula used is:

        vsh = (rt_shale / rt_log) * (rt_clean - rt_log) / (rt_clean - rt_shale)

        If rt_log is greater than 2 * rt_shale, the shale volume is calculated using the formula:

        vsh = 0.5 * (2 * vrt) ** (0.67 * (vrt + 1))

        Otherwise, the shale volume is simply vrt.

        Note:
            - The resistivity logs are assumed to be in ohmm.
            - The shale volume is returned in decimal units.
        """
        
        vrt = (rt_shale/rt_log)*(rt_clean-rt_log)/(rt_clean-rt_shale)
        if (rt_log > 2 * rt_shale):
            vsh = 0.5 * (2 * vrt) ** (0.67*(vrt+1))
        else:
            vsh = vrt
        return vsh

    def shale_vol_from_neut_den(neut_log, den_log, neut_clean1, den_clean1, neut_clean2, den_clean2, neut_clay, den_clay):
        """
        A function to calculate shale volume based on neutron-density logs.
        Parameters:
            neut_log (float): Neutron log value.
            den_log (float): Density log value.
            neut_clean1 (float): Neutron clean1 value.
            den_clean1 (float): Density clean1 value.
            neut_clean2 (float): Neutron clean2 value.
            den_clean2 (float): Density clean2 value.
            neut_clay (float): Neutron clay value.
            den_clay (float): Density clay value.
        Returns:
            float: The calculated shale volume.
        """
        term1 = (den_clean2-den_clean1)*(neut_log-neut_clean1)-(den_log-den_clean1)*(neut_clean2-neut_clean1)
        term2 = (den_clean2-den_clean1)*(neut_clay-neut_clean1)-(den_clay-den_clean1)*(neut_clean2-neut_clean1)
        return term1/term2
    

    def vclay_from_vshale(vshale, multiplier=1.2):
        """
        Converts a shale volume to clay volume using a multiplier.
        Parameters
        ----------
        vshale : float
            Shale volume
        multiplier : float
            Shale to clay multiplier (decimal)
        Returns
        -------
        float
            Returns a clay volume.
        References
        ----------
        Bhuyan, K. and Passey, Q. R. (1994) Clay estimation from GR and neutron-density porosity logs, SPWLA 35th Annual Logging Symposium, pp. 1–15.
        """
        return vshale * multiplier


    def porosity_from_density(rhob, rhob_matrix=2.65, rhob_fluid=1, limit_result=True, low_limit=0, high_limit=0.5):
        """
        Calculates the porosity from the bulk density using the petrophysics equation.

        Parameters:
            rhob (array): The bulk density value.
            rhob_matrix (float): The density value of the matrix.
            rhob_fluid (float): The density value of the fluid.

        Returns:
            float: The calculated porosity value.

        Raises:
            None
        """

        por = (rhob_matrix - rhob) / (rhob_matrix - rhob_fluid)

        if limit_result is True:
            return np.clip(por, low_limit, high_limit)
        else:
            return por
        return por
    
    
    def porosity_effective_from_neutron(neut_por_log, vsh_log, shale_porosity=0.3):
        """
        Calculate the effective porosity from a neutron porosity log and a shale log.

        Parameters:
            neut_por_log (float or array-like): The neutron porosity log values.
            vsh_log (float or array-like): The shale log values.
            shale_porosity (float, optional): The porosity value of the shale. Defaults to 0.3.

        Returns:
            float: The calculated effective porosity.
        """
        return neut_por_log - (vsh_log*shale_porosity)
    
    def porosity_effective_from_sonic(sonic_por_log, vsh_log, shale_porosity=0.3):
        """
        Calculate the effective porosity from a sonic porosity log and a shale log.

        Parameters:
            sonic_por_log (float or array-like): The sonic porosity log values.
            vsh_log (float or array-like): The shale log values.
            shale_porosity (float, optional): The porosity value of the shale. Defaults to 0.3.

        Returns:
            float: The calculated effective porosity.
        """
        return sonic_por_log - (vsh_log*shale_porosity)



    def sw_archie(porosity, rt, rw, archie_a = 1, archie_m = 2, archie_n = 2):
        """
        Calculate the water saturation using the Archie equation.

        Parameters:
            porosity (float): Porosity of the formation.
            rt (float): True resistivity of the formation (ohmm).
            rw (float): Formation water resistivity (ohmm).
            archie_a (float, optional): Archie's parameter A (tortuosity factor, unitless). Defaults to 1.
            archie_m (float, optional): Archie's parameter M (cementation constant, unitless). Defaults to 2.
            archie_n (float, optional): Archie's parameter N (water saturation exponent, unitless). Defaults to 2.

        Returns:
            float: Water saturation calculated using the Archie equation.

        References:
        Archie, G. E., (1942). The Electrical Resistivity Log as an Aid in Determining Some 
        Reservoir Characteristics. SPE Journal, 146 (1), pp. 54-62.
        """
        
        sw = ((archie_a / (porosity ** archie_m)) * (rw/rt))**(1/archie_n)
        return sw


    def sw_simandoux(phie, rt, rw, vsh, rt_shale, archie_a = 1, archie_m = 2, archie_n = 2):
        """
        Calculate the water saturation using the Simandoux equation.

        Parameters:
            phie (float): Effective Porosity (fraction)
            rt (float): True Resistivity (ohmm)
            rw (float): Formation Water Resistivity (ohmm)
            vsh (float): Bulk Volume Fraction of Shale (fraction)
            rt_shale (float): Shale Resistivity (ohmm)
            archie_a (float, optional): Archie's parameter A (tortuosity factor, unitless). Defaults to 1.
            archie_m (float, optional): Archie's parameter M (cementation constant, unitless). Defaults to 2.
            archie_n (float, optional): Archie's parameter N (water saturation exponent, unitless). Defaults to 2.

        Returns:
            float: Water saturation calculated using the Simandoux equation (fraction)

        References:
            Simandoux, P. (1963). Mesuresd ielectriques en milieu poreux, application a mesure des 
            saturations en eau, Etude du Comportment des massifs Argileux. Supplementary Issue, 
            Revue de I’Institut Francais du Petrol
        """

        A = (1 - vsh) * archie_a * rw / (phie ** archie_m)
        B = A * vsh / (2 * rt_shale)
        C = A / rt

        sw = ((B ** 2 + C)**0.5 - B) ** (2 / archie_n)
        return sw

 
    def sw_waxman_smits(phie, vsh, resd, rw_ft, ft, a=1, m=2, n=2, dens_matrix=2.65):
        """
        Calculate water saturation from CEC method using the Waxman-Smits equation.

        Parameters:
            phie (float): Effective porosity (fractional).
            vsh (float): Shale volume (fractional).
            resd (float): Deep resistivity log reading (ohmm).
            rw_ft (float): Water resistivity at formation temperature (ohmm).
            ft (float): Formation temperature (degrees Fahrenheit or Celsius).
            rw: Water resistivity at 25 degrees celsius (ohm-m)
            a (float, optional): Formation factor (unitless). Default is 1.
            m (float, optional): Cementation exponent (unitless). Default is 2.
            n (float, optional): Saturation exponent (unitless). Default is 2.
            dens_matrix (float, optional): Matrix density (gm/cc or kg/m3). Default is 2.65.

        Returns:
            sw_cec (float): Water saturation from CEC method (fractional).
        """
        if phie > 0:
            cec = 10 ** (1.9832 * vsh - 2.4473)
            rw = (rw_ft) * (ft + 21.5) / 46.5
            b = 4.6 * (1 - 0.6 * np.exp(-0.77 / rw))
            f = a / (phie ** m)
            qv = cec * (1 - phie) * dens_matrix / phie
            sw_cec = 0.5 * ((- b * qv * rw) + ((b * qv * rw) ** 2 + 4 * f * rw_ft / resd) ** 0.5) ** (2 / n)
        else:
            sw_cec = 1
        return sw_cec





