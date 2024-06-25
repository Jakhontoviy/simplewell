import numpy as np
import client.services.log_processing_service as lps


class petrophysics:
    '''
    Petrophysics basic equations
    '''

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

        if type(rhob) == 'client.objects.log.Log' and type(rhob_matrix) == 'client.objects.log.Log' and type(rhob_fluid) == 'client.objects.log.Log':
            rhob, rhob_matrix, rhob_fluid = lps.LogProcessingService.interpolate_to_common_reference([
                                                                                                     rhob, rhob_matrix, rhob_fluid])
        por = (rhob_matrix - rhob) / (rhob_matrix - rhob_fluid)

        if limit_result is True:
            return np.clip(por, low_limit, high_limit)
        else:
            return por
        return por

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

    def vsh_from_rt(rt_log, rt_clean, rt_shale):
        vrt = (rt_shale/rt_log)*(rt_clean-rt_log)/(rt_clean-rt_shale)
        if (rt_log > 2 * rt_shale):
            vclrt = 0.5 * (2 * vrt) ** (0.67*(vrt+1))
        else:
            vclrt = vrt
        return vclrt

    def vsh_nd(neut_log, den_log, neut_clean1, den_clean1, neut_clean2, den_clean2, neut_clay, den_clay):
        term1 = (den_clean2-den_clean1)*(neut_log-neut_clean1)-(den_log-den_clean1)*(neut_clean2-neut_clean1)
        term2 = (den_clean2-den_clean1)*(neut_clay-neut_clean1)-(den_clay-den_clean1)*(neut_clean2-neut_clean1)
        vclnd = term1/term2
        return vclnd

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

    def sw_archie(porosity, rt, rw, archie_a = 1, archie_m = 2, archie_n = 2):
        '''
        References:
        Archimedes, R. C. (1994) Shale fluid properties. 3rd ed. London: John Wiley and Sons.
        https://www.spec2000.net/01-quickmath.htm
        '''
        sw = ((archie_a / (porosity ** archie_m)) * (rw/rt))**(1/archie_n)
        return sw


    def sw_simandoux(phie, rt, rw, vsh, rt_shale, archie_a = 1, archie_m = 2, archie_n = 2):
        A = (1 - vsh) * archie_a * rw / (phie ** archie_m)
        B = A * vsh / (2 * rt_shale)
        C = A / rt

        sw = ((B ** 2 + C)**0.5 - B) ** (2 / archie_n)
        return sw

    def sw_waxman_smits(swb,fluid_density, salinity):
        '''
        Solve Qv from Swb based on Hill, Shirley and Klein technique (1979)
        Parameters
        ----------
        fluid_density : float
            Density of the fluid
        Salinity : float
            Salinity of the fluid in parts per million (ppm)
         
        References:
        Hill, H.J., Shirley, O.J., Klein, G.E.: “Bound Water in Shaley Sands - Its Relation to Qv and Other Formation Properties”, Log Analyst, May-June 1979.
        '''
        qv = swb/(0.6425/((fluid_density*salinity)**0.5) + 0.22)
        # m* apparent = log10(Rw /(Rt*(1 + Rw*B*Qv))) / log10(PHIT)




