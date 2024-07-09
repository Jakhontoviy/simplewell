import numpy as np
import client.services.log_processing_service as lps


class shale_volume:
    '''
    Petrophysics basic equations
    '''

    def shale_volume_from_gammaray(gr_log, gr_clean, gr_shale, method="linear", limit_result=True, low_limit=0, high_limit=1):
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
        gr_clean : float
            Value representing a 100% clean interval.
        gr_shale : float
            Value representing either 100% clay or 100% shale.
        gr_log : float
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

        vsh = (gr_log - gr_clean)/(gr_shale - gr_clean)

        if method == "linear":
            result = vsh
        elif method == "larionov-young":
            result = 0.083 * ((2**(3.71 * vsh)) - 1)
        elif method == "larionov-old":
            result = 0.33 * ((2**(2 * vsh)) - 1)
        elif method == "steiber":
            result = vsh / (3 - 2 * vsh)
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
        """
        Calculates the volume of shale from the ratio of resistivity logs.

        Parameters:
            rt_log (float): The resistivity log value.
            rt_clean (float): The resistivity log value of clean formation.
            rt_shale (float): The resistivity log value of shale formation.

        Returns:
            float: The volume of shale calculated from the ratio of resistivity logs.
        """
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
    

class porosity:

    def porosity_total_from_density(den, den_matrix=2.65, den_fluid=1, limit_result=True, low_limit=0, high_limit=0.5):
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

        phit_d = (den_matrix - den) / (den_matrix - den_fluid)

        if limit_result is True:
            return np.clip(phit_d, low_limit, high_limit)
        else:
            return phit_d
        
    def porosity_shale_from_density(den_shale, den_matrix, den_fl, limit_result=True, low_limit=0, high_limit=0.5):
        phit_d_shale = (den_shale - den_matrix) / (den_fl - den_matrix)
    
        if limit_result is True:
            return np.clip(phit_d_shale, low_limit, high_limit)
        else:
            return phit_d_shale
    
    def porosity_effective_from_density(den, den_matrix, den_fl, den_shale, vsh, limit_result=True, low_limit=0, high_limit=0.5):
        phit_d = (den - den_matrix) / (den_fl - den_matrix)
        phit_d_sh = (den_shale - den_matrix) / (den_fl - den_matrix)
        phie_d = phit_d - vsh * phit_d_sh
        
        if limit_result is True:
            return np.clip(phie_d, low_limit, high_limit)
        else:
            return phie_d


    def porosity_shale_from_sonic_willie(dt_shale, dt_matrix, dt_fluid, limit_result=True, low_limit=0, high_limit=0.5):
        '''
        Willie-TimeAverage
        '''
        phit_s_shale = (dt_shale-dt_matrix)/(dt_fluid-dt_matrix)
    
        if limit_result is True:
            return np.clip(phit_s_shale, low_limit, high_limit)
        else:
            return phit_s_shale


    def porosity_total_from_sonic_willie(dt_log, dt_matrix, dt_fluid, cp = 1, limit_result=True, low_limit=0, high_limit=0.5):
        '''
        Willie-TimeAverage
        '''
        phit_s=(1/cp)*(dt_log-dt_matrix)/(dt_fluid-dt_matrix)

        if limit_result is True:
            return np.clip(phit_s, low_limit, high_limit)
        else:
            return phit_s
    
    def porosity_effective_from_sonic_willie(dt_log, dt_matrix, dt_fluid, dt_shale, vsh, cp = 1, limit_result=True, low_limit=0, high_limit=0.5):
        '''
        Willie-TimeAverage
        '''
        phit_s=(1/cp)*(dt_log-dt_matrix)/(dt_fluid-dt_matrix)
        phit_s_sh = (dt_shale-dt_matrix)/(dt_fluid-dt_matrix)
        phie_s = phit_s - vsh * phit_s_sh

        if limit_result is True:
            return np.clip(phie_s, low_limit, high_limit)
        else:
            return phie_s
    
    
    def porosity_total_from_sonic_rhg(dt_log, dt_matrix, alpha = 0.67, limit_result=True, low_limit=0, high_limit=0.5):
        '''
        #Raymer-Hunt-Gardner (the alpha(5/8) ranges from 0.625-0.70, 0.67-most, 0.60-gas reservoirs)
        '''
        phit_s = (alpha)*(dt_log-dt_matrix)/(dt_log)
    
        if limit_result is True:
            return np.clip(phit_s, low_limit, high_limit)
        else:
            return phit_s
    
    def porosity_effective_from_acoustic_rhg(dt_log, dt_matrix, dt_fluid, dt_shale, vsh, alpha = 0.67, limit_result=True, low_limit=0, high_limit=0.5):
        '''
        #Raymer-Hunt-Gardner (the alpha(5/8) ranges from 0.625-0.70, 0.67-most, 0.60-gas reservoirs)
        '''
        phit_s = alpha * (dt_log-dt_matrix)/(dt_log)
        phit_s_shale = (dt_shale-dt_matrix)/(dt_fluid-dt_matrix)
        phie_s = phit_s - vsh * phit_s_shale
        
        if limit_result is True:
            return np.clip(phie_s, low_limit, high_limit)
        else:
            return phie_s
    
    def porosity_effective_from_neutron(neut, neut_sh, vsh, limit_result=True, low_limit=0, high_limit=0.5):
        phie_n = (neut-vsh*neut_sh)

        if limit_result is True:
            return np.clip(phie_n, low_limit, high_limit)
        else:
            return phie_n

        #Neutron-Density
    def porosity_total_from_neutron_density(phin, phid):
        phit_nd= (phin + phid) / 2
        return phit_nd

    def porosity_total_from_neutron_density_gas_corrected(phin, phid):
        phit_nd_gas_corr = ((phin**2 + phid**2)/2)**(0.5)    #for gas intervals (nphi<dphi = crossover)
        return phit_nd_gas_corr

    

class saturation:

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

    def sw_waxman_smits_temp(swb,fluid_density, salinity):
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

    def sw_waxmansmits(rw, T, RwT, rt, PHIT, PHIE, den_fl, Swb, Rw75, Qv, B, m_cem = 2, mslope = 3):
        '''
        waxmansmits(Rw, Rt, PhiT, aa, mm, CEC)
        **Waxman-Smits CEC method obtains Qv from Hill, Shirley and Klein
          Eq solved for n=2
        *Input parameters:
         - PHIT - total porosity
         - m_cem -  cementation exponent is adjusted for Swb
         - Rw - formation water resistivity ohmm
         - B - cation mobility (mho cm2 / meq)
         - Qv - concentration of exchange cations per volume unit (meq/ml pore space)
         - CEC - cation exchange capacity of shale(meq/100 gm of sample)
         - den_ma - mineral graind density (g/cc)
         - m_cem - best determined from SCAL      
        *Returns:
         - Sw_Total_WS - total water saturation from Waxman-Smits
         - Sw_WS =(  ((1/PHIT**mstar)*Rw)/Rt*(1+Rw*B*QV)/Sw )**(1/nstar)
         - RwT = T        # Temperature of Rw measument
         - T = 150.       # Reservoir temperature in DegF
         - TC=(T-32)/1.8  # Temp DegC
        '''
        
        #convert m_cem to mstar with increase in mstar with increase in Swb
        mstar = m_cem + mslope*Swb

        Rw75=((RwT+6.77)*rwa)/(75+6.77)

        # Salinity in KPPM
        SAL=(10**((3.562-math.log10(Rw75-0.0123))/0.955))/1000

        B = math.exp(7.3-28.242/math.log(T)-0.2266*math.log(rwa)) 

        Bdacy=(1-0.83*math.exp(-math.exp(-2.38+(42.17/TC))/rwa))*(-3.16+1.59*math.log(TC))**2 #SCA Paper SCA2006-29

        #Crain's Waxman-Smits in lieu of using iterative. 
        #Swc = 0.5 * ((- B * Qv * RW2) + ((B * Qv *  RW2)**2 + 4 * F * RW@FT / RESD) ^ 0.5)**(2 / N)
        swT = 0.5 * (    (- B*Qv*Rw75) + (  (B*Qv*Rw75)**2  +  4*(1/PHIT**mstar)*rw/rt)**0.5)**(2/2)
        return swT


    
    def sw_dualwater(Rw, T, RwT, Rt, PHIT, PHIE, Swb):
        '''
        dualwater(Rw, Rt, PHIT, por_shale, Vsh, Rsh)
        **Dual-Water (clavier, 1977) with later modifications/rearrangements.
          Formulas from Doveton "Principles of mathematical petrophysics"
        *Input parameters:
         - PHIT - total porosity
         - por_shale - shale porosity
         - Rw - formation water resistivity [ohmm]
         - Swb - clay-bound water saturation 
         - Sw_dw - total water saturation
         *Returns:
         - Sw_dw - Total water saturation (or water saturation in Total pore space)
         - CBVWT - Coates DW CBW in Total Porossity system

         1. Coates, G.R., Gardner, J.S., and Miller, D.L., 1994, 
            Applying pulse-echo NMR to shaly sand formation evaluation, 
            paper B, 35th Annual SPWLA Logging Symposium Transactions, 22 p.
        '''

        #--------------------------------------------------
        #
        #  BEGINNING OF MRIAN AS PROPOSED BY COATES, et al
        #
        #----- COMPUTE BASIC DATA -------------------------
        #RwT = T
        CW = (T+7.)/(Rw*(RwT+7.))
        #----  CCW FROM COATES   
        CCW = 0.000216*(T-16.7)*(T+504.4)
        Ct   = 1/Rt
        #--------------------------------------------------
        #Swb = 1.0 - (PHIE/PHIT)
        #Swia=BVI/PHIE, but no NMR BVI here
        Swia = Swb #estimate
        #--------------------------------------------------
        CBW = Swb * PHIT
        #----- COMPUTE DCW ------------------------------
        ALPHA = 1 #for now with Salinity > 40,000 ppm
        #DCWW=CW+ALPHA*Swb*(CCW-CW)
        #        DCWI=CW+ALPHA*SWBI*(CCW-CW)
        #----- W @ Sw = 1 -----------------------------------
        #WW=math.log10(Ct/DCWW)/(math.log10(PHIT))
        #----- W @ Sw AT BVI --------------------------------
        #        WI=LOG10(CT/DCWI)/(LOG10(BVIT))
        #----- THE ACTUAL W ---------------------------------
        #Wq = 0.4*Swia + 1.65
        Wq = 0.4*Swia + 1.9
        #----- WW AND WI CONSTRAN W ------------------------
        #----- COMPUTE CBVW TOTAL -----------------------
        #AA=CW
        #BB=ALPHA*(CBW)*(CW-CCW)
        #CC=Ct
        #CBVWA = (BB + math.sqrt(BB*BB + 4*AA*CC))/(2*AA)
        CBVWA = (ALPHA*(CBW)*(CW-CCW) + ((ALPHA*(CBW)*(CW-CCW))**2 + 4*CW*Ct)**(1/2))/(2*CW)
        CBVWT = CBVWA**(2/Wq)
        #---- COMPUTE Bulk Volume Water CBVWE in Effective System ----------
        #CBVWE = CBVWT-CBW    
        Sw_dw = CBVWT/PHIT
        return Sw_dw
        #------------------------------------------------------------
        #      END OF GEORGE COATES' MRIAN                                        
        #-------------------------------------------------------------

    def watervolume_from_watersaturation(sw, phit, limit_result=True, low_limit=0, high_limit=1):
        bvw = sw * phit

        if limit_result is True:
            return np.clip(bvw, low_limit, high_limit)
        else:
            return bvw
