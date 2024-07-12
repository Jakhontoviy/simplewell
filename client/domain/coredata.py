class capillary_data:

    '''
    capillary data calculation equations
    '''

    def calculate_j_ift(pc_lab, perm, poro, sigma_lab=72, costh_lab=1):
        """
        Calculates the value of J_ift using the given parameters.

        Parameters:
            pc_lab (float): The value of pressure in laboratory conditions, bar
            perm (float): The value of permeability, mD
            poro (float): The value of porosity, v/v
            sigma_lab (float): The value of sigma_lab, dyn/cm. default value is 72
            costh_lab (float): The value of costh_lab, default value is 1

        Returns:
            float: The calculated value of J_ift.
        """

        j_ift = 0.31832*(pc_lab/(sigma_lab*costh_lab))*((perm/poro)**0.5)
        return j_ift