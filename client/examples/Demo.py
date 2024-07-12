import client.domain.petrophysics as pp
import numpy as np

rhob_vals_list = np.array([2.4, 2.3, 2])


por = pp.porosity.porosity_total_from_density(rhob_vals_list, 2.65, 1.0)

por_shale = pp.porosity.porosity_shale_from_density(2.4, 2.65, 1.0)

vsh = pp.shale_volume.vsh_from_rt(6, 2, 4)

print(por.round(3))

print(por_shale)

print(vsh)

