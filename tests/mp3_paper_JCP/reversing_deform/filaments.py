import numpy as np
from helpers import nested_dict

def calc_filament_diags(geo_data, field_data, nys):
    filament_diags = nested_dict()

    taus = np.linspace(0.10, 1.0, 19)
    eps = 1e-12

    for ny in nys:
        g = geo_data[ny]['g']

        for opt in field_data[ny]:
            lfs = []
            for tau in taus:
                area_0 = np.sum(g[field_data[ny][opt]['0.0']['cb'] >= tau - eps])
                area_h = np.sum(g[field_data[ny][opt]['2.5']['cb'] >= tau - eps])
                if area_0 < eps:
                    lf = 0
                else:
                    lf = 100 * area_h / area_0
                lfs.append(lf)

            filament_diags[ny][opt] = list(zip(taus, lfs))

    return filament_diags
