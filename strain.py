from data import ternary_alloys_bowings as bow, binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()
ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]
temperature = np.linspace(0.0, 500.0, 1000)

parameter_keys = ["Eg", "VBO", "delta_so", "alc_temp", "alpha", "beta", "ac", "av", "b", "c11", "c12"]
var = np.linspace(0.001, 0.999, 1000)
fixed = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

ternary_params = dict()
for mat in ternary_mat:
    ternary_params[mat] = dict()

quaternary_params = {"fix_x": dict(), "fix_y": dict()}
for y in fixed:
    key = '{0:1.1f}'.format(y)
    quaternary_params["fix_x"][key] = dict()
    quaternary_params["fix_y"][key] = dict()
ternary_params_with_temperature = []
quaternary_params_with_temperature = []
for T in temperature:
    for par_key in parameter_keys:

        if par_key == "Eg":
            p_AlAs = par.AlAs[par_key] - par.AlAs["alpha"] * T ** 2 / (T + par.AlAs["beta"])
            p_GaAs = par.GaAs[par_key] - par.GaAs["alpha"] * T ** 2 / (T + par.GaAs["beta"])
            p_AlSb = par.AlSb[par_key] - par.AlSb["alpha"] * T ** 2 / (T + par.AlSb["beta"])
            p_GaSb = par.GaSb[par_key] - par.GaSb["alpha"] * T ** 2 / (T + par.GaSb["beta"])
        elif par_key == "alc_temp":
            p_AlAs = par.AlAs[par_key](T)
            p_GaAs = par.GaAs[par_key](T)
            p_AlSb = par.AlSb[par_key](T)
            p_GaSb = par.GaSb[par_key](T)
        else:
            p_AlAs = par.AlAs[par_key]
            p_GaAs = par.GaAs[par_key]
            p_AlSb = par.AlSb[par_key]
            p_GaSb = par.GaSb[par_key]

        if par_key == "Eg":
            b_AlGaAs = bow.AlGaAs[par_key](var)
            b_AlGaSb = bow.AlGaSb[par_key](var)
            b_AlAsSb = bow.AlAsSb[par_key]
            b_GaAsSb = bow.GaAsSb[par_key]
        else:
            b_AlGaAs = bow.AlGaAs[par_key]
            b_AlGaSb = bow.AlGaSb[par_key]
            b_AlAsSb = bow.AlAsSb[par_key]
            b_GaAsSb = bow.GaAsSb[par_key]

        ternary_params["AlGaAs"][par_key] = var * p_AlAs + (1 - var) * p_GaAs - var * (1 - var) * b_AlGaAs
        ternary_params["AlGaSb"][par_key] = var * p_AlSb + (1 - var) * p_GaSb - var * (1 - var) * b_AlGaSb
        ternary_params["AlAsSb"][par_key] = var * p_AlAs + (1 - var) * p_AlSb - var * (1 - var) * b_AlAsSb
        ternary_params["GaAsSb"][par_key] = var * p_GaAs + (1 - var) * p_GaSb - var * (1 - var) * b_GaAsSb

        for x in fixed:
            p_AlGaAs = x * p_AlAs + (1 - x) * p_GaAs - x * (1 - x) * b_AlGaAs
            p_AlGaSb = x * p_AlSb + (1 - x) * p_GaSb - x * (1 - x) * b_AlGaSb
            p = x * (1 - x) * (var * p_AlGaAs + (1 - var) * p_AlGaSb) + \
                var * (1 - var) * (x * ternary_params["AlAsSb"][par_key] + (1 - x) * ternary_params["GaAsSb"][par_key])
            quaternary_params["fix_x"]['{0:1.1f}'.format(x)][par_key] = p / (x * (1 - x) + var * (1 - var))
        for y in fixed:
            p_AlAsSb = y * p_AlAs + (1 - y) * p_AlSb - y * (1 - y) * b_AlAsSb
            p_GaAsSb = y * p_GaAs + (1 - y) * p_GaSb - y * (1 - y) * b_GaAsSb
            p = var * (1 - var) * (y * ternary_params["AlGaAs"][par_key] + (1 - y) * ternary_params["AlGaSb"][par_key]) \
                + y * (1 - y) * (var * p_AlAsSb + (1 - var) * p_GaAsSb)
            quaternary_params["fix_y"]['{0:1.1f}'.format(y)][par_key] = p / (var * (1 - var) + y * (1 - y))

    ternary_params_with_temperature.append(ternary_params)
    quaternary_params_with_temperature.append(quaternary_params)

deformation_keys = ["eps_par", "eps_orth", "dEc_hydro", "dEv_hydro", "dEv_biax", "dEv_biax_plus", "dEv_biax_minus"]
deformation_params = []
# GaAs is our substrate
for i, T in enumerate(temperature):
    deformation_params.append(dict())
    deformation_params[i] = {"fix_x": dict(), "fix_y": dict()}
    for fix_arg in ["fix_x", "fix_y"]:
        for j, fix in enumerate(fixed):
            fix_key = '{0:1.1f}'.format(fix)
            deformation_params[i][fix_arg][fix_key] = dict()

            deformation_params[i][fix_arg][fix_key]["eps_par"] = \
                (par.GaAs["alc_temp"](T) - quaternary_params_with_temperature[i][fix_arg][fix_key]["alc_temp"]) / \
                quaternary_params_with_temperature[i][fix_arg][fix_key]["alc_temp"]

            deformation_params[i][fix_arg][fix_key]["eps_orth"] = \
                -2.0 * quaternary_params_with_temperature[i][fix_arg][fix_key]["c12"] / \
                quaternary_params_with_temperature[i][fix_arg][fix_key]["c12"] * \
                deformation_params[i][fix_arg][fix_key][
                    "eps_par"]

            deformation_params[i][fix_arg][fix_key]["dEc_hydro"] = \
                quaternary_params_with_temperature[i][fix_arg][fix_key]["ac"] * \
                (deformation_params[i][fix_arg][fix_key]["eps_orth"] + 2 * deformation_params[i][fix_arg][fix_key][
                    "eps_par"])

            deformation_params[i][fix_arg][fix_key]["dEv_hydro"] = \
                quaternary_params_with_temperature[i][fix_arg][fix_key]["av"] * \
                (deformation_params[i][fix_arg][fix_key]["eps_orth"] + 2 * deformation_params[i][fix_arg][fix_key][
                    "eps_par"])

            deformation_params[i][fix_arg][fix_key]["dEv_biax"] = \
                quaternary_params_with_temperature[i][fix_arg][fix_key]["b"] * \
                (deformation_params[i][fix_arg][fix_key]["eps_orth"] - deformation_params[i][fix_arg][fix_key][
                    "eps_par"])

            delta_so = quaternary_params_with_temperature[i][fix_arg][fix_key]["delta_so"]
            dEv_biax = deformation_params[i][fix_arg][fix_key]["dEv_biax"]
            deformation_params[i][fix_arg][fix_key]["dEv_biax_plus"] = 1 / 2 * (dEv_biax - delta_so + np.sqrt(
                9 * dEv_biax ** 2 + 2 * dEv_biax * delta_so + delta_so ** 2))

print('Program finished in ', time.time() - start, ' s.')
