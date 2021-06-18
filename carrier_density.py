import random

from data import ternary_alloys_bowings as bow, binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
import scipy.integrate as integrate
from matplotlib.gridspec import GridSpec
from scipy.optimize import brentq as bq
from matplotlib import cm


def GaAs_Ec_DOS(E, T):
    Eg_T = Eg["GaAs"] - (alpha["GaAs"] * (T ** 2) / (T + beta["GaAs"]))
    tmp = E - VBO["GaAs"] - Eg_T

    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_e["GaAs"] ** (3 / 2) * tmp ** (1 / 2)
    return dos


def GaAs_Ehh_DOS(E, T):
    tmp = -(E - VBO["GaAs"])
    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_hh["GaAs"] ** (3 / 2) * tmp ** (1 / 2)
    return dos


def GaAs_Elh_DOS(E, T):
    tmp = -(E - VBO["GaAs"])
    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_lh["GaAs"] ** (3 / 2) * tmp ** (1 / 2)
    return dos


def GaAs_Esh_DOS(E, T):
    tmp = -(E - VBO["GaAs"] + delta_so["GaAs"])
    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_so["GaAs"] ** (3 / 2) * tmp ** (1 / 2)
    return dos


def GaAs_DOS(E, T):
    return np.array([GaAs_Ec_DOS(E, T), GaAs_Elh_DOS(E, T), GaAs_Ehh_DOS(E, T), GaAs_Esh_DOS(E, T)])


#  ---------------------------------------------------------------------------------------------------


def quat_Ec_DOS(E, T):
    m_e_tmp = np.zeros(len(temperature))
    e = np.zeros(len(temperature))
    for i in range(len(temperature)):
        m_e_tmp[i] = quaternary_params_with_temperature[i]["m_e"]
        e[i] = energies[i]["Ec"]
    m_e_val = np.interp(T, temperature, m_e_tmp)
    en = np.interp(T, temperature, e)
    tmp = E - en
    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_e_val ** (3 / 2) * tmp ** (1 / 2)
    return dos


def quat_Ehh_DOS(E, T):
    m_e_tmp = np.zeros(len(temperature))
    e = np.zeros(len(temperature))
    for i in range(len(temperature)):
        m_e_tmp[i] = quaternary_params_with_temperature[i]["m_hh"]
        e[i] = energies[i]["Ev_hh"]
    m_e_val = np.interp(T, temperature, m_e_tmp)
    en = np.interp(T, temperature, e)
    tmp = -(E - en)
    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_e_val ** (3 / 2) * tmp ** (1 / 2)
    return dos


def quat_Elh_DOS(E, T):
    m_e_tmp = np.zeros(len(temperature))
    e = np.zeros(len(temperature))
    for i in range(len(temperature)):
        m_e_tmp[i] = quaternary_params_with_temperature[i]["m_lh"]
        e[i] = energies[i]["Ev_lh"]
    m_e_val = np.interp(T, temperature, m_e_tmp)
    en = np.interp(T, temperature, e)
    tmp = -(E - en)
    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_e_val ** (3 / 2) * tmp ** (1 / 2)
    return dos


def quat_Esh_DOS(E, T):
    m_e_tmp = np.zeros(len(temperature))
    e = np.zeros(len(temperature))
    for i in range(len(temperature)):
        m_e_tmp[i] = quaternary_params_with_temperature[i]["m_so"]
        e[i] = energies[i]["Ev_sh"]
    m_e_val = np.interp(T, temperature, m_e_tmp)
    en = np.interp(T, temperature, e)
    tmp = -(E - en)
    if isinstance(tmp, (list, tuple, np.ndarray)):
        tmp[tmp < 0] = 0
    else:
        if tmp < 0:
            tmp = 0
    dos = constDOS * m_e_val ** (3 / 2) * tmp ** (1 / 2)
    return dos


def quat_DOS(E, T):
    return np.array([quat_Ec_DOS(E, T), quat_Elh_DOS(E, T), quat_Ehh_DOS(E, T), quat_Esh_DOS(E, T)])


def DOS(E, T):
    return (2 * GaAs_DOS(E, T) * GaAs_layer1 + quat_DOS(E, T) * AlGaAsSb_layer) / (2 * GaAs_layer1 + AlGaAsSb_layer)


def DOS_Ec(E, T):
    return (2 * GaAs_Ec_DOS(E, T) * GaAs_layer1 + quat_Ec_DOS(E, T) * AlGaAsSb_layer) / (
            2 * GaAs_layer1 + AlGaAsSb_layer)


def DOS_Ehh(E, T):
    return (2 * GaAs_Ehh_DOS(E, T) * GaAs_layer1 + quat_Ehh_DOS(E, T) * AlGaAsSb_layer) / (
            2 * GaAs_layer1 + AlGaAsSb_layer)


def DOS_Elh(E, T):
    return (2 * GaAs_Elh_DOS(E, T) * GaAs_layer1 + quat_Elh_DOS(E, T) * AlGaAsSb_layer) / (
            2 * GaAs_layer1 + AlGaAsSb_layer)


def DOS_Esh(E, T):
    return (2 * GaAs_Esh_DOS(E, T) * GaAs_layer1 + quat_Esh_DOS(E, T) * AlGaAsSb_layer) / (
            2 * GaAs_layer1 + AlGaAsSb_layer)


def FD(E, T, Ef):
    return 1 / (1 + np.exp((E - Ef) / (kb * T)))


def FDholes(E, T, Ef):
    return 1 / (1 + np.exp(-(E - Ef) / (kb * T)))


def gap_min_max(T):
    ec_tmp = np.zeros(len(temperature))
    ev_sh_tmp = np.zeros(len(temperature))
    ev_lh_tmp = np.zeros(len(temperature))
    ev_hh_tmp = np.zeros(len(temperature))
    for i in range(len(temperature)):
        ec_tmp[i] = energies[i]["Ec"]
        ev_lh_tmp[i] = energies[i]["Ev_lh"]
        ev_hh_tmp[i] = energies[i]["Ev_hh"]
        ev_sh_tmp[i] = energies[i]["Ev_sh"]

    ec = np.interp(T, temperature, ec_tmp)
    ev_lh = np.interp(T, temperature, ev_lh_tmp)
    ev_hh = np.interp(T, temperature, ev_hh_tmp)
    ev_sh = np.interp(T, temperature, ev_sh_tmp)

    eg_T_GaAs = Eg["GaAs"] - (alpha["GaAs"] * (T ** 2) / (T + beta["GaAs"]))
    ec_GaAs = VBO["GaAs"] + eg_T_GaAs
    ehh_GaAs = VBO["GaAs"]
    elh_GaAs = VBO["GaAs"]
    eso_GaAs = VBO["GaAs"]
    return [min(ec, ec_GaAs), max(ev_hh, ev_lh, ev_sh, ehh_GaAs, elh_GaAs, eso_GaAs)]


def int_n(efv, Ef_max, T):
    intval = integrate.quad(lambda xx: FD(xx, T, efv) * DOS_Ec(xx, T),
                            Ef_max, Ef_max + 10, limit=70)
    return intval[0]


def int_p(efv, Ef_min, T):
    intval = integrate.quad(lambda xx: (FDholes(xx, T, efv)) * (DOS_Ehh(xx, T) + DOS_Elh(xx, T)),
                            Ef_min - 5, Ef_min, limit=70)
    return intval[0]


def fun(efv, Ef_max, Ef_min, T):
    return int_n(efv, Ef_max, T) - int_p(efv, Ef_min, T)


start = time.time()

# --------------------- SETTINGS ------------------------------------
temperature = np.linspace(0, 650, 300)

# varied = [0.5, 0.1]  # value of x
varied = [0.0, 0.1, 0.2, 0.25, 0.5, 1.0]  # value of x
col = []
for index in range(len(varied)):
    col.append(list(plt.cm.tab10(index)))
fixed = 0.9  # value of y

# structure parameters
GaAs_layer1 = 500
AlGaAsSb_layer = 500  # from 100 to 500 nm
GaAs_layer2 = 500
lb = '(d)'
# --------------------- GaAs prelim --------------------------------------

VBO = {
    "GaAs": par.GaAs["VBO"],
    "AlAs": par.AlAs["VBO"],
    "AlSb": par.AlSb["VBO"],
    "GaSb": par.GaSb["VBO"],
}
Eg = {
    "GaAs": par.GaAs["Eg"],
    "AlAs": par.AlAs["Eg"],
    "AlSb": par.AlSb["Eg"],
    "GaSb": par.GaSb["Eg"],
}
gamma1 = {
    "GaAs": par.GaAs["gamma1"],
    "AlAs": par.AlAs["gamma1"],
    "AlSb": par.AlSb["gamma1"],
    "GaSb": par.GaSb["gamma1"],
}
gamma2 = {
    "GaAs": par.GaAs["gamma2"],
    "AlAs": par.AlAs["gamma2"],
    "AlSb": par.AlSb["gamma2"],
    "GaSb": par.GaSb["gamma2"],
}
m_e = {
    "GaAs": par.GaAs["m_e"],
    "AlAs": par.AlAs["m_e"],
    "AlSb": par.AlSb["m_e"],
    "GaSb": par.GaSb["m_e"],
}
m_hh = {
    "GaAs": par.GaAs["m_hh"],
    "AlAs": par.AlAs["m_hh"],
    "AlSb": par.AlSb["m_hh"],
    "GaSb": par.GaSb["m_hh"],
}
m_lh = {
    "GaAs": par.GaAs["m_lh"],
    "AlAs": par.AlAs["m_lh"],
    "AlSb": par.AlSb["m_lh"],
    "GaSb": par.GaSb["m_lh"],
}
m_so = {
    "GaAs": par.GaAs["m_so"],
    "AlAs": par.AlAs["m_so"],
    "AlSb": par.AlSb["m_so"],
    "GaSb": par.GaSb["m_so"],
}
alpha = {
    "GaAs": par.GaAs["alpha"],
    "AlAs": par.AlAs["alpha"],
    "AlSb": par.AlSb["alpha"],
    "GaSb": par.GaSb["alpha"],
}
beta = {
    "GaAs": par.GaAs["beta"],
    "AlAs": par.AlAs["beta"],
    "AlSb": par.AlSb["beta"],
    "GaSb": par.GaSb["beta"],
}
delta_so = {
    "GaAs": par.GaAs["delta_so"],
    "AlAs": par.AlAs["delta_so"],
    "AlSb": par.AlSb["delta_so"],
    "GaSb": par.GaSb["delta_so"],
}
kb = 8.617333262145 * 10 ** (-5)
const = 0.0380998212
constDOS = 1 / (2 * np.pi ** 2) * (1 / const) ** (3 / 2)
#
# Ec_GaAs = dict()
# Ehh_GaAs = dict()
# Elh_GaAs = dict()
# Eso_GaAs = dict()
# Ec_DOS_GaAs = dict()
# Ehh_DOS_GaAs = dict()
# Elh_DOS_GaAs = dict()
# Eso_DOS_GaAs = dict()
# E_GaAs = dict()
# Eg_T_GaAs = dict()
# kz = 0
# for material in materials:
#     Ec_GaAs[material] = np.zeros([len(temperature)], dtype=float)
#     Ehh_GaAs[material] = np.zeros([len(temperature)], dtype=float)
#     Elh_GaAs[material] = np.zeros([len(temperature)], dtype=float)
#     Eso_GaAs[material] = np.zeros([len(temperature)], dtype=float)
#     Eg_T_GaAs[material] = np.zeros(len(temperature))
#     E_GaAs[material] = np.zeros([len(temperature), 500], dtype=float)
#     for i, T in enumerate(temperature):
#         Eg_T_GaAs[material][i] = Eg[material] - (alpha[material] * (T ** 2) / (T + beta[material]))
#         Ec_GaAs[material][i] = VBO[material] + Eg_T_GaAs[material][i] + const * kz ** 2 * 1 / m_e[material]
#         # Ehh[material][i] = VBO[material] - const * kz ** 2 * (gamma1[material] - 2 * gamma2[material])
#         # Elh[material][i] = VBO[material] - const * kz ** 2 * (gamma1[material] + 2 * gamma2[material])
#         Eso_GaAs[material][i] = VBO[material] - const * kz ** 2 * gamma1[material] - delta_so[material]
#         E_GaAs[material][i, :] = np.linspace(Eso_GaAs[material][i], Ec_GaAs[material][i], 500)

# ------------------------------------------------------------------------------------------------
# --------------------------------  AlGaAsSb DOS -------------------------------------------------


ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]
parameter_keys = ["Eg", "VBO", "delta_so", "gamma1", "gamma2",
                  "m_e", "m_hh", "m_lh", "m_so", "ac", "av", "b", "c11", "c12",
                  "alc_temp"]
#
#

nq = []
pq = []
Efq = []
n_eq = []
Ef_eq = []
Ef_min_q = []
dosq = []
for q, var in enumerate(varied):
    ternary_params = dict()
    for mat in ternary_mat:
        ternary_params[mat] = dict()

    quaternary_params = dict()
    ternary_params_with_temperature = []
    quaternary_params_with_temperature = []
    for k, T in enumerate(temperature):
        for par_key in parameter_keys:

            if par_key == "Eg":
                p_AlAs = par.AlAs[par_key] - (par.AlAs["alpha"] * (T ** 2) / (T + par.AlAs["beta"]))
                p_GaAs = par.GaAs[par_key] - (par.GaAs["alpha"] * (T ** 2) / (T + par.GaAs["beta"]))
                p_AlSb = par.AlSb[par_key] - (par.AlSb["alpha"] * (T ** 2) / (T + par.AlSb["beta"]))
                p_GaSb = par.GaSb[par_key] - (par.GaSb["alpha"] * (T ** 2) / (T + par.GaSb["beta"]))
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

            y = fixed
            x = var
            # x = fixed
            # y = var
            p_AlAsSb = y * p_AlAs + (1 - y) * p_AlSb - y * (1 - y) * b_AlAsSb
            p_GaAsSb = y * p_GaAs + (1 - y) * p_GaSb - y * (1 - y) * b_GaAsSb
            p = x * (1 - x) * (y * ternary_params["AlGaAs"][par_key] + (1 - y) * ternary_params["AlGaSb"][par_key]) \
                + y * (1 - y) * (x * p_AlAsSb + (1 - x) * p_GaAsSb)
            quaternary_params[par_key] = p / (x * (1 - x) + y * (1 - y))

        ternary_params_with_temperature.append(copy.deepcopy(ternary_params))
        quaternary_params_with_temperature.append(copy.deepcopy(quaternary_params))

    deformation_keys = ["eps_par", "eps_orth", "dEc_hydro", "dEv_hydro", "dEv_biax", "dEv_biax_plus", "dEv_biax_minus"]
    deformation_params = []
    #
    energies_keys = ["Ec", "Ev_hh", "Ev_lh", "Ev_sh"]
    energies = []
    energies_no_strain = []
    # GaAs is our substrate
    for i, T in enumerate(temperature):
        deformation_params.append(dict())
        energies.append(dict())
        energies_no_strain.append(dict())

        deformation_params[i]["eps_par"] = \
            (par.GaAs["alc_temp"](T) - quaternary_params_with_temperature[i]["alc_temp"]) / \
            quaternary_params_with_temperature[i]["alc_temp"]

        deformation_params[i]["eps_orth"] = \
            -2.0 * (quaternary_params_with_temperature[i]["c12"] /
                    quaternary_params_with_temperature[i]["c11"]) * \
            deformation_params[i][
                "eps_par"]

        deformation_params[i]["dEc_hydro"] = \
            quaternary_params_with_temperature[i]["ac"] * \
            (deformation_params[i]["eps_orth"] + 2 * deformation_params[i][
                "eps_par"])

        deformation_params[i]["dEv_hydro"] = \
            quaternary_params_with_temperature[i]["av"] * \
            (deformation_params[i]["eps_orth"] + 2 * deformation_params[i][
                "eps_par"])

        deformation_params[i]["dEv_biax"] = \
            quaternary_params_with_temperature[i]["b"] * \
            (deformation_params[i]["eps_orth"] - deformation_params[i][
                "eps_par"])

        ds = quaternary_params_with_temperature[i]["delta_so"]
        dEv_biax = deformation_params[i]["dEv_biax"]
        deformation_params[i]["dEv_biax_plus"] = 0.5 * (dEv_biax - ds + np.sqrt(
            9 * dEv_biax ** 2 + 2 * dEv_biax * ds + ds ** 2))
        deformation_params[i]["dEv_biax_minus"] = 0.5 * (dEv_biax - ds - np.sqrt(
            9 * dEv_biax ** 2 + 2 * dEv_biax * ds + ds ** 2))

        energies[i]["Ec"] = quaternary_params_with_temperature[i]["VBO"] + \
                            quaternary_params_with_temperature[i]["Eg"] + \
                            deformation_params[i]["dEc_hydro"]
        energies_no_strain[i]["Ec"] = quaternary_params_with_temperature[i][
                                          "VBO"] + \
                                      quaternary_params_with_temperature[i][
                                          "Eg"]
        energies[i]["Ev_hh"] = quaternary_params_with_temperature[i]["VBO"] + \
                               deformation_params[i]["dEv_hydro"] - \
                               deformation_params[i]["dEv_biax"]
        energies_no_strain[i]["Ev_hh"] = quaternary_params_with_temperature[i][
            "VBO"]

        energies[i]["Ev_lh"] = quaternary_params_with_temperature[i]["VBO"] + \
                               deformation_params[i]["dEv_hydro"] + \
                               deformation_params[i]["dEv_biax_plus"]
        energies_no_strain[i]["Ev_lh"] = quaternary_params_with_temperature[i][
            "VBO"]

        energies[i]["Ev_sh"] = quaternary_params_with_temperature[i]["VBO"] + \
                               deformation_params[i]["dEv_hydro"] + \
                               deformation_params[i]["dEv_biax_minus"]
        energies_no_strain[i]["Ev_sh"] = quaternary_params_with_temperature[i][
                                             "VBO"] - ds


    # Ec = np.zeros([len(temperature)], dtype=float)
    # Ehh = np.zeros([len(temperature)], dtype=float)
    # Elh = np.zeros([len(temperature)], dtype=float)
    # Eso = np.zeros([len(temperature)], dtype=float)
    # Ec_ns = np.zeros([len(temperature)], dtype=float)
    # Ehh_ns = np.zeros([len(temperature)], dtype=float)
    # Elh_ns = np.zeros([len(temperature)], dtype=float)
    # Eso_ns = np.zeros([len(temperature)], dtype=float)
    # Eg_T = np.zeros(len(temperature))
    # E = np.zeros([len(temperature), 500], dtype=float)
    # kz = 0
    # for i, T in enumerate(temperature):
    #     Ec[i] = energies[i]["Ec"] + \
    #             const * kz ** 2 * 1 / quaternary_params_with_temperature[i]["m_e"]
    #     Ehh[i] = energies[i]["Ev_hh"] - \
    #              const * kz ** 2 * (quaternary_params_with_temperature[i]["gamma1"] -
    #                                 2 * quaternary_params_with_temperature[i]["gamma2"])
    #     Elh[i] = energies[i]["Ev_lh"] - \
    #              const * kz ** 2 * (quaternary_params_with_temperature[i]["gamma1"] +
    #                                 2 * quaternary_params_with_temperature[i]["gamma2"])
    #     Eso[i] = energies[i]["Ev_sh"] - \
    #              const * kz ** 2 * quaternary_params_with_temperature[i]["gamma1"]
    #     Ec_ns[i] = energies_no_strain[i]["Ec"] + \
    #                const * kz ** 2 * 1 / quaternary_params_with_temperature[i]["m_e"]
    #     Ehh_ns[i] = energies_no_strain[i]["Ev_hh"] - \
    #                 const * kz ** 2 * (quaternary_params_with_temperature[i]["gamma1"] -
    #                                    2 * quaternary_params_with_temperature[i]["gamma2"])
    #     Elh_ns[i] = energies_no_strain[i]["Ev_lh"] - \
    #                 const * kz ** 2 * (quaternary_params_with_temperature[i]["gamma1"] +
    #                                    2 * quaternary_params_with_temperature[i]["gamma2"])
    #     Eso_ns[i] = energies_no_strain[i]["Ev_sh"] - \
    #                 const * kz ** 2 * quaternary_params_with_temperature[i]["gamma1"]
    #     E[i, :] = np.linspace(Eso[i], Ec[i], 500)
    #  --------------- Calculation of n and p dependence of Ef ----------------------------------------
    # T = 300
    # [Ef_max, Ef_min] = gap_min_max(T)
    # Ef_values = np.linspace(Ef_min, Ef_max, 200)
    # # Ef_values = np.linspace(-0.2, 0.1, 200)
    # n = []
    # p = []
    # for ef in Ef_values:
    #     intval1 = integrate.quad(lambda x, efv=ef: FD(x, T, efv)*DOS_Ec(x, T),
    #                              Ef_max, Ef_max+10, limit=70)
    #     intval2 = integrate.quad(lambda x, efv=ef: (FDholes(x, T, efv))*(DOS_Ehh(x, T)+DOS_Elh(x, T)),
    #                              Ef_min-5, Ef_min, limit=70)
    #     n.append(intval1[0])
    #     p.append(intval2[0])
    # nq.append(np.array(n))
    # pq.append(np.array(p))
    # Efq.append(Ef_values)
    #  -----------------------------------------------------------------------------------
    # Ef_values = np.linspace(Ef_min-2, Ef_max+1, 300)
    # dosq.append(DOS_Ehh(Ef_values,T)+DOS_Elh(Ef_values,T))
    # dosq.append(DOS_Ec(Ef_values,T))
    # Efq.append(Ef_values)
    #  ------------------ Calculation of equilibrium Fermi levels -------------------------
    T_values = np.linspace(25, 600, 120)
    EfT = np.zeros(len(T_values))
    nT = np.zeros(len(T_values))
    for i, T in enumerate(T_values):
        [Ef_max, Ef_min] = gap_min_max(T)
        Ef_values = np.linspace(Ef_min, Ef_max, 50)
        dE = Ef_values[3] - Ef_values[0]
        n = []
        p = []
        for ef in Ef_values:
            intval1 = integrate.quad(lambda x, efv=ef: FD(x, T, efv) * DOS_Ec(x, T),
                                     Ef_max, Ef_max + 10, limit=70)
            intval2 = integrate.quad(lambda x, efv=ef: (FDholes(x, T, efv)) * (DOS_Ehh(x, T) + DOS_Elh(x, T)),
                                     Ef_min - 5, Ef_min, limit=70)
            n.append(intval1[0])
            p.append(intval2[0])
        narr = np.array(n)
        parr = np.array(p)
        diff_abs = np.abs(narr - parr)
        idx = np.argmin(diff_abs)
        Ef_close_to_zero = Ef_values[idx]
        # Ef_zeros = np.linspace(Ef_close_to_zero-dE, Ef_close_to_zero+dE, 200)
        EfT[i] = bq(
            fun, Ef_close_to_zero - dE, Ef_close_to_zero + dE,
            args=(Ef_max, Ef_min, T))
        nT[i] = int_n(EfT[i], Ef_max, T)
        p_GaAs = par.GaAs["Eg"] - (par.GaAs["alpha"] * (T ** 2) / (T + par.GaAs["beta"]))
        EfT[i] = EfT[i] - p_GaAs/2
    n_eq.append(nT)
    Ef_eq.append(EfT)

#  ------------------ Plot options -----------------------------------------
linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
linewidth = 0.8
fontsize = 14
tick_fontsize = 15
legend_fontsize = 12

# for q in range(len(varied)):
#     plt.plot(Efq[q], dosq[q],label='$x = '+str(varied[q])+'$',color=col[q])
# plt.legend(loc='best',fontsize=legend_fontsize+4)
# plt.xlabel("Energia $[eV]$",fontsize=fontsize+4)
# plt.ylabel("Gęstość stanów",fontsize=fontsize+4)
# plt.xticks(fontsize=tick_fontsize+2)
# plt.yticks(fontsize=tick_fontsize+2)
# plt.xlim([-1.5,-0.5])
# # plt.xlim([0.2,1.4])
# plt.ylim([0.0,3.0])
# # plt.ylim([0.0,0.12])
# ax = plt.gca()
# ax.text(0.1, 0.9, lb,
#         verticalalignment='bottom', horizontalalignment='right',
#         transform=ax.transAxes, fontsize=16)
# plt.tight_layout()
# folder = "report/Figures/carriers/"
# filename = 'dos_val_L_' + str(AlGaAsSb_layer)
# ext = '.pdf'
# plt.tight_layout()
# plt.savefig(folder+filename+ext)
# plt.savefig(folder+filename+'.png')
# plt.show()

#  ------------------ Plots of stationary concentration and Fermi energy as a function of T
fig = plt.figure()
# gs = GridSpec(nrows=1, ncols=2,
#               left=0.13, right=0.98, bottom=0.11, top=0.98,
#               wspace=0.45, hspace=0.05)

gs = GridSpec(nrows=1, ncols=1,
              left=0.15, right=0.98, bottom=0.11, top=0.98,
              wspace=0.01, hspace=0.05)
ax1 = fig.add_subplot(gs[0])
# ax2 = fig.add_subplot(gs[1])
ax1.text(0.15, 0.9, lb,
         verticalalignment='bottom', horizontalalignment='right',
         transform=ax1.transAxes, fontsize=13)
ax1.set_xlabel("Temperatura $[K]$", fontsize=fontsize)
# ax2.set_xlabel("Temperatura $[K]$", fontsize=fontsize)
ax1.set_ylabel("Poziom Fermiego $[eV]$", fontsize=fontsize)
# ax2.set_ylabel("Koncentracja nośników $[cm^{-3}]$", fontsize=fontsize)
# ax1.set_ylim([0.55, 0.78])
ax1.set_ylim([-0.8, -0.63])
# ax2.set_ylim([1e-18, 1e14])
for q in range(len(varied)):
    ax1.plot(T_values, Ef_eq[q],color=col[q],label='$x = ' + str(varied[q]) + '$')
    # ax2.semilogy(T_values, n_eq[q] * 1e21, color=col[q], label='$y = ' + str(varied[q]) + '$')
# plt.legend(loc='lower right', fontsize=legend_fontsize)
folder = "report/Figures/carriers/"
filename = 'fermi_Eg_x_L_' + str(AlGaAsSb_layer)
ext = '.pdf'
plt.tight_layout()
plt.savefig(folder+filename+ext)
plt.savefig(folder+filename+'.png')
plt.show()
#  ------------------ Plots of concentrations dependence on Ef ---------------
# for q in range(len(varied)):
    # plt.plot(Efq[q], nq[q], linestyle='-', color=col[q])
    # plt.semilogy(Efq[q]-Efq[q][0], nq[q]*1e21, label='$x = '+str(varied[q])+'$', linestyle='solid', color=col[q])
    # plt.plot(Efq[q], pq[q], label='$x = '+str(varied[q])+'$', linestyle='--', color=col[q])
    # plt.semilogy(Efq[q]-Efq[q][0], pq[q]*1e21, linestyle='dashed', color=col[q])

# plt.ylim([1e-15, 1e22])
# plt.xlim([-1.8, 1.8])
# ax = plt.gca()
# ax.text(0.1, 0.8, '(d)',
#         verticalalignment='bottom', horizontalalignment='right',
#         transform=ax.transAxes, fontsize=13)
# plt.xlabel('Poziom Fermiego $[eV]$', fontsize=fontsize)
# plt.ylabel('Koncentracja nośników $[cm^{-3}]$', fontsize=fontsize)
# plt.legend(loc='center right', fontsize=legend_fontsize)
#
# Ef = Efq[0]+1
# dos = DOS_Ec(Ef,T)
# fd = FD(Ef, T, Ef[100])
#
# fig, ax = plt.subplots()
# ax.plot(Ef, dos, label='$D_{c}(E)$')
# ax.plot(Ef, fd, label='$f(E)$')
# ax.plot(Ef, dos*fd, label='$D_{c}\cdot f$')
# plt.fill_between(Ef, dos*fd, color='tab:green', alpha=.2)
# plt.ylim([0,0.2])
# plt.xlim([0.2,1.3])
# plt.legend(loc='best', fontsize=legend_fontsize)
# # ax.set(yticklabels=[])
# plt.xlabel("Energia $[eV]$",fontsize=fontsize)
# plt.tight_layout()
# folder = "report/Figures/carriers/"
# filename = 'integrand'
# ext = '.pdf'
# plt.savefig(folder+filename+ext)
# plt.savefig(folder+filename+'.png')
# plt.show()
#
# folder = "report/Figures/carriers/"
# filename = 'concentration_L_' + str(AlGaAsSb_layer)
# ext = '.pdf'
# plt.tight_layout()
# plt.savefig(folder+filename+ext)
# plt.savefig(folder+filename+'.png')
# plt.show()
#  ----------------------------------------------------------------
print('Program finished in ', time.time() - start, ' s.')
