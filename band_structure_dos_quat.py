from brokenaxes import brokenaxes
from data import ternary_alloys_bowings as bow, binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec

start = time.time()

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
const = 0.0380998212
constDOS = 1 / (2 * np.pi ** 2) * (1 / const) ** (3 / 2)

ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]
parameter_keys = ["Eg", "VBO", "delta_so", "gamma1", "gamma2",
                  "m_e", "m_hh", "m_lh", "m_so", "ac", "av", "b", "c11", "c12",
                  "alc_temp"]
#
var = 0.9
fixed = np.array([0.5])
txt = '(d)'
temperature = np.array([0, 150, 300, 450])
# temperature = np.array([300])
kz_values = np.linspace(-1, 1, 200)

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
    deformation_params[i] = {"fix_x": dict(), "fix_y": dict()}
    energies.append(dict())
    energies_no_strain.append(dict())
    energies[i] = {"fix_x": dict(), "fix_y": dict()}
    energies_no_strain[i] = {"fix_x": dict(), "fix_y": dict()}

    for fix_arg in ["fix_x", "fix_y"]:
        for j, fix in enumerate(fixed):
            fix_key = '{0:1.1f}'.format(fix)
            deformation_params[i][fix_arg][fix_key] = dict()
            energies[i][fix_arg][fix_key] = dict()
            energies_no_strain[i][fix_arg][fix_key] = dict()

            deformation_params[i][fix_arg][fix_key]["eps_par"] = \
                (par.GaAs["alc_temp"](T) - quaternary_params_with_temperature[i][fix_arg][fix_key]["alc_temp"]) / \
                quaternary_params_with_temperature[i][fix_arg][fix_key]["alc_temp"]

            deformation_params[i][fix_arg][fix_key]["eps_orth"] = \
                -2.0 * (quaternary_params_with_temperature[i][fix_arg][fix_key]["c12"] /
                        quaternary_params_with_temperature[i][fix_arg][fix_key]["c11"]) * \
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
            deformation_params[i][fix_arg][fix_key]["dEv_biax_plus"] = 0.5 * (dEv_biax - delta_so + np.sqrt(
                9 * dEv_biax ** 2 + 2 * dEv_biax * delta_so + delta_so ** 2))
            deformation_params[i][fix_arg][fix_key]["dEv_biax_minus"] = 0.5 * (dEv_biax - delta_so - np.sqrt(
                9 * dEv_biax ** 2 + 2 * dEv_biax * delta_so + delta_so ** 2))

            energies[i][fix_arg][fix_key]["Ec"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                  quaternary_params_with_temperature[i][fix_arg][fix_key]["Eg"] + \
                                                  deformation_params[i][fix_arg][fix_key]["dEc_hydro"]
            energies_no_strain[i][fix_arg][fix_key]["Ec"] = quaternary_params_with_temperature[i][fix_arg][fix_key][
                                                                "VBO"] + \
                                                            quaternary_params_with_temperature[i][fix_arg][fix_key][
                                                                "Eg"]
            energies[i][fix_arg][fix_key]["Ev_hh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] - \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax"]
            energies_no_strain[i][fix_arg][fix_key]["Ev_hh"] = quaternary_params_with_temperature[i][fix_arg][fix_key][
                "VBO"]

            energies[i][fix_arg][fix_key]["Ev_lh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax_plus"]
            energies_no_strain[i][fix_arg][fix_key]["Ev_lh"] = quaternary_params_with_temperature[i][fix_arg][fix_key][
                "VBO"]

            energies[i][fix_arg][fix_key]["Ev_sh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax_minus"]
            energies_no_strain[i][fix_arg][fix_key]["Ev_sh"] = quaternary_params_with_temperature[i][fix_arg][fix_key][
                                                                   "VBO"] - delta_so

Ec = np.zeros([len(temperature), len(kz_values)], dtype=float)
Ehh = np.zeros([len(temperature), len(kz_values)], dtype=float)
Elh = np.zeros([len(temperature), len(kz_values)], dtype=float)
Eso = np.zeros([len(temperature), len(kz_values)], dtype=float)
Ec_ns = np.zeros([len(temperature), len(kz_values)], dtype=float)
Ehh_ns = np.zeros([len(temperature), len(kz_values)], dtype=float)
Elh_ns = np.zeros([len(temperature), len(kz_values)], dtype=float)
Eso_ns = np.zeros([len(temperature), len(kz_values)], dtype=float)
Eg_T = np.zeros(len(temperature))
E = np.zeros([len(temperature), 200], dtype=float)

filename = ''
path = 'report/Figures/band_str/'
linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
linewidth = 0.8
fontsize = 14
tick_fontsize = 15
legend_fontsize = 12
leg_no_T = [r'$E_c$', r'$E_{v,hh}$', r'$E_{v,lh}$', r'$E_{v,sh}$']
for fix_arg in ["fix_x", "fix_y"]:
    for n, fix in enumerate(fixed):
        fix_key = '{0:1.1f}'.format(fix)
        # ylims = ((np.min(Eso) - 0.15, np.max(Ehh) + 0.15),
        #          (np.min(Ec_ns) - 0.15, np.max(Ec) + 0.15))
        ylims = ((-1.25,0.2),(1.9,2.55))
        fig = plt.figure()
        sps = GridSpec(nrows=1, ncols=2, width_ratios=(7, 4),
                       left=0.13, right=0.95, bottom=0.13, top=0.9,
                       wspace=0.1, hspace=0.05)
        ax = brokenaxes(ylims=ylims, hspace=.1, subplot_spec=sps[0], d=0.01)
        ax_DOS = brokenaxes(ylims=ylims, hspace=.1, subplot_spec=sps[1], d=0.01)
        ax_DOS.tick_params(axis="y", labelleft=False)
        for i, T in enumerate(temperature):
            for j, kz in enumerate(kz_values):
                Ec[i, j] = energies[i][fix_arg][fix_key]["Ec"] + \
                           const * kz ** 2 * 1 / quaternary_params_with_temperature[i][fix_arg][fix_key]["m_e"]
                Ehh[i, j] = energies[i][fix_arg][fix_key]["Ev_hh"] - \
                            const * kz ** 2 * (quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma1"] -
                                               2 * quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma2"])
                Elh[i, j] = energies[i][fix_arg][fix_key]["Ev_lh"] - \
                            const * kz ** 2 * (quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma1"] +
                                               2 * quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma2"])
                Eso[i, j] = energies[i][fix_arg][fix_key]["Ev_sh"] - \
                            const * kz ** 2 * quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma1"]
                Ec_ns[i, j] = energies_no_strain[i][fix_arg][fix_key]["Ec"] + \
                              const * kz ** 2 * 1 / quaternary_params_with_temperature[i][fix_arg][fix_key]["m_e"]
                Ehh_ns[i, j] = energies_no_strain[i][fix_arg][fix_key]["Ev_hh"] - \
                               const * kz ** 2 * (quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma1"] -
                                                  2 * quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma2"])
                Elh_ns[i, j] = energies_no_strain[i][fix_arg][fix_key]["Ev_lh"] - \
                               const * kz ** 2 * (quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma1"] +
                                                  2 * quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma2"])
                Eso_ns[i, j] = energies_no_strain[i][fix_arg][fix_key]["Ev_sh"] - \
                               const * kz ** 2 * quaternary_params_with_temperature[i][fix_arg][fix_key]["gamma1"]
            E[i, :] = np.linspace(np.min(Eso[i, :]), np.max(Ec[i, :]), 200)

        Ec_DOS = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        Ehh_DOS = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        Elh_DOS = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        Eso_DOS = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        Ec_DOSns = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        Ehh_DOSns = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        Elh_DOSns = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        Eso_DOSns = np.zeros([len(temperature), len(E[0, :])], dtype=float)
        for i, T in enumerate(temperature):
            for j, E_val in enumerate(E[i, :]):
                Ec_DOS[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_e"] ** (3 / 2) * \
                               (E_val - energies[i][fix_arg][fix_key]["Ec"]) ** (1 / 2)
                Ehh_DOS[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_hh"] ** (3 / 2) * \
                                (-(E_val - energies[i][fix_arg][fix_key]["Ev_hh"])) ** (1 / 2)
                Elh_DOS[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_lh"] ** (3 / 2) * \
                                (-(E_val - energies[i][fix_arg][fix_key]["Ev_lh"])) ** (1 / 2)
                Eso_DOS[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_so"] ** (3 / 2) * \
                                (-(E_val - energies[i][fix_arg][fix_key]["Ev_sh"])) ** (1 / 2)

                Ec_DOSns[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_e"] ** (3 / 2) * \
                                 (E_val - energies_no_strain[i][fix_arg][fix_key]["Ec"]) ** (1 / 2)
                Ehh_DOSns[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_hh"] ** (
                            3 / 2) * \
                                  (-(E_val - energies_no_strain[i][fix_arg][fix_key]["Ev_hh"])) ** (1 / 2)
                Elh_DOSns[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_lh"] ** (
                            3 / 2) * \
                                  (-(E_val - energies_no_strain[i][fix_arg][fix_key]["Ev_lh"])) ** (1 / 2)
                Eso_DOSns[i, j] = constDOS * quaternary_params_with_temperature[i][fix_arg][fix_key]["m_so"] ** (
                            3 / 2) * \
                                  (-(E_val - energies_no_strain[i][fix_arg][fix_key]["Ev_sh"])) ** (1 / 2)

        # ============ Temperature dependent strain ==========================================================
        colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728']
        for i in range(len(temperature)):
            ax.plot(kz_values, Ec[i, :], color=colors[0], linestyle=linestyles[i])
            ax_DOS.plot(Ec_DOS[i, :], E[i], color=colors[0], linestyle=linestyles[i])
            ax.plot(kz_values, Ehh[i, :], color=colors[1], linestyle=linestyles[i])
            ax.plot(kz_values, Elh[i, :], color=colors[2], linestyle=linestyles[i])
            ax.plot(kz_values, Eso[i, :], color=colors[3], linestyle=linestyles[i])
            ax_DOS.plot(Ehh_DOS[i, :], E[i], color=colors[1], linestyle=linestyles[i])
            ax_DOS.plot(Elh_DOS[i, :], E[i], color=colors[2], linestyle=linestyles[i])
            ax_DOS.plot(Eso_DOS[i, :], E[i], color=colors[3], linestyle=linestyles[i])
        patches = []
        lines = []
        for col_ind in range(len(temperature)):  # len(temperature)
            patches.append(mpatches.Patch(color=colors[col_ind], label=leg_no_T[col_ind]))
            label = "T K".replace("T", str(temperature[col_ind]))
            lines.append(mlines.Line2D([], [], color='black', linestyle=linestyles[col_ind],
                                       label=label))

        hh = [patches, lines]
        flat_list = sum(hh, [])
        handles, labels = ax_DOS.get_legend_handles_labels()
        # fig.legend(handles=flat_list, loc='upper center', ncol=4, fontsize=legend_fontsize)

        filename = ''
        if fix_arg == "fix_x":
            filename = "x_" + str(fix) + "_y_" + str(var)
        else:
            filename = "x_" + str(var) + "_y_" + str(fix)
        # ax.set_ylabel(r'Energia [eV]')
        ax_DOS.set_xlim(left=0, right=None)
        fig.text(0.2, 0.83, txt, va='center')
        fig.text(0.0, 0.37, r'Energia [eV]', fontsize=fontsize, rotation='vertical')
        fig.text(0.32, 0.06, r'$k_z\;[1/nm]$', va='center',fontsize=fontsize)
        fig.text(0.73, 0.06, r'$DOS\;[stany/eV]$', va='center',fontsize=legend_fontsize)
        plt.savefig(path + filename + '.png')
        plt.savefig(path + filename + '.pdf')
        plt.clf()
        plt.gca().set_prop_cycle(None)
        # ============================================================================================
        #========================= Strain - no strain ================================================
        # colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728']
        #
        # ax.plot(kz_values, Ec[2, :], color=colors[0], linestyle=linestyles[1])
        # ax_DOS.plot(Ec_DOS[2, :], E[2], color=colors[0], linestyle=linestyles[1])
        # ax.plot(kz_values, Ehh[2, :], color=colors[1], linestyle=linestyles[1])
        # ax.plot(kz_values, Elh[2, :], color=colors[2], linestyle=linestyles[1])
        # ax.plot(kz_values, Eso[2, :], color=colors[3], linestyle=linestyles[1])
        # ax_DOS.plot(Ehh_DOS[2, :], E[2], color=colors[1], linestyle=linestyles[1])
        # ax_DOS.plot(Elh_DOS[2, :], E[2], color=colors[2], linestyle=linestyles[1])
        # ax_DOS.plot(Eso_DOS[2, :], E[2], color=colors[3], linestyle=linestyles[1])
        #
        # ax.plot(kz_values, Ec_ns[2, :], color=colors[0], linestyle=linestyles[0])
        # ax_DOS.plot(Ec_DOSns[2, :], E[2], color=colors[0], linestyle=linestyles[0])
        # ax.plot(kz_values, Ehh_ns[2, :], color=colors[1], linestyle=linestyles[0])
        # ax.plot(kz_values, Elh_ns[2, :], color=colors[2], linestyle=linestyles[0])
        # ax.plot(kz_values, Eso_ns[2, :], color=colors[3], linestyle=linestyles[0])
        # ax_DOS.plot(Ehh_DOSns[2, :], E[2], color=colors[1], linestyle=linestyles[0])
        # ax_DOS.plot(Elh_DOSns[2, :], E[2], color=colors[2], linestyle=linestyles[0])
        # ax_DOS.plot(Eso_DOSns[2, :], E[2], color=colors[3], linestyle=linestyles[0])
        #
        # patches = []
        # lines = []
        # for col_ind in range(4):  # len(temperature)
        #     patches.append(mpatches.Patch(color=colors[col_ind], label=leg_no_T[col_ind]))
        # label1 = "z odkształceniami"
        # lines.append(mlines.Line2D([], [], color='black', linestyle=linestyles[2],
        #                            label=label1))
        # label2 = "bez odkształceń"
        # lines.append(mlines.Line2D([], [], color='black', linestyle=linestyles[0],
        #                            label=label2))
        #
        # hh = [patches, lines]
        # flat_list = sum(hh, [])
        # handles, labels = ax_DOS.get_legend_handles_labels()
        # # fig.legend(handles=flat_list, loc='upper center', ncol=3, fontsize=legend_fontsize)
        #
        # filename = ''
        # if fix_arg == "fix_x":
        #     filename = "ns_x_" + str(fix) + "_y_" + str(var)
        # else:
        #     filename = "ns_x_" + str(var) + "_y_" + str(fix)
        # ax_DOS.set_xlim(left=0, right=None)
        # fig.text(0.2, 0.83, txt, va='center')
        # fig.text(0.0, 0.37, r'Energia [eV]', fontsize=fontsize, rotation='vertical')
        # fig.text(0.32, 0.06, r'$k_z\;[1/nm]$', va='center',fontsize=fontsize)
        # fig.text(0.73, 0.06, r'$DOS\;[stany/eV]$', va='center',fontsize=legend_fontsize)
        # plt.savefig(path + filename + '.png')
        # plt.savefig(path + filename + '.pdf')
        # plt.clf()
        # plt.gca().set_prop_cycle(None)
print('Program finished in ', time.time() - start, ' s.')
