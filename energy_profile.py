from sys import argv

from data import ternary_alloys_bowings as bow, binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec
from brokenaxes import brokenaxes
import sys

start = time.time()

# --------------------- SETTINGS ------------------------------------
temperature = np.array([0, 150, 300, 450])
var = np.array([0.5, 0.9])
fixed = np.array([0.0, 0.1, 0.2, 0.5])
# var = np.array([float(sys.argv[1])])
# fixed = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
titles = ['x = 0.1 y = 0.5', 'x = 0.0 y = 0.9', 'x = 0.2 y = 0.9', 'x = 0.5 y = 0.9']
boundary_AB = 225.0
boundary_BA = 275.0
xx = np.arange(0, 500, 0.1)
leg = [r'$E_c$ $T K$', r'$E_{v,hh}$ $T K$', r'$E_{v,lh}$ $T K$', r'$E_{v,sh}$ $T K$']
leg_no_T = [r'$E_c$', r'$E_{v,hh}$', r'$E_{v,lh}$', r'$E_{v,sh}$']

cmap = cm.get_cmap('viridis', 4)
colors = cmap(np.linspace(0, 1, 4))
linestyles = ['solid', 'dotted', 'dashed', 'dashdot']
linewidth = 0.8
fontsize = 11
tick_fontsize = 15
legend_fontsize = 10

folder = 'report/Figures/structure/'
ext = '.pdf'
# -------------------------------------------------------------------

ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]
parameter_keys = ["Eg", "VBO", "delta_so", "alc_temp", "alpha", "beta", "ac", "av", "b", "c11", "c12"]

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
            energies_no_strain[i][fix_arg][fix_key]["Ec"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                  quaternary_params_with_temperature[i][fix_arg][fix_key]["Eg"]

            energies[i][fix_arg][fix_key]["Ev_hh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] - \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax"]
            energies_no_strain[i][fix_arg][fix_key]["Ev_hh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"]

            energies[i][fix_arg][fix_key]["Ev_lh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax_plus"]
            energies_no_strain[i][fix_arg][fix_key]["Ev_lh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"]

            energies[i][fix_arg][fix_key]["Ev_sh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax_minus"]
            energies_no_strain[i][fix_arg][fix_key]["Ev_sh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] - delta_so

Ec = np.zeros((len(temperature), len(xx)))
Ev_hh = np.zeros((len(temperature), len(xx)))
Ev_lh = np.zeros((len(temperature), len(xx)))
Ev_sh = np.zeros((len(temperature), len(xx)))
Ec_ns = np.zeros((len(temperature), len(xx)))
Ev_hh_ns = np.zeros((len(temperature), len(xx)))
Ev_lh_ns = np.zeros((len(temperature), len(xx)))
Ev_sh_ns = np.zeros((len(temperature), len(xx)))

VBO_GaAS = par.GaAs["VBO"]

fig = plt.figure()
sps = GridSpec(nrows=2, ncols=2, hspace=0.2, wspace=0.25)
spss = [sps[0, 0], sps[0, 1], sps[1, 0], sps[1, 1]]
spsi = 0
filename = ''
baxes = []
ylims = [((-1.4, 0.1), (0.4, 1.2)), ((-1.25, -0.55), (0.5, 0.9)),
         ((-1.2, -0.7), (0.54, 0.9)), ((-1.5, -0.7), (0.4, 1.3))]
textt = ['(a)', '(b)', '(c)', '(d)']
textx = 50
texty = [1.0, 0.8, 0.9, 1.0]
for var_i in range(2):
    for fix_arg in ["fix_x", "fix_y"]:
        for j, fix in enumerate(fixed):
            fix_key = '{0:1.1f}'.format(fix)
            # textt = ''
            if fix_arg == 'fix_x':
                # textt = r'$Al_{' + '%.1f' % var[var_i] + r'}Ga_{' + '%.1f' % (
                #         1 - var[var_i]) + '}As_{' + '%.1f' % fix + '}Sb_{' + '%.1f' % (1 - fix) + '}$'
                title = 'x = ' + str(fix) + ' y = ' + str(var[var_i])
            else:
                # textt = r'$Al_{' + '%.1f' % fix + r'}Ga_{' + '%.1f' % (
                #         1 - fix) + '}As_{' + '%.1f' % var[var_i] + '}Sb_{' + '%.1f' % (1 - var[var_i]) + '}$'
                title = 'x = ' + str(var[var_i]) + ' y = ' + str(fix)
            if title not in titles:
                continue
            else:
                print(title)
            for k, T in enumerate(temperature):
                Eg_GaAs = par.GaAs["Eg"] - (par.GaAs["alpha"] * (T ** 2) / (T + par.GaAs["beta"]))
                Ec_quat = energies[k][fix_arg][fix_key]["Ec"][var_i]
                Ev_hh_quat = energies[k][fix_arg][fix_key]["Ev_hh"][var_i]
                Ev_lh_quat = energies[k][fix_arg][fix_key]["Ev_lh"][var_i]
                Ev_sh_quat = energies[k][fix_arg][fix_key]["Ev_sh"][var_i]

                Ec_quat_ns = energies_no_strain[k][fix_arg][fix_key]["Ec"][var_i]
                Ev_hh_quat_ns = energies_no_strain[k][fix_arg][fix_key]["Ev_hh"][var_i]
                Ev_lh_quat_ns = energies_no_strain[k][fix_arg][fix_key]["Ev_lh"][var_i]
                Ev_sh_quat_ns = energies_no_strain[k][fix_arg][fix_key]["Ev_sh"][var_i]
                for i, x in enumerate(xx):
                    if x < boundary_AB or x > boundary_BA:  # in GaAs
                        Ec[k, i] = Eg_GaAs + VBO_GaAS
                        Ev_hh[k, i] = VBO_GaAS
                        Ev_lh[k, i] = VBO_GaAS
                        Ev_sh[k, i] = VBO_GaAS - par.GaAs["delta_so"]
                        Ec_ns[k, i] = Eg_GaAs + VBO_GaAS
                        Ev_hh_ns[k, i] = VBO_GaAS
                        Ev_lh_ns[k, i] = VBO_GaAS
                        Ev_sh_ns[k, i] = VBO_GaAS - par.GaAs["delta_so"]
                    elif boundary_AB < x < boundary_BA:
                        Ec[k, i] = Ec_quat
                        Ev_hh[k, i] = Ev_hh_quat
                        Ev_lh[k, i] = Ev_lh_quat
                        Ev_sh[k, i] = Ev_sh_quat
                        Ec_ns[k, i] = Ec_quat_ns
                        Ev_hh_ns[k, i] = Ev_hh_quat_ns
                        Ev_lh_ns[k, i] = Ev_lh_quat_ns
                        Ev_sh_ns[k, i] = Ev_sh_quat_ns

                    elif x == boundary_AB or x == boundary_BA:
                        Ec[k, i] = (Ec_quat + Eg_GaAs + VBO_GaAS) / 2
                        Ev_hh[k, i] = (Ev_hh_quat + VBO_GaAS) / 2
                        Ev_lh[k, i] = (Ev_lh_quat + VBO_GaAS) / 2
                        Ev_sh[k, i] = (Ev_sh_quat + VBO_GaAS - par.GaAs["delta_so"]) / 2
                        Ec_ns[k, i] = (Ec_quat_ns + Eg_GaAs + VBO_GaAS) / 2
                        Ev_hh_ns[k, i] = (Ev_hh_quat_ns + VBO_GaAS) / 2
                        Ev_lh_ns[k, i] = (Ev_lh_quat_ns + VBO_GaAS) / 2
                        Ev_sh_ns[k, i] = (Ev_sh_quat_ns + VBO_GaAS - par.GaAs["delta_so"]) / 2
                    else:
                        print("Error")

            baxes.append(brokenaxes(ylims=ylims[spsi], hspace=.1, subplot_spec=spss[spsi], d=0.01))

            for k in range(len(temperature)):
                baxes[spsi].plot(xx, Ec[k, :], color=colors[0], linestyle=linestyles[k], linewidth=linewidth)
                baxes[spsi].plot(xx, Ev_hh[k, :], color=colors[1], linestyle=linestyles[k], linewidth=linewidth)
                baxes[spsi].plot(xx, Ev_lh[k, :], color=colors[2], linestyle=linestyles[k], linewidth=linewidth)
                baxes[spsi].plot(xx, Ev_sh[k, :], color=colors[3], linestyle=linestyles[k], linewidth=linewidth)
                # baxes[spsi].plot(xx, Ec_ns[k, :], color=colors[0], linestyle=linestyles[0], linewidth=linewidth)
                # baxes[spsi].plot(xx, Ev_hh_ns[k, :], color=colors[1], linestyle=linestyles[0], linewidth=linewidth)
                # baxes[spsi].plot(xx, Ev_lh_ns[k, :], color=colors[2], linestyle=linestyles[0], linewidth=linewidth)
                # baxes[spsi].plot(xx, Ev_sh_ns[k, :], color=colors[3], linestyle=linestyles[0], linewidth=linewidth)
            if spsi == 0:
                baxes[spsi].annotate(text='50 nm', xy=(boundary_AB - 5.0, 0.0), xytext=(boundary_BA + 11.0, -0.04),
                                     arrowprops={'arrowstyle': '<->', 'shrinkA': 0, 'shrinkB': 0,
                                                 'mutation_scale': 5})

            spsi += 1
            patches = []
            lines = []
            for col_ind in range(4): #len(temperature)
                patches.append(mpatches.Patch(color=colors[col_ind], label=leg_no_T[col_ind]))
                label = "T K".replace("T", str(temperature[col_ind]))
                lines.append(mlines.Line2D([], [], color='black', linestyle=linestyles[col_ind],
                                           label=label))
            # label1 = "z odkształceniami"
            # lines.append(mlines.Line2D([], [], color='black', linestyle=linestyles[2],
            #                            label=label1))
            # label2 = "bez odkształceń"
            # lines.append(mlines.Line2D([], [], color='black', linestyle=linestyles[0],
            #                            label=label2))
            if spsi == 4:
                hh = [patches, lines]
                flat_list = sum(hh, [])
                handles, labels = baxes[spsi - 1].get_legend_handles_labels()
                fig.legend(handles=flat_list, loc='upper center', ncol=4, fontsize=legend_fontsize)


filename = folder + 'temp_strain' + ext
fig.text(0.5, 0.03, 'Grubość warstwy [nm]', ha='center')
fig.text(0.02, 0.5, 'Energia [eV]', va='center', rotation='vertical')
plt.text(0.1, 0.9, textt[3], ha='center', fontsize=7)
plt.text(0.1, 2.1, textt[1], ha='center', fontsize=7)
plt.text(-1.15, 0.9, textt[2], ha='center', fontsize=7)
plt.text(-1.15, 2.1, textt[0], ha='center', fontsize=7)
plt.savefig(filename)
plt.show()

print('Program finished in ', time.time() - start, ' s.')
