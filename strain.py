from data import ternary_alloys_bowings as bow, binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time
import copy

start = time.time()
ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]
temperature = np.linspace(0.0, 500.0, 1000)

parameter_keys = ["Eg", "VBO", "delta_so", "alc_temp", "alpha", "beta", "ac", "av", "b", "c11", "c12"]
#
var = np.array([0.1, 0.3, 0.5])
fixed = np.array([0.1, 0.3, 0.5])
series_num = '1'
#
# var = np.array([0.2, 0.4, 0.6])
# fixed = np.array([0.2, 0.4, 0.6])
# series_num = '2'

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
# GaAs is our substrate
for i, T in enumerate(temperature):
    deformation_params.append(dict())
    deformation_params[i] = {"fix_x": dict(), "fix_y": dict()}
    energies.append(dict())
    energies[i] = {"fix_x": dict(), "fix_y": dict()}

    for fix_arg in ["fix_x", "fix_y"]:
        for j, fix in enumerate(fixed):
            fix_key = '{0:1.1f}'.format(fix)
            deformation_params[i][fix_arg][fix_key] = dict()
            energies[i][fix_arg][fix_key] = dict()

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
            energies[i][fix_arg][fix_key]["Ev_hh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] - \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax"]
            energies[i][fix_arg][fix_key]["Ev_lh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax_plus"]
            energies[i][fix_arg][fix_key]["Ev_sh"] = quaternary_params_with_temperature[i][fix_arg][fix_key]["VBO"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_hydro"] + \
                                                     deformation_params[i][fix_arg][fix_key]["dEv_biax_minus"]

fontsize = 13
tick_fontsize = 17
legend_fontsize = 10
labels = [r'$Al_{x}Ga_{1-x}As$', r'$Al_{x}Ga_{1-x}Sb$', r'$AlAs_{x}Sb_{1-x}$', r'$GaAs_{x}Sb_{1-x}$']
folder = 'report/Figures/strain/'
ext = '.pdf'

filenames = ['ec', 'ev_hh', 'ev_lh', 'ev_sh']
filenames = [folder + i + ext for i in filenames]
xlabel = 'Temperatura (K)'


filenames_x = [folder + i + '_x' + ext for i in filenames]
filenames_y = [folder + i + '_y' + ext for i in filenames]
filenames_fix = {"fix_x": filenames_x, "fix_y": filenames_y}
#
names = ["eg", "alc"]
labs = [r'Przerwa wzbroniona (eV)', r'Parametr sieci ($\AA$)']
# for k, par_key in enumerate(["Eg", "alc_temp"]):
#     for i, mat in enumerate(ternary_mat):
#         for j, v in enumerate(var):
#             ydata = []
#             for t in range(len(temperature)):
#                 pp = ternary_params_with_temperature[t][mat][par_key][j]
#                 ydata.append(pp)
#             label = labels[i].replace("1-x", str(1-v), 1)
#             label = label.replace("x", str(v), )
#             plt.plot(temperature, ydata, label=label)
#         plt.gca().set_ylabel(labs[k], fontsize=fontsize)
#         plt.gca().set_xlabel(xlabel, fontsize=fontsize)
#         if par_key == "alc_temp":
#             plt.legend(loc='best', fontsize=legend_fontsize)
#         plt.yticks(fontsize=tick_fontsize)
#         plt.xticks(fontsize=tick_fontsize)
#         plt.tight_layout()
#         # plt.xlim([0.0, 1.0])
#     plt.savefig(folder + r'ter_'+names[k] + series_num + ext)
#     plt.gca().clear()
#
for k, par_key in enumerate(["Eg", "alc_temp"]):
    for fix in fixed:
        key = '{0:1.1f}'.format(fix)
        for j, v in enumerate(var):
            ydata = []
            for t in range(len(temperature)):
                pp = quaternary_params_with_temperature[t]["fix_x"][key][par_key][j]
                ydata.append(pp)
            # label = r'$Al_{' + '%.1f' % fix + r'}Ga_{' + '%.1f' % (1 - fix) + '}As_{' + '%.1f' % v + \
            #         r'}Sb_{' + '%.1f' % (1 - v) + '}$'
            label = r'x = '+'%.1f' % fix + ', y = '+'%.1f' % v
            plt.plot(temperature, ydata, label=label)
        plt.gca().set_ylabel(labs[k], fontsize=fontsize)
        plt.gca().set_xlabel(xlabel, fontsize=fontsize)
        if par_key == "alc_temp":
            plt.legend(loc='best', fontsize=legend_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.tight_layout()
        # plt.xlim([0.0, 1.0])
    plt.savefig(folder + names[k] + series_num + ext)
    plt.gca().clear()

names = ["eps_par", "eps_orth"]
labs = ["Odkształcenie planarne (%)", "Odkształcenie prostopadłe (%)"]
for k, par_key in enumerate(["eps_par", "eps_orth"]):
    for fix in fixed:
        key = '{0:1.1f}'.format(fix)
        for j, v in enumerate(var):
            ydata = []
            for t in range(len(temperature)):
                pp = deformation_params[t]["fix_x"][key][par_key][j]*100
                ydata.append(pp)
            # label = r'$Al_{' + '%.1f' % fix + r'}Ga_{' + '%.1f' % (1 - fix) + '}As_{' + '%.1f' % v + \
            #         r'}Sb_{' + '%.1f' % (1 - v) + '}$'
            label = r'x = ' + '%.1f' % fix + ', y = ' + '%.1f' % v
            plt.plot(temperature, ydata, label=label)
        plt.gca().set_ylabel(labs[k], fontsize=fontsize)
        plt.gca().set_xlabel(xlabel, fontsize=fontsize)
        # plt.legend(loc='best', fontsize=legend_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.tight_layout()
        # plt.xlim([0.0, 1.0])
    plt.savefig(folder + names[k] + series_num + ext)
    plt.gca().clear()

names = energies_keys
labs = [r'Pasmo przewodnictwa (eV)', r'Pasma walencyjne ciężkodziurowe (eV)',
        r'Pasmo walencyjne lekkodziurowe (eV)',
        r'Pasmo walencyjne rozszczepione (eV)']
for k, par_key in enumerate(names):
    for fix in fixed:
        key = '{0:1.1f}'.format(fix)
        for j, v in enumerate(var):
            ydata = []
            for t in range(len(temperature)):
                pp = energies[t]["fix_x"][key][par_key][j]
                ydata.append(pp)
            # label = r'$Al_{' + '%.1f' % fix + r'}Ga_{' + '%.1f' % (1 - fix) + '}As_{' + '%.1f' % v + \
            #         r'}Sb_{' + '%.1f' % (1 - v) + '}$'
            label = r'x = ' + '%.1f' % fix + ', y = ' + '%.1f' % v
            plt.plot(temperature, ydata, label=label)
        plt.gca().set_ylabel(labs[k], fontsize=fontsize)
        plt.gca().set_xlabel(xlabel, fontsize=fontsize)
        if par_key == "Ev_hh":
            plt.legend(loc='best', fontsize=legend_fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.tight_layout()
        # plt.xlim([0.0, 1.0])
    plt.savefig(folder + names[k] + series_num + ext)
    plt.gca().clear()
#
# names = energies_keys
# labs = r'Energia (eV)'
# leg = [r"$E_c$", r'$E_{v,hh}$', r'$E_{v,lh}$', r'$E_{v,sh}$']
# for k, par_key in enumerate(names):
#     for fix in [0.4]:
#         key = '{0:1.1f}'.format(fix)
#         for j, v in enumerate([0.4]):
#             ydata = []
#             for t in range(len(temperature)):
#                 pp = energies[t]["fix_x"][key][par_key][j]
#                 ydata.append(pp)
#             plt.plot(temperature, ydata, label=leg[k])
#     plt.gca().set_ylabel(labs, fontsize=fontsize)
#     plt.gca().set_xlabel(xlabel, fontsize=fontsize)
#     plt.legend(loc='best', fontsize=legend_fontsize)
#     plt.yticks(fontsize=tick_fontsize)
#     plt.xticks(fontsize=tick_fontsize)
#     plt.tight_layout()
#         # plt.xlim([0.0, 1.0])
# plt.savefig(folder + 'one' + ext)
# plt.gca().clear()

print('Program finished in ', time.time() - start, ' s.')
