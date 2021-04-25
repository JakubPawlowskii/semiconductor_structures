from data import ternary_alloys_bowings as bow, binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
from scipy import optimize


def fun(hc, b, ni, f):
    return b / (2 * np.pi * f) * ((1 - 0.25 * ni) / (1 + ni)) * (np.log(hc / b) + 1) - hc


start = time.time()

legend_fontsize = 12
fontsize = 13
tick_fontsize = 14
ext = '.pdf'
folder = 'report/Figures/thickness/'
ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]
temperature = np.array([300])
parameter_keys = ["alc_temp", "c11", "c12"]
var = np.linspace(0.0, 1.0, 1000)
fixed = np.array([0.1, 0.3, 0.5, 0.8])

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
crit_thick = []
burg = []
poiss = []
deff = []
for i in range(len(temperature)):
    crit_thick.append(dict())
    for fix_type in ["fix_x", "fix_y"]:
        crit_thick[i][fix_type] = []
        for j, fix in enumerate(fixed):
            crit_thick[i][fix_type].append([])
            for k, v in enumerate(var):
                a = quaternary_params_with_temperature[i][fix_type]['{0:1.1f}'.format(fix)]["alc_temp"][k]
                a_sub = par.GaAs["alc_temp"](temperature[i])
                c12 = quaternary_params_with_temperature[i][fix_type]['{0:1.1f}'.format(fix)]["c12"][k]
                c11 = quaternary_params_with_temperature[i][fix_type]['{0:1.1f}'.format(fix)]["c11"][k]

                burgers = a / np.sqrt(2)
                poisson = c12 / (c11 + c12)
                deformation = np.abs((a_sub - a) / a)
                burg.append(burgers)
                poiss.append(poisson)
                deff.append(deformation)
                a = 0.001
                val_a = fun(a, burgers, poisson, deformation)
                b_range = np.linspace(0.001, 400, 2000)
                b = 0
                for n in b_range:
                    val_b = fun(n, burgers, poisson, deformation)
                    # print(val_b)
                    if val_a * val_b < 0:
                        b = n
                        break
                if b == 0:
                    raise Exception("Could not find proper b")
                hc = optimize.brentq(
                    fun, a, b,
                    args=(burgers, poisson, deformation)
                )
                crit_thick[i][fix_type][j].append(hc)

hh = np.linspace(0.01, 20, 1000)
ii = [7113, 6585, 4535, 3373, 2758]
for index in ii:
    funval = fun(hh, burg[index], poiss[index], deff[index])
    plt.plot(hh, funval)
plt.ylim([-75, 45])
plt.axhline(y=0, color='k', linewidth=0.5)
plt.xlabel("Grubość krytyczna $(\AA)$", fontsize=fontsize)
plt.ylabel("f$(h_c)$ $(\AA)$", fontsize=fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.savefig(folder+'fhc'+ext)
plt.gca().clear()
for i in range(len(temperature)):
    for fix_type in ['fix_x', 'fix_y']:
        for j in range(len(fixed)):
            hc = crit_thick[i][fix_type][j]
            if fix_type == 'fix_x':
                label = r'$Al_{y}Ga_{1-y}As_{' + '%.1f' % fixed[j] + r'}Sb_{' + '%.1f' % (1 - fixed[j]) + '}$'
            else:
                label = r'$Al_{' + '%.1f' % fixed[j] + r'}Ga_{' + '%.1f' % (1 - fixed[j]) + '}As_{x}Sb_{1-x}$'
            plt.plot(var, hc, label=label)
        if fix_type == 'fix_x':
            xlabel = 'Ułamek molowy y'
            plt.text(0.1, 2.3, "(a)", fontsize=fontsize)
        else:
            xlabel = 'Ułamek molowy x'
            plt.text(0.0, 2.18, "(b)", fontsize=fontsize)

        plt.legend(loc='best', fontsize=legend_fontsize)
        # plt.legend(bbox_to_anchor=(0, 1, 1, 0), loc="lower left", ncol=2, fontsize=fontsize-4)
        plt.gca().set_ylabel(r'Grubość krytyczna $(\AA)$', fontsize=fontsize)
        plt.gca().set_xlabel(xlabel, fontsize=fontsize)
        plt.yticks(fontsize=tick_fontsize)
        plt.xticks(fontsize=tick_fontsize)
        plt.tight_layout()
        # if fix_type == 'fix_x':
        #     plt.savefig(folder + 'fix_x' + ext)
        # else:
        #     plt.savefig(folder + 'fix_y' + ext)
        plt.gca().clear()


print('Program finished in ', time.time() - start, ' s.')
