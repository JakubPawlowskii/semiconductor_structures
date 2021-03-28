from data import ternary_alloys_bowings as bow, binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time

start = time.time()
ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]

parameter_keys = ["Eg", "VBO", "delta_so", "alc", "m_e", "m_hh", "m_lh"]
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

# # Al_x Ga_1-x As = AlAs + GaAs
# # Al_x Ga_1-x Sb = AlSb + GaSb
# # Al As_x Sb_1-x = AlAs + AlSb
# # Ga As_x Sb_1-x = GaAs + GaSb

# Quaternary material
# Al_x Ga_1-x As_y Sb_1-y
# W -> AlAs (x), X -> GaAs (x), Y -> AlSb (y), Z -> GaSb (y)


for par_key in parameter_keys:
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
        p = var * (1 - var) * (y * ternary_params["AlGaAs"][par_key] + (1 - y) * ternary_params["AlGaSb"][par_key]) + \
            y * (1 - y) * (var * p_AlAsSb + (1 - var) * p_GaAsSb)
        quaternary_params["fix_y"]['{0:1.1f}'.format(y)][par_key] = p / (var * (1 - var) + y * (1 - y))

fontsize = 18
tick_fontsize = 14
legend_fontsize = 12
labels = [r'$Al_{x}Ga_{1-x}As$', r'$Al_{x}Ga_{1-x}Sb$', r'$AlAs_{x}Sb_{1-x}$', r'$GaAs_{x}Sb_{1-x}$']
folder = 'report/Figures/ternary/'
ext = '.pdf'
filenames = ['eg', 'vbo', 'delta_so', 'alc', 'm_e', 'm_hh', 'm_lh']
xlabel = 'Ułamek molowy x'
ylabels = ['Przerwa energetyczna (eV)', 'Położenie pasma walencyjnego (eV)', 'Rozszczepienie spin-orbita (eV)',
           r'Parametr sieci ($\AA$)',
           r'Masa efektywna elektronu ($\frac{m_e^{\ast}}{m_0}$)',
           r'Masa efektywna dziury ciężkiej ($\frac{m_{hh}^{\ast}}{m_0}$)',
           r'Masa efektywna dziury lekkiej ($\frac{m_{lh}^{\ast}}{m_0}$)']
filenames = [folder + i + ext for i in filenames]


for i, param_key in enumerate(parameter_keys):
    for j, mat in enumerate(ternary_mat):
        plt.plot(var, ternary_params[mat][param_key], label=labels[j])
    if param_key == "Eg":
        plt.legend(loc='best', fontsize=legend_fontsize + 2)
    plt.gca().set_xlabel(xlabel, fontsize=fontsize)
    plt.gca().set_ylabel(ylabels[i], fontsize=fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.tight_layout()
    plt.xlim([0, 1])
    plt.savefig(filenames[i])
    plt.gca().clear()

for j, mat in enumerate(ternary_mat):
    plt.plot(ternary_params[mat]["alc"], ternary_params[mat]["Eg"], label=labels[j])
plt.legend(loc='best', fontsize=legend_fontsize)
plt.gca().set_xlabel(ylabels[-1], fontsize=fontsize-2)
plt.gca().set_ylabel(ylabels[0], fontsize=fontsize-2)
plt.yticks(fontsize=tick_fontsize-2)
plt.xticks(fontsize=tick_fontsize-2)
plt.savefig(folder + 'Eg_alc' + ext)
plt.gca().clear()

xlabel = 'Ułamek molowy '
folder = 'report/Figures/quaternary/'
filenames = ['quat_eg', 'quat_vbo', 'quat_delta_so', 'quat_alc', 'quat_m_e', 'quat_m_hh', 'quat_m_lh']
filenames_x = [folder + i + '_x' + ext for i in filenames]
filenames_y = [folder + i + '_y' + ext for i in filenames]
filenames_fix = {"fix_x": filenames_x, "fix_y": filenames_y}

for i, par_key in enumerate(parameter_keys):
    for ind, fix_val in enumerate(fixed):
        label = ''
        if fix_val == 0.0:
            label = r'$GaAs_{y}Sb_{1-y}$'
        elif fix_val == 1.0:
            label = r'$AlAs_{y}Sb_{1-y}$'
        else:
            label = r'$Al_{' + '%.1f' % fix_val + r'}Ga_{' + '%.1f' % (1 - fix_val) + '}As_{y}Sb_{1-y}$'
        plt.plot(var, quaternary_params["fix_x"]['{0:1.1f}'.format(fix_val)][par_key], label=label)
    plt.gca().set_ylabel(ylabels[i], fontsize=fontsize)
    plt.gca().set_xlabel(xlabel + 'y', fontsize=fontsize)
    if par_key == "VBO":
        plt.legend(loc='best', fontsize=legend_fontsize + 2)
    plt.yticks(fontsize=tick_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.tight_layout()
    plt.xlim([0.0, 1.0])
    plt.savefig(filenames_fix["fix_x"][i])
    plt.gca().clear()

    for ind, fix_val in enumerate(fixed):
        label = ''
        if fix_val == 0.0:
            label = r'$GaAs_{x}Sb_{1-x}$'
        elif fix_val == 1.0:
            label = r'$AlAs_{x}Sb_{1-x}$'
        else:
            label = r'$Al_{' + '%.1f' % fix_val + r'}Ga_{' + '%.1f' % (1 - fix_val) + '}As_{x}Sb_{1-x}$'
        plt.plot(var, quaternary_params["fix_y"]['{0:1.1f}'.format(fix_val)][par_key], label=label)
    plt.gca().set_ylabel(ylabels[i], fontsize=fontsize)
    plt.gca().set_xlabel(xlabel + 'x', fontsize=fontsize)
    if par_key == "Eg":
        plt.legend(loc='best', fontsize=legend_fontsize + 2)
    plt.yticks(fontsize=tick_fontsize)
    plt.xticks(fontsize=tick_fontsize)
    plt.tight_layout()
    plt.xlim([0.0, 1.0])
    plt.savefig(filenames_fix["fix_y"][i])
    plt.gca().clear()

print('Program finished in ', time.time() - start, ' s.')
