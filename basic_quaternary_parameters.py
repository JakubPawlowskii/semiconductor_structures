import binary_materials_parameters as par
import ternary_alloys_bowings as bow
import numpy as np
import matplotlib.pyplot as plt
import time


# def ternary(molar_frac, first_material, second_material, bowing):
#     return molar_frac * first_material + (1 - molar_frac) * second_material - (1 - molar_frac) * molar_frac * bowing
#
#
# def quaternary(molar_frac1, molar_frac2, wx, yz, wy, xz):
#     den = molar_frac1 * (1 - molar_frac1) + molar_frac2 * (1 - molar_frac2)
#     if abs(den) < 1e-6:
#         return np.NaN
#     else:
#         return (molar_frac1 * (1 - molar_frac1) * (molar_frac2 * wx + (1 - molar_frac2) * yz) +
#                 molar_frac2 * (1 - molar_frac2) * (molar_frac1 * wy + (1 - molar_frac1) * xz)) / den
#
#
# # Al_x Ga_1-x As = AlAs + GaAs
# # Al_x Ga_1-x Sb = AlSb + GaSb
# # Al As_x Sb_1-x = AlAs + AlSb
# # Ga As_x Sb_1-x = GaAs + GaSb
#
# start = time.time()
#
# xx = np.linspace(0.0, 1.0, 1001)
# len_x = len(xx)
#
# AlGaAs_Eg = np.empty(len_x)
# AlGaAs_VBO = np.empty(len_x)
# AlGaAs_delta_SO = np.empty(len_x)
# AlGaAs_m_e = np.empty(len_x)
# AlGaAs_m_hh = np.empty(len_x)
# AlGaAs_m_lh = np.empty(len_x)
# AlGaAs_alc = np.empty(len_x)
#
# AlGaSb_Eg = np.empty(len_x)
# AlGaSb_VBO = np.empty(len_x)
# AlGaSb_delta_SO = np.empty(len_x)
# AlGaSb_m_e = np.empty(len_x)
# AlGaSb_m_hh = np.empty(len_x)
# AlGaSb_m_lh = np.empty(len_x)
# AlGaSb_alc = np.empty(len_x)
#
# AlAsSb_Eg = np.empty(len_x)
# AlAsSb_VBO = np.empty(len_x)
# AlAsSb_delta_SO = np.empty(len_x)
# AlAsSb_m_e = np.empty(len_x)
# AlAsSb_m_hh = np.empty(len_x)
# AlAsSb_m_lh = np.empty(len_x)
# AlAsSb_alc = np.empty(len_x)
#
# GaAsSb_Eg = np.empty(len_x)
# GaAsSb_VBO = np.empty(len_x)
# GaAsSb_delta_SO = np.empty(len_x)
# GaAsSb_m_e = np.empty(len_x)
# GaAsSb_m_hh = np.empty(len_x)
# GaAsSb_m_lh = np.empty(len_x)
# GaAsSb_alc = np.empty(len_x)

# Al_x Ga_1-x As = AlAs + GaAs
# Al_x Ga_1-x Sb = AlSb + GaSb
# Al As_x Sb_1-x = AlAs + AlSb
# Ga As_x Sb_1-x = GaAs + GaSb

# # Ternary materials
# for ind, x in enumerate(xx):
#     AlGaAs_Eg[ind] = ternary(x, par.AlAs_Eg, par.GaAs_Eg, bow.AlGaAs_b_Eg(x))
#     AlGaAs_VBO[ind] = ternary(x, par.AlAs_VBO, par.GaAs_VBO, bow.AlGaAs_b_VBO)
#     AlGaAs_delta_SO[ind] = ternary(x, par.AlAs_delta_SO, par.GaAs_delta_SO, bow.AlGaAs_b_delta_SO)
#     AlGaAs_m_e[ind] = ternary(x, par.AlAs_m_e, par.GaAs_m_e, bow.AlGaAs_b_m_e)
#     AlGaAs_m_hh[ind] = ternary(x, par.AlAs_mhh, par.GaAs_mhh, bow.AlGaAs_b_mhh)
#     AlGaAs_m_lh[ind] = ternary(x, par.AlAs_mlh, par.GaAs_mlh, bow.AlGaAs_b_mlh)
#     AlGaAs_alc[ind] = ternary(x, par.AlAs_alc, par.GaAs_alc, bow.AlGaAs_b_alc)
#
#     AlGaSb_Eg[ind] = ternary(x, par.AlSb_Eg, par.GaSb_Eg, bow.AlGaSb_b_Eg(x))
#     AlGaSb_VBO[ind] = ternary(x, par.AlSb_VBO, par.GaSb_VBO, bow.AlGaSb_b_VBO)
#     AlGaSb_delta_SO[ind] = ternary(x, par.AlSb_delta_SO, par.GaSb_delta_SO, bow.AlGaSb_b_delta_SO)
#     AlGaSb_m_e[ind] = ternary(x, par.AlSb_m_e, par.GaSb_m_e, bow.AlGaSb_b_m_e)
#     AlGaSb_m_hh[ind] = ternary(x, par.AlSb_mhh, par.GaSb_mhh, bow.AlGaSb_b_mhh)
#     AlGaSb_m_lh[ind] = ternary(x, par.AlSb_mlh, par.GaSb_mlh, bow.AlGaSb_b_mlh)
#     AlGaSb_alc[ind] = ternary(x, par.AlSb_alc, par.GaSb_alc, bow.AlGaSb_b_alc)
#
#     AlAsSb_Eg[ind] = ternary(x, par.AlAs_Eg, par.AlSb_Eg, bow.AlAsSb_b_Eg)
#     AlAsSb_VBO[ind] = ternary(x, par.AlAs_VBO, par.AlSb_VBO, bow.AlAsSb_b_VBO)
#     AlAsSb_delta_SO[ind] = ternary(x, par.AlAs_delta_SO, par.AlSb_delta_SO, bow.AlAsSb_b_delta_SO)
#     AlAsSb_m_e[ind] = ternary(x, par.AlAs_m_e, par.AlSb_m_e, bow.AlAsSb_b_m_e)
#     AlAsSb_m_hh[ind] = ternary(x, par.AlAs_mhh, par.AlSb_mhh, bow.AlAsSb_b_mhh)
#     AlAsSb_m_lh[ind] = ternary(x, par.AlAs_mlh, par.AlSb_mlh, bow.AlAsSb_b_mlh)
#     AlAsSb_alc[ind] = ternary(x, par.AlAs_alc, par.AlSb_alc, bow.AlAsSb_b_alc)
#
#     GaAsSb_Eg[ind] = ternary(x, par.GaAs_Eg, par.GaSb_Eg, bow.GaAsSb_b_Eg)
#     GaAsSb_VBO[ind] = ternary(x, par.GaAs_VBO, par.GaSb_VBO, bow.GaAsSb_b_VBO)
#     GaAsSb_delta_SO[ind] = ternary(x, par.GaAs_delta_SO, par.GaSb_delta_SO, bow.GaAsSb_b_delta_SO)
#     GaAsSb_m_e[ind] = bow.GaAsSb_m_e(x)
#     GaAsSb_m_hh[ind] = ternary(x, par.GaAs_mhh, par.GaSb_mhh, bow.GaAsSb_b_mhh)
#     GaAsSb_m_lh[ind] = ternary(x, par.GaAs_mlh, par.GaSb_mlh, bow.GaAsSb_b_mlh)
#     GaAsSb_alc[ind] = ternary(x, par.GaAs_alc, par.GaSb_alc, bow.GaAsSb_b_alc)

# Eg = [AlAsSb_Eg, GaAsSb_Eg, AlGaSb_Eg, AlGaAs_Eg]
# vbo = [AlAsSb_VBO, GaAsSb_VBO, AlGaSb_VBO, AlGaAs_VBO]
# delta_so = [AlAsSb_delta_SO, GaAsSb_delta_SO, AlGaSb_delta_SO, AlGaAs_delta_SO]
# m_e = [AlAsSb_m_e, GaAsSb_m_e, AlGaSb_m_e, AlGaAs_m_e]
# m_hh = [AlAsSb_m_hh, GaAsSb_m_hh, AlGaSb_m_hh, AlGaAs_m_hh]
# m_lh = [AlAsSb_m_lh, GaAsSb_m_lh, AlGaSb_m_lh, AlGaAs_m_lh]
# alc = [AlAsSb_alc, GaAsSb_alc, AlGaSb_alc, AlGaAs_alc]

# data = [Eg, vbo, delta_so, m_e, m_hh, m_lh, alc]
#
# labels = [r'$AlAs_{x}Sb_{1-x}$', r'$GaAs_{x}Sb_{1-x}$',
#           r'$Al_{x}Ga_{1-x}Sb$', r'$Al_{x}Ga_{1-x}As$']
#
# folder = 'plots/ternary/'
# ext = '.png'
# filenames = ['eg', 'vbo', 'delta_so', 'm_e', 'm_hh', 'm_lh', 'alc']
# xlabel = 'Ułamek molowy x'
# ylabels = ['Przerwa energetyczna (eV)', 'Położenie pasma walencyjnego (eV)', 'Rozszczepienie spin-orbita (eV)',
#            r'Masa efektywna elektronu ($\frac{m_e^{\ast}}{m_0}$)',
#            r'Masa efektywna dziury ciężkiej ($\frac{m_{hh}^{\ast}}{m_0}$)',
#            r'Masa efektywna dziury lekkiej ($\frac{m_{lh}^{\ast}}{m_0}$)',
#            r'Parametr sieci ($\AA$)']
# filenames = [folder + i + ext for i in filenames]
#
# for i, param in enumerate(data):
#     for j, material in enumerate(param):
#         plt.plot(xx, material, label=labels[j])
#     plt.legend(loc='best')
#     plt.gca().set_xlabel(xlabel)
#     plt.gca().set_ylabel(ylabels[i])
#     plt.tight_layout()
#     plt.xlim([0, 1])
#     plt.savefig(filenames[i])
#     plt.gca().clear()
#
# for i in range(len(alc)):
#     plt.plot(alc[i], Eg[i], label=labels[i])
# plt.legend(loc='best')
# plt.gca().set_xlabel(ylabels[-1])
# plt.gca().set_ylabel(ylabels[0])
# plt.savefig(folder + 'Eg_alc' + ext)
# plt.gca().clear()

# Quaternary material
# Al_x Ga_1-x As_y Sb_1-y
# W -> AlAs (x), X -> GaAs (x), Y -> AlSb (y), Z -> GaSb (y)
#
# composition = np.linspace(0.0, 1.0, 1001)
# fixed_comp = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
#
# fixed_y_Eg = np.empty([len(fixed_comp), len(composition)])
# fixed_x_Eg = np.empty([len(fixed_comp), len(composition)])
#
# fixed_y_vbo = np.empty([len(fixed_comp), len(composition)])
# fixed_x_vbo = np.empty([len(fixed_comp), len(composition)])
#
# fixed_y_delta_so = np.empty([len(fixed_comp), len(composition)])
# fixed_x_delta_so = np.empty([len(fixed_comp), len(composition)])
#
# fixed_y_alc = np.empty([len(fixed_comp), len(composition)])
# fixed_x_alc = np.empty([len(fixed_comp), len(composition)])
#
# fixed_y_m_e = np.empty([len(fixed_comp), len(composition)])
# fixed_x_m_e = np.empty([len(fixed_comp), len(composition)])
#
# fixed_y_m_hh = np.empty([len(fixed_comp), len(composition)])
# fixed_x_m_hh = np.empty([len(fixed_comp), len(composition)])
#
# fixed_y_m_lh = np.empty([len(fixed_comp), len(composition)])
# fixed_x_m_lh = np.empty([len(fixed_comp), len(composition)])
#
# for i, y in enumerate(fixed_comp):
#     for j, x in enumerate(composition):
#         fixed_y_Eg[i, j] = quaternary(x, y, AlGaAs_Eg[j], AlGaSb_Eg[j],
#                                       ternary(y, par.AlAs_Eg, par.AlSb_Eg, bow.AlAsSb_b_Eg),
#                                       ternary(y, par.GaAs_Eg, par.GaSb_Eg, bow.GaAsSb_b_Eg))
#         fixed_y_vbo[i, j] = quaternary(x, y, AlGaAs_VBO[j], AlGaSb_VBO[j],
#                                        ternary(y, par.AlAs_VBO, par.AlSb_VBO, bow.AlAsSb_b_VBO),
#                                        ternary(y, par.GaAs_VBO, par.GaSb_VBO, bow.GaAsSb_b_VBO))
#         fixed_y_delta_so[i, j] = quaternary(x, y, AlGaAs_delta_SO[j], AlGaSb_delta_SO[j],
#                                             ternary(y, par.AlAs_delta_SO, par.AlSb_delta_SO, bow.AlAsSb_b_delta_SO),
#                                             ternary(y, par.GaAs_delta_SO, par.GaSb_delta_SO, bow.GaAsSb_b_delta_SO))
#         fixed_y_alc[i, j] = quaternary(x, y, AlGaAs_alc[j], AlGaSb_alc[j],
#                                        ternary(y, par.AlAs_alc, par.AlSb_alc, bow.AlAsSb_b_alc),
#                                        ternary(y, par.GaAs_alc, par.GaSb_alc, bow.GaAsSb_b_alc))
#         fixed_y_m_e[i, j] = quaternary(x, y, AlGaAs_m_e[j], AlGaSb_m_e[j],
#                                        ternary(y, par.AlAs_m_e, par.AlSb_m_e, bow.AlAsSb_b_m_e),
#                                        bow.GaAsSb_m_e(y))
#         fixed_y_m_hh[i, j] = quaternary(x, y, AlGaAs_m_hh[j], AlGaSb_m_hh[j],
#                                         ternary(y, par.AlAs_mhh, par.AlSb_mhh, bow.AlAsSb_b_mhh),
#                                         ternary(y, par.GaAs_mhh, par.GaSb_mhh, bow.GaAsSb_b_mhh))
#         fixed_y_m_lh[i, j] = quaternary(x, y, AlGaAs_m_lh[j], AlGaSb_m_lh[j],
#                                         ternary(y, par.AlAs_mlh, par.AlSb_mlh, bow.AlAsSb_b_mlh),
#                                         ternary(y, par.GaAs_mlh, par.GaSb_mlh, bow.GaAsSb_b_mlh))
#
# for i, x in enumerate(fixed_comp):
#     for j, y in enumerate(composition):
#         fixed_x_Eg[i, j] = quaternary(x, y, ternary(x, par.AlAs_Eg, par.GaAs_Eg, bow.AlGaAs_b_Eg(x)),
#                                       ternary(x, par.AlSb_Eg, par.GaSb_Eg, bow.AlGaSb_b_Eg(x)),
#                                       AlAsSb_Eg[j], GaAsSb_Eg[j])
#         fixed_x_vbo[i, j] = quaternary(x, y, ternary(x, par.AlAs_VBO, par.GaAs_VBO, bow.AlGaAs_b_VBO),
#                                        ternary(x, par.AlSb_VBO, par.GaSb_VBO, bow.AlGaSb_b_VBO),
#                                        AlAsSb_VBO[j], GaAsSb_VBO[j])
#         fixed_x_delta_so[i, j] = quaternary(x, y,
#                                             ternary(x, par.AlAs_delta_SO, par.GaAs_delta_SO, bow.AlGaAs_b_delta_SO),
#                                             ternary(x, par.AlSb_delta_SO, par.GaSb_delta_SO, bow.AlGaSb_b_delta_SO),
#                                             AlAsSb_delta_SO[j], GaAsSb_delta_SO[j])
#         fixed_x_alc[i, j] = quaternary(x, y, ternary(x, par.AlAs_alc, par.GaAs_alc, bow.AlGaAs_b_alc),
#                                        ternary(x, par.AlSb_alc, par.GaSb_alc, bow.AlGaSb_b_alc),
#                                        AlAsSb_alc[j], GaAsSb_alc[j])
#         fixed_x_m_e[i, j] = quaternary(x, y, ternary(x, par.AlAs_m_e, par.GaAs_m_e, bow.AlGaAs_b_m_e),
#                                        ternary(x, par.AlSb_m_e, par.GaSb_m_e, bow.AlGaSb_b_m_e),
#                                        AlAsSb_m_e[j], GaAsSb_m_e[j])
#         fixed_x_m_hh[i, j] = quaternary(x, y, ternary(x, par.AlAs_mhh, par.GaAs_mhh, bow.AlGaAs_b_mhh),
#                                         ternary(x, par.AlSb_mhh, par.GaSb_mhh, bow.AlGaSb_b_mhh),
#                                         AlAsSb_m_hh[j], GaAsSb_m_hh[j])
#         fixed_x_m_lh[i, j] = quaternary(x, y, ternary(x, par.AlAs_mlh, par.GaAs_mlh, bow.AlGaAs_b_mlh),
#                                         ternary(x, par.AlSb_mlh, par.GaSb_mlh, bow.AlGaSb_b_mlh),
#                                         AlAsSb_m_lh[j], GaAsSb_m_lh[j])

# fixed_x = [fixed_x_Eg, fixed_x_vbo, fixed_x_delta_so, fixed_x_alc, fixed_x_m_e, fixed_x_m_hh, fixed_x_m_lh]
# fixed_y = [fixed_y_Eg, fixed_y_vbo, fixed_y_delta_so, fixed_y_alc, fixed_y_m_e, fixed_y_m_hh, fixed_y_m_lh]
#
# xlabel = 'Ułamek molowy '
# ylabels = ['Przerwa energetyczna (eV)', 'Położenie pasma walencyjnego (eV)', 'Rozszczepienie spin-orbita (eV)',
#            r'Parametr sieci ($\AA$)',
#            r'Masa efektywna elektronu ($\frac{m_e^{\ast}}{m_0}$)',
#            r'Masa efektywna dziury ciężkiej ($\frac{m_{hh}^{\ast}}{m_0}$)',
#            r'Masa efektywna dziury lekkiej ($\frac{m_{lh}^{\ast}}{m_0}$)']
# folder = 'plots/'
# ext = '.png'
# filenames = ['quat_eg', 'quat_vbo', 'quat_delta_so', 'quat_alc', 'quat_m_e', 'quat_m_hh', 'quat_m_lh']
# filenames_x = [folder + i + '_x' + ext for i in filenames]
# filenames_y = [folder + i + '_y' + ext for i in filenames]
# for ind, fix in enumerate(fixed_x):
#     for i in range(len(fixed_comp)):
#         label = ''
#         if fixed_comp[i] == 0.0:
#             label = r'$GaAs_{y}Sb_{1-y}$'
#         elif fixed_comp[i] == 1.0:
#             label = r'$AlAs_{y}Sb_{1-y}$'
#         else:
#             label = r'$Al_{' + '%.1f' % fixed_comp[i] + r'}Ga_{' + '%.1f' % (1 - fixed_comp[i]) + '}As_{y}Sb_{1-y}$'
#         plt.plot(composition, fix[i, :], label=label)
#     plt.gca().set_ylabel(ylabels[ind])
#     plt.gca().set_xlabel(xlabel + 'y')
#     plt.legend(loc='best', fontsize=9)
#     plt.tight_layout()
#     plt.xlim([0.0, 1.0])
#     plt.savefig(filenames_x[ind])
#     plt.gca().clear()
#
# for ind, fix in enumerate(fixed_y):
#     for i in range(len(fixed_comp)):
#         label = ''
#         if fixed_comp[i] == 0.0:
#             label = r'$Al_{x}Ga_{1-x}Sb$'
#         elif fixed_comp[i] == 1.0:
#             label = r'$Al_{x}Ga_{1-x}As$'
#         else:
#             label = r'$Al_{x}Ga_{1-x}As_{' + '%.1f' % fixed_comp[i] + '}Sb_{' + '%.1f' % (1 - fixed_comp[i]) + '}$'
#         plt.plot(composition, fix[i, :], label=label)
#     plt.gca().set_ylabel(ylabels[ind])
#     plt.gca().set_xlabel(xlabel + 'x')
#     plt.legend(loc='best', fontsize=9)
#     plt.tight_layout()
#     plt.xlim([0.0, 1.0])
#     plt.savefig(filenames_y[ind])
#     plt.gca().clear()

start = time.time()
ternary_mat = ["AlGaAs", "AlGaSb", "AlAsSb", "GaAsSb"]
parameter_keys = ["Eg", "VBO", "delta_so", "alc", "m_e", "m_hh", "m_lh"]
var = np.linspace(0.001, 0.999, 1000)

ternary_params = dict()
for mat in ternary_mat:
    ternary_params[mat] = dict()

for par_key in parameter_keys:
    p_AlAs = par.AlAs[par_key]
    p_GaAs = par.GaAs[par_key]
    p_AlSb = par.AlSb[par_key]
    p_GaSb = par.GaSb[par_key]

    if par_key == "Eg":
        b_WX = bow.AlGaAs[par_key](var)
        b_YZ = bow.AlGaSb[par_key](var)
        b_WY = bow.AlAsSb[par_key]
        b_XZ = bow.GaAsSb[par_key]
    else:
        b_WX = bow.AlGaAs[par_key]
        b_YZ = bow.AlGaSb[par_key]
        b_WY = bow.AlAsSb[par_key]
        b_XZ = bow.GaAsSb[par_key]

    ternary_params["AlGaAs"][par_key] = var * p_AlAs + (1 - var) * p_GaAs - var * (1 - var) * b_WX
    ternary_params["AlGaSb"][par_key] = var * p_AlSb + (1 - var) * p_GaSb - var * (1 - var) * b_YZ
    ternary_params["AlAsSb"][par_key] = var * p_AlAs + (1 - var) * p_AlSb - var * (1 - var) * b_WY
    ternary_params["GaAsSb"][par_key] = var * p_GaAs + (1 - var) * p_GaSb - var * (1 - var) * b_XZ


labels = [r'$AlAs_{x}Sb_{1-x}$', r'$GaAs_{x}Sb_{1-x}$',
          r'$Al_{x}Ga_{1-x}Sb$', r'$Al_{x}Ga_{1-x}As$']

folder = 'plots/ternary/'
ext = '.png'
filenames = ['eg', 'vbo', 'delta_so', 'm_e', 'm_hh', 'm_lh', 'alc']
xlabel = 'Ułamek molowy x'
ylabels = ['Przerwa energetyczna (eV)', 'Położenie pasma walencyjnego (eV)', 'Rozszczepienie spin-orbita (eV)',
           r'Masa efektywna elektronu ($\frac{m_e^{\ast}}{m_0}$)',
           r'Masa efektywna dziury ciężkiej ($\frac{m_{hh}^{\ast}}{m_0}$)',
           r'Masa efektywna dziury lekkiej ($\frac{m_{lh}^{\ast}}{m_0}$)',
           r'Parametr sieci ($\AA$)']
filenames = [folder + i + ext for i in filenames]

for i, param_key in enumerate(parameter_keys):
    for j, mat in enumerate(ternary_mat):
        plt.plot(var, ternary_params[mat][param_key], label=labels[j])
    plt.legend(loc='best')
    plt.gca().set_xlabel(xlabel)
    plt.gca().set_ylabel(ylabels[i])
    plt.tight_layout()
    plt.xlim([0, 1])
    plt.savefig(filenames[i])
    plt.gca().clear()

for j, mat in enumerate(ternary_mat):
    plt.plot(ternary_params[mat]["alc"], ternary_params[mat]["Eg"], label=labels[j])
plt.legend(loc='best')
plt.gca().set_xlabel(ylabels[-1])
plt.gca().set_ylabel(ylabels[0])
plt.savefig(folder + 'Eg_alc' + ext)
plt.gca().clear()


print('Program finished in ', time.time() - start, ' s.')
