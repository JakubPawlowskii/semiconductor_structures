from brokenaxes import brokenaxes
from data import binary_materials_parameters as par
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib import cm
from matplotlib.gridspec import GridSpec

start = time.time()
temperature = np.array([0, 150, 300, 450])
kz_values = np.linspace(-1, 1, 200)

Ec = dict()
Ehh = dict()
Elh = dict()
Eso = dict()

Ec_DOS = dict()
Ehh_DOS = dict()
Elh_DOS = dict()
Eso_DOS = dict()
E = dict()
Eg_T = dict()

materials = ["GaAs"]
# ylims = ((-1.8, -1.0), (1.5, 2.1)) #AlAs
# ylims = ((-1.4, 0.0), (0.6, 1.8)) #GaSb
ylims = ((-1.45, -0.5), (0.4, 1.4)) #GaAs
# ylims = ((-1.4, -0.3), (1.8, 2.4)) #AlSb

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
constDOS = 1 /(2 * np.pi ** 2) * (1 / const) ** (3 / 2)
for material in materials:
    Ec[material] = np.zeros([len(temperature), len(kz_values)], dtype=float)
    Ehh[material] = np.zeros([len(temperature), len(kz_values)], dtype=float)
    Elh[material] = np.zeros([len(temperature), len(kz_values)], dtype=float)
    Eso[material] = np.zeros([len(temperature), len(kz_values)], dtype=float)
    Eg_T[material] = np.zeros(len(temperature))
    E[material] = np.zeros([len(temperature), 200], dtype=float)
    for i, T in enumerate(temperature):

        Eg_T[material][i] = Eg[material] - (alpha[material] * (T ** 2) / (T + beta[material]))

        for j, kz in enumerate(kz_values):
            Ec[material][i, j] = VBO[material] + Eg_T[material][i] + const * kz ** 2 * 1 / m_e[material]
            Ehh[material][i, j] = VBO[material] - const * kz ** 2 * (gamma1[material] - 2 * gamma2[material])
            Elh[material][i, j] = VBO[material] - const * kz ** 2 * (gamma1[material] + 2 * gamma2[material])
            Eso[material][i, j] = VBO[material] - const * kz ** 2 * gamma1[material] - delta_so[material]

        E[material][i, :] = np.linspace(np.min(Eso[material][i, :]), np.max(Ec[material][i, :]), 200)

for material in materials:
    Ec_DOS[material] = np.zeros([len(temperature), len(E[material][0, :])], dtype=float)
    Ehh_DOS[material] = np.zeros([len(temperature), len(E[material][0, :])], dtype=float)
    Elh_DOS[material] = np.zeros([len(temperature), len(E[material][0, :])], dtype=float)
    Eso_DOS[material] = np.zeros([len(temperature), len(E[material][0, :])], dtype=float)
    for i, T in enumerate(temperature):
        for j, E_val in enumerate(E[material][i, :]):
            Ec_DOS[material][i, j] = constDOS * m_e[material] ** (3 / 2) * (
                        E_val - VBO[material] - Eg_T[material][i]) ** (1 / 2)
            if np.isnan(Ec_DOS[material][i,j]):
                Ec_DOS[material][i, j] = 1e-20
            Ehh_DOS[material][i, j] = constDOS * m_hh[material] ** (3 / 2) * (-(E_val - VBO[material])) ** (1 / 2)
            if np.isnan(Ehh_DOS[material][i,j]):
                Ehh_DOS[material][i, j] = 1e-20
            Elh_DOS[material][i, j] = constDOS * m_lh[material] ** (3 / 2) * (-(E_val - VBO[material])) ** (1 / 2)
            if np.isnan(Elh_DOS[material][i,j]):
                Elh_DOS[material][i, j] = 1e-20
            Eso_DOS[material][i, j] = constDOS * m_so[material] ** (3 / 2) * (
                -(E_val - VBO[material] + delta_so[material])) ** (1 / 2)
            if np.isnan(Eso_DOS[material][i,j]):
                Eso_DOS[material][i, j] = 1e-20

# fig = plt.figure(figsize=(8, 8))
# gs = fig.add_gridspec(1, 2, width_ratios=(7, 4),
#                       left=0.1, right=0.9, bottom=0.1, top=0.9,
#                       wspace=0.05, hspace=0.05)
# ax = fig.add_subplot(gs[0])
# ax_DOS = fig.add_subplot(gs[1], sharey=ax)

fig = plt.figure()
sps = GridSpec(nrows=1, ncols=2, width_ratios=(7, 4),
               left=0.12, right=0.95, bottom=0.13, top=0.9,
               wspace=0.1, hspace=0.05)


ax = brokenaxes(ylims=ylims, hspace=.1, subplot_spec=sps[0], d=0.01)
ax_DOS = brokenaxes(ylims=ylims, hspace=.1, subplot_spec=sps[1], d=0.01)
ax_DOS.tick_params(axis="y", labelleft=False)
ax_DOS.set_xlim([0,0.2])
filename = ''
path = 'report/Figures/band_str/'
linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
linewidth = 0.8
fontsize = 14
tick_fontsize = 15
legend_fontsize = 12
leg_no_T = [r'$E_c$', r'$E_{v,hh}$', r'$E_{v,lh}$', r'$E_{v,sh}$']
for material in materials:
    # for i in range(len(temperature)):
    #     l1 = ax.plot(kz_values, Ec[material][i, :], linestyle=linestyles[i])
    #     ax_DOS.plot(Ec_DOS[material][i, :], E[material][i], color=l1[0][0].get_color(), linestyle=linestyles[i],
    #                 label="$E_c,\;"+ "T = " + str(temperature[i])+"\;K$")
    #
    # l2 = ax.plot(kz_values, Ehh[material][0, :])
    # l3 = ax.plot(kz_values, Elh[material][0, :])
    # l4 = ax.plot(kz_values, Eso[material][0, :])
    # ax_DOS.plot(Ehh_DOS[material][0, :], E[material][0], color=l2[0][0].get_color(), label="$E_{v,hh}$")
    # ax_DOS.plot(Elh_DOS[material][0, :], E[material][0], color=l3[0][0].get_color(), label="$E_{v,lh}$")
    # ax_DOS.plot(Eso_DOS[material][0, :], E[material][0], color=l4[0][0].get_color(), label="$E_{v,sh}$")

    colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728']
    for i in range(len(temperature)):
        ax.plot(kz_values, Ec[material][i, :], color=colors[0], linestyle=linestyles[i])
        ax_DOS.plot(Ec_DOS[material][i, :], E[material][i], color=colors[0], linestyle=linestyles[i])
        ax.plot(kz_values, Ehh[material][i, :], color=colors[1], linestyle=linestyles[i])
        ax.plot(kz_values, Elh[material][i, :], color=colors[2], linestyle=linestyles[i])
        ax.plot(kz_values, Eso[material][i, :], color=colors[3], linestyle=linestyles[i])
        ax_DOS.plot(Ehh_DOS[material][i, :], E[material][i], color=colors[1], linestyle=linestyles[i])
        ax_DOS.plot(Elh_DOS[material][i, :], E[material][i], color=colors[2], linestyle=linestyles[i])
        ax_DOS.plot(Eso_DOS[material][i, :], E[material][i], color=colors[3], linestyle=linestyles[i])
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

    filename = material+'.pdf'
    # ax_DOS.legend(loc='best')
    # ax.set_ylabel(r'Energia [eV]', fontsize)
    ax_DOS.set_xlim(left=0, right=None)
    txt = ''
    if material == "AlAs":
        txt = '(a)'
    elif material == 'GaAs':
        txt = '(b)'
    elif material == 'AlSb':
        txt = '(c)'
    elif material == 'GaSb':
        txt = '(d)'
    fig.text(0.2,0.85,txt,va='center')
    fig.text(0.0,0.37,r'Energia [eV]', fontsize=fontsize, rotation='vertical' )
    fig.text(0.32, 0.06, r'$k_z\;[1/nm]$', va='center',fontsize=fontsize)
    fig.text(0.7, 0.06, r'$DOS\;[stany/eV]$', va='center',fontsize=fontsize)
    # plt.savefig(path+filename)
    # ax.fig.gca().clear()
    # ax_DOS.fig.gca().clear()
    # plt.gca().set_prop_cycle(None)
    plt.show()
print('Program finished in ', time.time() - start, ' s.')
