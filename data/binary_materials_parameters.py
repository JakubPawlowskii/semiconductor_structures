# File storing binary materials parameters

# GaAs
def alc_GaAs_temp(temperature):
    return 5.65325 + 3.88*10**(-5)*(temperature-300)


GaAs = {
    "Eg": 1.519,
    "VBO": -0.80,
    "delta_so": 0.341,
    "alc": 5.65325,
    "m_e": 0.067,
    "m_hh": 0.55,
    "m_lh": 0.083,
    "alpha": 0.5405*(10**(-3)),
    "beta": 204,
    "alc_temp": alc_GaAs_temp,
    "ac": -7.17,
    "av": -1.16,
    "b": -2.0,
    "c11": 1221,
    "c12": 566
}

# GaSb


def alc_GaSb_temp(temperature):
    return 6.0959 + 4.72*10**(-5)*(temperature-300)


GaSb = {
    "Eg": 0.812,
    "VBO": -0.03,
    "delta_so": 0.76,
    "alc": 6.0959,
    "m_e": 0.039,
    "m_hh": 0.37,
    "m_lh": 0.043,
    "alpha": 0.417*(10**(-3)),
    "beta": 140,
    "alc_temp": alc_GaSb_temp,
    "ac": -7.5,
    "av": -0.8,
    "b": -2.0,
    "c11": 884.2,
    "c12": 402.6
}

# AlAs


def alc_AlAs_temp(temperature):
    return 5.6611+2.90*10**(-5)*(temperature-300)


AlAs = {
    "Eg": 3.099,
    "VBO": -1.33,
    "delta_so": 0.28,
    "alc": 5.6611,
    "m_e": 0.15,
    "m_hh": 0.81,
    "m_lh": 0.16,
    "alpha": 0.885*(10**(-3)),
    "beta": 530,
    "alc_temp": alc_AlAs_temp,
    "ac": -5.64,
    "av": -2.47,
    "b": -2.3,
    "c11": 1250,
    "c12": 534
}

# AlSb


def alc_AlSb_temp(temperature):
    return 6.1355 + 2.60*10**(-5)*(temperature-300)


AlSb = {
    "Eg": 2.386,
    "VBO": -0.41,
    "delta_so": 0.676,
    "alc": 6.1355,
    "m_e": 0.14,
    "m_hh": 0.9,
    "m_lh": 0.13,
    "alpha": 0.42*(10**(-3)),
    "beta": 140,
    "alc_temp": alc_AlSb_temp,
    "ac": -4.5,
    "av": -1.4,
    "b": -1.35,
    "c11": 876.9,
    "c12": 434.1
}



