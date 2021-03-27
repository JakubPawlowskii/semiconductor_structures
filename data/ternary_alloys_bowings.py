# Al_x Ga_1-x As
def AlGaAs_b_Eg(x):
    return -0.127 + 1.310*x


def AlGaSb_b_Eg(x):
    return -0.044 + 1.22*x


AlGaAs = {
    "Eg": AlGaAs_b_Eg,
    "VBO": 0.0,
    "delta_so": 0.0,
    "alc": 0.0,
    "m_e": 0.0,
    "m_hh": 0.0,
    "m_lh": 0.0
}


# Al_x Ga_1-x Sb

AlGaSb = {
    "Eg": AlGaSb_b_Eg,
    "VBO": 0.0,
    "delta_so": 0.3,
    "alc": 0.0,
    "m_e": 0.0,
    "m_hh": 0.0,
    "m_lh": 0.0
}


# Al As_x Sb_1-x

AlAsSb = {
    "Eg": 0.8,
    "VBO": -1.71,
    "delta_so": 0.15,
    "alc": 0.0,
    "m_e": 0.0,
    "m_hh": 0.0,
    "m_lh": 0.0
}


# Ga As_x Sb_1-x

GaAsSb = {
    "Eg": 1.43,
    "VBO": -1.06,
    "delta_so": 0.6,
    "alc": 0.0,
    "m_e": 0.014,
    "m_hh": 0.0,
    "m_lh": 0.0
}
