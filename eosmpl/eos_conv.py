import math
import numpy as np

# Metric prefixes (log exponents)
dict_pref = {
    'y'  : -24,
    'z'  : -21,
    'a'  : -18,
    'f'  : -15, 
    'p'  : -12,
    'n'  : -9,
    'mu' : -6,
    'm'  : -3,
    'c'  : -2,
    'd'  : -1,
    'x'  : 0,
    'da' : 1,
    'h'  : 2,
    'k'  : 3,
    'M'  : 6,
    'G'  : 9,
    'T'  : 12,
    'P'  : 15,
    'E'  : 18,
    'Z'  : 21,
    'Y'  : 24
}

# Unit conversion multipliers (to SI base units)
dict_conv = {
    'l' : # Length (to meters)
    {
        'm'  : 1.0,
        'Mm' : 1e6,
        'km' : 1000.0,
        'cm' : 0.01,
        'mm' : 0.001,
        'ml' : 5280.0*0.3048,
        'yd' : 0.9144,
        'ft' : 0.3048,
        'in' : 0.0254,
        'au' : 149597870700.0,
        'AU' : 149597870700.0,
        'RE' : 6.3781e6,
        'RJ' : 7.1492e7,
        'nm' : 1852.0
    },
    'A' : # Area (to square meters)
    {
        'm2'   : 1.0,
        'km2'  : 1e6,
        'cm2'  : 1e-4,
        'mm2'  : 1e-6,
        'ml2'  : (5280.0*0.3048)**2,
        'yd2'  : 0.83612736,
        'ft2'  : 0.09290304,
        'in2'  : 6.4516e-4,
        'barn' : 1e-28
    },
    'V' : # Volume (to cubic meters)
    {
        'm3'  : 1.0,
        'dm3' : 1e-3,
        'l'   : 1e-3,
        'cm3' : 1e-6,
        'mm3' : 1e-9,
        'ml3' : (5280.0*0.3048)**3,
        'yd3' : 0.9144**3,
        'ft3' : 0.3048**3,
        'in3' : 0.0254**3,
        'gal' : 3.785412e-3
    },
    'v' : # Velocity (to m/s)
    {
        'm/s'  : 1.0,
        'Mm/s' : 1e6,
        'km/s' : 1e3,
        'km/h' : 1.0/3.6,
        'ml/h' : 5280.0*0.3048/3600,
        'ml/s' : 5280.0*0.3048,
        'ft/s' : 0.3048,
        'kn'   : 1852/3.6,
        'a'    : np.sqrt(1.4*287.058*288.15), #ISA SL conditions
        'c'    : 299792458
    },
    'm' : # Mass (to grams)
    {
        'g'    : 1.0,
        't'    : 1e6,
        'kg'   : 1e3,
        'slug' : 14593.90294,
        'lb'   : 453.59237,
        'oz'   : 28.349523125,
        'Mo'   : 1.98855e33,
        'ME'   : 5.9722e24,
        'u'    : 1.660539040e-24
    },
    'D' : # Density (to g/m3)
    {
        'g/m3'     : 1.0,
        'kg/m3'    : 1e9,
        'g/cm3'    : 1e6,
        'kg/dm3'   : 1e6,    
        'kg/m3'    : 1e3,
        'g/dm3'    : 1e3,
        'slug/ft3' : 14593.90294/(0.3048**3),
        'lb/ft3'   : 453.59237/(0.3048**3)
    },
    'F' : # Force (to Newtons)
    {
        'N'   : 1.0,
        'kN'  : 1e3,
        'g_E' : 9.80665,
        'lbf' : 4.448222,
        'dyn' : 1e-5
    },
    'p' : # Pressure (to Pascales)
    {
        'Pa'   : 1.0,
        'atm'  : 101325,
        'at'   : 98066.5,
        'bar'  : 1e5,
        'kPa'  : 1e3,
        'mbar' : 1e2,
        'psi'  : 6.894757e3
    },
    'E' : # Energy (to Joules)
    {
        'J'     : 1.0,
        'kJ'    : 1e3,
        'Wh'    : 3600,
        'kWh'   : 3600000,
        'ft.lb' : 1.355818,
        'cal'   : 4.184,
        'kcal'  : 4184,
        'eV'    : 1.6021766208e-19
    },
    'P' : # Power (to Joules)
    {
        'W'  : 1.0,
        'kW' : 1e3,
        'hp' : 745.69987158227022
    },
    't' : # Time (to seconds)
    {
        's'    : 1.0,
        'sec'    : 1.0,
        'yr'   : 31557600,
        'yr_s' : 365.25636*86400,
        'mo'   : 30.0*86400.0,
        'mo_s' : 27.321661*86400,
        'd'    : 86400.0,
        'dy'    : 86400.0,
        'd_s'  : 86164.090530833,
        'dy_s'  : 86164.090530833,
        'hr'   : 3600.0,
        'min'  : 60.0
    }
}

# Temperature units
dict_temp = ['F', 'K', 'C', 'R']

# Metric prefix to log
def mprf(i):
    return dict_pref[i]

# Temperature conversion
def ctmp(i,o,iv):
    if i != o:
        if i == 'C':
            x = float(iv) + 273.15
        elif i == 'F':
            x =(float(iv) + 459.67) * (5.0/9.0)
        elif i == 'R':
            x = float(iv)*(5.0/9.0)
        elif i == 'K':
            x = float(iv)
        if o == 'C':
            return x - 273.15
        elif o == 'F':
            return x * (9.0/5.0) - 459.67
        elif o == 'R':
            return x * (9.0/5.0)
        elif o == 'K':
            return x
    else:
        return float(iv)

# Unit conversion
def unit(iv,i,o,d=3):
    try:
        d = int(d)
    except TypeError:
        raise TypeError('Only integers allowed as decimal places')
    d = int(d)
    if d < 0:
        raise ValueError('Only positive decimal places allowed')
    try:
        iv = float(iv)
    except TypeError:
        raise TypeError('Invalid value input')
    iv = float(iv)
    isplit = i.split('-')
    if len(isplit) == 2:
        ilog = float(mprf(isplit[0]))
        ib = isplit[1]
    elif len(isplit) == 1:
        ilog = float(0)
        ib = isplit[0]
    else:
        raise SyntaxError('Invalid unit input')

    osplit = o.split('-')
    if len(osplit) == 2:
        olog = float(mprf(osplit[0]))
        ob = osplit[1]
    elif len(osplit) == 1:
        olog = float(0)
        ob = osplit[0]
    else:
        raise SyntaxError('Invalid unit input')

    if ib in dict_temp and ob in dict_temp:
        return ctmp(ib,ob,iv)
    else:
        for unit_parents, unit_children in dict_conv.items():
            if ib in unit_children and ob in unit_children:
                if unit_parents == 'A':
                    ilog = ilog*2.0
                    olog = olog*2.0
                elif unit_parents == 'V':
                    ilog = ilog*3.0
                    olog = olog*3.0
                output = float(iv*10**ilog)*float(unit_children[ib]) * (1.0 / (float(unit_children[ob])*float(10**olog)))
                return output
        raise SyntaxError('Illegal unit combination')

# Geometric altitude to geopotential altitude and vice versa (def. inp. meters; Earth cond.)
# Geometric to geopotential
def hg_to_h(hg):
    hg = float(hg)
    return (hg * 6.3781e6) / (6.3781e6 + hg)
# Geopotential to geometric
def h_to_hg(h):
    h = float(h)
    return (h * 6.3781e6) / (6.3781e6 - h)