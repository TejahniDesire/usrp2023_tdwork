import numpy as np
from astropy import constants as const
from astropy import units as u

# TODO Check repo by jason (B field values)| V2 check to V3| Put Plots in Overleaf with defined parameters

# All r in units of Rg

G = const.G.cgs
c = const.c.cgs
me = const.m_e.cgs
mp = const.m_p.cgs
e = const.e.esu
kB = const.k_B.cgs

# Units
cm = u.cm
cmcubed = 1 / (cm ** 3)
kelv = 1 * u.K
grams = 1 * u.g
gauss = 1 * u.cm ** (-1 / 2) * u.g ** (1 / 2) * 1 / (1 * u.s)  # Gauss units in cgs
rads = 1 * u.rad
Hz = 1 * u.Hz

# KWARGS

kw_n_th0 = 1.23 * 10 ** 4 * cmcubed
kw_t_e0 = 8.1 * 10 ** 9 * kelv
kw_beta = 1
kw_mass = 1.989 * 10 ** 42 * grams  # m87 mass
#kw_mass = 9 * 10 ** 39 * grams
kw_p_dens = -.7
kw_p_temp = -.84
kw_theta_b = 60 * (np.pi / 180) * rads
kw_nu = 230 * 10 ** 9 * Hz
kw_scale_height = .5


def rg_func(mass=kw_mass):
    return G * mass / c ** 2


def rb_func(mass=kw_mass):
    return 20 * G * mass / (c ** 2)


def te_func(r, mass=kw_mass, t_e0=kw_t_e0, p_temp=kw_p_temp):
    if t_e0.unit != u.K:
        raise ValueError('n_th0 needs units of Kelvin')
    rb = rb_func(mass)
    rg = rg_func(mass)
    return t_e0 * (r * rg / rb) ** p_temp


def theta_e_func(r, mass=kw_mass, t_e0=kw_t_e0, p_temp=kw_p_temp):
    temp = te_func(r, mass, t_e0, p_temp)
    #return (kB * temp / (me * c ** 2)+1e-10).to(u.dimensionless_unscaled)
    return (kB * temp / (me * c ** 2)).to(u.dimensionless_unscaled)


def nth_func(r, mass=kw_mass, n_th0=kw_n_th0, p_dens=kw_p_dens):
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    rg = rg_func(mass)
    rb = rb_func(mass)
    return n_th0 * (r * rg / rb) ** p_dens


def b_func(r, mass=kw_mass, n_th0=kw_n_th0):
    rg = rg_func(mass)
    rb = rb_func(mass)
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    return 1.96789 * gauss * (r * rg / rb) ** -0.85  # From parameter fitting for KWARGS
    # return np.sqrt(nth * 8 * np.pi * mp * c ** 2 / (6 * beta * r))


def nu_c_func(r, mass=kw_mass, n_th0=kw_n_th0, t_e0=kw_t_e0, p_temp=kw_p_temp, theta_b=kw_theta_b):
    b_field = b_func(r, mass, n_th0)
    nu_b = (e * b_field / (2 * np.pi * me * c)).to(u.Hz)

    theta_e = theta_e_func(r, mass, t_e0, p_temp)
    return 3 / 2 * nu_b * np.sin(theta_b) * theta_e ** 2


def synchrotron_func(x):
    return 2.5651 * (1 + 1.92 * (x ** (-1 / 3)) + (0.9977 * x ** (-2 / 3))) * np.exp(-1.8899 * x ** (1 / 3))


def j_emission_nu_func(r, nu=kw_nu, mass=kw_mass, n_th0=kw_n_th0, p_dens=kw_p_dens, t_e0=kw_t_e0, p_temp=kw_p_temp,
                       theta_b=kw_theta_b):
    n = nth_func(r, mass, n_th0, p_dens)
    theta_e = theta_e_func(r, mass, t_e0, p_temp)
    nu_c = nu_c_func(r, mass, n_th0, t_e0, p_temp, theta_b)
    x = nu / nu_c
    return (n * e ** 2 * nu * synchrotron_func(x) / (2 * np.sqrt(3) * c * theta_e ** 2)).to(
        u.erg / (u.cm ** 3 * u.s * u.Hz))


def specific_intensity(r, nu=kw_nu, mass=kw_mass, n_th0=kw_n_th0, p_dens=kw_p_dens, t_e0=kw_t_e0, p_temp=kw_p_temp,
                       theta_b=kw_theta_b, scale_height=kw_scale_height):
    j_coeff = j_emission_nu_func(r, nu, mass, n_th0, p_dens, t_e0, p_temp, theta_b)
    rg = rg_func(mass)
    size = r * scale_height * rg
    return (size * j_coeff).to(u.erg / u.cm ** 2)


def bright_temp(r, nu=kw_nu, mass=kw_mass, n_th0=kw_n_th0, p_dens=kw_p_dens, t_e0=kw_t_e0, p_temp=kw_p_temp,
                theta_b=kw_theta_b, scale_height=kw_scale_height):
    s_intensity = specific_intensity(r, nu, mass, n_th0, p_dens, t_e0, p_temp, theta_b, scale_height)
    return (c ** 2 / (2 * nu ** 2 * kB) * s_intensity).to(u.K)


'''
def size_func(r, mass, scale_height):
    rg = rg_func(mass)
    return r * scale_height * rg
    
def nu_b_func(r, mass, n_th0, p_dens, beta):
    b_field = b_func(r, mass, n_th0, p_dens, beta)
    return (e * b_field / (2*np.pi * me * c)).to(u.Hz)
    
def nu_func(r, x, mass, n_th0, p_dens, beta, theta_b, gamma):
    #return x * nu_c_func(r, mass, n_th0, p_dens, beta, theta_b, gamma)
    
def x_func(r, nu, mass, n_th0, t_e0, p_dens, p_temp, beta, theta_b):
    return nu / nu_c_func(r, mass, n_th0, t_e0, p_dens, p_temp, beta, theta_b)

'''
