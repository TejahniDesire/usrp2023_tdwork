import numpy as np
from astropy import constants as const
from astropy import units as u

# All r in units of Rg

G = const.G.cgs
c = const.c.cgs
me = const.m_e.cgs
mp = const.m_p.cgs
e = const.e.esu
kB = const.k_B.cgs


def rg_func(mass):
    return G * mass / c ** 2


def rb_func(mass):
    return 20 * G * mass / (c ** 2)


def size_func(r, mass, scale_height):
    rg = rg_func(mass)
    return r * scale_height * rg


def te_func(r, mass, t_e0, p_temp):
    if t_e0.unit != u.K:
        raise ValueError('n_th0 needs units of Kelvin')
    rb = rb_func(mass)
    rg = rg_func(mass)
    return t_e0 * (r*rg/rb) ** p_temp


def theta_e_func(r, mass, t_e0, p_temp):
    temp = te_func(r, mass, t_e0, p_temp)
    return (kB * temp / (me * c ** 2)).to(u.dimensionless_unscaled)


def nth_func(r, mass, n_th0, p_dens):
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    rg = rg_func(mass)
    rb = rb_func(mass)
    return n_th0 * (r * rg / rb) ** p_dens


def b_func(r, mass, n_th0, p_dens, beta):
    rg = rg_func(mass)
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    nth = nth_func(r, mass, n_th0, p_dens)
    return np.sqrt(nth * 8 * np.pi * mp * c ** 2 / (6 * beta * r))


def nu_b_func(r, mass, n_th0, p_dens, beta):
    b_field = b_func(r, mass, n_th0, p_dens, beta)
    return (e * b_field / (2*np.pi * me * c)).to(u.Hz)


def nu_c_func(r, mass, n_th0, p_dens, beta, theta_b, gamma):
    nu_b = nu_b_func(r, mass, n_th0, p_dens, beta)
    return 3/2*nu_b*np.sin(theta_b)*gamma**2


def x_func(r, nu, mass, n_th0, p_dens, beta, theta_b, gamma):
    return nu / nu_c_func(r, mass, n_th0, p_dens, beta, theta_b, gamma)


def nu_func(r, x, mass, n_th0, p_dens, beta, theta_b, gamma):
    return x * nu_c_func(r, mass, n_th0, p_dens, beta, theta_b, gamma)


def synchrotron_func(x):
    return 2.5651 * (1 + 1.92 * (x ** (-1/3)) + (0.9977*x**(-2/3)))*np.exp(-1.8899*x**(1/3))


def j_emission_x_func(r, x, mass, n_th0, p_dens, t_e0, p_temp, beta, theta_b, gamma):
    n = nth_func(r, mass, n_th0, p_dens)
    theta_e = theta_e_func(r, mass, t_e0, p_temp)
    nu = nu_func(r, x, mass, n_th0, p_dens, beta, theta_b, gamma)
    return (n * e ** 2 * nu * synchrotron_func(x) / (2 * np.sqrt(3) * c * theta_e ** 2)).to(u.erg / (u.cm ** 3 * u.s * u.Hz))


def specific_intensity(r, x, mass, n_th0, p_dens, t_e0, p_temp, beta, theta_b, gamma, scale_height):  
    j_coeff = j_emission_x_func(r, x, mass, n_th0, p_dens, t_e0, p_temp, beta, theta_b, gamma)
    size = size_func(r, mass, scale_height)
    return (size * j_coeff).to(u.erg / u.cm **2)


def bright_temp(r, x, mass, n_th0, p_dens, t_e0, p_temp, beta, theta_b, gamma, scale_height):
    s_intensity = specific_intensity(r, x, mass, n_th0, p_dens, t_e0, p_temp, beta, theta_b, gamma, scale_height)
    nu = nu_func(r, x, mass, n_th0, p_dens, beta, theta_b, gamma)
    return (const.c.cgs ** 2 / (2 * nu ** 2 * const.k_B.cgs) * s_intensity).to(u.K)

