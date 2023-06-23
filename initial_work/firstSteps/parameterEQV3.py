import numpy as np
from astropy import constants as const
from astropy import units as u


def nth_func(r, n_th0, p, mass):
    rb = rb_func(mass)
    return n_th0 * (r/rb) ** p


def te_func(r, t_e0, p, mass):
    rb = rb_func(mass)
    return t_e0 * (r/rb) ** p


def b_func(r, n_th0, p, mass, beta):
    if n_th0.unit != u.cm ** -3:
        raise ValueError('n_th0 needs units of number density')
    nth = nth_func(r*u.cm, n_th0, p, mass)
    return np.sqrt(nth * 8 * np.pi * const.m_p.cgs * const.c.cgs ** 2 / (6 * beta * r))


def rb_func(mass):
    return 20 * const.G.cgs * mass / (const.c.cgs ** 2)