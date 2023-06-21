import numpy as np
from astropy import constants as const
from astropy import units as u
# All units are in cgs


def bright_temp(nu, intensity_nu):
    """Calculate Brightness temperature given photon frequency and specific intensity recieved in far field

    Args:
        nu (Float): Photon frequency
        intensity_nu (Float): Photon intensity

    Returns:
        Float: Brightness temperature
    """  
    return (const.c.cgs ** 2 / (2 * nu ** 2 * const.k_B.cgs) * intensity_nu).to(u.K)


def specific_intensity(j_coeff, mass):
    """Calculate photon specific intensity at far field given emission coefficient and encountered black hole mass

    Args:
        j_coeff (Float): Emission coefficient
        mass (Float): Encountered black hole mass

    Returns:
        Float: Photon specific intensity
    """    
    return (5 * (2*const.G.cgs*mass / (const.c.cgs ** 2)) * j_coeff).to(u.erg / u.cm **2)


def j_emission_nu(nu, n_dens, theta_e, theta_b, b_field, gamma):
    """Calculate photon emission coefficient given frequency and parameters

    Args:
        nu (Float): Photon frequency
        n_dens (Float): Number density of photons
        theta_e (Float): Dimensionless photon energy number
        theta_b (Float): Angle between magnetic field and wave vector
        b_field (Float): Magnetic field strength
        gamma (Float): _description_

    Returns:
        FLoat: Photon emission coefficient
    """      
    nu_b = (const.e.esu * b_field / (2*np.pi * const.m_e.cgs * const.c.cgs)).to(u.Hz)
    nu_c = 3/2*nu_b*np.sin(theta_b)*gamma**2
    x = nu / nu_c
    return ((n_dens * (const.e.esu ** 2) * nu / (2*np.sqrt(3) * const.c.cgs * theta_e ** 2)) \
            * intensity(x)).to(u.erg / (u.cm ** 3 * u.s * u.Hz))


def j_emission_x(x, n_dens, theta_e, theta_b, b_field, gamma):
    """""Calculate photon emission coefficient with dimensionless frequency and parameters

    Args:
        x (Float): Dimensionless frequency variable
        n_dens (Float): Number density of photons
        theta_e (Float): Dimensionless photon energy number
        theta_b (Float): Angle between magnetic field and wave vector
        b_field (Float): Magnetic field strength
        gamma (Float): _description_

    Returns:
        Float: Photon emission coefficient
    """    
    nu_b = (const.e.esu * b_field / (2 * np.pi * const.m_e.cgs * const.c.cgs)).to(u.Hz)
    nu_c = 3 / 2 * nu_b * np.sin(theta_b) * gamma ** 2
    nu = x * nu_c
    return ((n_dens * (const.e.esu ** 2) * nu / (2*np.sqrt(3) * const.c.cgs * theta_e ** 2)) \
            * intensity(x)).to(u.erg / (u.cm ** 3 * u.s * u.Hz))


def intensity(x):
    return 2.5651 * (1 + 1.92 * (x ** (-1/3)) + (0.9977*x**(-2/3)))*np.exp(-1.8899*x**(1/3))
