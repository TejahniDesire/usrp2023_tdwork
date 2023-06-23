import matplotlib.pyplot as plt
from astropy import units as u

def plot_with_units(x_values, y_values, label=None, xlabel='' ,ylabel=''):
    plt.plot(x_values, y_values,  label=label)
    if str(x_values[0].unit) != u.dimensionless_unscaled:
        plt.xlabel(str(xlabel)+ " (" + str(x_values[0].unit) + ")")
    else:
        plt.xlabel(str(xlabel))
        
    if str(y_values[0].unit) != u.dimensionless_unscaled:
        plt.ylabel(str(ylabel)+ " (" + str(y_values[0].unit) + ")")
    else:
        plt.ylabel(str(ylabel))


