import astropy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure


def plot_with_units(x_values, y_values, label=None):
    plt.plot(x_values, y_values,  label=label)
    plt.xlabel("(" + str(x_values[0].unit) + ")")
    plt.ylabel("(" + str(y_values[0].unit) + ")")

