import matplotlib.pyplot as plt


def plot_with_units(x_values, y_values, label=None, xlabel=None ,ylabel=None):
    plt.plot(x_values, y_values,  label=label)
    plt.xlabel(str(xlabel)+ " (" + str(x_values[0].unit) + ")")
    plt.ylabel(str(ylabel)+ " (" + str(y_values[0].unit) + ")")


