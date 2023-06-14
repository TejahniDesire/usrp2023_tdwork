import tools
import numpy as np
x = np.linspace(1, 2, 10)
x2 = x
plot = tools.StdPlot()

plot.plot(x, x2)
# plot.plot(x, np.exp(x))
# plot.set_axes("x", "y")
