from matplotlib import pyplot
from numpy import arange, mean, std

x = [99.4183578491, 199.506011963, 299.507415771, 399.524291992, 499.418365479,
     599.378479004, 699.528015137, 799.580505371, 899.462768555, 999.611877441,
     1099.47290039, 1199.52978516]
y = [599.74609375, 599.810424805, 599.979187012, 599.985656738, 600.090393066,
     600.268676758, 600.105102539, 600.241516113, 600.316040039, 600.132446289,
     600.706481934, 600.829528809]

step = 0.02
n_estimated_l = arange(599, 600, 0.02)

estimated_means = []
estimated_stds = []
for idx_n, n_ in enumerate(n_estimated_l):
    y_estimated_l = []
    for idx_x, x_ in enumerate(x):
        y_estimated = x_ * 0.001 + n_
        y_difference = abs(y[idx_x] - y_estimated)
        y_estimated_l.append(y_difference)
    estimated_means.append(mean(y_estimated_l))
    estimated_stds.append(std(y_estimated_l))

plot_size = [16.53, 11.69]
plot_dpi = 100

fig = pyplot.figure(figsize=plot_size, dpi=plot_dpi)
ax_1 = fig.add_subplot(1, 1, 1)
ax_1.plot(n_estimated_l, estimated_means)
ax_1.plot(n_estimated_l, estimated_stds)
pyplot.show()

