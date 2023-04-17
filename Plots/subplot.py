import matplotlib.pyplot as plt

# generate your plots here
plot_a = plt.plot([1,2,3],[4,5,6])
plot_b = plt.plot([1,2,3],[3,2,1])
plot_c = plt.plot([1,2,3],[6,5,4])
plot_d = plt.plot([1,2,3],[1,2,3])
plot_e = plt.plot([1,2,3],[2,2,2])

# create a figure and subplots
fig, axs = plt.subplots(1, 5, figsize=(15,3))

# add each plot to the corresponding subplot
axs[0].plot(plot_a)
axs[1].plot(plot_b)
axs[2].plot(plot_c)
axs[3].plot(plot_d)
axs[4].plot(plot_e)

# add labels to the figure
fig.text(0.05, 0.95, "a", fontsize=14, fontweight='bold')
fig.text(0.25, 0.95, "b", fontsize=14, fontweight='bold')
fig.text(0.45, 0.95, "c", fontsize=14, fontweight='bold')
fig.text(0.65, 0.95, "d", fontsize=14, fontweight='bold')
fig.text(0.85, 0.95, "e", fontsize=14, fontweight='bold')

# display the figure
plt.show()