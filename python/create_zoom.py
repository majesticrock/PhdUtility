import numpy as np

# ax                    -   Axes to plot on, needs to be an object; no plt wrapper
# inset_x/ypos, 
#   inset_width/height  -   Position and size of the inset in relative coordinates (0 to 1)
# x/ylim                -   x and y limit of the zoomed region, expects tuples (lower limit, upper limit)
# y_funcs               -   Functions generating the y_data from x_data. May only take 1 argument. Need to be in order
#                               If it is None, the data is extracted from the plot
# skip_lines            -   If some lines should be skipped, you can pass the corresponding indizes in an array
# x/yticklabels         -   Default behaviour: "None" matplotlib default behaviour
# **kwargs              -   Is forwarded to the creation of the inset axes
def create_zoom(ax, inset_xpos, inset_ypos, inset_width, inset_height, xlim=(0, .1), ylim=(0, .1), y_funcs=None, skip_lines=[], mark_inset=True, xticks=None, yticks=None, xticklabels=None, yticklabels=None, **kwargs):
    axins = ax.inset_axes([inset_xpos, inset_ypos, inset_width, inset_height], xlim=xlim, ylim=ylim, **kwargs)
    if xticks is not None:
        axins.set_xticks(xticks, xticklabels)
    if yticks is not None:
        axins.set_yticks(yticks, yticklabels)
        
    if mark_inset: 
        ax.indicate_inset_zoom(axins, edgecolor="black")
    
    i=0
    for j, line in enumerate(ax.lines):
        if j in skip_lines:
            continue
        
        if y_funcs is None:
            # Extract the data from the line
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            new_line, = axins.plot(x_data, y_data)
        else:
            x_data = np.linspace(xlim[0], xlim[1], 100)
            new_line, = axins.plot(x_data, (y_funcs[i])(x_data))
        i+=1
        
        # Create a new object for the inset plot
        # I do not know how this works, but it does
        new_line.update_from(line)
        new_line.remove()
        new_line.axes = axins
        new_line.set_transform(axins.transData)
        axins.add_line(new_line)
    return axins
        
## Usage example
#import matplotlib.pyplot as plt
#import numpy as np
#x = np.linspace(0, 10, 100)
#fig, ax = plt.subplots()
#
#ax.plot(x, np.sin(x),           linestyle="-.", label="sin", color="orange", marker="v", markevery=4)
#ax.plot(x, np.cos(x),           linestyle="--", label="cos", color="purple", marker="o", markevery=3)
#ax.plot(x, np.sin(np.sin(x)),   linestyle="-", label="sin(sin)", color="green", marker="X", markevery=5)
#create_zoom(ax, inset_xpos=0.45, inset_ypos=0.25, inset_width=0.25, inset_height=0.35, xlim=(0.65, 0.95), ylim=(0.6, 0.8))
#
#ax.legend()
#ax.set_xlabel("$x$")
#ax.set_ylabel("$y$")
#
#fig.tight_layout()
#plt.show()