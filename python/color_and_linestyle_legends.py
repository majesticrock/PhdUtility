from collections import OrderedDict

def color_labels_dict(labels, ax):
    lines = ax.get_lines()
    # Get unique colors used in the plot
    unique_colors = list(OrderedDict.fromkeys([line.get_color() for line in lines]))
    # Create a dictionary mapping unique colors to custom labels
    color_labels = {unique_colors[i]: labels[i] if i < len(labels) else '' for i in range(len(unique_colors))}
    return color_labels

def linestyle_labels_dict(labels, ax):
    lines = ax.get_lines()
    # Get unique linestyles used in the plot
    unique_linestyles = list(OrderedDict.fromkeys([line.get_linestyle() for line in lines]))
    # Create a dictionary mapping unique linestyles to custom labels
    linestyle_labels = {unique_linestyles[i]: labels[i] if i < len(labels) else '' for i in range(len(unique_linestyles))}
    return linestyle_labels

def color_and_linestyle_legends(ax, color_legend_loc='upper right', linestyle_legend_loc='lower right', linestyle_legend_color='black', color_legend_title='', linestyle_legend_title='', color_labels=None, linestyle_labels=None):
    lines = ax.get_lines()
    # Create dictionaries to store unique colors and linestyles
    unique_colors = OrderedDict()
    unique_linestyles = OrderedDict()

    for line in lines:
        # Extract color and linestyle from each line
        color = line.get_color()
        linestyle = line.get_linestyle()
    
        # Store unique colors and linestyles in dictionaries
        if color not in unique_colors:
            unique_colors[color] = line
        if linestyle not in unique_linestyles:
            unique_linestyles[linestyle] = line

    # Create legends for colors and linestyles with specified locations and default linestyle color
    if color_labels:
        if isinstance(color_labels, dict):
            color_legend = ax.legend(unique_colors.values(), [color_labels.get(color, '') for color in unique_colors.keys()], title=color_legend_title, loc=color_legend_loc)
        else:
            color_legend = ax.legend(unique_colors.values(), [color_labels_dict(color_labels, ax).get(color, '') for color in unique_colors.keys()], title=color_legend_title, loc=color_legend_loc)
    else:
        color_legend = ax.legend(unique_colors.values(), unique_colors.keys(), title=color_legend_title, loc=color_legend_loc)
    
    if linestyle_labels:
        if isinstance(linestyle_labels, dict):
            linestyle_legend = ax.legend(unique_linestyles.values(), [linestyle_labels.get(linestyle, '') for linestyle in unique_linestyles.keys()], title=linestyle_legend_title, loc=linestyle_legend_loc)
        else:
            linestyle_legend = ax.legend(unique_linestyles.values(), [linestyle_labels_dict(linestyle_labels, ax).get(linestyle, '') for linestyle in unique_linestyles.keys()], title=linestyle_legend_title, loc=linestyle_legend_loc)
    else:
        linestyle_legend = ax.legend(unique_linestyles.values(), unique_linestyles.keys(), title=linestyle_legend_title, loc=linestyle_legend_loc)

    # Set the line colors in the linestyle legend
    if hasattr(linestyle_legend, 'legend_handles'):
        for handle in linestyle_legend.legend_handles:
            handle.set_color(linestyle_legend_color)
    else:
        for handle in linestyle_legend.legendHandles:
            handle.set_color(linestyle_legend_color)

    # Add both legends to the axes
    ax.add_artist(color_legend)  # Add color legend back to the axes

    return color_legend, linestyle_legend

## Example usage:
#import matplotlib.pyplot as plt
#plt.logging.getLogger('matplotlib.font_manager').disabled = True
## Create a figure and axes
#fig, ax = plt.subplots()
#
## Plot example lines with different colors and linestyles
#ax.plot([1, 2, 3], [1, 2, 3], color='red', linestyle='-' )
#ax.plot([1, 2, 3], [2, 3, 4], color='red', linestyle='--')
#ax.plot([1, 2, 3], [3, 4, 5], color='red', linestyle=':' )
#
#ax.plot([1, 2, 3], [5, 4, 3], color='green', linestyle='-' )
#ax.plot([1, 2, 3], [4, 3, 2], color='green', linestyle='--')
#ax.plot([1, 2, 3], [3, 2, 1], color='green', linestyle=':' )
#
#ax.plot([1, 2, 3], [1, 1, 1], color='blue', linestyle='-' )
#ax.plot([1, 2, 3], [2, 2, 2], color='blue', linestyle='--')
#ax.plot([1, 2, 3], [3, 3, 3], color='blue', linestyle=':' )
#
## Define custom color and linestyle labels
## Both presented variants are possible
## The first one relies on dicts and the other one on arrays (but creates the dicts internally)
##color_labels = color_labels_dict(['Red Label', 'Green Label', 'Blue Label'], ax)
##linestyle_labels = linestyle_labels_dict(['Solid Line', 'Dashed Line', 'Dotted Line'], ax)
#color_labels = ['Red Label', 'Green Label', 'Blue Label']
#linestyle_labels = ['Solid Line', 'Dashed Line', 'Dotted Line']
#
## Call the function to create legends with optional settings, including custom labels for colors and linestyles
#color_and_linestyle_legends(ax, color_legend_loc='upper right', linestyle_legend_loc='lower right', linestyle_legend_color='black', 
#                            color_legend_title='Color title', linestyle_legend_title='Linestyle title', color_labels=color_labels, linestyle_labels=linestyle_labels)
#
## Show the plot
#plt.show()