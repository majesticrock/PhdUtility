import numpy as np
from matplotlib.colors import ListedColormap
import colorspacious as cs

def perceptual_colormap(rgb_colors, name="custom_lab", n=256):
    """Create a perceptually uniform colormap interpolated in LAB space."""
    lab_colors = np.array([cs.cspace_convert(c, "sRGB1", "CIELab") for c in rgb_colors])
    t = np.linspace(0, 1, n)
    lab_interp = np.array([
        np.interp(t, np.linspace(0, 1, len(lab_colors)), lab_colors[:, i])
        for i in range(3)
    ]).T
    rgb = cs.cspace_convert(lab_interp, "CIELab", "sRGB1")
    rgb = np.clip(rgb, 0, 1)
    return ListedColormap(rgb, name=name)

# ------------------------------
# Define base colormaps
# ------------------------------

# Red colormap
mrock_red = perceptual_colormap(
    [
        (1, 1, 1),       # white
        (1.0, 0.65, 0.1),# slightly warmer orange
        (1.0, 0.2, 0.0), # scarlet
        (0.55, 0.0, 0.0),# dark red
        (0, 0, 0)        # black
    ],
    name="mrock_red"
)

# Blue colormap
mrock_blue = perceptual_colormap(
    [
        (1, 1, 1),       # white
        (0.35, 0.75, 1.0), # light cyan
        (0.15, 0.45, 0.85),  # true blue
        (0.05, 0.1, 0.7),  # rich dark blue
        (0, 0, 0)          # black
    ],
    name="mrock_blue"
)

mrock_green = perceptual_colormap(
    [
        (1.0, 1.0, 1.0),         # white
        (0.9, 0.97, 0.55),       # pale yellow-green
        (0.70, 0.88, 0.10),      # vivid yellow-green (lighter than #84b819)
        (0.25, 0.70, 0.15),      # bright true green
        (0.0, 0.50, 0.35),       # deep turquoise-green (cooler)
        (0.0, 0.25, 0.2),        # dark teal/forest
        (0.0, 0.0, 0.0)          # black
    ],
    name="mrock_green"
)

mrock_purple = perceptual_colormap(
    [
        (1.0, 1.0, 1.0),        # white
        (0.95, 0.80, 1.0),      # light lavender
        (0.75, 0.45, 0.95),     # bright magenta-violet
        (0.45, 0.15, 0.80),     # vivid purple (#4D19B8 region)
        (0.25, 0.05, 0.55),     # indigo-purple
        (0.08, 0.0, 0.25),      # deep blue-violet (almost black)
        (0.0, 0.0, 0.0)         # black
    ],
    name="mrock_purple"
)

mrock_orange = perceptual_colormap(
    [
        (1.0, 1.0, 1.0),       # white
        (1.0, 0.85, 0.55),     # pale orange-yellow
        (1.0, 0.5, 0.06),      # #ff7f0e bright orange
        (0.45, 0.2, 0.0),      # dark burnt orange
        (0.0, 0.0, 0.0)        # black
    ],
    name="mrock_orange"
)


# ------------------------------
# Reverse versions
# ------------------------------
mrock_red_r = mrock_red.reversed()
mrock_blue_r = mrock_blue.reversed()
mrock_green_r = mrock_green.reversed()
mrock_purple_r = mrock_purple.reversed()
mrock_orange_r = mrock_orange.reversed()

# ------------------------------
# Diverging colormap: Red -> White -> Blue
# ------------------------------
def create_diverging_from_existing(red_cmap, blue_cmap, n=256, name="mrock_diverging"):
    red_colors = red_cmap.colors[::2]
    blue_colors = blue_cmap.colors[::2]
    return ListedColormap(np.vstack([red_colors, blue_colors]), name=name)

mrock_diverging = create_diverging_from_existing(mrock_red_r, mrock_blue)
mrock_diverging_r = mrock_diverging.reversed()

mrock_dark_diverging = create_diverging_from_existing(mrock_red, mrock_blue_r)
mrock_dark_diverging_r = mrock_dark_diverging.reversed()

mrock_tu_diverging = create_diverging_from_existing(mrock_purple_r, mrock_green)
mrock_tu_diverging_r = mrock_tu_diverging.reversed()

mrock_tu_dark_diverging = create_diverging_from_existing(mrock_purple, mrock_green_r)
mrock_tu_dark_diverging_r = mrock_tu_dark_diverging.reversed()

# ------------------------------
# Quick test plot
# ------------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack([gradient, gradient])
    
    fig, axes = plt.subplots(9, 1, figsize=(8,8))
    
    axes[0].imshow(gradient, aspect='auto', cmap=mrock_red)
    axes[0].set_title("mrock_red")
    axes[0].axis("off")
    
    axes[1].imshow(gradient, aspect='auto', cmap=mrock_blue)
    axes[1].set_title("mrock_blue")
    axes[1].axis("off")
     
    axes[2].imshow(gradient, aspect='auto', cmap=mrock_diverging)
    axes[2].set_title("mrock_diverging")
    axes[2].axis("off")

    axes[3].imshow(gradient, aspect='auto', cmap=mrock_dark_diverging)
    axes[3].set_title("mrock_dark_diverging")
    axes[3].axis("off")
    
    axes[4].imshow(gradient, aspect='auto', cmap=mrock_green)
    axes[4].set_title("mrock_green")
    axes[4].axis("off")
    
    axes[5].imshow(gradient, aspect='auto', cmap=mrock_purple)
    axes[5].set_title("mrock_purple")
    axes[5].axis("off")
    
    axes[6].imshow(gradient, aspect='auto', cmap=mrock_orange)
    axes[6].set_title("mrock_orange")
    axes[6].axis("off")
    
    axes[7].imshow(gradient, aspect='auto', cmap=mrock_tu_diverging)
    axes[7].set_title("mrock_tu_diverging")
    axes[7].axis("off")

    axes[8].imshow(gradient, aspect='auto', cmap=mrock_tu_dark_diverging)
    axes[8].set_title("mrock_tu_dark_diverging")
    axes[8].axis("off")
    
    plt.tight_layout()
    plt.show()