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




# ------------------------------
# Reverse versions
# ------------------------------
mrock_red_r = ListedColormap(mrock_red.colors[::-1], name="mrock_red_r")
mrock_blue_r = ListedColormap(mrock_blue.colors[::-1], name="mrock_blue_r")

# ------------------------------
# Diverging colormap: Red -> White -> Blue
# ------------------------------
def create_diverging_from_existing(red_cmap, blue_cmap, n=256, name="mrock_diverging"):
    red_colors = red_cmap.colors[::2]
    blue_colors = blue_cmap.colors[::2]
    return ListedColormap(np.vstack([red_colors, blue_colors]), name=name)

mrock_diverging = create_diverging_from_existing(mrock_red_r, mrock_blue)
mrock_diverging_r = ListedColormap(mrock_diverging.colors[::-1], name="mrock_diverging_r")

# ------------------------------
# Quick test plot
# ------------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack([gradient, gradient])
    
    fig, axes = plt.subplots(5, 1, figsize=(8,5))
    
    axes[0].imshow(gradient, aspect='auto', cmap=mrock_red)
    axes[0].set_title("mrock_red")
    axes[0].axis("off")
    
    axes[1].imshow(gradient, aspect='auto', cmap=mrock_red_r)
    axes[1].set_title("mrock_red_r")
    axes[1].axis("off")
    
    axes[2].imshow(gradient, aspect='auto', cmap=mrock_blue)
    axes[2].set_title("mrock_blue")
    axes[2].axis("off")
    
    axes[3].imshow(gradient, aspect='auto', cmap=mrock_blue_r)
    axes[3].set_title("mrock_blue_r")
    axes[3].axis("off")
    
    axes[4].imshow(gradient, aspect='auto', cmap=mrock_diverging)
    axes[4].set_title("mrock_diverging")
    axes[4].axis("off")
    
    plt.tight_layout()
    plt.show()