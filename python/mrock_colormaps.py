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

def colormap_dL_dh(begin, dL, dh, name="dLdh", n=256, fade_to_white=16, fade_to_black=16, chroma_curve=None):
    if chroma_curve is None:
        colors = [ np.array([begin[0] + t * dL, begin[1], begin[2] + t * dh]) 
                    for t in np.linspace(0, 1, n - fade_to_white - fade_to_black)]
    else:
        colors = [ np.array([begin[0] + t * dL, begin[1] * chroma, begin[2] + t * dh]) 
                    for t, chroma in zip(np.linspace(0, 1, n - fade_to_white - fade_to_black), chroma_curve)]
    if fade_to_white > 0:
        lab_white = np.array([100, 0, colors[-1][2] + dh * fade_to_white / n])
        transition_colors = np.array([np.interp(
                np.linspace(0, 1, fade_to_white),
                np.array([0., 1.]),
                np.array([colors[-1][i], lab_white[i]]))
            for i in range(3)]).T
        colors.extend(transition_colors)
        
    if fade_to_black > 0:
        lab_black = np.array([0, 0, colors[0][2] - dh * fade_to_black / n])
        transition_colors = np.array([np.interp(
                np.linspace(0, 1, fade_to_black),
                np.array([0., 1.]),
                np.array([colors[0][i], lab_black[i]]))
             for i in range(3)]).T
        colors.reverse()
        colors.extend(transition_colors)
        colors.reverse()
        
    rgb = cs.cspace_convert(colors, "CIELCh", "sRGB1")
    rgb = np.clip(rgb, 0, 1)
    return ListedColormap(rgb, name=name)

def create_diverging_from_existing(first_cmap, second_cmap, name="mrock_diverging", skip_first=1):
    red_colors = first_cmap.colors[::skip_first]
    blue_colors = second_cmap.colors
    return ListedColormap(np.vstack([red_colors, blue_colors]), name=name)

# ------------------------------
# Define base colormaps
# ------------------------------


# Red colormap
mrock_red = perceptual_colormap(
    [
        (1, 1, 1),       # white
        cs.cspace_convert((75, 80, 80), "CIELCh", "sRGB1"),# slightly warmer orange
        cs.cspace_convert((50, 100, 45), "CIELCh", "sRGB1"), # scarlet
        cs.cspace_convert((25, 80, 10), "CIELCh", "sRGB1"),# dark red
        (0, 0, 0)        # black
    ],
    name="mrock_red"
)

# Blue colormap
mrock_blue = perceptual_colormap(
    [
        (1, 1, 1),       # white
        cs.cspace_convert((75, 70, 200), "CIELCh", "sRGB1"), # light cyan
        cs.cspace_convert((50, 80, 250), "CIELCh", "sRGB1"),   # true blue
        cs.cspace_convert((25, 90, 300), "CIELCh", "sRGB1"), # rich dark blue
        (0, 0, 0)          # black
    ],
    name="mrock_blue"
)

mrock_green = perceptual_colormap([
        (1, 1, 1),
        cs.cspace_convert((75, 80, 120), "CIELCh", "sRGB1"),
        cs.cspace_convert((50, 90, 170), "CIELCh", "sRGB1"), 
        cs.cspace_convert((25, 80, 220), "CIELCh", "sRGB1"), 
        (0,0,0)
    ],
    name="mrock_green"
)



# ------------------------------
# Reverse versions
# ------------------------------
mrock_red_r = mrock_red.reversed()
mrock_blue_r = mrock_blue.reversed()
mrock_green_r = mrock_green.reversed()

mrock_diverging = create_diverging_from_existing(mrock_red_r, mrock_blue)
mrock_diverging_r = mrock_diverging.reversed()

mrock_dark_diverging = create_diverging_from_existing(mrock_red, mrock_blue_r)
mrock_dark_diverging_r = mrock_dark_diverging.reversed()

mrock_diverging_low_center = create_diverging_from_existing(
    mrock_blue, 
    perceptual_colormap([
        (0,0,0),
        cs.cspace_convert((22, 100, 350), "CIELCh", "sRGB1"),
        cs.cspace_convert((44, 100, 22), "CIELCh", "sRGB1"),
        cs.cspace_convert((66, 100, 54), "CIELCh", "sRGB1"),
        cs.cspace_convert((88, 100, 86), "CIELCh", "sRGB1"),
    ], n=256), skip_first=2)
mrock_diverging_low_center_r = mrock_diverging_low_center.reversed()

mrock_diverging_low_center_green = create_diverging_from_existing(
    mrock_green, 
    perceptual_colormap([
        (0,0,0),
        cs.cspace_convert((25, 80, 350), "CIELCh", "sRGB1"),
        cs.cspace_convert((50, 75, 22), "CIELCh", "sRGB1"),
        cs.cspace_convert((75, 75, 54), "CIELCh", "sRGB1"),
        cs.cspace_convert((100, 60, 86), "CIELCh", "sRGB1"),
    ], n=256), skip_first=2)
mrock_diverging_low_center_green_r = mrock_diverging_low_center.reversed()

# ------------------------------
# Quick test plot
# ------------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack([gradient, gradient])
    
    cmaps = [ 
             ("mrock_red", mrock_red),
             ("mrock_blue", mrock_blue),
             ("mrock_green", mrock_green),
             ("mrock_diverging", mrock_diverging),
             ("mrock_dark_diverging", mrock_dark_diverging),
             ("mrock_diverging_low_center", mrock_diverging_low_center),
             ("mrock_diverging_low_center_green", mrock_diverging_low_center_green)
        ]
    
    fig, axes = plt.subplots(len(cmaps), 1, figsize=(8,8))
    
    for ax, cmap in zip(axes, cmaps):
        ax.imshow(gradient, aspect='auto', cmap=cmap[1])
        ax.set_title(cmap[0])
        ax.axis("off")

    plt.tight_layout()
    plt.show()