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

def create_diverging_from_existing(first_cmap, second_cmap, name="diverging", skip_first=1):
    red_colors = first_cmap.colors[::skip_first]
    blue_colors = second_cmap.colors
    return ListedColormap(np.vstack([red_colors, blue_colors]), name=name)

# ------------------------------
# Define base colormaps
# ------------------------------


# Red colormap
red = perceptual_colormap(
    [
        (1, 1, 1),       # white
        cs.cspace_convert((75, 78, 80), "CIELCh", "sRGB1"),# slightly warmer orange
        cs.cspace_convert((50, 89, 45), "CIELCh", "sRGB1"), # scarlet
        cs.cspace_convert((25, 48, 10), "CIELCh", "sRGB1"),# dark red
        (0, 0, 0)        # black
    ],
    name="red"
)

# Blue colormap
blue = perceptual_colormap(
    [
        (1, 1, 1),       # white
        cs.cspace_convert((75, 44, 200), "CIELCh", "sRGB1"), # light cyan
        cs.cspace_convert((50, 40, 250), "CIELCh", "sRGB1"),   # true blue
        cs.cspace_convert((25, 108, 300), "CIELCh", "sRGB1"), # rich dark blue
        (0, 0, 0)          # black
    ],
    name="blue"
)

green = perceptual_colormap([
        (1, 1, 1),
        cs.cspace_convert((75, 83, 120), "CIELCh", "sRGB1"),
        cs.cspace_convert((50, 38, 170), "CIELCh", "sRGB1"), 
        cs.cspace_convert((25, 20, 220), "CIELCh", "sRGB1"), 
        (0,0,0)
    ],
    name="green"
)



# ------------------------------
# Reverse versions
# ------------------------------
red_r = red.reversed()
blue_r = blue.reversed()
green_r = green.reversed()

diverging = create_diverging_from_existing(red_r, blue)
diverging_r = diverging.reversed()

dark_diverging = create_diverging_from_existing(red, blue_r)
dark_diverging_r = dark_diverging.reversed()

diverging_low_center = create_diverging_from_existing(
    blue, 
    perceptual_colormap([
        (0,0,0),
        cs.cspace_convert((22, 45, 350), "CIELCh", "sRGB1"),
        cs.cspace_convert((44, 74, 22), "CIELCh", "sRGB1"),
        cs.cspace_convert((66, 85, 54), "CIELCh", "sRGB1"),
        cs.cspace_convert((88, 63, 86), "CIELCh", "sRGB1"),
    ], n=256), skip_first=2)
diverging_low_center_r = diverging_low_center.reversed()

diverging_low_center_green = create_diverging_from_existing(
    green, 
    perceptual_colormap([
        (0,0,0),
        cs.cspace_convert((25, 49, 350), "CIELCh", "sRGB1"),
        cs.cspace_convert((50, 82, 22), "CIELCh", "sRGB1"),
        cs.cspace_convert((75, 55, 54), "CIELCh", "sRGB1"),
        cs.cspace_convert((100, 40, 86), "CIELCh", "sRGB1"),
    ], n=256), skip_first=2)
diverging_low_center_green_r = diverging_low_center.reversed()

blackidis = perceptual_colormap([
        cs.cspace_convert((100, 40, 60), "CIELCh", "sRGB1"),
        cs.cspace_convert((80, 89, 120), "CIELCh", "sRGB1"),
        cs.cspace_convert((60, 40, 180), "CIELCh", "sRGB1"),
        cs.cspace_convert((40, 31, 240), "CIELCh", "sRGB1"), 
        cs.cspace_convert((20, 80, 300), "CIELCh", "sRGB1"), 
        (0,0,0)
    ], n=256)
blackidis_r = blackidis.reversed()

# ------------------------------
# Quick test plot
# ------------------------------
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack([gradient, gradient])
    
    cmaps = [ 
             ("red", red),
             ("blue", blue),
             ("green", green),
             ("diverging", diverging),
             ("dark_diverging", dark_diverging),
             ("diverging_low_center", diverging_low_center),
             ("diverging_low_center_green", diverging_low_center_green),
             ("blackidis", blackidis)
        ]
    
    fig, axes = plt.subplots(len(cmaps), 1, figsize=(8,8))
    
    for ax, cmap in zip(axes, cmaps):
        ax.imshow(gradient, aspect='auto', cmap=cmap[1])
        ax.set_title(cmap[0])
        ax.axis("off")

    plt.tight_layout()
    plt.show()