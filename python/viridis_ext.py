import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import colorspacious as cs

def extended_colormap(cmap_name, n=256, extension=0.2):
    n_base = int(n / (1 + extension))
    base = plt.get_cmap(cmap_name, n_base)
    colors_base_rgb = base(np.linspace(0, 1, n_base))[:, :3]

    last_rgb = colors_base_rgb[-1]
    last_lab = cs.cspace_convert(last_rgb, "sRGB1", "CIELab")

    white_lab = np.array([100.0, 0.0, 0.0])

    n_ext = n - n_base
    labs_ext = np.array([
        (1 - t) * last_lab + t * white_lab
        for t in np.linspace(0, 1, n_ext)
    ])
    rgb_ext = np.array([cs.cspace_convert(lab, "CIELab", "sRGB1") for lab in labs_ext])
    rgb_ext = np.clip(rgb_ext, 0, 1)

    colors = np.vstack([colors_base_rgb, rgb_ext])
    alphas = np.ones((colors.shape[0], 1))
    colors = np.hstack([colors, alphas])

    cmap = LinearSegmentedColormap.from_list(f"{cmap_name}_ext", colors, N=n)
    cmap_r = cmap.reversed()

    return cmap, cmap_r

viridis_ext, viridis_ext_r = extended_colormap("viridis")
inferno_ext, inferno_ext_r = extended_colormap("inferno")

if __name__ == "__main__":
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack([gradient, gradient])

    fig, axes = plt.subplots(nrows=4, figsize=(6,4))
    axes[0].imshow(gradient, aspect='auto', cmap=viridis_ext)
    axes[0].set_title("viridis_ext")
    axes[0].axis("off")

    axes[1].imshow(gradient, aspect='auto', cmap=viridis_ext_r)
    axes[1].set_title("viridis_ext_r")
    axes[1].axis("off")

    axes[2].imshow(gradient, aspect='auto', cmap=inferno_ext)
    axes[2].set_title("inferno_ext")
    axes[2].axis("off")

    axes[3].imshow(gradient, aspect='auto', cmap=inferno_ext_r)
    axes[3].set_title("inferno_ext_r")
    axes[3].axis("off")

    plt.tight_layout()
    plt.show()
