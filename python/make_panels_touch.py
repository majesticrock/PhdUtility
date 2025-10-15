import matplotlib.pyplot as plt
import numpy as np

def iter_except(seq, n):
    n = n % len(seq)
    for i, x in enumerate(seq):
        if i != n:
            yield x

def make_panels_touch(fig, axes, touch_x=True, touch_y=True, x_for_constrained=0, y_for_constrained=-1, w_pad=0, h_pad=0, wspace=0, hspace=0):
    """
    Adjusts the spacing between subplots in a matplotlib figure so that panels touch each other
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure containing the axes to adjust.
    axes : numpy.ndarray or list of matplotlib.axes.Axes
        The collection of axes (as returned by plt.subplots) to adjust.
    touch_x : bool, optional
        If True, remove vertical spacing between panels (columns touch). Default is True.
    touch_y : bool, optional
        If True, remove horizontal spacing between panels (rows touch). Default is True.
    x_for_constrained : int, optional
        Index of the column to keep x-axis tick labels visible when touch_y is True. Default is 0.
    y_for_constrained : int, optional
        Index of the row to keep y-axis tick labels visible when touch_x is True. Default is -1 (last row).
    Notes
    -----
    - Temporarily sets the figure's layout engine to 'constrained' to remove spacing, then restores it.
    - Restores the original visibility state of tick labels after adjusting layout.
    - Designed for use with grid-like axes (e.g., from plt.subplots).
    """
    if not hasattr(axes, '__len__'):
        return
    
    if isinstance(axes, list):
        axes = np.array(axes)
    
    tick_state = {}
    for ax in axes.flat:
        tick_state[ax] = {
            "labelleft": ax.yaxis.get_ticklabels()[0].get_visible() if ax.yaxis.get_ticklabels() else False,
            "labelbottom": ax.xaxis.get_ticklabels()[0].get_visible() if ax.xaxis.get_ticklabels() else False,
        }
    
    if hasattr(axes[0], '__len__'):
        xstates = []
        ystates = []
        if touch_x:
            for axs in iter_except(axes, y_for_constrained):
                for ax in axs:
                    ax.tick_params(labelleft=False)
                     
            ax0 = axes[y_for_constrained, 0]
            sharey = all(ax0.get_shared_y_axes().joined(ax0, ax) for ax in axes[y_for_constrained, 1:])
            if sharey:
                ax = ax0
                ystate = {
                        "yticks": ax.get_yticks().copy(),
                        "yticklabels": [lbl.get_text() for lbl in ax.get_yticklabels()],
                        "ylim" : ax.get_ylim()
                    }
                if len(ystate["yticks"]) > 0:
                    idx = len(ystate["yticks"]) // 2
                    if ystate["yticklabels"]:
                        ax.set_yticks([ystate["yticks"][idx]], [ystate["yticklabels"][idx]])
                    else:
                        ax.set_yticks([ystate["yticks"][idx]])
                else:
                    ax.set_yticks([])
            else:
                for ax in axes[y_for_constrained, :]:
                    ystates.append({
                        "yticks": ax.get_yticks().copy(),
                        "yticklabels": [lbl.get_text() for lbl in ax.get_yticklabels()],
                        "ylim" : ax.get_ylim()
                    })
                    ystate = ystates[-1]
                    if len(ystate["yticks"]) > 0:
                        idx = len(ystate["yticks"]) // 2
                        if ystate["yticklabels"]:
                            ax.set_yticks([ystate["yticks"][idx]], [ystate["yticklabels"][idx]])
                        else:
                            ax.set_yticks([ystate["yticks"][idx]])
                    else:
                        ax.set_yticks([])
        ################################################################################
        if touch_y:
            for axs in iter_except(axes.T, x_for_constrained):
                for ax in axs:
                    ax.tick_params(labelbottom=False)
            
            
            ax0 = axes[0, x_for_constrained]
            sharex = all(ax0.get_shared_x_axes().joined(ax0, ax) for ax in axes[1:, x_for_constrained])
            if sharex:
                ax = ax0
                xstate = {
                        "xticks": ax.get_xticks().copy(),
                        "xticklabels": [lbl.get_text() for lbl in ax.get_xticklabels()],
                        "xlim" : ax.get_xlim()
                    }
                if len(xstate["xticks"]) > 0:
                    idx = len(xstate["xticks"]) // 2
                    if xstate["xticklabels"]:
                        ax.set_xticks([xstate["xticks"][idx]], [xstate["xticklabels"][idx]])
                    else:
                        ax.set_xticks([xstate["xticks"][idx]])
                else:
                    ax.set_xticks([])
            else:
                for ax in axes[:, x_for_constrained]:
                    xstates.append({
                        "xticks": ax.get_xticks().copy(),
                        "xticklabels": [lbl.get_text() for lbl in ax.get_xticklabels()],
                        "xlim" : ax.get_xlim()
                    })
                    xstate = xstates[-1]
                    if len(xstate["xticks"]) > 0:
                        idx = len(xstate["xticks"]) // 2
                        if xstate["xticklabels"]:
                            ax.set_xticks([xstate["xticks"][idx]], [xstate["xticklabels"][idx]])
                        else:
                            ax.set_xticks([xstate["xticks"][idx]])
                    else:
                        ax.set_xticks([])
    ################################################################################
    else:
        skip = x_for_constrained if touch_y else y_for_constrained
        for ax in iter_except(axes, skip):
            if touch_y:
                ax.tick_params(labelbottom=False)
            else:
                ax.tick_params(labelleft=False)
        ################################################################################
        if touch_x:
            ax = axes[skip]
            ystate = {
                "yticks": ax.get_yticks().copy(),
                "yticklabels": [lbl.get_text() for lbl in ax.get_yticklabels()],
                "ylim" : ax.get_ylim()
            }
            if len(ystate["yticks"]) > 0:
                idx = len(ystate["yticks"]) // 2
                if ystate["yticklabels"]:
                    ax.set_yticks([ystate["yticks"][idx]], [ystate["yticklabels"][idx]])
                else:
                    ax.set_yticks([ystate["yticks"][idx]])
            else:
                ax.set_yticks([])
        ################################################################################
        if touch_y:
            ax = axes[skip]
            xstate = {
                "xticks": ax.get_xticks().copy(),
                "xticklabels": [lbl.get_text() for lbl in ax.get_xticklabels()],
                "xlim" : ax.get_xlim()
            }
            if len(xstate["xticks"]) > 0:
                idx = len(xstate["xticks"]) // 2
                if xstate["xticklabels"]:
                    ax.set_xticks([xstate["xticks"][idx]], [xstate["xticklabels"][idx]])
                else:
                    ax.set_xticks([xstate["xticks"][idx]])
            else:
                ax.set_xticks([])
    ################################################################################
    
    fig.set_layout_engine('constrained')
    fig.get_layout_engine().set(w_pad=w_pad, h_pad=h_pad, hspace=hspace, wspace=wspace)
    fig.canvas.draw()
    fig.set_layout_engine('none')
    
    # restore previous state
    for ax, state in tick_state.items():
        ax.tick_params(labelleft=state["labelleft"], labelbottom=state["labelbottom"])
    
    if touch_y: 
        if hasattr(axes[0], '__len__'):
            if sharex:
                if xstate["xticklabels"]:
                    ax0.set_xticks(xstate["xticks"], xstate["xticklabels"])
                else:
                    ax0.set_xticks(xstate["xticks"])
                ax0.set_xlim(xstate["xlim"])
            else:
                for xstate, ax in zip(reversed(xstates), reversed(axes[:, x_for_constrained])):
                    if xstate["xticklabels"]:
                        ax.set_xticks(xstate["xticks"], xstate["xticklabels"])
                    else:
                        ax.set_xticks(xstate["xticks"])
                    ax.set_xlim(xstate["xlim"])
        else:
            axes[skip].set_xticks(xstate["xticks"], xstate["xticklabels"])
            axes[skip].set_xlim(xstate["xlim"])
    if touch_x: 
        if hasattr(axes[0], '__len__'):
            if sharey:
                if ystate["yticklabels"]:
                    ax0.set_yticks(ystate["yticks"], ystate["yticklabels"])
                else:
                    ax0.set_yticks(ystate["yticks"])
                ax0.set_ylim(ystate["ylim"])
            else:
                for ystate, ax in zip(reversed(ystates), reversed(axes[y_for_constrained, :])):
                    if ystate["yticklabels"]:
                        ax.set_yticks(ystate["yticks"], ystate["yticklabels"])
                    else:
                        ax.set_yticks(ystate["yticks"])
                    ax.set_ylim(ystate["ylim"])
        else:
            axes[skip].set_yticks(ystate["yticks"], ystate["yticklabels"])
            axes[skip].set_ylim(ystate["ylim"])


if __name__ == '__main__':
    plt.rcParams['font.size'] = 20
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['ytick.minor.visible'] = True

    # =========================================================
    # Original figure (2×1)
    # =========================================================
    X, Y = np.meshgrid(np.linspace(-np.pi, np.pi, 50, endpoint=False),
                       np.linspace(0, 1, 50, endpoint=False))
    Z_up  = np.sin(X + Y)
    Z_low = np.sin(X * Y)

    fig1, axes1 = plt.subplots(nrows=2, sharex=True)
    pc_up  = axes1[0].pcolormesh(X, Y, Z_up,  cmap="seismic", shading="auto")
    pc_low = axes1[1].pcolormesh(X, Y, Z_low, cmap="seismic", shading="auto")

    cbar = fig1.colorbar(pc_up, ax=axes1)
    cbar.set_label("$f(X, Y)$")

    axes1[-1].set_xlabel("X")
    axes1[0].set_ylabel("Addition\n$Y$")
    axes1[1].set_ylabel("Multiplication\n$Y$")

    make_panels_touch(fig1, axes1, touch_x=True, touch_y=False)
    #fig1.savefig("test1.pdf")

    # =========================================================
    # X and Y swapped (1×2 layout)
    # =========================================================
    Y, X = np.meshgrid(np.linspace(-np.pi, np.pi, 50, endpoint=False),
                       np.linspace(0, 1, 50, endpoint=False))
    Z_left  = np.sin(X + Y)
    Z_right = np.sin(X * Y)

    fig2, axes2 = plt.subplots(ncols=2, sharey=True)
    pc_left  = axes2[0].pcolormesh(X, Y, Z_left,  cmap="seismic", shading="auto")
    pc_right = axes2[1].pcolormesh(X, Y, Z_right, cmap="seismic", shading="auto")

    cbar = fig2.colorbar(pc_left, ax=axes2)
    cbar.set_label("$f(X, Y)$")

    axes2[0].set_xlabel("Addition\n$X$")
    axes2[1].set_xlabel("Multiplication\n$X$")
    axes2[0].set_ylabel("Y")

    make_panels_touch(fig2, axes2, touch_x=False, touch_y=True)
    #fig2.savefig("test2.pdf")

    # =========================================================
    # 3×2 grid of different functions
    # =========================================================
    X, Y = np.meshgrid(np.linspace(0, 2 * np.pi, 80, endpoint=False),
                       np.linspace(0, 1.5, 80, endpoint=False))

    Z_funcs = [
        np.sin(X + Y),
        np.cos(X * Y),
        np.sin(X) * np.cos(Y),
        np.sin(X * Y**2),
        np.exp(-X**2 - Y**2),
        np.sin(3*X) * np.sin(3*Y)
    ]

    fig3, axes3 = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8, 6))
    for ax, Z in zip(axes3.flat, Z_funcs):
        pc = ax.pcolormesh(X, Y, Z, cmap="seismic", shading="auto")

    cbar = fig3.colorbar(pc, ax=axes3, shrink=0.8)
    cbar.set_label("$f(X, Y)$")

    axes3[-1, 0].set_xlabel("X")
    axes3[-1, 1].set_xlabel("X")
    axes3[0, 0].set_ylabel("Y")
    axes3[1, 0].set_ylabel("Y")
    axes3[2, 0].set_ylabel("Y")

    axes3[0, -1].set_xticks([0, 2, 3.2, 5.1])

    make_panels_touch(fig3, axes3, touch_x=True, touch_y=True)
    #fig3.savefig("test3.pdf")
    
    plt.show()