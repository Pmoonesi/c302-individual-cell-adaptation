import numpy as np
import pandas as pd
import re

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#took this from c302

def get_muscle_position(muscle):
    if muscle == "MANAL" or muscle == "MVULVA":
        return 0, 0, 0

    pat1 = r"M([VD])([LR])(\d+)"
    pat2 = r"([VD])([LR])(\d+)"
    md = re.fullmatch(pat1, muscle)
    if not md:
        md = re.fullmatch(pat2, muscle)

    if md:
        dv = md.group(1)
        lr = md.group(2)
        idx = md.group(3)
        x = 80 * (1 if lr == "L" else -1)
        z = 80 * (-1 if dv == "V" else 1)
        y = -300 + 30 * int(idx) # bruh this is based on index already
        return x, y, z

    raise Exception("Unrecognized muscle name format %s" % muscle)



def plot_muscle_kymograph(
    dorsal_matrix,
    ventral_matrix,
    dorsal_names,
    ventral_names,
    time_axis=None,
    title="C. elegans Muscle Activity",
    convert_to_mV=True,
    save_path=None,
):
    """
    Plot kymographs for dorsal and ventral muscle activity.
 
    Parameters
    ----------
    dorsal_matrix  : np.ndarray, shape (n_timesteps, n_dorsal_muscles)
    ventral_matrix : np.ndarray, shape (n_timesteps, n_ventral_muscles)
    dorsal_names   : list of str, muscle names in anterior->posterior order
    ventral_names  : list of str, muscle names in anterior->posterior order
    time_axis      : np.ndarray or None — if None uses sample index
    convert_to_mV  : bool — multiply by 1000 if data is in SI volts
    save_path      : str or None
    """
 
    # --- transpose to (n_muscles, n_timesteps) for imshow ---
    D = dorsal_matrix.T
    V = ventral_matrix.T
 
    if convert_to_mV:
        D = D * 1000
        V = V * 1000
 
    if time_axis is None:
        time_axis = np.arange(dorsal_matrix.shape[0])
 
    # --- shared colour scale across both panels ---
    vmin = np.percentile(np.concatenate([D, V]), 2)
    vmax = np.percentile(np.concatenate([D, V]), 98)
 
    # --- layout ---
    fig = plt.figure(figsize=(16, 8), facecolor="#111118")
    fig.suptitle(title, color="white", fontsize=13, y=0.98)
 
    gs = gridspec.GridSpec(
        2, 2,
        width_ratios=[1, 0.025],
        height_ratios=[1, 1],
        hspace=0.08,
        wspace=0.03,
        left=0.08, right=0.92,
        top=0.93, bottom=0.08,
    )
 
    ax_d  = fig.add_subplot(gs[0, 0])   # dorsal heatmap
    ax_v  = fig.add_subplot(gs[1, 0])   # ventral heatmap
    ax_cb = fig.add_subplot(gs[:, 1])   # shared colorbar
 
    # --- helper to draw one heatmap ---
    def draw_kymograph(ax, matrix, names, xlabel=False):
        im = ax.imshow(
            matrix,
            aspect="auto",
            cmap="RdBu_r",
            vmin=vmin, vmax=vmax,
            interpolation="nearest",
            origin="upper",
            extent=[time_axis[0], time_axis[-1], len(names) - 0.5, -0.5],
        )
        ax.set_facecolor("#0a0a12")
        ax.tick_params(colors="white", labelsize=8)
        ax.spines[:].set_color("#333")
 
        # y ticks — muscle names
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, color="#cccccc", fontsize=7)
 
        if xlabel:
            ax.set_xlabel("Time (0.1ms)" if time_axis[-1] > 10 else "Timesteps",
                          color="white", fontsize=10)
        else:
            ax.set_xticklabels([])
 
        return im
 
    im = draw_kymograph(ax_d, D, dorsal_names,  xlabel=False)
    draw_kymograph(ax_v, V, ventral_names, xlabel=True)
 
    # --- labels ---
    ax_d.set_ylabel("Dorsal muscles\n(ant → post)",
                    color="white", fontsize=9)
    ax_v.set_ylabel("Ventral muscles\n(ant → post)",
                    color="white", fontsize=9)
 
    # --- colorbar ---
    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("Membrane Potential (mV)", color="white", fontsize=9)
    cb.ax.yaxis.set_tick_params(color="white", labelsize=8)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")
    cb.outline.set_edgecolor("#444")
 
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"Saved: {save_path}")
 
    plt.show()
    return fig
 
 
# ── example usage ────────────────────────────────────────────────────────────
if __name__ == "__main__":

    # DL, DR, VL (23), VR
    columns = ["timestep"]

    columns.extend([f"MDL{i:02d}" for i in range(1,25)])
    columns.extend([f"MDR{i:02d}" for i in range(1,25)])
    columns.extend([f"MVL{i:02d}" for i in range(1,24)])
    columns.extend([f"MVR{i:02d}" for i in range(1,25)])

    df = pd.read_table("c302_I0_Full_Interactome.muscles.dat", sep='\t', names=columns, index_col=False, usecols=range(96))

    dorsal_muscles = columns[1:49]
    ventral_muscles = columns[49:]

    sorted_dorsal_muscles = sorted(dorsal_muscles, key=lambda x: get_muscle_position(x)[1])
    sorted_ventral_muscles = sorted(ventral_muscles, key=lambda x: get_muscle_position(x)[1])

    dorsal_muscles_matrix = df.loc[:, sorted_dorsal_muscles].to_numpy()
    ventral_muscles_matrix = df.loc[:, sorted_ventral_muscles].to_numpy()
 
    # replace these with your actual data
    # dorsal_matrix  shape: (n_timesteps, n_dorsal_muscles)
    # ventral_matrix shape: (n_timesteps, n_ventral_muscles)
 
    plot_muscle_kymograph(
        dorsal_matrix=dorsal_muscles_matrix,    # your variable
        ventral_matrix=ventral_muscles_matrix,  # your variable
        dorsal_names=sorted_dorsal_muscles,      # your list
        ventral_names=sorted_ventral_muscles,    # your list
        time_axis=None,              # your time array in ms, or None
        convert_to_mV=True,             # set False if already in mV
        save_path="muscle_kymograph.png",
    )
 
