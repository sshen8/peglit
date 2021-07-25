import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import linkage_tree
from .utils import fill_bpp

def plot_bpp(seq, basepair_probs, seq_idx, seq_names=("Spacer", "Scaffold", "Template", "PBS", "Linker")):
    bpp = fill_bpp(basepair_probs)
    ticklabels_minor = tuple(seq)
    ticks_minor = np.arange(len(ticklabels_minor)) + 0.5
    fig = plt.figure(figsize=(8, 6.5))
    gs = fig.add_gridspec(ncols=2, nrows=2, width_ratios=[1, 0.03], height_ratios=[0.3, 0.7])
    ax = fig.add_subplot(gs[:, 0])
    subplt_cbar = fig.add_subplot(gs[0, 1])
    ax.set_aspect("equal")
    c = ax.pcolormesh(bpp, cmap="Reds", vmin=0., vmax=1.)
    ax.invert_yaxis()
    ax.tick_params(axis="both", which="major", length=0)
    ax.tick_params(axis="y", which="major", labelrotation=90)
    ax.set_xticks(ticks_minor, minor=True)
    ax.set_xticklabels(ticklabels_minor, minor=True, fontdict={'fontsize': 6})
    ax.set_yticks(ticks_minor, minor=True)
    ax.set_yticklabels(ticklabels_minor, minor=True, fontdict={'fontsize': 6})
    ticks_major = []
    ticklabels_major = []
    for subseq_idx, subseq_name in zip(seq_idx, seq_names):
        if subseq_idx.stop - subseq_idx.start == 0:
            continue
        ax.axvline(subseq_idx.stop, color="k", linewidth=0.5)
        ax.axhline(subseq_idx.stop, color="k", linewidth=0.5)
        ticks_major.append((subseq_idx.start + subseq_idx.stop) / 2. + 1e-2)
        ticklabels_major.append(subseq_name)
    ax.set_xticks(ticks_major)
    ax.set_xticklabels(["\n{:s}".format(tk) for tk in ticklabels_major], fontdict={'fontsize': 10, 'horizontalalignment': 'center'})
    ax.set_yticks(ticks_major)
    ax.set_yticklabels(["{:s}\n".format(tk) for tk in ticklabels_major], fontdict={'fontsize': 10, 'verticalalignment': 'center'})
    fig.colorbar(c, cax=subplt_cbar, label="Base pair probability")
    plt.tight_layout()

def plot_clusters(linker_good, linker_feats, n_clusters):
    # Create linkage matrix and then plot the dendrogram
    children, _, _, _, distances = linkage_tree(linker_feats, n_clusters=n_clusters, linkage="complete", return_distance=True)
    counts = np.zeros(children.shape[0])
    for i, merge in enumerate(children):
        current_count = 0
        for child_idx in merge:
            if child_idx < len(linker_good):
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - len(linker_good)]
        counts[i] = current_count
    linkage_matrix = np.column_stack([children, distances, counts]).astype(float)
    fig = plt.figure(figsize=(10, 7.5))
    gs = fig.add_gridspec(ncols=3, nrows=2, width_ratios=[0.1, 0.9, 0.03], height_ratios=[0.3, 0.7])
    subplt1 = fig.add_subplot(gs[:, 0])
    subplt2 = fig.add_subplot(gs[:, 1])
    subplt_cbar = fig.add_subplot(gs[0, 2])
    plt.sca(subplt1)
    dn = dendrogram(linkage_matrix, labels=linker_good, orientation="left", color_threshold=distances[-n_clusters + 1])
    subplt1.set_xticks([])
    subplt1.invert_yaxis()
    plt.box(on=None)
    cmap = plt.cm.get_cmap("Blues", np.max(linker_feats))
    c = subplt2.pcolormesh(linker_feats[dn["leaves"], :][:, dn["leaves"]],
                           vmin=-0.5, # offset so cbar ticks centered
                           vmax=np.max(linker_feats) - 0.5,
                           cmap=cmap)
    subplt2.axis("off")
    subplt2.invert_yaxis()
    plt.colorbar(c, cax=subplt_cbar, label="Edit distance")
    plt.tight_layout()
