"""Data-backed teaching views of GRETA's Literature tasks.

These views reconstruct the native gretapy 0.0.3 task definitions. They do not
change benchmark scores. A missing reference annotation is not a biological
negative, and TF pairs inferred from target overlap are not directed TF edges.
"""

import gretapy as gt
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control, hypergeom


COLORS = {"Pearson": "#2A6F97", "CollecTRI": "#D97732", "Random": "#777777"}
THRESHOLD = 0.01


def _color(name):
    return COLORS.get(name, "#6B5B95")


def _reference_markers(result, database):
    genes = set(result["adata"].var_names)
    db = gt.ds.read_db("hg38", database, verbose=False)
    db = db.loc[db[0].astype(str).isin(genes)]
    terms = result.get("run_info", {}).get("manifest", {}).get("terms", {})
    cats = terms.get(database)
    if cats is not None:
        # This is the same regular-expression filtering used by native _tfm.
        db = db.loc[db[1].str.contains("|".join(cats))]
    return set(db[0].astype(str))


def _tidy_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def plot_marker_example(result):
    """Show source-TF annotation counts and a small actual membership table.

    Returns ``(figure, tables)``. ``counts`` contains the native TP/FP/FN
    bookkeeping, using descriptive labels instead of interpreting reference
    absence as a verified negative. ``membership`` contains the full union of
    measured reference markers and network regulators; ``displayed_membership``
    records the deterministic illustrative subset shown in the figure.
    """
    networks = result["networks"]
    names = list(networks)
    databases = ["HPA", "TF-Marker"]
    markers = {db: _reference_markers(result, db) for db in databases}
    sources = {name: set(frame.source.astype(str)) for name, frame in networks.items()}
    counts = []
    for database in databases:
        reference = markers[database]
        for name in names:
            regulators = sources[name]
            n_listed = len(regulators & reference)
            counts.append(dict(
                network=name, database=database, network_tfs=len(regulators),
                reference_tfs=len(reference), listed_tfs=n_listed,
                not_in_reference=len(regulators - reference),
                reference_tfs_not_in_network=len(reference - regulators),
                precision=n_listed / len(regulators) if regulators else np.nan,
                recall=n_listed / len(reference) if reference else np.nan,
            ))
    counts = pd.DataFrame(counts)
    universe = sorted(set().union(*sources.values(), *markers.values()))
    membership = pd.DataFrame({"TF": universe})
    for name in names:
        membership[f"GRN: {name}"] = membership.TF.isin(sources[name])
    for database in databases:
        membership[f"Reference: {database}"] = membership.TF.isin(markers[database])
    network_columns = [f"GRN: {name}" for name in names]
    reference_columns = [f"Reference: {db}" for db in databases]
    # Include both annotation agreements and omissions, rather than displaying
    # only positive examples. Selection uses names, never network task scores.
    signature = membership[network_columns + reference_columns].astype(int).astype(str).agg("".join, axis=1)
    selected = (membership.assign(_signature=signature)
                .groupby("_signature", sort=True, group_keys=False)
                .head(2).sort_values(["_signature", "TF"], ascending=[False, True])
                .drop(columns="_signature").head(18))
    columns = network_columns + reference_columns
    fig, (ax_counts, ax_membership) = plt.subplots(
        1, 2, figsize=(14.8, 6.8), gridspec_kw={"width_ratios": [1.2, 1.25]},
        layout="constrained",
    )
    positions, labels = [], []
    y = 0
    for database in databases:
        for name in names:
            row = counts.loc[counts.database.eq(database) & counts.network.eq(name)].iloc[0]
            ax_counts.barh(y, row.listed_tfs, color=_color(name), height=0.72)
            ax_counts.barh(y, row.not_in_reference, left=row.listed_tfs,
                           color="#e3e3e3", height=0.72)
            ax_counts.text(row.network_tfs + 3, y, f"{row.listed_tfs}/{row.network_tfs}", va="center", fontsize=10)
            positions.append(y)
            labels.append(f"{database} · {name}")
            y += 1
        y += 0.7
    max_count = max(counts.network_tfs.max(), 1)
    ax_counts.set_xlim(0, max_count * 1.28)
    ax_counts.set_yticks(positions, labels)
    ax_counts.invert_yaxis()
    ax_counts.set_xlabel("Number of TFs used as GRN sources")
    ax_counts.set_title("A  TF identity determines this task", loc="left", weight="bold", pad=15)
    _tidy_axis(ax_counts)
    ax_counts.legend(handles=[Patch(color="#777777", label="Listed in reference"),
                              Patch(color="#e3e3e3", label="Not listed in reference")],
                     loc="lower left", bbox_to_anchor=(0, -0.22), frameon=False, fontsize=9)

    values = selected[columns].to_numpy(dtype=bool)
    rgba = np.ones((*values.shape, 4))
    for j, column in enumerate(columns):
        color = _color(names[j]) if j < len(names) else "#52666D"
        rgba[values[:, j], j] = to_rgba(color)
        rgba[~values[:, j], j] = to_rgba("#f0f0f0")
    ax_membership.imshow(rgba, aspect="auto", interpolation="nearest")
    for i, j in zip(*np.where(values), strict=True):
        ax_membership.text(j, i, "●", color="white", ha="center", va="center", fontsize=9)
    ax_membership.set_yticks(range(len(selected)), selected.TF)
    ax_membership.set_xticks(range(len(columns)), names + databases, rotation=35, ha="right")
    ax_membership.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax_membership.set_yticks(np.arange(-0.5, len(selected), 1), minor=True)
    ax_membership.grid(which="minor", color="white", linewidth=2)
    ax_membership.tick_params(which="both", length=0)
    ax_membership.axvline(len(names) - 0.5, color="#777777", linewidth=1.2)
    ax_membership.set_title("B  Actual TF membership examples", loc="left", weight="bold", pad=15)
    ax_membership.set_xlabel("GRN source lists                 Reference marker lists", labelpad=12)
    fig.suptitle("TF Markers checks the identities of included regulators", fontsize=15, weight="bold")
    fig.text(0.50, -0.025,
             "Filled cell = present in that list. Reference absence is unannotated, not a verified biological negative.\n"
             "Examples: at most two alphabetic TFs per membership pattern; counts use every TF.",
             ha="center", fontsize=9, color="#444444")
    return fig, {"counts": counts, "membership": membership,
                 "displayed_membership": selected.reset_index(drop=True)}


def _pair_statistics(frame, reference_tfs):
    """Reconstruct native _find_pairs after _tfp's reference-TF restriction.

    A hypergeometric upper tail is algebraically the one-sided Fisher exact
    p value. BH is performed over *all* eligible TF pairs, including pairs with
    no overlap (p=1), exactly as native gretapy 0.0.3 does.
    """
    graph = frame[["source", "target"]].drop_duplicates()
    graph = graph.loc[graph.source.isin(reference_tfs)]
    tfs = sorted(graph.source.unique())
    targets = sorted(graph.target.unique())
    columns = ["tf_a", "tf_b", "targets_a", "targets_b", "shared_targets",
               "background_targets", "pval", "padj", "significant"]
    if len(tfs) < 2:
        return pd.DataFrame(columns=columns)
    row = pd.Categorical(graph.source, categories=tfs).codes
    col = pd.Categorical(graph.target, categories=targets).codes
    incidence = sparse.csr_matrix((np.ones(len(graph), dtype=np.int64), (row, col)),
                                  shape=(len(tfs), len(targets)))
    sizes = np.asarray(incidence.sum(axis=1)).ravel()
    overlap = (incidence @ incidence.T).toarray()
    a, b = np.triu_indices(len(tfs), k=1)
    shared = overlap[a, b]
    pvals = hypergeom.sf(shared - 1, len(targets), sizes[a], sizes[b])
    qvals = false_discovery_control(pvals, method="bh")
    return pd.DataFrame(dict(
        tf_a=np.asarray(tfs)[a], tf_b=np.asarray(tfs)[b],
        targets_a=sizes[a], targets_b=sizes[b], shared_targets=shared,
        background_targets=len(targets), pval=pvals, padj=qvals,
        significant=qvals < THRESHOLD,
    ))


def _q_label(value):
    if pd.isna(value):
        return "not tested"
    if value == 0:
        return "< 1e-308"
    return f"{value:.2g}"


def plot_tf_pair_example(result, pair=None, database="IntAct"):
    """Plot actual shared targets, a native pair test, and reference agreement.

    By default select a reference-supported pair significant in the most GRNs;
    break ties by total shared-target count, then alphabetic pair name. This is
    an illustrative positive example, not a representative performance claim.
    ``pair=(tf_a, tf_b)`` selects any explicit pair. The returned tables expose
    all pair tests (including BH's full family), the selected pair, all its
    target memberships, and whole-network reference-overlap counts.
    """
    if database not in {"IntAct", "Europe PMC"}:
        raise ValueError("TF Pairs examples support IntAct or Europe PMC.")
    networks = result["networks"]
    names = list(networks)
    genes = set(result["adata"].var_names)
    db = gt.ds.read_db("hg38", database, verbose=False)
    db = db.loc[db[0].astype(str).isin(genes) & db[1].astype(str).isin(genes)]
    reference = {tuple(sorted((str(a), str(b)))) for a, b in zip(db[0], db[1], strict=True)}
    reference_tfs = set(db[0].astype(str)) | set(db[1].astype(str))
    tests = []
    for name, frame in networks.items():
        table = _pair_statistics(frame, reference_tfs)
        table.insert(0, "network", name)
        table["in_reference"] = [tuple(p) in reference for p in table[["tf_a", "tf_b"]].to_numpy()]
        tests.append(table)
    all_tests = pd.concat(tests, ignore_index=True)
    selection_rule = "Explicit user-selected pair"
    if pair is None:
        candidates = (all_tests.loc[all_tests.in_reference]
                      .groupby(["tf_a", "tf_b"], as_index=False)
                      .agg(n_significant=("significant", "sum"),
                           shared_total=("shared_targets", "sum")))
        candidates = candidates.loc[candidates.n_significant > 0]
        if candidates.empty:
            raise ValueError("No reference-supported significant pair; supply pair=(tf_a, tf_b) to inspect a specific pair.")
        chosen = candidates.sort_values(["n_significant", "shared_total", "tf_a", "tf_b"],
                                         ascending=[False, False, True, True]).iloc[0]
        pair = (chosen.tf_a, chosen.tf_b)
        selection_rule = ("Reference-supported pair with most significant GRNs, then largest total shared-target count, "
                          "then alphabetic order; illustrative positive example")
    if len(pair) != 2 or pair[0] == pair[1]:
        raise ValueError("Supply two different TF identifiers.")
    tf_a, tf_b = sorted(str(tf) for tf in pair)
    selected = all_tests.loc[all_tests.tf_a.eq(tf_a) & all_tests.tf_b.eq(tf_b)].copy()
    selected = selected.set_index("network").reindex(names).reset_index()
    selected["tf_a"], selected["tf_b"] = tf_a, tf_b
    selected["in_reference"] = (tf_a, tf_b) in reference
    selected["significant"] = selected.significant.fillna(False).astype(bool)
    selected["database"] = database
    selected["selection_rule"] = selection_rule

    target_sets = {}
    for name, frame in networks.items():
        target_sets[(name, tf_a)] = set(frame.loc[frame.source.eq(tf_a), "target"])
        target_sets[(name, tf_b)] = set(frame.loc[frame.source.eq(tf_b), "target"])
    union = sorted(set().union(*target_sets.values()))
    if not union:
        raise ValueError("Neither selected TF has targets in the supplied networks.")
    memberships = pd.DataFrame({"target": union})
    columns = []
    shared_columns = []
    for name in names:
        for tf in (tf_a, tf_b):
            col = f"{name}: {tf}"
            columns.append(col)
            memberships[col] = memberships.target.isin(target_sets[(name, tf)])
        col = f"shared: {name}"
        shared_columns.append(col)
        memberships[col] = memberships[columns[-2:]].all(axis=1)
    memberships["n_networks_shared"] = memberships[shared_columns].sum(axis=1)
    memberships["n_memberships"] = memberships[columns].sum(axis=1)
    displayed = (memberships.sort_values(["n_networks_shared", "n_memberships", "target"],
                                         ascending=[False, False, True]).head(18))
    counts = []
    for name in names:
        table = all_tests.loc[all_tests.network.eq(name)]
        significant = table.loc[table.significant]
        supported = int(significant.in_reference.sum())
        counts.append(dict(network=name, database=database,
                           tested_pairs=len(table), significant_pairs=len(significant),
                           in_reference=supported, not_in_reference=len(significant) - supported,
                           reference_pairs=len(reference), reference_not_recovered=len(reference) - supported,
                           precision=supported / len(significant) if len(significant) else 0.,
                           recall=supported / len(reference) if len(reference) else 0.))
    counts = pd.DataFrame(counts)

    fig = plt.figure(figsize=(15, 7.4), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.5, 1.2], height_ratios=[1, 1])
    ax_tile = fig.add_subplot(grid[:, 0])
    ax_pair = fig.add_subplot(grid[0, 1])
    ax_all = fig.add_subplot(grid[1, 1])
    values = displayed[columns].to_numpy(dtype=bool)
    rgba = np.ones((*values.shape, 4))
    for j in range(len(columns)):
        rgba[values[:, j], j] = to_rgba(_color(names[j // 2]))
        rgba[~values[:, j], j] = to_rgba("#f0f0f0")
    ax_tile.imshow(rgba, aspect="auto", interpolation="nearest")
    ax_tile.set_yticks(range(len(displayed)), displayed.target)
    ax_tile.set_xticks(range(len(columns)), [tf for _ in names for tf in (tf_a, tf_b)], rotation=40, ha="right")
    ax_tile.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax_tile.set_yticks(np.arange(-0.5, len(displayed), 1), minor=True)
    ax_tile.grid(which="minor", color="white", linewidth=2)
    ax_tile.tick_params(which="both", length=0)
    for j, name in enumerate(names):
        ax_tile.text(2 * j + 0.5, -0.9, name, ha="center", va="bottom", color=_color(name), weight="bold")
        if j:
            ax_tile.axvline(2 * j - 0.5, color="#888888", linewidth=1)
    ax_tile.set_title("A  Which targets are shared?", loc="left", pad=35, weight="bold")
    ax_tile.set_xlabel(f"Filled cell = TF–target edge; {len(displayed)} of {len(union)} targets shown\n"
                       "Rows favor shared targets, then alphabetic order", labelpad=10)

    positions = np.arange(len(names))
    ax_pair.barh(positions, selected.shared_targets.fillna(0), color=[_color(n) for n in names], height=0.65)
    ax_pair.set_yticks(positions, names)
    ax_pair.invert_yaxis()
    max_shared = max(selected.shared_targets.fillna(0).max(), 1)
    ax_pair.set_xlim(0, max_shared * 1.9)
    for i, row in selected.iterrows():
        shared = int(row.shared_targets) if pd.notna(row.shared_targets) else 0
        text = f"{shared} shared; q={_q_label(row.padj)}"
        if pd.notna(row.background_targets):
            text += f"\nTF degrees {int(row.targets_a)}/{int(row.targets_b)}; background {int(row.background_targets)}"
        ax_pair.text(shared + max_shared * .025, i, text, va="center", fontsize=8.5)
    ax_pair.set_xlabel("Number of shared targets")
    ax_pair.set_title(f"B  {tf_a} + {tf_b}: native pair test", loc="left", weight="bold", pad=12)
    _tidy_axis(ax_pair)

    ax_all.barh(positions, counts.in_reference, color=[_color(n) for n in names], height=0.65)
    ax_all.barh(positions, counts.not_in_reference, left=counts.in_reference, color="#dedede", height=0.65)
    ax_all.set_yticks(positions, names)
    ax_all.invert_yaxis()
    limit = max(counts.significant_pairs.max(), 1)
    ax_all.set_xlim(0, limit * 1.28)
    for i, row in counts.iterrows():
        ax_all.text(row.significant_pairs + limit * .025, i,
                    f"{row.in_reference}/{row.significant_pairs}", va="center", fontsize=9)
    ax_all.set_xlabel("All pairs with BH q < 0.01")
    ax_all.set_title(f"C  Which significant pairs are listed in {database}?", loc="left", weight="bold", pad=12)
    _tidy_axis(ax_all)
    ax_all.legend(handles=[Patch(color="#777777", label="Listed in reference"),
                           Patch(color="#dedede", label="Not listed in reference")],
                  loc="lower left", bbox_to_anchor=(0, -0.44), fontsize=8.5, frameon=False)
    ref_label = "listed" if (tf_a, tf_b) in reference else "not listed"
    fig.suptitle(f"TF Pairs: shared targets → significant overlap → reference comparison\n"
                 f"Illustrative pair {tf_a} + {tf_b} is {ref_label} in {database}; this is not a directed TF→TF edge",
                 fontsize=13.5, weight="bold")
    return fig, {"selected_pair": selected, "pair_tests": all_tests, "counts": counts,
                 "target_membership": memberships, "displayed_targets": displayed.reset_index(drop=True)}
