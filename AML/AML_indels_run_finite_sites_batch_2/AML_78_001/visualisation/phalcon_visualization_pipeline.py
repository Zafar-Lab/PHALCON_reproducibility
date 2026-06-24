import argparse
import os
import re
import graphviz as gv_lib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.cluster.hierarchy as sch
import seaborn as sns

from ete3 import Tree
import sys

sys.setrecursionlimit(10000)

def parse_args():
    p = argparse.ArgumentParser(
        prog="phalcon_visualization_pipeline",
        description="PHALCON heatmap + UMAP + tree with genotype-merged cluster colors.",
    )
    p.add_argument("-d", "--data",       type=str,   default="AML_78_001")
    p.add_argument("-s", "--seed",       type=int,   default=17215)
    p.add_argument("--min_dist",         type=float, default=0.2)
    p.add_argument("--n_neighbors",      type=int,   default=15)
    p.add_argument("--outdir",           type=str,   default="/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/AML_78_001/visualisation/")
    p.add_argument("--cluster_labels_file",                type=str, default='AML_78_001_indels_finite_sites_inferred_cluster_labels.txt')
    p.add_argument("--vcf_file",                           type=str, default='AML_78_001_indels_nonzero_dbsnp_finite_sites_outputInference.vcf')
    p.add_argument("--variant_function_file",              type=str, default='AML_78_001_indels_dbsnp_non_zero_finite_sites.avinput.variant_function')
    p.add_argument("--lklhd_file",                         type=str, default='AML_78_001_indels_finite_sites_final_lklhds.tsv')
    p.add_argument("--newick_file",                        type=str, default='AML_78_001_indels_finite_sites_inferred_tree.nw')
    p.add_argument("--final_df_file",         type=str, default='AML_78_001_indels_finite_sites_final_df.tsv')
    p.add_argument("--genotype_config_file",               type=str, default='Genotype configuration.tsv')
    return p.parse_args()


# step 0 – genotype matrix loading
def build_genotype_matrix(vcf_file: str, variant_function_file: str) -> pd.DataFrame:
    final_vcf        = np.loadtxt(vcf_file, dtype="object")
    final_vcf        = pd.DataFrame(final_vcf).to_numpy()
    variant_function = pd.read_csv(variant_function_file, sep="\t", header=None)
    variant_function.loc[variant_function[1].str.contains("EZH2",  case=False, na=False), 1] = 'EZH2'
    variant_function.loc[variant_function[1].str.contains("U2AF1",  case=False, na=False), 1] = 'U2AF1'
    variant_function = variant_function.to_numpy()


    n_variants = final_vcf.shape[0]
    n_cells    = final_vcf.shape[1] - 9

    genotype_arr = np.zeros((n_variants, n_cells + 3), dtype="object")
    for i in range(n_variants):
        genotype_arr[i][0] = final_vcf[i][0]
        genotype_arr[i][1] = variant_function[i][1] + "_" + str(final_vcf[i][1])
        genotype_arr[i][2] = variant_function[i][0]
        for j in range(n_cells):
            gt = final_vcf[i][9 + j].split(":")[0]
            if   gt == "0/1": genotype_arr[i][j + 3] = 1
            elif gt == "0/0": genotype_arr[i][j + 3] = 0
            else: print(f"  Warning: unexpected GT '{gt}' at variant {i}, cell {j}")

    df = pd.DataFrame(genotype_arr)
    df = df[~df[2].isin(["intronic", "intergenic", "ncRNA_intronic"])]
    df.drop([0, 2], axis=1, inplace=True)
    df.columns = range(df.shape[1])
    df.set_index([0], inplace=True)
    df.index.name = None
    df.columns = range(df.shape[1])
    df = df.astype("int")
    return df.T   # cells × variants


# step 1 – cluster merging

def compute_merged_labels(df_genotype: pd.DataFrame, original_labels: np.ndarray) -> tuple:
    unique_orig = sorted(set(original_labels.tolist()))
    cluster_vectors = {}
    for cid in unique_orig:
        cell_idx = np.where(original_labels == cid)[0]
        cluster_vectors[cid] = tuple(df_genotype.iloc[cell_idx[0]].tolist())

    vector_to_clusters: dict = {}
    for cid, vec in cluster_vectors.items():
        vector_to_clusters.setdefault(vec, []).append(cid)

    original_to_merged = {}
    merge_report       = {}
    for vec, cids in vector_to_clusters.items():
        merged_id = min(cids)
        merge_report[merged_id] = sorted(cids)
        for cid in cids:
            original_to_merged[cid] = merged_id

    merged_labels = np.array([original_to_merged[c] for c in original_labels])
    return merged_labels, original_to_merged, merge_report


def print_merge_report(merge_report: dict, data: str, outdir: str) -> None:
    print("\n[MERGE] Genotype-based cluster merging:")
    rows = []
    for merged_id, orig_ids in sorted(merge_report.items()):
        status = "unchanged" if len(orig_ids) == 1 else f"merged from {orig_ids}"
        print(f"  Merged cluster {merged_id:3d}  ←  original clusters {orig_ids}  [{status}]")
        rows.append({"merged_cluster": merged_id, "original_clusters": str(orig_ids), "status": status})
    path = os.path.join(outdir, f"{data}_merge_report.tsv")
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    print(f"  Report saved: {path}")


# step 2 – color palette

_DISTINCT_COLORS_HEX = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
    "#dcbeff", "#9A6324", "#fffac8", "#800000", "#aaffc3",
    "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#ffffff",
]

def _hex_to_rgb255(h: str) -> tuple:
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

def build_palette(merged_labels: np.ndarray) -> dict:
    unique_merged = sorted(set(merged_labels.tolist()))
    n = len(unique_merged)
    if n <= len(_DISTINCT_COLORS_HEX):
        colors_255 = [_hex_to_rgb255(h) for h in _DISTINCT_COLORS_HEX[:n]]
    else:
        tab20 = sns.color_palette("tab20", n)
        colors_255 = [(int(r*255), int(g*255), int(b*255)) for r, g, b in tab20]
    return {cid: colors_255[idx] for idx, cid in enumerate(unique_merged)}

def rgb255_to_hex(rgb: tuple) -> str:
    return "#{:02x}{:02x}{:02x}".format(*rgb)

def rgb255_to_01(rgb: tuple) -> tuple:
    return tuple(v / 255.0 for v in rgb)

# heatmap

def plot_heatmap(data, df_genotype, merged_labels, palette, outdir):
    print("\n[HEATMAP] Plotting …")

    df = df_genotype.copy()
    df["Cluster"] = merged_labels

    mutation_sum  = df.drop(columns=["Cluster"]).sum(axis=1)
    no_mut_cells  = df[mutation_sum == 0]
    mutated_cells = df[mutation_sum > 0]

    print(f"  Cells with no mutations : {len(no_mut_cells)}")
    print(f"  Cells with mutations    : {len(mutated_cells)}")

    df_clustered  = mutated_cells.groupby("Cluster").first()
    linkage       = sch.linkage(df_clustered, method="ward")
    dendro        = sch.dendrogram(linkage, no_plot=True)
    cluster_order = [df_clustered.index[i] for i in dendro["leaves"]]

    mutated_sorted = mutated_cells.set_index("Cluster").loc[cluster_order].reset_index()

    # reset_index so strip_colors[i] aligns with clustermap column i
    combined      = pd.concat([no_mut_cells, mutated_sorted], axis=0).reset_index(drop=True)
    df_final      = combined.drop(columns=["Cluster"]).astype("int")

    
    cluster_strip = combined["Cluster"]

    cmap = mcolors.LinearSegmentedColormap.from_list("", ["#D8D3CD", "#792A03"])
    ax = sns.clustermap(
    df_final.T,
    method="ward",
    col_cluster=True,
    row_cluster=True,
    cmap=cmap,
    figsize=(18, 6),
    xticklabels=False,
    cbar_pos=(0.35, -0.01, 0.30, 0.03)
)
    ax.cax.clear()

    ax.fig.colorbar(
        ax.ax_heatmap.collections[0],
        cax=ax.cax,
        orientation="horizontal"
    )
# make colorbar horizontal
    

    ax.cax.xaxis.set_ticks_position("bottom")
    ax.ax_heatmap.yaxis.set_ticks_position("left")

    strip_colors = [rgb255_to_01(palette[c]) for c in cluster_strip]
    reordered    = ax.dendrogram_col.reordered_ind
    strip_colors = [strip_colors[i] for i in reordered]

    ax1 = ax.ax_heatmap
    for i, color in enumerate(strip_colors):
        ax1.add_patch(plt.Rectangle((i, -0.5), 1, 0.5, color=color, clip_on=False))

    unique_in_plot = sorted(set(cluster_strip.tolist()))
    patches = [
        mpatches.Patch(color=rgb255_to_01(palette[c]), label=f"Cluster {c}")
        for c in unique_in_plot
    ]
    ax.ax_heatmap.legend(
        handles=patches, title="Cluster",
        bbox_to_anchor=(1.01, 0.7), loc="upper left", borderaxespad=0, fontsize=8,
    )
    ax.ax_row_dendrogram.set_visible(False)
    ax.ax_col_dendrogram.set_visible(False)

    for ext in ("svg", "png"):
        path = os.path.join(outdir, f"{data}_heatmap.{ext}")
        ax.fig.savefig(path, bbox_inches="tight", dpi=150)
        print(f"  Saved: {path}")
    plt.close("all")


# umap

def plot_umap(data, lklhd_file, merged_labels, palette, seed, min_dist, n_neighbors, outdir):
    print("\n[UMAP] Computing embedding …")

    data_df = pd.read_csv(lklhd_file, sep="\t", header=None)
    cluster = merged_labels.astype(int)

    adata = sc.AnnData(X=data_df)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.umap(adata, min_dist=min_dist, random_state=seed)

    unique_merged = sorted(set(cluster.tolist()))
    adata.obs["batch"] = pd.Categorical(
        cluster.astype(str), categories=[str(c) for c in unique_merged],
    )
    adata.uns["batch_colors"] = [rgb255_to_hex(palette[c]) for c in unique_merged]

    fig = sc.pl.umap(adata, color="batch", show=False, return_fig=True)
    path = os.path.join(outdir, f"{data}_umap.svg")
    fig.savefig(path, bbox_inches="tight")
    print(f"  Saved: {path}")
    plt.close("all")


# clonal tree

def _decrease_numbers(text: str) -> str:
    # Decrement only the number after 'c' prefix (c1→0, c2→1 …)
    # then strip the 'c' so leaf names become plain 0-based integers.
    # Branch lengths, chromosome names, genomic positions are untouched.
    text = re.sub(r"(?<=c)(\d+)", lambda m: str(int(m.group(0)) - 1), text)
    text = re.sub(r"c(\d+)", r"\1", text)
    return text


def _findPar(tree_nw, split_list):
    for node1 in tree_nw.traverse("preorder"):
        des = list(node1.iter_descendants())
        ans = list(node1.get_ancestors())
        for node2 in tree_nw.traverse("preorder"):
            if len(ans) == 0: break
            if node1 != node2 and node2 not in ans:
                if node2 not in des and node1.up != node2.up:
                    m1 = {node1.name} if node1.is_leaf() else set(node1.get_leaf_names())
                    m2 = {node2.name} if node2.is_leaf() else set(node2.get_leaf_names())
                    m  = m1 | m2
                    if set(np.array(list(m), dtype="int")) == {i for i,v in enumerate(split_list) if v==1}:
                        return node1, node2


def _findLoss(tree_nw, split_list):
    for node1 in tree_nw.traverse("preorder"):
        des = list(node1.iter_descendants())
        ans = list(node1.get_ancestors())
        for node2 in tree_nw.traverse("preorder"):
            if len(des) == 0: break
            if node1 != node2 and node2 not in ans:
                if node2 in des:
                    m1  = {node1.name} if node1.is_leaf() else set(node1.get_leaf_names())
                    nm2 = {node2.name} if node2.is_leaf() else set(node2.get_leaf_names())
                    m   = m1 - nm2
                    if set(np.array(list(m), dtype="int")) == {i for i,v in enumerate(split_list) if v==1}:
                        return node1, node2


def _findNode(tree_nw, split_list, config):
    if config in ("Het", "Clonal"):
        mutated = [i for i,v in enumerate(split_list) if v==1]
        if len(mutated) == 1: return tree_nw & str(mutated[0])
        return tree_nw.get_common_ancestor([str(i) for i in mutated])
    elif config == "Par":
        return _findPar(tree_nw, split_list)
    elif config in ("Loss", "Back"):
        return _findLoss(tree_nw, split_list)
    else:
        print("Problem: unknown config", config); return 0


def _load_newick(newick_file: str) -> Tree:
    with open(newick_file) as f:
        raw = f.read()
    tmp = "_tmp_phalcon.nw"
    with open(tmp, "w") as f:
        f.write(_decrease_numbers(raw))
        f.write(";")
    t = Tree(tmp)
    os.remove(tmp)
    return t


def _build_clonal_config(fs_vf_file, fs_df_file, geno_cfg_file, original_labels_list):
    vf = pd.read_csv(fs_vf_file, sep="\t", header=None)
    vf = vf[~vf[0].isin(["intronic", "intergenic", "ncRNA_intronic"])]
    vf.loc[vf[1].str.contains("EZH2",  case=False, na=False), 1] = "EZH2"
    vf.loc[vf[1].str.contains("U2AF1", case=False, na=False), 1] = "U2AF1"

    final_variants = []
    for row in vf.to_numpy():
        final_variants.append(row[3]-1 if row[6]=="-" else row[3])

    orig_pos = pd.read_csv(fs_df_file, sep="\t", header=None)
    orig_pos = orig_pos[orig_pos[1].isin(final_variants)]

    cfg    = pd.read_csv(geno_cfg_file, sep="\t", header=None)
    before = pd.read_csv(fs_df_file,   sep="\t", header=None)

    numcells = cfg.shape[1] - 1
    cfg.insert(0, "chr",  before[0])
    cfg.insert(1, "site", before[1])
    cfg = cfg[cfg["site"].isin(orig_pos[1])]

    int_cols = list(range(numcells))
    cfg["sum"] = cfg[int_cols].sum(axis=1)
    cfg = cfg[cfg["sum"] != 0].copy()
    cfg.reset_index(drop=True, inplace=True)
    cfg.drop("sum", axis=1, inplace=True)

    cfg2 = cfg.T.drop_duplicates().T
    columns = {}
    for i in cfg2.columns[2:-1]:
        columns[i] = original_labels_list[i]

    n_orig = len(set(original_labels_list))
    for c in range(n_orig):
        if c not in columns.values():
            idx = original_labels_list.index(c)
            cfg2.insert(cfg2.shape[1]-1, idx, cfg[idx])
            columns[idx] = c

    cfg2.rename(columns=columns, inplace=True)
    cols = ["chr", "site"] + list(range(n_orig)) + [numcells]
    cfg2 = cfg2[cols]
    return np.array(cfg2, dtype="object"), n_orig


def _assign_labels(tree_nw, clonal_info, n_orig):
    for i, node in enumerate(tree_nw.traverse("preorder")):
        node.temp  = i + n_orig
        node.arrow = []
        node.label = []

    n_rows = clonal_info.shape[0]
    print(f"  Assigning labels for {n_rows} mutations …")
    for i in range(n_rows):
        split_list = clonal_info[i][2:-1]
        config     = clonal_info[i][-1]
        chr_       = clonal_info[i][0][3:]
        site       = clonal_info[i][1]
        lbl        = f"{chr_}:{site}"
        if config in ("Clonal", "Het"):
            node = _findNode(tree_nw, split_list, config)
            if node: node.label.append(lbl)
        elif config in ("Par", "Back", "Loss"):
            result = _findNode(tree_nw, split_list, config)
            if result:
                n1, n2 = result
                n1.label.append(lbl)
                n2.label.append(lbl + ("-" if config in ("Back", "Loss") else ""))
        else:
            print(f"  Warning: unknown config '{config}' row {i}")
    print("  Label assignment complete.")


def _remove_empty_nodes(tree_nw):
    """
    Snapshot empty internal nodes, delete them (ete3 rewires children to
    parent automatically), repeat until none remain.
    Never mutates the tree while iterating — avoids iterator corruption.
    """
    while True:
        to_delete = [
            n for n in tree_nw.traverse("preorder")
            if not n.is_leaf() and not n.label and n.up is not None
        ]
        if not to_delete:
            break
        for n in to_delete:
            if n.up is not None:   # may have been detached in this same pass
                n.delete(prevent_nondicotomic=False)


def _rebuild_arrows(tree_nw):
    """Reassign node.temp (fresh sequential IDs) and rebuild arrow lists."""
    for i, node in enumerate(tree_nw.traverse("preorder")):
        node.temp  = i + 1
        node.arrow = []
    for node in tree_nw.traverse("preorder"):
        if not node.is_leaf():
            node.arrow.extend(node.get_children())


def _gene_name_map(fs_vf_file):
    arr  = pd.read_csv(fs_vf_file, sep="\t", header=None).to_numpy(dtype="object")
    gmap = {}
    for row in arr:
        site = row[3]-1 if row[6]=="-" else row[3]
        key  = f"{row[2][3:]}:{site}"
        # EZH2 check first, then U2AF1, then default gene name
        if "EZH2"  in str(row[1]): gmap[key] = "EZH2"
        elif "U2AF1,U2AF1L5" in str(row[1]): gmap[key] = "U2AF1"
        else:                         gmap[key] = row[1]
    return gmap


def _translate_label(raw, gmap):
    for pos, gene in gmap.items():
        if pos in raw:
            m    = re.search(r":(.+)", pos)
            site = m.group(1) if m else pos
            raw  = raw.replace(pos, f"{gene}_{site}")
    return raw


def _write_gv(tree_nw, gmap, original_to_merged, palette,
              cells_per_original, out_path, with_boxes):
    """
    Write a graphviz .gv file.
    with_boxes=True  → clean tree + colored box nodes on leaves
    with_boxes=False → reference tree (all mutation nodes, no boxes)

    Node IDs:
      tree nodes → node.temp (1-based sequential, reassigned after deletion)
      box nodes  → 90000 + orig_cluster_id  (safely outside tree node range)
    """
    BOX_OFFSET = 90000
    live_temps = {n.temp for n in tree_nw.traverse("preorder")}

    with open(out_path, "w") as f:
        f.write("digraph CellPhyTree {\n")
        f.write("    graph [rankdir=TB];\n")
        f.write("    node  [fontname=Helvetica, fontsize=11];\n\n")

        # Mutation / clone nodes
        for node in tree_nw.traverse("preorder"):
            parts   = [_translate_label(lbl, gmap) for lbl in node.label]
            lbl_str = "\\n".join(parts) if parts else ""
            f.write(f'    {node.temp} [label="{lbl_str}", '
                    f'style=filled, fillcolor="#f5f5f5", shape=ellipse];\n')

        # Colored box nodes
        if with_boxes:
            f.write("\n    // Box nodes — one per original PHALCON cluster\n")
            for node in tree_nw.traverse("preorder"):
                if not node.is_leaf():
                    continue
                try:
                    orig_cid = int(node.name)
                except ValueError:
                    continue
                merged_cid = original_to_merged.get(orig_cid, orig_cid)
                hex_color  = rgb255_to_hex(palette[merged_cid])
                n_cells    = cells_per_original.get(orig_cid, 0)
                box_id     = BOX_OFFSET + orig_cid
                f.write(f'    {box_id} [label="{n_cells} cells", '
                        f'shape=box, style="filled,rounded", '
                        f'fillcolor="{hex_color}", fontsize=10];\n')

        # Tree edges — use get_children() (live topology, not stale arrows)
        f.write("\n    // Tree edges\n")
        for node in tree_nw.traverse("preorder"):
            if node.is_leaf():
                continue
            for child in node.get_children():
                if child.temp in live_temps:
                    f.write(f"    {node.temp} -> {child.temp};\n")

        # Leaf → box edges
        if with_boxes:
            f.write("\n    // Leaf → box edges\n")
            for node in tree_nw.traverse("preorder"):
                if not node.is_leaf():
                    continue
                try:
                    orig_cid = int(node.name)
                except ValueError:
                    continue
                f.write(f"    {node.temp} -> {BOX_OFFSET + orig_cid};\n")

        f.write("}\n")


def plot_tree(data, newick_file, fs_vf_file, fs_df_file, geno_cfg_file,
              original_labels, original_to_merged, palette, outdir):
    print("\n[TREE] Building clonal tree …")

    tree_nw   = _load_newick(newick_file)
    orig_list = original_labels.astype(int).tolist()

    clonal_info, n_orig = _build_clonal_config(fs_vf_file, fs_df_file, geno_cfg_file, orig_list)

    print("  [TREE] Step 1: assigning mutation labels …", flush=True)
    _assign_labels(tree_nw, clonal_info, n_orig)

    print("  [TREE] Step 2: building gene name map …", flush=True)
    gmap = _gene_name_map(fs_vf_file)

    orig_arr           = original_labels.astype(int)
    cells_per_original = {c: int(np.sum(orig_arr == c)) for c in sorted(set(orig_arr.tolist()))}

    # Reference tree (with empty nodes, no boxes)
    print("  [TREE] Step 3: writing reference tree …", flush=True)
    _rebuild_arrows(tree_nw)
    out_gv_full = os.path.join(outdir, f"{data}_clonal_tree_with_empty.gv")
    _write_gv(tree_nw, gmap, original_to_merged, palette,
              cells_per_original, out_gv_full, with_boxes=False)
    print(f"  Saved: {out_gv_full}")

    # Remove empty nodes → rebuild → write clean tree
    print("  [TREE] Step 4: removing empty nodes …", flush=True)
    _remove_empty_nodes(tree_nw)
    print("  [TREE] Step 5: rebuilding arrows …", flush=True)
    _rebuild_arrows(tree_nw)
    print("  [TREE] Step 6: writing clean tree …", flush=True)
    out_gv = os.path.join(outdir, f"{data}_clonal_tree.gv")
    _write_gv(tree_nw, gmap, original_to_merged, palette,
              cells_per_original, out_gv, with_boxes=True)
    print(f"  Saved: {out_gv}")

    # Render both .gv files to PNG using graphviz.render()
    # graphviz.render() reads the .gv and writes <name>.png — no duplicate file.
    for gv_path, label in [(out_gv_full, "full"), (out_gv, "clean")]:
        try:
            out_png = gv_path.replace(".gv", ".png")
            gv_lib.render("dot", format="png", filepath=gv_path, outfile=out_png)
            print(f"  Rendered ({label}): {out_png}")
        except Exception as e:
            print(f"  Warning: graphviz render failed for {label} — {e}")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════
def main():
    args   = parse_args()
    data   = args.data
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    cluster_labels_file   = args.cluster_labels_file
    vcf_file              = args.vcf_file
    variant_function_file = args.variant_function_file
    lklhd_file            = args.lklhd_file
    newick_file           = args.newick_file
    fs_vf_file            = args.variant_function_file
    fs_df_file            = args.final_df_file
    geno_cfg_file         = args.genotype_config_file

    print(f"\n[INIT] Loading cluster labels: {cluster_labels_file}")
    original_labels = np.loadtxt(cluster_labels_file).astype(int)
    n_orig = len(set(original_labels.tolist()))
    print(f"  {len(original_labels)} cells  |  {n_orig} original PHALCON clusters")

    print("\n[INIT] Building exonic genotype matrix …")
    df_genotype = build_genotype_matrix(vcf_file, variant_function_file)
    print(f"  Matrix shape: {df_genotype.shape}  (cells × exonic variants)")

    merged_labels, original_to_merged, merge_report = compute_merged_labels(
        df_genotype, original_labels)
    n_merged = len(set(merged_labels.tolist()))
    print_merge_report(merge_report, data, outdir)
    print(f"\n  Original clusters : {n_orig}")
    print(f"  Merged clusters   : {n_merged}  ← drives the color palette")

    palette = build_palette(merged_labels)
    print("\n  Merged cluster → color:")
    for cid, rgb in palette.items():
        print(f"    Merged cluster {cid:3d}  →  {rgb255_to_hex(rgb)}")
    print("\n  Original → merged mapping:")
    for orig, merged in sorted(original_to_merged.items()):
        print(f"    orig {orig:3d}  →  merged {merged:3d}  →  {rgb255_to_hex(palette[merged])}")

    plot_heatmap(data=data, df_genotype=df_genotype,
                 merged_labels=merged_labels, palette=palette, outdir=outdir)

    plot_umap(data=data, lklhd_file=lklhd_file, merged_labels=merged_labels,
              palette=palette, seed=args.seed, min_dist=args.min_dist,
              n_neighbors=args.n_neighbors, outdir=outdir)

    plot_tree(data=data, newick_file=newick_file,
              fs_vf_file=fs_vf_file, fs_df_file=fs_df_file,
              geno_cfg_file=geno_cfg_file, original_labels=original_labels,
              original_to_merged=original_to_merged, palette=palette, outdir=outdir)

    print(f"\n✓ Pipeline complete. Outputs in: {outdir}")


if __name__ == "__main__":
    main()