# Script for scVelo analysis on Seurat-processed samples
# Author: Serafiina Jaatinen

from pathlib import Path
import anndata
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import dataframe_image as dfi

DATA_DIR = Path(".")
FIGURE_DIR = Path("figures")

FIGURE_DIR.mkdir(exist_ok=True)

samples = {"L1": {"loom": "L1.loom", "loom_prefix": "L1:"},
    "L2": {"loom": "GB2-L2.loom", "loom_prefix": "GB2-L2:"},
    "L3": {"loom": "GB2-L3.loom", "loom_prefix": "GB2-L3:"}
}

scv.settings.verbosity = 3

for sample_name, cfg in samples.items():
    print(f"Processing {sample_name}")

    adata = anndata.read_loom(DATA_DIR / cfg["loom"])

    cell_ids = pd.read_csv(DATA_DIR / f"{sample_name}_norm_cells.csv")
    umap = pd.read_csv(DATA_DIR / f"{sample_name}_norm_embeddings.csv", index_col=0)
    metadata = pd.read_csv(DATA_DIR / f"{sample_name}_norm_metadata.csv", index_col=0)

    # Match loom cell IDs to Seurat cell IDs
    adata.obs_names = adata.obs_names.str.replace(rf"^{cfg['loom_prefix']}", "", regex=True)
    adata.obs_names = adata.obs_names.str.replace(r"x$", "-1", regex=True)

    # Filter common cells
    common_cells = np.intersect1d(adata.obs_names, cell_ids.iloc[:, 0])
    adata = adata[common_cells, :].copy()

    # Add Seurat metadata and embeddings
    adata.obs = metadata.loc[common_cells]
    adata.obsm["X_umap"] = umap.loc[common_cells].to_numpy()

    # Ensure unique gene names
    adata.var_names_make_unique()

    # Save processed AnnData object
    adata.write(DATA_DIR / f"{sample_name}_filtered_for_scvelo_res0.8.h5ad")

    # scVelo preprocessing
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)

    # Velocity analysis
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)

    # Velocity stream plot
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="neftel_state",
        legend_loc="right margin",
        figsize=(8, 8),
        min_mass=0,
        save=f"{sample_name}_scvelo_stream.png"
    )

    # Cell cycle analysis
    scv.tl.score_genes_cell_cycle(adata)
    # Seurat cell cycle scores
    scv.pl.scatter(
        adata,
        color_gradients=["S.Score", "G2M.Score"],
        smooth=True,
        perc=[5, 95],
        figsize=(6, 6),
        save=f"{sample_name}_scvelo_cellcycle_seurat.png"
    )
    # scVelo cell cycle scores
    scv.pl.scatter(
        adata,
        color_gradients=["S_score", "G2M_score"],
        smooth=True,
        perc=[5, 95],
        figsize=(6, 6),
        save=f"{sample_name}_scvelo_cellcycle.png"
    )

    # Velocity confidence and pseudotime
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_pseudotime(adata)
    keys = ["velocity_length", "velocity_confidence", "velocity_pseudotime"]
    scv.pl.scatter(
        adata,
        c=keys,
        cmap="coolwarm",
        perc=[5, 95],
        figsize=(6, 6),
        save=f"{sample_name}_scvelo_velocity.png"
    )
    # Mean velocity metrics by Neftel state
    df = (adata.obs.groupby("neftel_state")[keys].mean().T)
    styled_df = df.style.background_gradient(cmap="coolwarm", axis=1)
    dfi.export(styled_df, FIGURE_DIR / f"{sample_name}_scvelo_neftel_velocities.png")

    # Transition graph
    scv.pl.velocity_graph(
        adata,
        threshold=0.1,
        figsize=(6, 6),
        save=f"{sample_name}_scvelo_transitions.png"
    )

    print(f"Finished {sample_name}")

