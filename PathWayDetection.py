"""
=============================================================================
POST-HOC BIOLOGICAL DISCOVERY: TCGA BREAST CANCER NETWORK REWIRING
=============================================================================

This script performs the secondary biological analysis presented in Section VI
of the manuscript. It computes pairwise Frobenius dissimilarities among breast
cancer molecular subtypes and identifies genes with the largest covariance
rewiring scores.

Features
--------
1. Pairwise subtype covariance dissimilarity.
2. Gene-level covariance rewiring score.
3. Z-score normalization for scale-invariant covariance estimation.

Input
-----
tcga_brca_5000genes_4subtypes.csv

Output
------
- Pairwise Frobenius distance matrix.
- Ranked list of top rewired genes.

=============================================================================
"""

import pandas as pd
import numpy as np
from itertools import combinations

# Input dataset
path = "tcga_brca_5000genes_4subtypes.csv"


def get_final_biological_results():

    print("--- Biological Discovery Analysis ---")

    df = pd.read_csv(path)

    # Identify candidate gene columns
    clinical_blacklist = [
        "age", "count", "weight", "score",
        "status", "load", "id"
    ]

    gene_candidates = [
        col for col in df.columns
        if not any(x in col.lower() for x in clinical_blacklist)
    ]

    gene_candidates = [
        col for col in gene_candidates
        if "unnamed" not in col.lower() and col != "subtype"
    ]

    # Select molecular subtypes
    subtypes = [
        "BRCA_LumA",
        "BRCA_LumB",
        "BRCA_Her2",
        "BRCA_Basal",
    ]

    df = df[df["subtype"].isin(subtypes)].copy()

    # Data cleaning
    data = df[gene_candidates].apply(pd.to_numeric, errors="coerce")
    data = data.dropna(axis=1, how="all")
    data = data.fillna(data.median())
    data = data.loc[:, data.std() > 1e-6]

    final_genes = data.columns.tolist()

    print(f"Number of valid genes: {len(final_genes)}")

    # Z-score normalization
    scaled_data = (data - data.mean()) / data.std()

    print("\n--- Computing Covariance Differences ---")

    cov_dict = {}

    for st in subtypes:
        subset = scaled_data[df["subtype"] == st].values
        cov_dict[st] = np.cov(subset, rowvar=False)

    # Pairwise Frobenius distances
    results_dist = []

    for g1, g2 in combinations(subtypes, 2):

        dist = np.linalg.norm(
            cov_dict[g1] - cov_dict[g2],
            ord="fro"
        )

        results_dist.append(
            {
                "Comparison": f"{g1} vs {g2}",
                "Distance": dist,
            }
        )

    dist_df = (
        pd.DataFrame(results_dist)
        .sort_values(by="Distance", ascending=False)
    )

    print("\n[Pairwise Structural Dissimilarity]")
    print(dist_df)

    # Gene-level rewiring score
    avg_cov = np.mean(list(cov_dict.values()), axis=0)

    gene_scores = []

    for i, gene in enumerate(final_genes):

        rewiring_magnitude = sum(
            np.sum(np.square(cov_dict[st][i, :] - avg_cov[i, :]))
            for st in subtypes
        )

        gene_scores.append(
            {
                "Gene": gene,
                "Score": rewiring_magnitude,
            }
        )

    top_genes = (
        pd.DataFrame(gene_scores)
        .sort_values(by="Score", ascending=False)
    )

    print("\n[Top 10 Rewired Genes]")
    print(top_genes.head(10))


if __name__ == "__main__":
    get_final_biological_results()
