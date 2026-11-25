"""
OSCC GSE275870: QC, PCA and paired differential expression
Metastasis (Me) versus primary tumor (Tm) in patients with both sample types.

Run:  python 01_qc_and_paired_de.py
Needs: pandas numpy matplotlib scikit-learn pydeseq2
Put GSE275870_counts_csv.gz in a folder named data next to this script.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

os.makedirs("results", exist_ok=True)

# 1. Load counts and keep protein coding genes
raw = pd.read_csv("data/GSE275870_counts_csv.gz", low_memory=False)
raw = raw[raw["bioType"] == "protein_coding"].drop_duplicates("gene")
counts = raw.set_index("gene").iloc[:, 4:]          # genes x samples

# 2. Sample table from column names: P001Tm1 = patient P001, primary tumor, replicate 1
meta = pd.DataFrame({"sample": counts.columns})
meta["patient"] = meta["sample"].str[:4]
meta["tissue"] = meta["sample"].str[4:6].map({"Tm": "primary", "Me": "metastasis"})
meta = meta.set_index("sample")

# 3. Filter low expression: at least 10 counts in at least 20 samples
keep = (counts >= 10).sum(axis=1) >= 20
counts = counts.loc[keep]
print(f"Genes kept: {counts.shape[0]}")

# 4. PCA on log2 CPM, all 108 samples
cpm = counts / counts.sum() * 1e6
logcpm = np.log2(cpm + 1)
pca = PCA(n_components=2)
pcs = pca.fit_transform(logcpm.T)
fig, ax = plt.subplots(figsize=(6, 5))
for tissue, color in [("primary", "#2a6f97"), ("metastasis", "#d62828")]:
    idx = (meta["tissue"] == tissue).values
    ax.scatter(pcs[idx, 0], pcs[idx, 1], s=22, c=color, label=tissue, alpha=0.8)
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.legend()
fig.tight_layout()
fig.savefig("results/pca_all_samples.png", dpi=200)
plt.close(fig)

# 5. Paired set: first primary and first metastasis sample for each patient with both
paired_samples = []
for patient, grp in meta.groupby("patient"):
    if {"primary", "metastasis"} <= set(grp["tissue"]):
        paired_samples.append(grp[grp["tissue"] == "primary"].index[0])
        paired_samples.append(grp[grp["tissue"] == "metastasis"].index[0])
paired_meta = meta.loc[paired_samples]
paired_counts = counts[paired_samples].T             # samples x genes
print(f"Paired patients: {paired_meta['patient'].nunique()}")

# 6. Differential expression with patient as covariate
dds = DeseqDataSet(counts=paired_counts, metadata=paired_meta,
                   design="~patient + tissue", quiet=True)
dds.deseq2()
stats = DeseqStats(dds, contrast=["tissue", "metastasis", "primary"], quiet=True)
stats.summary()
res = stats.results_df.sort_values("padj")
res.to_csv("results/de_metastasis_vs_primary_paired.csv")

sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)]
print(f"Significant genes: {len(sig)} (up {int((sig.log2FoldChange > 0).sum())}, "
      f"down {int((sig.log2FoldChange < 0).sum())})")

# 7. Volcano plot
res = res.dropna(subset=["padj"])
x, y = res["log2FoldChange"], -np.log10(res["padj"])
colors = np.where((res["padj"] < 0.05) & (x > 1), "#d62828",
                  np.where((res["padj"] < 0.05) & (x < -1), "#2a6f97", "#bbbbbb"))
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(x, y, s=6, c=colors)
ax.axvline(-1, ls="--", lw=0.6, c="k")
ax.axvline(1, ls="--", lw=0.6, c="k")
ax.axhline(-np.log10(0.05), ls="--", lw=0.6, c="k")
ax.set_xlabel("log2 fold change (metastasis vs primary)")
ax.set_ylabel("-log10 adjusted p")
for gene in res.head(10).index:
    ax.annotate(gene, (res.loc[gene, "log2FoldChange"], -np.log10(res.loc[gene, "padj"])), fontsize=7)
fig.tight_layout()
fig.savefig("results/volcano_metastasis_vs_primary.png", dpi=200)
plt.close(fig)

# 8. Heatmap of the top 30 genes by adjusted p, paired samples, z-scored log2 CPM
top = res.head(30).index
lc = np.log2(paired_counts.T.loc[top] / paired_counts.T.sum() * 1e6 + 1)
z = lc.sub(lc.mean(axis=1), axis=0).div(lc.std(axis=1), axis=0)
order = list(paired_meta.sort_values("tissue", ascending=False).index)
z = z[order]
fig, ax = plt.subplots(figsize=(7, 7))
im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax.set_yticks(range(len(top)))
ax.set_yticklabels(top, fontsize=7)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([f"{s[:4]} {meta.loc[s, 'tissue'][:3]}" for s in order], rotation=90, fontsize=6)
fig.colorbar(im, ax=ax, label="z-score")
fig.tight_layout()
fig.savefig("results/heatmap_top30.png", dpi=200)
plt.close(fig)
