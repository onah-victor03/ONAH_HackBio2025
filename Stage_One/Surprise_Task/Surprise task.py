"""
HackBio Internship - Stage 1 Surprise Task
Team: Glycine
# Author: Onah Victor
# Task a:  Use the normalized gene expression dataset to plot a clustered heatmap of the top differentially expressed genes between _HBR_ and _UHR_ samples.

"""

import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt

# dataset URL (as you provided)
NORMALIZED_URL = ("https://raw.githubusercontent.com/HackBio-Internship/"
                  "2025_project_collection/refs/heads/main/Python/Dataset/"
                  "hbr_uhr_top_deg_normalized_counts.csv")

sb.set(style="white", context="notebook")

def plot_heatmap(normalized_url: str = NORMALIZED_URL,
                 top_n: int = 25,
                 cmap: str = "Blues",
                 out_file: str = "fig_1A_heatmap.png"):
    """
    Load normalized counts, pick top_n genes (first rows in file assumed top),
    z-score rows and produce a clustered heatmap.

    Parameters
    ----------
    normalized_url : str
        URL to normalized counts. Expected format: 'gene' column + sample columns.
    top_n : int
        Number of genes to display (for readability).
    cmap : str
        Matplotlib colormap (Blues recommended).
    out_file : str
        Output PNG filename.
    """
    df = pd.read_csv(normalized_url)
    if "gene" in df.columns:
        df = df.set_index("gene")
    # choose top_n genes by order (dataset already top DEGs)
    df_top = df.iloc[:top_n].copy()
    # ensure numeric
    df_top = df_top.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # z-score each row (gene)
    mat = df_top.values
    row_means = np.mean(mat, axis=1, keepdims=True)
    row_stds = np.std(mat, axis=1, ddof=1, keepdims=True)
    row_stds[row_stds == 0] = 1.0
    zmat = (mat - row_means) / row_stds
    zdf = pd.DataFrame(zmat, index=df_top.index, columns=df_top.columns)

    # clustered heatmap
    g = sb.clustermap(zdf, cmap=cmap, linewidths=0.5, figsize=(10, 10))
    g.ax_heatmap.set_xlabel("Samples")
    g.ax_heatmap.set_ylabel("Genes")
    plt.suptitle("Clustered Heatmap of Top DEGs (HBR vs UHR)", y=1.02)
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.show()
    print("Saved heatmap to", out_file)


if __name__ == "__main__":
    plot_heatmap()
    
#Task b: Plot 'log2FoldChange' vs 'log10(Padj)' the DEG results.
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DEG_URL = ("https://raw.githubusercontent.com/HackBio-Internship/"
           "2025_project_collection/refs/heads/main/Python/Dataset/"
           "hbr_uhr_deg_chr22_with_significance.csv")

sb.set(style="whitegrid", context="notebook")

def plot_volcano(deg_url: str = DEG_URL,
                 lfc_col: str = "log2FoldChange",
                 padj_col: str = "padj",
                 out_file: str = "fig_1B_volcano.png",
                 lfc_cutoff: float = 1.0):
    """
    Create a volcano plot from DEG results.

    Parameters
    ----------
    deg_url : str
        URL to DEG CSV. Expected columns: gene, log2FoldChange, padj, (optional Significance).
    lfc_col : str
        Column name for log2 fold change.
    padj_col : str
        Column name for adjusted p-value.
    out_file : str
        Output filename.
    lfc_cutoff : float
        Threshold for considering up/down regulation.
    """
    deg = pd.read_csv(deg_url)
    # ensure padj exists
    if padj_col not in deg.columns:
        # try to detect padj-like column
        possible = [c for c in deg.columns if 'pad' in c.lower() or 'adj' in c.lower()]
        if possible:
            padj_col = possible[0]
            print("Using detected padj column:", padj_col)
        else:
            raise ValueError("Adjusted p-value column (padj) not found.")

    # replace zeros (avoid -log10(0))
    deg[padj_col] = deg[padj_col].apply(lambda x: np.nextafter(0, 1) if x <= 0 else x)

    # determine significance (use provided column if exists)
    sig_col = None
    for c in ['Significance', 'significance', 'sig']:
        if c in deg.columns:
            sig_col = c
            break

    if sig_col:
        deg['significance'] = deg[sig_col].astype(str).str.lower().replace({
            'upregulated': 'up', 'downregulated': 'down'
        }).fillna('ns')
    else:
        deg['significance'] = 'ns'
        deg.loc[(deg[lfc_col] >= lfc_cutoff) & (deg[padj_col] < 0.05), 'significance'] = 'up'
        deg.loc[(deg[lfc_col] <= -lfc_cutoff) & (deg[padj_col] < 0.05), 'significance'] = 'down'

    color_map = {'up': 'green', 'down': 'orange', 'ns': 'grey'}
    colors = deg['significance'].map(color_map).fillna('grey')

    plt.figure(figsize=(9, 6))
    plt.scatter(deg[lfc_col], -np.log10(deg[padj_col]), c=colors, s=20, alpha=0.7)
    plt.axvline(x=lfc_cutoff, color='black', linestyle='--')
    plt.axvline(x=-lfc_cutoff, color='black', linestyle='--')
    plt.xlabel("log2FoldChange")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano Plot (Chr22)")

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green', label='Upregulated', markersize=6),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', label='Downregulated', markersize=6),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='grey', label='Not significant', markersize=6)
    ]
    plt.legend(handles=legend_elements, title="Significance", loc="upper right")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.show()
    print("Saved volcano plot to", out_file)


if __name__ == "__main__":
    plot_volcano()
   

# Task c: Plot 'texture_mean' vs 'radius_mean'
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

DATA_URL = ("https://raw.githubusercontent.com/HackBio-Internship/"
            "2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv")

sb.set(style="whitegrid", context="notebook")

def plot_radius_vs_texture(data_url: str = DATA_URL,
                           out_file: str = "fig_1C_radius_texture.png"):
    """
    Plot texture_mean vs radius_mean colored by diagnosis (M/B).
    """
    df = pd.read_csv(data_url)
    required = {'radius_mean', 'texture_mean', 'diagnosis'}
    if not required.issubset(df.columns):
        raise ValueError("Missing required columns for scatter plot: " + str(required - set(df.columns)))

    plt.figure(figsize=(8, 6))
    sb.scatterplot(data=df, x="radius_mean", y="texture_mean", hue="diagnosis", s=70, alpha=0.8)
    plt.xlabel("radius_mean")
    plt.ylabel("texture_mean")
    plt.title("radius_mean vs texture_mean by diagnosis")
    plt.legend(title="diagnosis")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.show()
    print("Saved scatter plot to", out_file)


if __name__ == "__main__":
    plot_radius_vs_texture()
# Task d: Compute the correlation matrix of six key features: 'radius_mean', 'texture_mean'. 'perimeter_mean', 'area_mean', 'smoothness_mean', 'compactness_mean'.
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

DATA_URL = ("https://raw.githubusercontent.com/HackBio-Internship/"
            "2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv")

sb.set(style="white", context="notebook")

def plot_correlation_heatmap(data_url: str = DATA_URL,
                             features=None,
                             out_file: str = "fig_1D_corrheatmap.png"):
    """
    Compute correlation matrix and plot as annotated heatmap.
    """
    if features is None:
        features = ['radius_mean', 'texture_mean', 'perimeter_mean',
                    'area_mean', 'smoothness_mean', 'compactness_mean']

    df = pd.read_csv(data_url)
    missing = set(features) - set(df.columns)
    if missing:
        raise ValueError("Missing expected features: " + str(missing))

    corr_mat = df[features].corr()
    plt.figure(figsize=(7, 6))
    sb.heatmap(corr_mat, annot=True, fmt=".2f", cmap="Blues", linewidths=0.5, square=True,
               cbar_kws={"shrink": 0.7})
    plt.title("Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.show()
    print("Saved correlation heatmap to", out_file)


if __name__ == "__main__":
    plot_correlation_heatmap()
# Task e: Plot 'compactness_mean' vs 'smoothness_mean'
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

DATA_URL = ("https://raw.githubusercontent.com/HackBio-Internship/"
            "2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv")

sb.set(style="whitegrid", context="notebook")

def plot_smoothness_vs_compactness(data_url: str = DATA_URL,
                                   out_file: str = "fig_1E_smooth_compact.png"):
    df = pd.read_csv(data_url)
    required = {'smoothness_mean', 'compactness_mean', 'diagnosis'}
    if not required.issubset(df.columns):
        raise ValueError("Missing required columns: " + str(required - set(df.columns)))

    plt.figure(figsize=(8, 6))
    sb.scatterplot(data=df, x="smoothness_mean", y="compactness_mean", hue="diagnosis", s=70, alpha=0.85)
    plt.xlabel("smoothness_mean")
    plt.ylabel("compactness_mean")
    plt.title("compactness_mean vs smoothness_mean by diagnosis")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(title="diagnosis")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.show()
    print("Saved scatter to", out_file)


if __name__ == "__main__":
    plot_smoothness_vs_compactness()
# Task f: Plot kernel density estimates (KDE) 'area_mean'
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

DATA_URL = ("https://raw.githubusercontent.com/HackBio-Internship/"
            "2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv")

sb.set(style="whitegrid", context="notebook")

def plot_area_kde(data_url: str = DATA_URL, out_file: str = "fig_1F_area_kde.png"):
    df = pd.read_csv(data_url)
    required = {'area_mean', 'diagnosis'}
    if not required.issubset(df.columns):
        raise ValueError("Missing required columns: " + str(required - set(df.columns)))

    plt.figure(figsize=(9, 6))
    # KDE for Malignant
    sb.kdeplot(data=df[df['diagnosis'] == 'M'], x='area_mean', fill=True, alpha=0.45, label='Malignant (M)', color='red')
    # KDE for Benign
    sb.kdeplot(data=df[df['diagnosis'] == 'B'], x='area_mean', fill=True, alpha=0.45, label='Benign (B)', color='skyblue')
    plt.xlabel("area_mean")
    plt.ylabel("Density")
    plt.title("KDE of area_mean by diagnosis")
    plt.legend(title="diagnosis")
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.show()
    print("Saved KDE to", out_file)


if __name__ == "__main__":
    plot_area_kde()