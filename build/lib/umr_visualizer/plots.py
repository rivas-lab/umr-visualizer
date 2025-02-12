import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.express as px

# import stats  # Import your stats.py module
from . import stats 

def plot_proportion_variance_explained(table, fitted_values, gene_filter=None):
    """
    Create a bar plot of the proportion of variance explained per moderator 
    and the full model (R²) for each gene separately.
    
    Parameters:
        table (pd.DataFrame): The dataset containing 'Beta', 'SE', 'Gene', and moderator columns.
        fitted_values (np.array): The predicted y-values (fitted values) using all moderators.
        gene_filter (list, optional): List of gene names to filter the analysis.
    """
    # Identify moderator columns (assuming they begin with "Mod")
    moderators = [col for col in table.columns if col.startswith("Mod")]

    # Get unique gene names (or use the provided filter)
    unique_genes = table['Gene'].unique() if gene_filter is None else gene_filter

    # For each moderator, compute the per-gene variance explained.
    # Each value in the dictionary is a DataFrame with columns ['Gene', 'Proportion_Variance_Explained'].
    variance_explained = {
        mod: stats.proportion_variance_explained_per_gene(table, fitted_values, mod, gene_filter)
        for mod in moderators
    }
    
    # Compute the full model R² per gene (using the full data for each gene)
    full_model_values = {}
    for gene in unique_genes:
        gene_data = table[table['Gene'] == gene]
        gene_fitted_values = fitted_values[table['Gene'] == gene]
        ssr_full = stats.ssr_mod_all(gene_data, gene_fitted_values)
        sst = stats.sst_all(gene_data)
        full_model_values[gene] = 1 - (ssr_full / sst)

    # For each gene, build and show a bar plot.
    for gene in unique_genes:
        # Gather moderator-specific variance explained for this gene
        gene_variance = {}
        for mod in moderators:
            # Get the row for this gene from the DataFrame returned by proportion_variance_explained_per_gene
            df_mod = variance_explained[mod]
            gene_variance[mod] = df_mod[df_mod['Gene'] == gene]['Proportion_Variance_Explained'].values[0]

        # Add the full model R² for the gene
        gene_variance["Full Model"] = full_model_values[gene]
        
        # print variance correctly
        print(f"Variance Explained for Gene: {gene}")
        for mod in moderators:
            print(f"Reduced for Moderator: {mod}")
            print(variance_explained[mod])
        
        print(f"Full Model for Gene: {gene} is {full_model_values[gene]} and computed as 1 - ({ssr_full} / {sst}) on the gene_data alone")
            
            
        # Create the bar plot for this gene
        plt.figure(figsize=(10, 6))
        x_labels = list(gene_variance.keys())
        y_values = list(gene_variance.values())
        plt.bar(x_labels, y_values, color="skyblue")
        plt.xlabel("Moderators")
        plt.ylabel("Proportion of Variance Explained (Rj² / R_total²)")
        plt.title(f"Proportion of Variance Explained for Gene: {gene}")
        plt.xticks(rotation=45)
        #plt.ylim(0, 1)
        plt.tight_layout()
        plt.show()

    

def sns_proportion_variance_explained(table, fitted_values):
    """
    Create a barplot of the proportion of variance explained per moderator and the full model.
    """
    # Identify moderator columns
    moderators = [col for col in table.columns if col.startswith("Mod")]

    # Compute R² for the full model
    ssr_full = stats.ssr_mod_all(table, fitted_values)
    sst_value = stats.sst_all(table)
    r_squared_full = 1 - (ssr_full / sst_value)

    # Store variance explained in dictionary
    variance_explained = {mod: stats.proportion_variance_explained(table, fitted_values, mod) for mod in moderators}
    variance_explained["Full Model"] = r_squared_full  # Add full model R²

    # Convert to DataFrame
    df = pd.DataFrame(variance_explained.items(), columns=["Moderator", "Proportion Explained"])

    # Set Seaborn theme
    sns.set_theme(style="whitegrid")

    # Create bar plot
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x="Moderator", y="Proportion Explained", data=df, palette="coolwarm")

    plt.xlabel("Moderators", fontsize=14)
    plt.ylabel(r"Proportion of Variance Explained ($R_j^2 / R_{\text{total}}^2$)", fontsize=14)
    plt.title("Proportion of Variance Explained per Moderator and Full Model", fontsize=16, pad=20)  # Add padding to avoid overlap

    plt.xticks(rotation=45, fontsize=12)
    plt.ylim(0, 1)

    # Add value labels inside bars
    for i, v in enumerate(df["Proportion Explained"]):
        color = "white" if v > 0.5 else "black"  # White text if bar is dark
        ax.text(i, v - 0.05, f"{v:.2f}", ha='center', fontsize=12, fontweight='bold', color=color)

    plt.show()

def sns2_proportion_variance_explained(table, fitted_values):
    """
    Create a horizontal bar plot of the proportion of variance explained per moderator and the full model,
    with labels displayed inside the bars.
    """
    # Identify moderator columns
    moderators = [col for col in table.columns if col.startswith("Mod")]

    # Compute R² for the full model
    ssr_full = stats.ssr_mod_all(table, fitted_values)
    sst_value = stats.sst_all(table)
    r_squared_full = 1 - (ssr_full / sst_value)

    # Store variance explained in dictionary
    variance_explained = {mod: stats.proportion_variance_explained(table, fitted_values, mod) for mod in moderators}
    variance_explained = {"Full Model": r_squared_full, **variance_explained}  # Ensure Full Model is first

    # Convert to DataFrame
    df = pd.DataFrame(variance_explained.items(), columns=["Moderator", "Proportion Explained"])
    df["Proportion Explained"] *= 100  # Convert to percentage

    # Set Seaborn theme
    sns.set_theme(style="whitegrid")

    # Create horizontal bar plot
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(y="Moderator", x="Proportion Explained", data=df, palette="coolwarm")

    plt.xlabel("Proportion of Variance Explained (%)", fontsize=14)
    plt.ylabel("Moderators", fontsize=14)
    plt.title("Proportion of Variance Explained per Moderator and Full Model", fontsize=16)
    
    # Add value labels inside the bars
    for i, v in enumerate(df["Proportion Explained"]):
        ax.text(v - 10, i, f"{v:.1f}%", va='center', ha='right', fontsize=12, fontweight='bold', color="white")

    plt.xlim(0, 100)  # Set x-axis to range from 0 to 100%
    plt.show()


def plotly_variance_explained(table, fitted_values):
    """
    Create an interactive horizontal bar chart using Plotly with the Full Model placed first and visually separated.
    """
    # Identify moderator columns
    moderators = [col for col in table.columns if col.startswith("Mod")]

    # Compute R² for the full model
    ssr_full = stats.ssr_mod_all(table, fitted_values)
    sst_value = stats.sst_all(table)
    r_squared_full = 1 - (ssr_full / sst_value)

    # Store variance explained in dictionary
    variance_explained = {mod: stats.proportion_variance_explained(table, fitted_values, mod) for mod in moderators}
    
    # Ensure Full Model is first and add a visual separator (empty string as a fake category)
    variance_explained = {"Full Model": r_squared_full, "": None, **variance_explained}  # Add separator

    # Convert to DataFrame
    df = pd.DataFrame(variance_explained.items(), columns=["Moderator", "Proportion Explained"])
    df["Proportion Explained"] *= 100  # Convert to percentage

    # Assign colors: Full Model different color, separator empty, others use gradient
    df["Color"] = df["Moderator"].apply(lambda x: "red" if x == "Full Model" else "blue" if x != "" else "white")

    # Create interactive horizontal bar plot
    fig = px.bar(df, x="Proportion Explained", y="Moderator", text="Proportion Explained",
                 title="Proportion of Variance Explained per Moderator and Full Model",
                 labels={"Proportion Explained": "Proportion of Variance Explained (%)"},
                 orientation='h', color="Proportion Explained", color_continuous_scale="blues")

    # Ensure separator row has no bars
    fig.update_traces(texttemplate='%{text:.1f}%', textposition='inside')
    fig.update_layout(
        xaxis=dict(range=[0, 100]), 
        yaxis=dict(categoryorder="array", categoryarray=df["Moderator"].tolist()),  # Ensure order with separator
        template="plotly_white",
        coloraxis_showscale=False  # Hide color bar if not needed
    )

    fig.show()

