import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.io as pio


# import stats  # Import your stats.py module
from . import stats 

def plot_proportion_variance_explained(table, fitted_values, style="matplotlib", gene_filter=None):
    """
    Create a bar plot of the proportion of variance explained per moderator 
    and the full model (R²) for each gene separately.
    
    Parameters:
        table (pd.DataFrame): The dataset containing 'BETA', 'SE', 'Gene', and moderator columns.
        fitted_values (np.array): The predicted y-values (fitted values) using all moderators.
        gene_filter (list, optional): List of gene names to filter the analysis.
        style (sttr): The style of the plot. Can be 'matplotlib', 'seaborn', or 'plotly'.
    """
    # Find the index of "SE" and select all columns that come after it
    se_index = table.columns.get_loc("SE") 
    moderators = table.columns[se_index + 1 :].tolist()  
    print(f"Moderators identified: {moderators}")  # Debug: Print moderator columns

    # Get unique gene names (or use the provided filter)
    unique_genes = table['gene'].unique() if gene_filter is None else gene_filter
    print(f"Unique genes found: {unique_genes}") 

    # For each moderator, compute the per-gene variance explained.
    # Each value in the dictionary is a DataFrame with columns ['Gene', 'Proportion_Variance_Explained', 'R_squared_explained'].
    variance_explained = {}
    for mod in moderators:
        print(f"Computing variance explained for moderator: {mod}")  # Debug progress
        variance_explained[mod] = stats.proportion_variance_explained_per_gene(table, fitted_values, mod, gene_filter)


    # For each gene, build and show a bar plot.
    for gene in unique_genes:
        print(f"Generating plot for gene: {gene}")  # Debug progress
        # Extract variance explained for this gene
        gene_variance = {}
        for mod in moderators:
            df_mod = variance_explained[mod]
            print(df_mod.head())   # Debug output
            value = df_mod[df_mod['gene'] == gene]['Proportion_Variance_Explained'].values
            
            if len(value) == 0:
                print(f"WARNING: No variance explained data for moderator {mod} and gene {gene}")
                gene_variance[mod] = np.nan
            else:
                gene_variance[mod] = value[0]

        print(f"Variance data for gene {gene}: {gene_variance}")  # Debug output
            
        if style == "matplotlib":    
            # Create the bar plot for this gene
            plt.figure(figsize=(10, 6))
            x_labels = list(gene_variance.keys())
            y_values = list(gene_variance.values())
            plt.bar(x_labels, y_values, color="skyblue")
            plt.xlabel("Moderators")
            plt.ylabel("Proportion of Variance Explained (Rj² / R_total²)")
            plt.title(f"Proportion of Variance Explained for Gene: {gene}")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.show()

        elif style == "seaborn":
            # Create the bar plot for this gene
            df = pd.DataFrame(gene_variance.items(), columns=["Moderator", "Proportion Explained"])
            sns.set_theme(style="whitegrid")
            plt.figure(figsize=(10, 6))
            ax = sns.barplot(x="Moderator", y="Proportion Explained", data=df, hue="Moderator", palette="coolwarm", legend=False)
            plt.xlabel("Moderators", fontsize=14)
            plt.ylabel(r"Proportion of Variance Explained ($R_j^2 / R_{\text{total}}^2$)", fontsize=14)
            plt.title(f"Proportion of Variance Explained for Gene: {gene}", fontsize=16, pad=20)
            plt.xticks(rotation=45, fontsize=12)
            plt.ylim(0, 1)
            for i, v in enumerate(df["Proportion Explained"]):
                color = "white" if v > 0.5 else "black"
                ax.text(i, v - 0.05, f"{v:.2f}", ha='center', fontsize=12, fontweight='bold', color=color)
            plt.show()

        elif style == "plotly":
            # Create the interactive bar plot for this gene
            df = pd.DataFrame(gene_variance.items(), columns=["Moderator", "Proportion Explained"])
            df["Proportion Explained"] *= 100
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

        else:
            raise ValueError("Invalid style argument. Must be 'matplotlib', 'seaborn', or 'plotly'.")


def plot_r2(table, fitted_values, style="matplotlib", gene_filter=None):
    """
    Create a bar plot of the r_squared for the reduced models and full models
    Parameters:
        table (pd.DataFrame): The dataset containing 'BETA', 'SE', 'Gene', and moderator columns.
        fitted_values (np.array): The predicted y-values (fitted values) using all moderators.
        gene_filter (list, optional): List of gene names to filter the analysis.
        style (sttr): The style of the plot. Can be 'matplotlib', 'seaborn', or 'plotly'.
    """
    # Find the index of "SE" and select all columns that come after it
    se_index = table.columns.get_loc("SE") 
    moderators = table.columns[se_index + 1 :].tolist()  
    print(f"Moderators identified: {moderators}")  # Debug: Print moderator columns

    # Get unique gene names (or use the provided filter)
    unique_genes = table['gene'].unique() if gene_filter is None else gene_filter
    print(f"Unique genes found: {unique_genes}") 

    full_model_values = {}
    variance_explained = {}

    for mod in moderators:
        variance_explained[mod] = stats.proportion_variance_explained_per_gene(table, fitted_values, mod, gene_filter)
    # Compute the full model R² per gene (using the full data for each gene)
    full_model_values = {}
    for gene in unique_genes:
        print(f"Processing full model R² for gene: {gene}")
        gene_data = table[table['gene'] == gene]
        gene_fitted_values = fitted_values[table['gene'] == gene]
        ssr_full = stats.ssr_mod_all(gene_data, gene_fitted_values)
        sst = stats.sst_all(gene_data)
        print(f"For gene {gene} : SSR Full= {ssr_full}, SST= {sst}")  # Debug output
        if sst == 0:
            print(f"WARNING: SST is zero for gene {gene}, skipping R² computation.")
            full_model_values[gene] = np.nan
        else:
            full_model_values[gene] = 1 - (ssr_full / sst)

        # Debug: Print computed R² values
        print("Computed Full Model R² values:", full_model_values)

    for gene in unique_genes:
        gene_r2_values = {
            mod: variance_explained[mod].loc[variance_explained[mod]['gene'] == gene, 'R_squared_reduced'].values[0]
            if gene in variance_explained[mod]['gene'].values else np.nan
            for mod in moderators
        }
        gene_r2_values["Full Model"] = full_model_values[gene]
        
        # plots for each gene the r2 values
        x_labels = list(gene_r2_values.keys())
        y_values = list(gene_r2_values.values())
        
        if style == "matplotlib":    
            plt.figure(figsize=(10, 6))
            plt.bar(x_labels, y_values, color="skyblue")
            plt.xlabel("Moderators")
            plt.ylabel("R^2 Values")
            plt.title(f"R^2 Values for Gene: {gene}")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.show()
        elif style == "seaborn":
            df = pd.DataFrame(gene_r2_values.items(), columns=["Moderator", "R2 Values"])
            plt.figure(figsize=(10, 6))
            sns.barplot(x="Moderator", y="R2 Values", data=df, palette="coolwarm")
            plt.title(f"R^2 Values for Gene: {gene}")
            plt.xticks(rotation=45)
            plt.show()
        elif style == "plotly":
            df = pd.DataFrame(gene_r2_values.items(), columns=["Moderator", "R2 Values"])
            fig = px.bar(df, x="R2 Values", y="Moderator", orientation='h' ,color="R2 Values")
            fig.update_layout(title=f"R^2 Values for Gene: {gene}")
            fig.show()
        else:
            raise ValueError("Invalid style argument.")
