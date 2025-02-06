import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import plotly.express as px

# import stats  # Import your stats.py module
from . import stats 

def plot_proportion_variance_explained(table, fitted_values):
    """
    Create a barplot of the proportion of variance explained per moderator and the full model.
    
    Parameters:
        table (pd.DataFrame): The dataset containing 'Beta', 'SE', and moderator columns.
        fitted_values (np.array): The predicted y-values (fitted values) using all moderators.
    """
    # Identify moderator columns
    moderators = [col for col in table.columns if col.startswith("Mod")]

    # Compute proportion of variance explained for each moderator
    variance_explained = {mod: stats.proportion_variance_explained(table, fitted_values, mod) for mod in moderators}

    # Compute R² for the full model
    ssr_full = stats.ssr_mod_all(table, fitted_values)
    sst_value = stats.sst_all(table)
    r_squared_full = 1 - (ssr_full / sst_value)
    variance_explained["Full Model"] = r_squared_full  # Add full model R²

    # Create bar plot
    plt.figure(figsize=(10, 6))
    plt.bar(variance_explained.keys(), variance_explained.values(), color="skyblue")
    plt.xlabel("Moderators")
    plt.ylabel(r"Proportion of Variance Explained ($R_j^2 / R_{\text{total}}^2$)")
    plt.title("Proportion of Variance Explained per Moderator and Full Model")
    plt.xticks(rotation=45)
    plt.ylim(0, 1)  # R² values range from 0 to 1
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

