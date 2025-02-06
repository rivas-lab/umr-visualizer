import numpy as np
import pandas as pd

def sst_all(table):
    """
    Compute the total sum of squares (SST).
    
    Parameters:
        table (pd.DataFrame): The dataset containing 'Beta' and 'SE' columns.
        fitted_values (np.array): The predicted y-values (fitted values).
    
    Returns:
        float: The SST value.
    """
    # Extract y values (Beta) and weights (wi = 1 / SE^2)
    y = table['Beta'].values
    w = 1 / (table['SE'].values ** 2)

    # Compute the weighted mean of y
    y_w = np.sum(w * y) / np.sum(w)

    # Compute SST
    sst_value = np.sum(w * (y - y_w) ** 2)
    
    return sst_value


def ssr_mod_all(table, fitted_values):
    """
    Compute the sum of squares for regression (SSR) using all moderators.
    
    Parameters:
        table (pd.DataFrame): The dataset containing 'Beta' and 'SE' columns.
        fitted_values (np.array): The predicted y-values (fitted values).
    
    Returns:
        float: The SSR value.
    """
    
    # Extract y values (Beta) and weights
    y = table['Beta'].values
    w = 1 / (table['SE'].values ** 2)

    # # Extract full design matrix (excluding 'Beta' and 'SE')
    # X = table[[col for col in table.columns if col.startswith("Mod")]].values

    # # Compute weighted least squares estimate
    # beta_full = np.linalg.inv(X.T @ w @ X) @ (X.T @ w @ y)

    # # Compute fitted values
    # fitted_values = X @ beta_full

    # Compute SSR
    ssr_value = np.sum(w * (y - fitted_values) ** 2)
    
    return ssr_value



def ssr_mod_reduced(table, fitted_values, mod):
    """
    Compute the SSR after removing a specific moderator.
    
    Parameters:
        table (pd.DataFrame): The dataset containing 'Beta', 'SE', and moderator columns.
        fitted_values (np.array): The predicted y-values (fitted values) using all moderators.
        mod (str): The column name of the moderator to remove.
    
    Returns:
        float: The SSR value after removing the moderator.
    """
    # Create a reduced model by removing the specified moderator
    predictors = [col for col in table.columns if col.startswith("Mod") and col != mod]
    
    # Fit a new regression model without the removed moderator
    X_reduced = table[predictors].values  # Extract remaining moderators
    y = table['Beta'].values
    w = 1 / (table['SE'].values ** 2)

    # Compute new fitted values (β*X_reduced)
    W = np.diag(w)
    beta_reduced = np.linalg.inv(X_reduced.T @ W @ X_reduced) @ X_reduced.T @ W @ y
    fitted_values_reduced = np.dot(X_reduced, beta_reduced)

    # Compute SSR for the reduced model
    ssr_reduced_value = np.sum(w * (y - fitted_values_reduced) ** 2)
    
    return ssr_reduced_value

def proportion_variance_explained(table, fitted_values, mod):
    """
    Compute the proportion of variance explained by a given moderator.
    
    Parameters:
        table (pd.DataFrame): The dataset containing 'Beta', 'SE', and moderator columns.
        mod (str): The column name of the moderator for which to compute R².
    
    Returns:
        float: The proportion of variance explained (R²) by the given moderator.
    """
    # Compute SSR for full model
    ssr_full = ssr_mod_all(table, fitted_values)

    # Compute SST
    sst = sst_all(table)

    # Compute R_squared for the full model
    r_squared_full = 1 - (ssr_full / sst)

    # Compute SSR for reduced model (excluding the moderator)
    ssr_reduced = ssr_mod_reduced(table, fitted_values, mod)

    # Compute proportion of variance explained
    r_squared_j = (ssr_reduced - ssr_full) / ssr_reduced

    return r_squared_j/r_squared_full
