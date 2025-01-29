# Allows users to import functions directly from umr_visualizer like this:
# from umr_visualizer import sst, plot_proportion_variance_explained

# Import functions from stats.py
from .stats import sst_all, ssr_mod_all, ssr_mod_reduced, proportion_variance_explained

# Import functions from plots.py
from .plots import plot_proportion_variance_explained, plot2_proportion_variance_explained, plot3_proportion_variance_explained, plotly_variance_explained