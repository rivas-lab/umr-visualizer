# This file defines how the package should be installed.
from setuptools import setup, find_packages

setup(
    name="umr-visualizer",  # Package name
    version="0.1.2",  # Initial version
    author="Salma Zainana",
    author_email="szainana@stanford.edu",
    description="A Python package for visualizing and analyzing a meta-regression models",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/umr-visualizer",  # Your GitHub repo URL
    packages=find_packages(),  # Automatically discover all modules
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
