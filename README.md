# RDKit Examples

A template project for working with RDKit, the open-source cheminformatics toolkit.

## Setup

### Option 1: Using Conda Environment File (Recommended)

```bash
# Create environment from file
conda env create -f environment.yml

# Activate environment
conda activate rdkit-env
```

### Option 2: Manual Conda Installation

```bash
# Create a new conda environment
conda create -n rdkit-env python=3.12

# Activate environment
conda activate rdkit-env

# Install RDKit from conda-forge
conda install -c conda-forge rdkit

# Install additional dependencies
conda install -c conda-forge numpy pandas matplotlib jupyter pytest
```

### Option 3: Using pip

```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
# On Linux/Mac:
source venv/bin/activate
# On Windows:
# venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

**Note:** Conda installation is strongly recommended as it provides pre-compiled binaries and better dependency management for RDKit.

## Project Structure

```
RDKit-examples/
├── environment.yml              # Conda environment specification
├── requirements.txt             # Python package dependencies
├── example.py                  # Basic RDKit example script
├── complete_analysis.py        # Interactive RDKit analysis script
├── analysis.ipynb              # Comprehensive RDKit analysis notebook
├── my_molecule.sdf             # Sample molecule in SDF format
├── my_molecule_formats.csv     # Molecule data in various formats
├── my_molecule_properties.csv  # Calculated molecular properties
├── my_molecule_structure.png   # Molecule structure image
└── README.md                   # This file
```

## Quick Start

Run the basic example script:

```bash
python example.py
```

Run the interactive analysis script:

```bash
python complete_analysis.py
```

Or explore the comprehensive analysis notebook:

```bash
jupyter notebook analysis.ipynb
```

## Analysis Tools

### Interactive Script (`complete_analysis.py`)

A command-line tool for interactive chemical analysis. Run it to analyze molecules with a menu-driven interface covering all major RDKit functionalities.

### Analysis Notebook (`analysis.ipynb`)

The `analysis.ipynb` notebook provides a comprehensive guide to using RDKit for chemical analysis, including:

- Molecule parsing and validation from various formats (SMILES, InChI, MOL)
- Format conversion and export
- Molecular property calculation (descriptors, fingerprints)
- Substructure searching and analysis
- Stereochemistry handling
- Charge neutralization and standardization
- Visualization and plotting
- Data export to CSV and image formats

The notebook uses sample molecules and generates the output files found in the project directory.

## Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [RDKit GitHub](https://github.com/rdkit/rdkit)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)

## License

This project is a template for educational and development purposes.
