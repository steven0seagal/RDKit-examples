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
├── environment.yml      # Conda environment specification
├── requirements.txt     # Python package dependencies
├── example.py          # Basic RDKit example script
└── README.md           # This file
```

## Quick Start

Run the example script:

```bash
python example.py
```

## Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [RDKit GitHub](https://github.com/rdkit/rdkit)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)

## License

This project is a template for educational and development purposes.
