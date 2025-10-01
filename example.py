#!/usr/bin/env python
"""
Basic RDKit example demonstrating common operations.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd


def main():
    print("RDKit Example - Basic Molecular Operations\n")

    # Create molecules from SMILES
    molecules = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ethanol": "CCO",
        "Benzene": "c1ccccc1"
    }

    results = []

    for name, smiles in molecules.items():
        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            # Calculate molecular properties
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hba = Descriptors.NumHAcceptors(mol)
            hbd = Descriptors.NumHDonors(mol)

            results.append({
                "Name": name,
                "SMILES": smiles,
                "Molecular Weight": f"{mw:.2f}",
                "LogP": f"{logp:.2f}",
                "H-Bond Acceptors": hba,
                "H-Bond Donors": hbd
            })

            print(f"{name}:")
            print(f"  SMILES: {smiles}")
            print(f"  Molecular Weight: {mw:.2f}")
            print(f"  LogP: {logp:.2f}")
            print(f"  H-Bond Acceptors: {hba}")
            print(f"  H-Bond Donors: {hbd}")
            print()

    # Create a DataFrame with results
    df = pd.DataFrame(results)
    print("\nSummary Table:")
    print(df.to_string(index=False))

    # Generate 2D coordinates for visualization
    print("\n✓ Molecules processed successfully!")
    print("✓ You can visualize these molecules using Draw.MolToImage(mol)")


if __name__ == "__main__":
    main()
