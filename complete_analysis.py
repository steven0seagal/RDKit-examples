"""
RDKit Chemical Compound Analyzer - Interactive Script
Comprehensive tool for chemical compound analysis, format conversion, and property calculation
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, Fragments, Lipinski
from rdkit.Chem import rdMolDescriptors, rdmolops
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator
import sys

class ChemicalAnalyzer:
    """Comprehensive chemical compound analyzer using RDKit"""
    
    def __init__(self, input_str, input_format='smiles'):
        """
        Initialize with a chemical structure
        
        Args:
            input_str: Chemical structure string
            input_format: Format type ('smiles', 'inchi', 'mol', 'sdf')
        """
        self.input_format = input_format.lower()
        self.mol = self._parse_molecule(input_str)
        
        if self.mol is None:
            raise ValueError(f"Failed to parse molecule from {input_format}")
    
    def _parse_molecule(self, input_str):
        """Parse molecule from various formats"""
        if self.input_format == 'smiles':
            return Chem.MolFromSmiles(input_str)
        elif self.input_format == 'inchi':
            return Chem.MolFromInchi(input_str)
        elif self.input_format in ['mol', 'sdf']:
            return Chem.MolFromMolBlock(input_str)
        else:
            return Chem.MolFromSmiles(input_str)
    
    def convert_formats(self):
        """Convert molecule to multiple formats"""
        if self.mol is None:
            return None
        
        formats = {
            'SMILES': Chem.MolToSmiles(self.mol),
            'Canonical_SMILES': Chem.MolToSmiles(self.mol, canonical=True),
            'Isomeric_SMILES': Chem.MolToSmiles(self.mol, isomericSmiles=True),
            'InChI': Chem.MolToInchi(self.mol),
            'InChIKey': Chem.MolToInchiKey(self.mol),
            'Mol_Formula': rdMolDescriptors.CalcMolFormula(self.mol),
            'MOL_Block': Chem.MolToMolBlock(self.mol)
        }
        return formats
    
    def calculate_properties(self):
        """Calculate various molecular properties"""
        if self.mol is None:
            return None
        
        props = {
            # Basic properties
            'Molecular_Weight': Descriptors.MolWt(self.mol),
            'Exact_Mass': Descriptors.ExactMolWt(self.mol),
            'Heavy_Atom_Count': Descriptors.HeavyAtomCount(self.mol),
            'Num_Atoms': self.mol.GetNumAtoms(),
            'Num_Bonds': self.mol.GetNumBonds(),
            'Num_Rings': Descriptors.RingCount(self.mol),
            'Aromatic_Rings': Descriptors.NumAromaticRings(self.mol),
            
            # Lipinski's Rule of Five
            'LogP': Descriptors.MolLogP(self.mol),
            'HBD': Descriptors.NumHDonors(self.mol),
            'HBA': Descriptors.NumHAcceptors(self.mol),
            'Rotatable_Bonds': Descriptors.NumRotatableBonds(self.mol),
            
            # Additional descriptors
            'TPSA': Descriptors.TPSA(self.mol),
            'Molar_Refractivity': Descriptors.MolMR(self.mol),
            'Fraction_CSP3': Lipinski.FractionCSP3(self.mol),
            'Num_Stereocenters': len(Chem.FindMolChiralCenters(self.mol, includeUnassigned=True)),
            'Num_Defined_Stereocenters': len(Chem.FindMolChiralCenters(self.mol)),
        }
        
        # Lipinski violations
        violations = 0
        if props['Molecular_Weight'] > 500: violations += 1
        if props['LogP'] > 5: violations += 1
        if props['HBD'] > 5: violations += 1
        if props['HBA'] > 10: violations += 1
        props['Lipinski_Violations'] = violations
        
        return props
    
    def generate_fingerprints(self):
        """Generate various molecular fingerprints"""
        if self.mol is None:
            return None
        morgan_generator_r2 = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        morgan_generator_r3 = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
        fps = {
            'Morgan (R=2, 2048 bits)': morgan_generator_r2.GetFingerprint(self.mol),
            'Morgan (R=3, 2048 bits)': morgan_generator_r3.GetFingerprint(self.mol),
            'RDKit Topological': Chem.RDKFingerprint(self.mol),
            'MACCS Keys (166 bits)': AllChem.GetMACCSKeysFingerprint(self.mol),
            'Atom Pair': rdFingerprintGenerator.GetAtomPairGenerator().GetFingerprint(self.mol),
            # 'Topological Torsion': AllChem.GetTopologicalTorsionFingerprint(mol),
            "Topological Torsion": rdFingerprintGenerator.GetTopologicalTorsionGenerator().GetFingerprint(self.mol)
        }
        return fps
    
    def get_substructures(self):
        """Extract important substructures for similarity searching"""
        if self.mol is None:
            return None
        
        substructures = {}
        
        # Murcko scaffold
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(self.mol)
            substructures['Murcko_Scaffold'] = Chem.MolToSmiles(scaffold)
        except:
            substructures['Murcko_Scaffold'] = None
        
        # Generic scaffold
        try:
            generic_scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
            substructures['Generic_Scaffold'] = Chem.MolToSmiles(generic_scaffold)
        except:
            substructures['Generic_Scaffold'] = None
        
        # Functional groups
        substructures['Functional_Groups'] = self._identify_functional_groups()
        
        # Ring systems
        ring_info = self.mol.GetRingInfo()
        substructures['Ring_Systems'] = {
            'Num_Rings': ring_info.NumRings(),
            'Ring_Sizes': [len(ring) for ring in ring_info.AtomRings()]
        }
        
        return substructures
    
    def _identify_functional_groups(self):
        """Identify common functional groups"""
        groups = {
            'Alcohols': Fragments.fr_Al_OH(self.mol),
            'Aldehydes': Fragments.fr_aldehyde(self.mol),
            'Ketones': Fragments.fr_ketone(self.mol),
            'Carboxylic_Acids': Fragments.fr_COO(self.mol),
            'Esters': Fragments.fr_ester(self.mol),
            'Amines': Fragments.fr_NH2(self.mol) + Fragments.fr_NH1(self.mol) + Fragments.fr_NH0(self.mol),
            'Amides': Fragments.fr_amide(self.mol),
            'Aromatic_Rings': Fragments.fr_benzene(self.mol),
            'Halides': Fragments.fr_halogen(self.mol),
        }
        return {k: v for k, v in groups.items() if v > 0}
    
    def neutralize_charges(self):
        """Neutralize charged molecules"""
        if self.mol is None:
            return None
        
        uncharger = rdMolStandardize.Uncharger()
        neutral_mol = uncharger.uncharge(self.mol)
        
        return {
            'Original_SMILES': Chem.MolToSmiles(self.mol),
            'Neutralized_SMILES': Chem.MolToSmiles(neutral_mol),
            'Original_Charge': Chem.GetFormalCharge(self.mol),
            'Neutralized_Charge': Chem.GetFormalCharge(neutral_mol),
            'Neutralized_Mol': neutral_mol
        }
    
    def analyze_stereochemistry(self):
        """Analyze stereochemistry information"""
        if self.mol is None:
            return None
        
        stereo_info = {
            'Has_Stereochemistry': False,
            'Chiral_Centers': [],
            'Double_Bond_Stereo': [],
            'Num_Unspecified_Stereo': 0
        }
        
        # Find chiral centers
        chiral_centers = Chem.FindMolChiralCenters(self.mol, includeUnassigned=True)
        if chiral_centers:
            stereo_info['Has_Stereochemistry'] = True
            stereo_info['Chiral_Centers'] = [
                {'Atom_Idx': idx, 'Chirality': chirality} 
                for idx, chirality in chiral_centers
            ]
        
        # Check double bond stereochemistry
        for bond in self.mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                stereo = bond.GetStereo()
                if stereo != Chem.BondStereo.STEREONONE:
                    stereo_info['Has_Stereochemistry'] = True
                    stereo_info['Double_Bond_Stereo'].append({
                        'Bond_Idx': bond.GetIdx(),
                        'Atoms': (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                        'Stereo': str(stereo)
                    })
        
        # Count unspecified stereocenters
        unspecified = len([c for c in chiral_centers if c[1] == '?'])
        stereo_info['Num_Unspecified_Stereo'] = unspecified
        
        return stereo_info
    
    def draw_molecule(self, filename='molecule.png', size=(400, 400)):
        """Draw molecule to file"""
        if self.mol is None:
            return False
        
        try:
            img = Draw.MolToImage(self.mol, size=size)
            img.save(filename)
            return True
        except Exception as e:
            print(f"Error drawing molecule: {e}")
            return False
    
    def calculate_similarity(self, other_smiles):
        """Calculate Tanimoto similarity with another molecule"""
        if self.mol is None:
            return None
        
        other_mol = Chem.MolFromSmiles(other_smiles)
        if other_mol is None:
            return None
        
        fp1 = AllChem.GetMorganFingerprintAsBitVect(self.mol, 2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(other_mol, 2, nBits=2048)
        
        return TanimotoSimilarity(fp1, fp2)
    
    def full_analysis(self):
        """Run complete analysis"""
        results = {
            'Formats': self.convert_formats(),
            'Properties': self.calculate_properties(),
            'Fingerprints_Info': {k: f"{len(v)} bits" for k, v in self.generate_fingerprints().items()},
            'Substructures': self.get_substructures(),
            'Neutralization': self.neutralize_charges(),
            'Stereochemistry': self.analyze_stereochemistry()
        }
        return results


def interactive_mode():
    """Interactive mode for chemical analysis"""
    print("=" * 70)
    print("RDKit Chemical Compound Analyzer - Interactive Mode")
    print("=" * 70)
    
    while True:
        print("\nEnter a chemical structure (or 'quit' to exit):")
        user_input = input("SMILES/InChI: ").strip()
        
        if user_input.lower() in ['quit', 'exit', 'q']:
            print("Exiting...")
            break
        
        if not user_input:
            continue
        
        try:
            # Detect format
            input_format = 'inchi' if user_input.startswith('InChI=') else 'smiles'
            
            # Create analyzer
            analyzer = ChemicalAnalyzer(user_input, input_format)
            
            # Menu
            while True:
                print("\n" + "=" * 70)
                print("Analysis Options:")
                print("1. Convert to multiple formats")
                print("2. Calculate properties")
                print("3. Generate fingerprints")
                print("4. Extract substructures")
                print("5. Neutralize charges")
                print("6. Analyze stereochemistry")
                print("7. Draw molecule")
                print("8. Calculate similarity with another molecule")
                print("9. Full analysis")
                print("0. New molecule")
                print("=" * 70)
                
                choice = input("Select option: ").strip()
                
                if choice == '1':
                    formats = analyzer.convert_formats()
                    print("\n--- Multiple Formats ---")
                    for key, value in formats.items():
                        if key != 'MOL_Block':
                            print(f"{key}: {value}")
                
                elif choice == '2':
                    props = analyzer.calculate_properties()
                    print("\n--- Molecular Properties ---")
                    for key, value in props.items():
                        print(f"{key}: {value}")
                
                elif choice == '3':
                    fps = analyzer.generate_fingerprints()
                    print("\n--- Fingerprints Generated ---")
                    for key, fp in fps.items():
                        print(f"{key}: {len(fp)} bits, OnBits: {fp.GetNumOnBits()}")
                
                elif choice == '4':
                    subs = analyzer.get_substructures()
                    print("\n--- Substructures ---")
                    for key, value in subs.items():
                        print(f"{key}: {value}")
                
                elif choice == '5':
                    neutral = analyzer.neutralize_charges()
                    print("\n--- Charge Neutralization ---")
                    for key, value in neutral.items():
                        if key != 'Neutralized_Mol':
                            print(f"{key}: {value}")
                
                elif choice == '6':
                    stereo = analyzer.analyze_stereochemistry()
                    print("\n--- Stereochemistry Analysis ---")
                    for key, value in stereo.items():
                        print(f"{key}: {value}")
                
                elif choice == '7':
                    filename = input("Enter filename (default: molecule.png): ").strip()
                    if not filename:
                        filename = "molecule.png"
                    if analyzer.draw_molecule(filename):
                        print(f"Molecule saved to {filename}")
                    else:
                        print("Failed to draw molecule")
                
                elif choice == '8':
                    other = input("Enter SMILES of molecule to compare: ").strip()
                    similarity = analyzer.calculate_similarity(other)
                    if similarity is not None:
                        print(f"\nTanimoto Similarity: {similarity:.4f}")
                    else:
                        print("Failed to calculate similarity")
                
                elif choice == '9':
                    print("\n--- Full Analysis ---")
                    results = analyzer.full_analysis()
                    for section, data in results.items():
                        print(f"\n{section}:")
                        if isinstance(data, dict):
                            for key, value in data.items():
                                print(f"  {key}: {value}")
                        else:
                            print(f"  {data}")
                
                elif choice == '0':
                    break
                
                else:
                    print("Invalid option")
        
        except Exception as e:
            print(f"Error: {e}")
            print("Please check your input and try again.")


if __name__ == "__main__":
    # Example usage
    if len(sys.argv) > 1:
        # Command line mode
        smiles = sys.argv[1]
        analyzer = ChemicalAnalyzer(smiles)
        results = analyzer.full_analysis()
        
        print("\n=== Full Chemical Analysis ===\n")
        for section, data in results.items():
            print(f"{section}:")
            if isinstance(data, dict):
                for key, value in data.items():
                    print(f"  {key}: {value}")
            print()
    else:
        # Interactive mode
        interactive_mode()