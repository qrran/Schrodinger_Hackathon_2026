"""
Test script to verify fingerprint comparison matches user's example.
This demonstrates the exact implementation pattern requested.
"""

from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import rdFingerprintGenerator as rdFpGen

# Radius 2, 1024 bits
MORGAN_FP_SIZE = 1024
MORGAN_GENERATOR = rdFpGen.GetMorganGenerator(radius=2, fpSize=MORGAN_FP_SIZE)

# Example: Find molecules similar to aspirin in a list of molecules
aspirin_smiles = 'CC(=O)Oc1ccccc1C(=O)O'
query_mol = Chem.MolFromSmiles(aspirin_smiles)
query_fp = MORGAN_GENERATOR.GetFingerprint(query_mol)

other_smiles = ['O=C(O)c1ccccc1C(=O)O', 'COc1ccccc1C(=O)O', 'CC(=O)Oc1ccc(Cl)cc1']
other_mols = [Chem.MolFromSmiles(s) for s in other_smiles]

print("=" * 70)
print("Chemical Similarity Test - Morgan Fingerprints")
print("=" * 70)
print(f"\nQuery Molecule (Aspirin): {aspirin_smiles}")
print(f"Fingerprint Size: {MORGAN_FP_SIZE} bits")
print(f"Morgan Radius: 2")
print("\n" + "-" * 70)
print("Similarity Scores (Tanimoto Coefficient):")
print("-" * 70)

for i, mol in enumerate(other_mols):
    fp = MORGAN_GENERATOR.GetFingerprint(mol)
    similarity = TanimotoSimilarity(query_fp, fp)
    canonical_smiles = Chem.MolToSmiles(mol)
    print(f"{i+1}. {canonical_smiles:30s} → {similarity:.3f}")

print("\n" + "=" * 70)
print("Test completed successfully!")
print("=" * 70)

