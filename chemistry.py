"""
Chemistry utilities for molecular similarity search using RDKit and ChEMBL.
"""

import base64
from io import BytesIO
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import rdFingerprintGenerator as rdFpGen
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from chembl_webresource_client.new_client import new_client
from typing import List, Dict, Optional

# Configure Morgan fingerprint generator
# Radius 2, 1024 bits 
MORGAN_FP_SIZE = 1024
MORGAN_GENERATOR = rdFpGen.GetMorganGenerator(radius=2, fpSize=MORGAN_FP_SIZE)


def validate_smiles(smiles: str) -> Optional[Chem.Mol]:
    """
    Validate a SMILES string and return the RDKit molecule object.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        RDKit Mol object if valid, None if invalid
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return mol
    except Exception:
        return None


def calculate_fingerprint(mol: Chem.Mol):
    """
    Calculate Morgan fingerprint for a molecule.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Morgan fingerprint
    """
    return MORGAN_GENERATOR.GetFingerprint(mol)


def calculate_molecular_descriptors(mol: Chem.Mol) -> Dict[str, float]:
    """
    Calculate molecular descriptors for a molecule using RDKit.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Dictionary of descriptor names and values
    """
    try:
        properties = rdMolDescriptors.Properties()
        names = properties.GetPropertyNames()
        values = properties.ComputeProperties(mol)
        
        # Create dictionary of descriptors
        descriptors = {}
        for name, value in zip(names, values):
            descriptors[name] = round(value, 3) if isinstance(value, float) else value
        
        return descriptors
    except Exception as e:
        return {"error": f"Failed to calculate descriptors: {str(e)}"}


def generate_molecule_image(mol: Chem.Mol, size=(300, 300)) -> str:
    """
    Generate a 2D structure image of a molecule as a base64-encoded PNG.
    
    Args:
        mol: RDKit Mol object
        size: Tuple of (width, height) for the image
        
    Returns:
        Base64-encoded data URL string for direct HTML embedding
    """
    try:
        # Generate 2D coordinates if not present
        if mol.GetNumConformers() == 0:
            from rdkit.Chem import AllChem
            AllChem.Compute2DCoords(mol)
        
        # Generate image
        img = Draw.MolToImage(mol, size=size)
        
        # Convert to base64
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        
        # Return as data URL
        return f"data:image/png;base64,{img_str}"
    except Exception as e:
        return None


def search_similar_compounds_chembl(query_smiles: str, similarity_threshold: int = 70, limit: int = 20) -> List[Dict]:
    """
    Search ChEMBL for compounds similar to the query SMILES.
    
    Args:
        query_smiles: Query SMILES string
        similarity_threshold: Tanimoto similarity threshold (0-100)
        limit: Maximum number of results to return
        
    Returns:
        List of dictionaries containing compound information
    """
    try:
        # Validate query SMILES
        query_mol = validate_smiles(query_smiles)
        if query_mol is None:
            return {"error": "Invalid SMILES string"}
        
        # Query ChEMBL similarity endpoint
        similarity = new_client.similarity
        results = similarity.filter(smiles=query_smiles, similarity=similarity_threshold).only(
            'molecule_chembl_id', 
            'pref_name',
            'molecule_structures'
        )[:limit]
        
        # Calculate accurate Tanimoto similarity for each result
        query_fp = calculate_fingerprint(query_mol)
        compounds = []
        
        for compound in results:
            try:
                # Extract SMILES from molecule structures
                if compound.get('molecule_structures') and compound['molecule_structures'].get('canonical_smiles'):
                    compound_smiles = compound['molecule_structures']['canonical_smiles']
                    compound_mol = validate_smiles(compound_smiles)
                    
                    if compound_mol:
                        # Calculate precise Tanimoto similarity
                        compound_fp = calculate_fingerprint(compound_mol)
                        similarity_score = TanimotoSimilarity(query_fp, compound_fp)
                        
                        # Calculate molecular descriptors
                        descriptors = calculate_molecular_descriptors(compound_mol)
                        
                        # Generate 2D structure image
                        image = generate_molecule_image(compound_mol)
                        
                        compounds.append({
                            'chembl_id': compound.get('molecule_chembl_id', 'N/A'),
                            'name': compound.get('pref_name', 'Unknown'),
                            'smiles': compound_smiles,
                            'similarity': round(similarity_score, 3),
                            'descriptors': descriptors,
                            'image': image
                        })
            except Exception as e:
                # Skip compounds that cause errors
                continue
        
        # Sort by similarity (highest first)
        compounds.sort(key=lambda x: x['similarity'], reverse=True)
        
        return compounds
        
    except Exception as e:
        return {"error": f"ChEMBL API error: {str(e)}"}


def compare_smiles_list(query_smiles: str, smiles_list: List[str]) -> List[Dict]:
    """
    Compare a query SMILES against a list of SMILES strings.
    Calculate Tanimoto similarity for each.
    
    Args:
        query_smiles: Query SMILES string
        smiles_list: List of SMILES strings to compare against
        
    Returns:
        List of dictionaries with SMILES and similarity scores
    """
    query_mol = validate_smiles(query_smiles)
    if query_mol is None:
        return {"error": "Invalid query SMILES string"}
    
    query_fp = calculate_fingerprint(query_mol)
    results = []
    
    for smiles in smiles_list:
        mol = validate_smiles(smiles)
        if mol:
            fp = calculate_fingerprint(mol)
            similarity = TanimotoSimilarity(query_fp, fp)
            results.append({
                'smiles': smiles,
                'similarity': round(similarity, 3)
            })
    
    # Sort by similarity (highest first)
    results.sort(key=lambda x: x['similarity'], reverse=True)
    
    return results

