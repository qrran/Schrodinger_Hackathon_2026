# Chemical Similarity Search Application
## Stack the Future 2026

A Flask application that accepts SMILES strings for chemicals/drugs and returns chemically similar compounds using RDKit fingerprints and the ChEMBL database.

## 🚀 Quick Start

### 1. Activate Virtual Environment & Install Dependencies

```bash
cd Schrodinger_Hackathon_2026
source venv/bin/activate
pip install -r requirements.txt
```

### 2. Start the Flask Server

```bash
flask run
```

The application will be available at: **http://127.0.0.1:5001**

**Note:** If port 5001 is in use, specify a different port:
```bash
FLASK_RUN_PORT=5003 flask run
```
Then access at: **http://127.0.0.1:5003**

## What Was Built

### Files Created/Modified

1. **`requirements.txt`** - Added `chembl-webresource-client==0.10.8`
2. **`chemistry.py`** (NEW) - Core chemical similarity logic with:
   - Morgan fingerprint generation
   - Molecular descriptor calculation (43 properties)
   - 2D structure image generation
   - ChEMBL API integration
3. **`app.py`** - Added `/api/search-similar` endpoint
4. **`templates/my-template.html`** - Complete search interface with:
   - Structure image display
   - Expandable descriptor panels
   - Interactive results table

### Key Features

**Morgan Fingerprint Similarity** (radius=2, 1024 bits)
- Uses RDKit's Morgan fingerprint generator
- Tanimoto similarity coefficient for comparison

**ChEMBL Database Integration**
- Searches millions of compounds
- Configurable similarity threshold (default 70%)
- Returns ChEMBL IDs, names, and SMILES

**2D Structure Images**
- Visual representation of each compound
- Automatically generated PNG images (300x300 pixels)
- Base64-encoded for direct HTML embedding
- No external files needed

**Molecular Descriptors**
- 43 calculated properties for each compound
- Includes MW, LogP, TPSA, H-bond donors/acceptors, rotatable bonds, etc.
- Expandable details panel in the UI
- Uses RDKit's `rdMolDescriptors.Properties()`

**Beautiful Web Interface**
- Input form with validation
- Example SMILES buttons (Aspirin, Caffeine, Ibuprofen, Sumatriptan)
- Results table with structure images
- Color-coded similarity scores
- Expandable descriptor panels for detailed molecular properties
- Direct links to ChEMBL compound pages

**Drug Similarity Support**
- Searches based on active ingredient structure
- Works for any valid SMILES string

##  Example Usage

### Via Web Interface

1. Open http://127.0.0.1:5001 in your browser (or http://127.0.0.1:5003 if using alternate port)
2. Click "Aspirin" example button (or enter your own SMILES)
3. Adjust threshold and limit if desired
4. Click "Search Similar Compounds"
5. View results with:
   - 2D structure images for each compound
   - Similarity scores (color-coded)
   - ChEMBL IDs with direct links
6. Click "Details" button on any compound to see all 43 molecular descriptors

### Via API

```bash
curl -X POST http://127.0.0.1:5001/api/search-similar \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "limit": 10,
    "threshold": 70
  }'
```

### Example SMILES Strings

- **Aspirin**: `CC(=O)Oc1ccccc1C(=O)O`
- **Caffeine**: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
- **Ibuprofen**: `CC(C)Cc1ccc(cc1)C(C)C(=O)O`
- **Sumatriptan**: `CN(C)CCc1c[nH]c2ccc(cc12)CS(=O)(=O)N`

## 🔬 How It Works

### 1. User Submits SMILES

```
Input: CC(=O)Oc1ccccc1C(=O)O (Aspirin)
```

### 2. RDKit Validates & Generates Fingerprint

```python
MORGAN_GENERATOR = rdFpGen.GetMorganGenerator(radius=2, fpSize=1024)
query_mol = Chem.MolFromSmiles(smiles)
query_fp = MORGAN_GENERATOR.GetFingerprint(query_mol)
```

### 3. ChEMBL API Search

```python
similarity.filter(smiles=query_smiles, similarity=70)
```

Returns compounds with ≥70% Tanimoto similarity

### 4. Calculate Precise Similarity

```python
compound_fp = MORGAN_GENERATOR.GetFingerprint(compound_mol)
similarity_score = TanimotoSimilarity(query_fp, compound_fp)
```

### 5. Calculate Molecular Descriptors & Generate Images

For each compound:
```python
# Calculate 43 molecular descriptors
properties = rdMolDescriptors.Properties()
names = properties.GetPropertyNames()
values = properties.ComputeProperties(mol)

# Generate 2D structure image
img = Draw.MolToImage(mol, size=(300, 300))
img_str = base64.b64encode(buffered.getvalue()).decode()
```

### 6. Return Sorted Results

Results sorted by similarity (highest first) with:
- 2D structure image
- ChEMBL ID
- Compound name
- SMILES structure
- Tanimoto score (0.000-1.000)
- 43 molecular descriptors (expandable)

## Molecular Descriptors

Each search result includes 43 calculated molecular properties:

**Sample Properties:**
- `exactmw` - Exact molecular weight
- `amw` - Average molecular weight
- `lipinskiHBA` - Lipinski H-bond acceptors
- `lipinskiHBD` - Lipinski H-bond donors
- `NumRotatableBonds` - Number of rotatable bonds
- `NumHBD` - Number of hydrogen bond donors
- `NumHBA` - Number of hydrogen bond acceptors
- `NumHeavyAtoms` - Number of heavy atoms
- `NumAromaticRings` - Number of aromatic rings
- And 34 more...

**Example for Aspirin:**
```
exactmw             : 180.042
amw                 : 180.159
lipinskiHBA         : 4.0
lipinskiHBD         : 1.0
NumRotatableBonds   : 3.0
NumHeavyAtoms       : 13.0
NumAromaticRings    : 1.0
```

**How to View:**
1. Search for a compound
2. Click the "Details" button on any result
3. View all 43 descriptors in an organized grid

##  Test Results

### Aspirin Search (Threshold: 70%)

Found 5 similar compounds:
1. **ASPIRIN** - 1.000 (exact match)
2. **CARBASPIRIN** - 0.889
3. **CARBASPIRIN CALCIUM** - 0.606
4. Other aspirin derivatives

### Caffeine Search (Threshold: 70%)

Found 4 similar compounds:
1. **CAFFEINE** - 1.000 (exact match)
2. Thioxanthine derivative - 0.724
3. **BISDIONIN B** - 0.700
4. **CAFFEINE CITRATE** - 0.632

### Fingerprint Comparison Test

Verified with example code:
- `O=C(O)c1ccccc1C(=O)O` → 0.520 similarity to aspirin
- `COc1ccccc1C(=O)O` → 0.643 similarity to aspirin
- `CC(=O)Oc1ccc(Cl)cc1` → 0.364 similarity to aspirin

##  Drug Similarity

For drugs (which may contain multiple chemical components), the system:
1. Extracts the canonical SMILES of the active ingredient
2. Generates fingerprint for the active component
3. Searches ChEMBL for structurally similar compounds
4. Returns drugs with similar active ingredients

##  Technical Stack

- **Flask 3.1.0** - Web framework
- **RDKit 2024.3.6** - Cheminformatics library
- **ChEMBL Web Resource Client 0.10.8** - API access
- **Bootstrap 5.3.2** - Frontend UI
- **Morgan Fingerprints** - Similarity algorithm (radius=2, 1024 bits)
- **Tanimoto Coefficient** - Similarity metric

##  API Reference

### POST /api/search-similar

**Request:**
```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "limit": 20,        // optional, default: 20
  "threshold": 70     // optional, default: 70 (percent)
}
```

**Response:**
```json
{
  "query_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "count": 5,
  "results": [
    {
      "chembl_id": "CHEMBL25",
      "name": "ASPIRIN",
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "similarity": 1.0,
      "image": "data:image/png;base64,iVBORw0KGgo...",
      "descriptors": {
        "exactmw": 180.042,
        "amw": 180.159,
        "lipinskiHBA": 4.0,
        "lipinskiHBD": 1.0,
        "NumRotatableBonds": 3.0,
        "NumHBD": 1.0,
        "NumHBA": 2.0,
        "NumHeavyAtoms": 13.0
        // ... 35 more properties
      }
    }
  ]
}
```

**Error Response:**
```json
{
  "error": "Invalid SMILES string"
}
```

## 🔍 Tips for Best Results

1. **Adjust Threshold**: Lower threshold (e.g., 60%) finds more diverse compounds
2. **Increase Limit**: Get more results by increasing the limit parameter
3. **Use Canonical SMILES**: RDKit automatically canonicalizes, but consistent input helps
4. **Check ChEMBL Links**: Click ChEMBL IDs in results to see detailed compound information

##  Troubleshooting

**"Invalid SMILES string"**
- Verify your SMILES syntax
- Try using RDKit to canonicalize: `Chem.MolToSmiles(Chem.MolFromSmiles(your_smiles))`

**"ChEMBL API error"**
- Check internet connection
- ChEMBL API may be temporarily unavailable
- Try again after a moment

**No results found**
- Lower the similarity threshold
- Verify the compound exists in ChEMBL
- Try a more common chemical structure

**Flask won't run**
- Make sure to activate the virtual environment first: `source venv/bin/activate`
- Then run: `flask run`

##  Additional Resources

- [Flask Documentation](https://flask.palletsprojects.com/en/stable/)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [ChEMBL Database](https://www.ebi.ac.uk/chembl/)
- [Morgan Fingerprints Paper](https://pubs.acs.org/doi/10.1021/ci100050t)
- [Tanimoto Similarity](https://en.wikipedia.org/wiki/Jaccard_index)
