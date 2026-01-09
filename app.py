from flask import Flask, render_template, request, jsonify
from chemistry import search_similar_compounds_chembl, validate_smiles

app = Flask(__name__)


@app.route("/", methods=["GET"])
def root():
    return render_template("my-template.html")


@app.route("/api/search-similar", methods=["POST"])
def search_similar():
    """
    API endpoint to search for chemically similar compounds.
    
    Request JSON:
        {
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "limit": 10,  # optional, default 20
            "threshold": 70  # optional, default 70
        }
    
    Response JSON:
        {
            "query_smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "results": [
                {
                    "chembl_id": "CHEMBL25",
                    "name": "Aspirin",
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "similarity": 1.000
                },
                ...
            ]
        }
    """
    try:
        data = request.get_json()
        
        if not data or 'smiles' not in data:
            return jsonify({"error": "Missing 'smiles' parameter"}), 400
        
        query_smiles = data['smiles'].strip()
        limit = data.get('limit', 20)
        threshold = data.get('threshold', 70)
        
        # Validate SMILES
        if not validate_smiles(query_smiles):
            return jsonify({"error": "Invalid SMILES string"}), 400
        
        # Search for similar compounds
        results = search_similar_compounds_chembl(
            query_smiles, 
            similarity_threshold=threshold, 
            limit=limit
        )
        
        # Check for errors
        if isinstance(results, dict) and 'error' in results:
            return jsonify(results), 500
        
        return jsonify({
            "query_smiles": query_smiles,
            "count": len(results),
            "results": results
        })
        
    except Exception as e:
        return jsonify({"error": f"Server error: {str(e)}"}), 500
