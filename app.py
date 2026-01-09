from flask import Flask, render_template, request, jsonify
from chemistry import (
    search_similar_compounds_chembl,
    validate_smiles,
    get_smiles_from_name,
)

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

        if not data:
            return jsonify({"error": "Missing request body"}), 400

        # Allow searching by SMILES or by compound name
        query_smiles = None
        query_name = None

        smiles_in = data.get("smiles", "")
        name_in = data.get("compound_name", "")

        limit = data.get("limit", 20)
        threshold = data.get("threshold", 70)

        if smiles_in and smiles_in.strip():
            query_smiles = smiles_in.strip()
        elif name_in and name_in.strip():
            query_name = name_in.strip()
            # Resolve name to SMILES via ChEMBL
            resolved = get_smiles_from_name(query_name)
            if not resolved:
                return (
                    jsonify(
                        {
                            "error": "Invalid compound name. Please input a valid compound name."
                        }
                    ),
                    400,
                )
            query_smiles = resolved
        else:
            return (
                jsonify(
                    {"error": "Please provide a SMILES string or a compound name."}
                ),
                400,
            )

        # Validate SMILES
        if not validate_smiles(query_smiles):
            return jsonify({"error": "Invalid SMILES string"}), 400

        # Search for similar compounds
        results = search_similar_compounds_chembl(
            query_smiles, similarity_threshold=threshold, limit=limit
        )

        # Check for errors
        if isinstance(results, dict) and "error" in results:
            return jsonify(results), 500

        resp = {"query_smiles": query_smiles, "count": len(results), "results": results}

        if query_name:
            resp["query_name"] = query_name

        return jsonify(resp)

    except Exception as e:
        return jsonify({"error": f"Server error: {str(e)}"}), 500
