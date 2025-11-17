import streamlit as st
from Bio import Entrez, SeqIO
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import requests

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
Entrez.email = "rafekatchadorian27@gmail.com" 

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

def clean_description(raw_desc):
    """
    Cleans up messy NCBI descriptions.
    Input: "RecName: Full=Adrenergic toxin rho-elapitoxin-Dp1a; Short=..."
    Output: "Adrenergic toxin rho-elapitoxin-Dp1a"
    """
    # 1. Remove "RecName: Full=" pattern
    if "RecName: Full=" in raw_desc:
        # Split at "Full=" and take the part after it
        raw_desc = raw_desc.split("Full=")[1]
        # Remove anything after a semicolon (like "; Short=...")
        raw_desc = raw_desc.split(";")[0]
        
    # 2. Remove "Chain X, " pattern (common in PDB structures)
    elif "Chain " in raw_desc:
        # Split at the first comma and take the part after it
        parts = raw_desc.split(", ")
        if len(parts) > 1:
            raw_desc = parts[1]
            
    # 3. Remove generic "[Snake Name]" from the end if it exists
    raw_desc = raw_desc.split("[")[0].strip()
    
    return raw_desc

# -----------------------------------------------------------------------------
# CORE FUNCTIONS
# -----------------------------------------------------------------------------

def search_toxins(snake_name):
    """
    Searches NCBI Protein database for toxins.
    Returns a dictionary: { "Clean Name": "NCBI_ID" }
    """
    try:
        search_term = f"{snake_name} AND toxin NOT partial"
        
        # Search for IDs
        handle = Entrez.esearch(db="protein", term=search_term, retmax=20, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        
        if not id_list:
            return {}

        # Fetch details 
        fetch_handle = Entrez.efetch(db="protein", id=id_list, rettype="gb", retmode="text")
        
        results = {}
        for seq_record in SeqIO.parse(fetch_handle, "gb"):
            # Filter huge or tiny sequences
            if len(seq_record.seq) > 5000 or len(seq_record.seq) < 10:
                continue

            # CLEAN THE NAME
            pretty_name = clean_description(seq_record.description)
            
            # Create a label with ID for uniqueness
            # Format: "Mambalgin-1 (P0DKR6.1)"
            label = f"{pretty_name} ({seq_record.id})"
            results[label] = seq_record.id
            
        fetch_handle.close()
        return results

    except Exception as e:
        st.error(f"Error connecting to NCBI: {e}")
        return {}

def get_toxin_details(toxin_id):
    """
    Fetches the full sequence record for a specific Protein ID.
    Returns a Biopython SeqRecord object.
    """
    try:
        handle = Entrez.efetch(db="protein", id=toxin_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "gb")
        handle.close()
        return record
    except Exception as e:
        st.error(f"Error fetching toxin details: {e}")
        return None

def analyze_protein(sequence):
    """
    Takes a sequence string and returns a dictionary of physicochemical properties.
    Uses Bio.SeqUtils.ProtParam.
    """
    # ProteinAnalysis requires a string, not a Seq object
    analyser = ProteinAnalysis(str(sequence))
    
    # Calculate properties
    return {
        "molecular_weight": analyser.molecular_weight(),
        "instability_index": analyser.instability_index(),
        "isoelectric_point": analyser.isoelectric_point(),
        "amino_acid_percent": analyser.get_amino_acids_percent()
    }

def get_pdb_structure(uniprot_id):
    """
    Attempts to find and fetch PDB structure data for a protein.
    First tries to extract UniProt ID from the record ID,
    then queries the PDB API for associated structures.
    
    Returns: (pdb_id, pdb_data) tuple or (None, None) if not found
    """
    try:
        # Clean the UniProt ID (remove version numbers like .1, .2)
        clean_id = uniprot_id.split('.')[0]
        
        # Query PDB API to find structures associated with this UniProt ID
        search_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{clean_id}"
        response = requests.get(search_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # Get the first PDB ID from the results
            if clean_id in data and data[clean_id]:
                pdb_ids = list(data[clean_id].keys())
                if pdb_ids:
                    pdb_id = pdb_ids[0]  # Take the first structure
                    
                    # Now fetch the actual PDB file
                    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    pdb_response = requests.get(pdb_url, timeout=10)
                    
                    if pdb_response.status_code == 200:
                        return pdb_id, pdb_response.text
        
        # If PDBe search fails, try AlphaFold database
        alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{clean_id}-F1-model_v4.pdb"
        af_response = requests.get(alphafold_url, timeout=10)
        
        if af_response.status_code == 200:
            return f"AlphaFold-{clean_id}", af_response.text
            
    except Exception as e:
        st.warning(f"Could not fetch 3D structure: {e}")
    
    return None, None

# -----------------------------------------------------------------------------
# TEST BLOCK
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    print("--- RUNNING CLEANER TEST ---")
    test_snake = "King Cobra"
    hits = search_toxins(test_snake)
    
    if hits:
        print(f"\n✅ FOUND {len(hits)} CLEAN RESULTS:")
        for name, id in hits.items():
            print(f" - {name}")
    else:
        print("No results found.")
