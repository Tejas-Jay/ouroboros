import streamlit as st
from Bio import Entrez, SeqIO
import re # Added for text cleaning

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
Entrez.email = "student_project@bmsce.ac.in" 

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

# -----------------------------------------------------------------------------
# TEST BLOCK
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    print("--- RUNNING CLEANER TEST ---")
    test_snake = "Green Anaconda"
    hits = search_toxins(test_snake)
    
    if hits:
        print(f"\n✅ FOUND {len(hits)} CLEAN RESULTS:")
        for name, id in hits.items():
            print(f" - {name}")
    else:
        print("No results found.")