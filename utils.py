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
        raw_desc = raw_desc.split("Full=")[1]
        raw_desc = raw_desc.split(";")[0]
        
    # 2. Remove "Chain X, " pattern (common in PDB structures)
    elif "Chain " in raw_desc:
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
    analyser = ProteinAnalysis(str(sequence))
    
    return {
        "molecular_weight": analyser.molecular_weight(),
        "instability_index": analyser.instability_index(),
        "isoelectric_point": analyser.isoelectric_point(),
        "amino_acid_percent": analyser.get_amino_acids_percent()
    }

def get_uniprot_from_ncbi(ncbi_id):
    """
    Extracts UniProt ID from NCBI protein record.
    Returns UniProt ID or None.
    """
    try:
        # Clean the ID (remove version numbers)
        clean_id = ncbi_id.split('.')[0]
        
        # Check if it's already a UniProt ID pattern
        if clean_id.startswith(('P', 'Q', 'O', 'A', 'B', 'C')):
            # UniProt IDs typically start with these letters and are 6-10 chars
            if len(clean_id) >= 6 and len(clean_id) <= 10:
                return clean_id
        
        # Try to fetch cross-references from NCBI
        handle = Entrez.efetch(db="protein", id=ncbi_id, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # Look for UniProt cross-references
        if records and len(records) > 0:
            for feature in records[0].get('GBSeq_feature-table', []):
                for qual in feature.get('GBFeature_quals', []):
                    if qual.get('GBQualifier_name') == 'db_xref':
                        xref = qual.get('GBQualifier_value', '')
                        if xref.startswith('UniProtKB/'):
                            return xref.split('/')[-1].split(':')[-1]
        
        return None
    except:
        return None

def get_pdb_structure(ncbi_id, sequence):
    """
    Multi-source structure fetcher with fallback strategy:
    1. AlphaFold Database (AI-predicted, 200M+ structures)
    2. ESMFold (on-demand AI prediction from sequence)
    3. PDB (experimental structures)
    
    Returns: (source_name, pdb_data) tuple or (None, None) if all fail
    """
    
    # --- METHOD 1: AlphaFold Database (BEST COVERAGE) ---
    try:
        uniprot_id = get_uniprot_from_ncbi(ncbi_id)
        
        if uniprot_id:
            clean_id = uniprot_id.split('.')[0]
            st.info(f"🔍 Checking AlphaFold for UniProt: {clean_id}")
            
            # Try AlphaFold v4
            alphafold_url = f"https://alphafold.ebi.ac.uk/files/AF-{clean_id}-F1-model_v4.pdb"
            af_response = requests.get(alphafold_url, timeout=15)
            
            if af_response.status_code == 200:
                st.success("✅ Found AlphaFold structure!")
                return f"AlphaFold-{clean_id}", af_response.text
            
            # Try PDBe for experimental structures
            st.info("🔍 Checking PDB for experimental structure...")
            search_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{clean_id}"
            pdb_response = requests.get(search_url, timeout=10)
            
            if pdb_response.status_code == 200:
                data = pdb_response.json()
                if clean_id in data and data[clean_id]:
                    pdb_ids = list(data[clean_id].keys())
                    if pdb_ids:
                        pdb_id = pdb_ids[0]
                        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                        pdb_file = requests.get(pdb_url, timeout=10)
                        
                        if pdb_file.status_code == 200:
                            st.success("✅ Found experimental PDB structure!")
                            return f"PDB-{pdb_id}", pdb_file.text
    except Exception as e:
        st.warning(f"AlphaFold/PDB lookup issue: {e}")
    
    # --- METHOD 2: ESMFold (ON-DEMAND PREDICTION) ---
    try:
        st.info("🧠 Generating structure with ESMFold AI...")
        
        # ESMFold API endpoint
        esm_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        
        # Clean sequence (remove any whitespace/newlines)
        clean_seq = str(sequence).replace('\n', '').replace(' ', '').strip()
        
        # ESMFold has length limits (typically ~400 residues work well)
        if len(clean_seq) > 400:
            st.warning("⚠️ Sequence is long, ESMFold may take time or fail")
        
        # Make request
        response = requests.post(
            esm_url,
            data=clean_seq,
            headers={'Content-Type': 'text/plain'},
            timeout=60  # ESMFold can take time
        )
        
        if response.status_code == 200:
            st.success("✅ ESMFold prediction complete!")
            return "ESMFold-Prediction", response.text
        else:
            st.warning(f"ESMFold failed with status: {response.status_code}")
            
    except requests.Timeout:
        st.warning("⚠️ ESMFold timed out (sequence may be too long)")
    except Exception as e:
        st.warning(f"ESMFold error: {e}")
    
    # --- ALL METHODS FAILED ---
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
