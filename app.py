import streamlit as st
import utils
import matplotlib.pyplot as plt
import py3Dmol
from stmol import showmol

# -----------------------------------------------------------------------------
# PAGE CONFIGURATION
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Venom Profiler",
    page_icon="🐍",
    layout="wide"  # Changed to wide for better 3D visualization
)

# -----------------------------------------------------------------------------
# HEADER
# -----------------------------------------------------------------------------
st.title("🐍 Venom Profiler")
st.markdown("### Explore the molecular weapons of the snake world.")
st.markdown("---")

# -----------------------------------------------------------------------------
# SIDEBAR (Search)
# -----------------------------------------------------------------------------
with st.sidebar:
    st.header("🔍 Search Settings")
    snake_query = st.text_input("Enter Snake Name:", placeholder="e.g., Black Mamba")
    search_button = st.button("Find Toxins")
    
    st.markdown("---")
    st.markdown("### 🎨 3D Visualization Style")
    
    viz_style = st.selectbox(
        "Choose display style:",
        ["Cartoon", "Stick", "Sphere", "Line", "Cross"]
    )
    
    color_scheme = st.selectbox(
        "Color scheme:",
        ["spectrum", "chain", "secondary structure", "atom type"]
    )

# -----------------------------------------------------------------------------
# MAIN LOGIC
# -----------------------------------------------------------------------------

# Initialize session state to store results between clicks
if "toxin_results" not in st.session_state:
    st.session_state.toxin_results = {}

# When the user clicks "Find Toxins"
if search_button and snake_query:
    with st.spinner(f"Searching NCBI for '{snake_query}'..."):
        # Call our backend function
        results = utils.search_toxins(snake_query)
        
        if results:
            st.session_state.toxin_results = results
            st.success(f"Found {len(results)} toxins!")
        else:
            st.error("No toxins found. Try a different name (e.g., 'Dendroaspis polylepis').")

# -----------------------------------------------------------------------------
# RESULTS DISPLAY
# -----------------------------------------------------------------------------
if st.session_state.toxin_results:
    st.subheader("Select a Toxin to Analyze")
    
    # create the dropdown
    selected_label = st.selectbox(
        "Choose from the list:", 
        options=st.session_state.toxin_results.keys()
    )
    
    # Get the actual ID from the dictionary using the label
    selected_id = st.session_state.toxin_results[selected_label]
    
    st.info(f"You selected ID: **{selected_id}**")
    
    # Analysis Button
    if st.button("Analyze Toxin 🧬"):
        with st.spinner("Fetching sequence and calculating properties..."):
            # 1. Fetch the full record
            record = utils.get_toxin_details(selected_id)
            
            if record:
                # 2. Analyze properties
                stats = utils.analyze_protein(record.seq)
                
                # --- DISPLAY: BASIC INFO ---
                st.markdown("---")
                st.subheader(f"🧬 Analysis: {selected_label}")
                
                # Create two columns for layout
                info_col, struct_col = st.columns([1, 1])
                
                with info_col:
                    st.text_area("Sequence:", record.seq, height=150)
                    
                    # --- DISPLAY: METRICS ---
                    st.markdown("### 📊 Physicochemical Properties")
                    
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Molecular Weight", f"{stats['molecular_weight']:.2f} Da")
                    with col2:
                        st.metric("Isoelectric Point (pI)", f"{stats['isoelectric_point']:.2f}")
                    with col3:
                        stability = "Stable" if stats['instability_index'] < 40 else "Unstable"
                        st.metric("Stability", f"{stats['instability_index']:.2f} ({stability})")
                
                # --- 3D STRUCTURE VISUALIZATION ---
                with struct_col:
                    st.markdown("### 🔬 3D Structure Viewer")
                    
                    with st.spinner("Searching for 3D structure..."):
                        pdb_id, pdb_data = utils.get_pdb_structure(selected_id)
                        
                        if pdb_data:
                            st.success(f"✅ Structure found: {pdb_id}")
                            
                            # Create 3D visualization
                            view = py3Dmol.view(width=500, height=400)
                            view.addModel(pdb_data, 'pdb')
                            
                            # Apply selected style
                            style_map = {
                                "Cartoon": "cartoon",
                                "Stick": "stick",
                                "Sphere": "sphere",
                                "Line": "line",
                                "Cross": "cross"
                            }
                            
                            color_map = {
                                "spectrum": {"spectrum": {}},
                                "chain": {"chain": {}},
                                "secondary structure": {"ss": {}},
                                "atom type": {}
                            }
                            
                            view.setStyle({style_map[viz_style]: color_map[color_scheme]})
                            view.zoomTo()
                            view.spin(True)
                            
                            # Display the 3D structure
                            showmol(view, height=400, width=500)
                            
                            st.caption(f"🔄 Interactive 3D model - Use mouse to rotate, zoom, and explore")
                            
                            # Download button for PDB file
                            st.download_button(
                                label="📥 Download PDB File",
                                data=pdb_data,
                                file_name=f"{pdb_id}.pdb",
                                mime="text/plain"
                            )
                        else:
                            st.warning("⚠️ No 3D structure available for this protein")
                            st.info("3D structures are sourced from PDB and AlphaFold databases. Not all proteins have resolved structures.")
                
                # --- AMINO ACID COMPOSITION CHART (Full Width) ---
                st.markdown("---")
                st.markdown("### 🧪 Amino Acid Composition")
                
                # Prepare data for the chart
                aa_data = stats['amino_acid_percent']
                
                # Sort amino acids by percentage for better readability
                sorted_aa = sorted([(k, v * 100) for k, v in aa_data.items() if v > 0], 
                                   key=lambda item: item[1], reverse=True)
                
                labels = [item[0] for item in sorted_aa]
                percentages = [item[1] for item in sorted_aa]
                
                # Create the bar plot using Matplotlib
                fig, ax = plt.subplots(figsize=(12, 6))
                
                # Create horizontal bar chart
                ax.barh(labels, percentages, color='skyblue') 
                
                # Add labels and title
                ax.set_xlabel("Percentage (%)")
                ax.set_title("Amino Acid Composition of Toxin")
                ax.invert_yaxis()
                
                # Add percentage values on the bars
                for index, value in enumerate(percentages):
                    ax.text(value + 0.5, index, f'{value:.1f}%', va='center')
                
                # Show it in Streamlit
                st.pyplot(fig, use_container_width=True)
                
            else:
                st.error("Failed to fetch sequence details.")
else:
    st.info("👈 Start by searching for a snake in the sidebar.")
