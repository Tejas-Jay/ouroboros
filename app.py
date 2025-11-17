import streamlit as st
import utils
import matplotlib.pyplot as plt
import py3Dmol
import streamlit.components.v1 as components

# -----------------------------------------------------------------------------
# PAGE CONFIGURATION
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Venom Profiler",
    page_icon="🐍",
    layout="wide"
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
        ["Cartoon", "Stick", "Sphere", "Cross", "Line"]
    )
    
    color_scheme = st.selectbox(
        "Color scheme:",
        ["Rainbow Spectrum", "Secondary Structure", "Hydrophobicity", "Charge", "Atom Type"]
    )
    
    st.markdown("---")
    st.markdown("### 🏷️ Interactive Labels & Highlights")
    
    show_labels = st.checkbox("Show Residue Labels", value=True)
    label_frequency = st.slider("Label every Nth residue:", 5, 50, 15)
    
    show_key_residues = st.checkbox("Highlight Key Residues", value=True)
    
    st.markdown("---")
    show_surface = st.checkbox("Show Molecular Surface", value=False)
    if show_surface:
        surface_opacity = st.slider("Surface Opacity:", 0.1, 1.0, 0.5, 0.1)
    else:
        surface_opacity = 0.5

# -----------------------------------------------------------------------------
# HELPER FUNCTION FOR 3D VIEWER
# -----------------------------------------------------------------------------
def create_enhanced_3d_viewer(pdb_data, sequence, viz_style, color_scheme, 
                              show_labels, label_frequency, show_key_residues, 
                              show_surface, surface_opacity):
    """Creates an enhanced 3D molecular viewer with labels and highlighting"""
    
    # Create viewer
    view = py3Dmol.view(width=700, height=500)
    view.addModel(pdb_data, 'pdb')
    
    # Style mapping
    style_map = {
        "Cartoon": "cartoon",
        "Stick": "stick",
        "Sphere": "sphere",
        "Line": "line",
        "Cross": "cross"
    }
    
    # Apply base coloring based on scheme
    if color_scheme == "Rainbow Spectrum":
        view.setStyle({style_map[viz_style]: {'color': 'spectrum'}})
    elif color_scheme == "Secondary Structure":
        view.setStyle({style_map[viz_style]: {'color': 'ss'}})
    elif color_scheme == "Atom Type":
        view.setStyle({style_map[viz_style]: {}})
    elif color_scheme == "Hydrophobicity":
        # Color by residue hydrophobicity
        hydrophobic = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
        view.setStyle({}, {style_map[viz_style]: {'color': 'lightgray'}})
        view.setStyle({'resn': hydrophobic}, {style_map[viz_style]: {'color': 'orange'}})
    elif color_scheme == "Charge":
        # Base color
        view.setStyle({}, {style_map[viz_style]: {'color': 'white'}})
        # Positive (blue)
        view.setStyle({'resn': ['ARG', 'LYS', 'HIS']}, {style_map[viz_style]: {'color': 'blue'}})
        # Negative (red)
        view.setStyle({'resn': ['ASP', 'GLU']}, {style_map[viz_style]: {'color': 'red'}})
    
    # Add molecular surface
    if show_surface:
        view.addSurface(py3Dmol.VDW, {
            'opacity': surface_opacity,
            'color': 'lightblue'
        })
    
    # Highlight key residues with special styling
    if show_key_residues:
        # Cysteines (yellow spheres) - disulfide bonds
        view.addStyle(
            {'resn': 'CYS'},
            {'sphere': {'color': 'yellow', 'radius': 0.8}}
        )
        
        # Aromatic residues (green) - binding pockets
        view.addStyle(
            {'resn': ['PHE', 'TYR', 'TRP']},
            {'stick': {'color': 'green', 'radius': 0.3}}
        )
        
        # Prolines (magenta) - structural kinks
        view.addStyle(
            {'resn': 'PRO'},
            {'sphere': {'color': 'magenta', 'radius': 0.6}}
        )
    
    # Add residue labels
    if show_labels:
        seq_len = len(sequence)
        for i in range(0, seq_len, label_frequency):
            if i < seq_len:
                aa_code = sequence[i] if i < len(sequence) else 'X'
                view.addLabel(
                    f"{aa_code}{i+1}",
                    {
                        'position': {'resi': i+1},
                        'backgroundColor': 'black',
                        'backgroundOpacity': 0.8,
                        'fontColor': 'white',
                        'fontSize': 12,
                        'showBackground': True,
                        'alignment': 'center'
                    }
                )
    
    view.zoomTo()
    view.rotate(90, {'x': 0, 'y': 1, 'z': 0})  # Better initial angle
    view.spin(True)
    
    return view

# -----------------------------------------------------------------------------
# MAIN LOGIC
# -----------------------------------------------------------------------------

# Initialize session state
if "toxin_results" not in st.session_state:
    st.session_state.toxin_results = {}

# When the user clicks "Find Toxins"
if search_button and snake_query:
    with st.spinner(f"Searching NCBI for '{snake_query}'..."):
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
    
    selected_label = st.selectbox(
        "Choose from the list:", 
        options=st.session_state.toxin_results.keys()
    )
    
    selected_id = st.session_state.toxin_results[selected_label]
    
    st.info(f"You selected ID: **{selected_id}**")
    
    # Analysis Button
    if st.button("Analyze Toxin 🧬"):
        with st.spinner("Fetching sequence and calculating properties..."):
            record = utils.get_toxin_details(selected_id)
            
            if record:
                stats = utils.analyze_protein(record.seq)
                
                # --- DISPLAY: BASIC INFO ---
                st.markdown("---")
                st.subheader(f"🧬 Analysis: {selected_label}")
                
                # Create two columns
                info_col, struct_col = st.columns([1, 1])
                
                with info_col:
                    st.text_area("Sequence:", record.seq, height=150)
                    
                    # --- METRICS ---
                    st.markdown("### 📊 Physicochemical Properties")
                    
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Molecular Weight", f"{stats['molecular_weight']:.2f} Da")
                    with col2:
                        st.metric("Isoelectric Point (pI)", f"{stats['isoelectric_point']:.2f}")
                    with col3:
                        stability = "Stable" if stats['instability_index'] < 40 else "Unstable"
                        st.metric("Stability", f"{stats['instability_index']:.2f} ({stability})")
                
                # --- 3D STRUCTURE ---
                with struct_col:
                    st.markdown("### 🔬 3D Structure Viewer")
                    st.caption("📡 Sources: AlphaFold → ESMFold → PDB")
                    
                    with st.spinner("Fetching 3D structure from multiple sources..."):
                        source_name, pdb_data = utils.get_pdb_structure(selected_id, str(record.seq))
                        
                        if pdb_data:
                            # Display source badge
                            if "AlphaFold" in source_name:
                                st.success(f"🤖 AI Predicted: {source_name}")
                            elif "ESMFold" in source_name:
                                st.success(f"🧠 ESMFold Prediction")
                            elif "PDB" in source_name:
                                st.success(f"🔬 Experimental: {source_name}")
                            
                            # Create enhanced 3D viewer
                            view = create_enhanced_3d_viewer(
                                pdb_data, 
                                str(record.seq),
                                viz_style,
                                color_scheme,
                                show_labels,
                                label_frequency,
                                show_key_residues,
                                show_surface,
                                surface_opacity
                            )
                            
                            # Display
                            components.html(view._make_html(), height=520, scrolling=False)
                            
                            st.caption("🖱️ Click & drag to rotate • Scroll to zoom • Double-click to center")
                            
                            # Legend
                            if show_key_residues:
                                st.markdown("**🎯 Key Residues Highlighted:**")
                                col_a, col_b, col_c = st.columns(3)
                                with col_a:
                                    st.markdown("🟡 **Cysteines** - Disulfide bonds")
                                with col_b:
                                    st.markdown("🟢 **Aromatics** - Binding sites")
                                with col_c:
                                    st.markdown("🟣 **Prolines** - Structural kinks")
                            
                            # Structure info
                            with st.expander("📋 Detailed Structure Information"):
                                cys_count = str(record.seq).count('C')
                                pro_count = str(record.seq).count('P')
                                aromatic = str(record.seq).count('F') + str(record.seq).count('Y') + str(record.seq).count('W')
                                positive = str(record.seq).count('R') + str(record.seq).count('K') + str(record.seq).count('H')
                                negative = str(record.seq).count('D') + str(record.seq).count('E')
                                
                                st.markdown(f"""
                                **Sequence Properties:**
                                - Total Length: {len(str(record.seq))} amino acids
                                
                                **Structural Features:**
                                - 🟡 Cysteines (C): {cys_count} → Potential {cys_count//2} disulfide bridges
                                - 🟣 Prolines (P): {pro_count} → Rigid structural elements
                                - 🟢 Aromatics (F/Y/W): {aromatic} → Hydrophobic core/binding
                                
                                **Charge Distribution:**
                                - Positive (+): {positive} residues (Arg, Lys, His)
                                - Negative (-): {negative} residues (Asp, Glu)
                                - Net Charge: {positive - negative:+d}
                                - Charge Ratio: {positive/negative if negative > 0 else '∞':.2f}
                                
                                **Toxin Insights:**
                                - High cysteine content suggests a tightly folded, stable structure
                                - Net charge affects target binding and membrane interaction
                                - Aromatic residues often form toxin-receptor binding interfaces
                                """)
                            
                            # Download
                            st.download_button(
                                label="📥 Download PDB File",
                                data=pdb_data,
                                file_name=f"{source_name}.pdb",
                                mime="text/plain"
                            )
                        else:
                            st.error("❌ No structure available from any source")
                            st.info("""
                            **Tried sources:**
                            1. ✗ AlphaFold Database (200M+ structures)
                            2. ✗ ESMFold AI Prediction
                            3. ✗ PDB Experimental Structures
                            
                            This protein may be too new or not yet characterized.
                            """)
                
                # --- AMINO ACID COMPOSITION ---
                st.markdown("---")
                st.markdown("### 🧪 Amino Acid Composition")
                
                aa_data = stats['amino_acid_percent']
                sorted_aa = sorted([(k, v * 100) for k, v in aa_data.items() if v > 0], 
                                   key=lambda item: item[1], reverse=True)
                
                labels = [item[0] for item in sorted_aa]
                percentages = [item[1] for item in sorted_aa]
                
                fig, ax = plt.subplots(figsize=(12, 6))
                ax.barh(labels, percentages, color='skyblue') 
                ax.set_xlabel("Percentage (%)")
                ax.set_title("Amino Acid Composition of Toxin")
                ax.invert_yaxis()
                
                for index, value in enumerate(percentages):
                    ax.text(value + 0.5, index, f'{value:.1f}%', va='center')
                
                st.pyplot(fig, use_container_width=True)
                
            else:
                st.error("Failed to fetch sequence details.")
else:
    st.info("👈 Start by searching for a snake in the sidebar.")
