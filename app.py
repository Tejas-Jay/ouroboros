import streamlit as st
import utils  # We import the file we just finished!
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# PAGE CONFIGURATION
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Venom Profiler",
    page_icon="🐍",
    layout="centered"
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
    # .keys() gives the clean names we made: "Mambalgin-1 (P0DKR6.1)"
    selected_label = st.selectbox(
        "Choose from the list:", 
        options=st.session_state.toxin_results.keys()
    )
    
    # Get the actual ID from the dictionary using the label
    selected_id = st.session_state.toxin_results[selected_label]
    
    st.info(f"You selected ID: **{selected_id}**")
    
    # Placeholder for Sprint 3 (The Analysis)
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
                st.text_area("Sequence:", record.seq, height=100)
                
                # --- DISPLAY: METRICS ---
                # Create 3 columns for a cool layout
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.metric("Molecular Weight", f"{stats['molecular_weight']:.2f} Da")
                with col2:
                    st.metric("Isoelectric Point (pI)", f"{stats['isoelectric_point']:.2f}")
                with col3:
                    # Instability < 40 is stable
                    stability = "Stable" if stats['instability_index'] < 40 else "Unstable"
                    st.metric("Stability", f"{stats['instability_index']:.2f} ({stability})")

                # --- DISPLAY: VISUALIZATION (Pie Chart) ---
                # --- DISPLAY: VISUALIZATION (Bar Chart) ---
                st.markdown("### 🧪 Amino Acid Composition")
                
                # Prepare data for the chart
                aa_data = stats['amino_acid_percent']
                
                # Sort amino acids by percentage for better readability
                sorted_aa = sorted([(k, v * 100) for k, v in aa_data.items() if v > 0], 
                                   key=lambda item: item[1], reverse=True)
                
                labels = [item[0] for item in sorted_aa]
                percentages = [item[1] for item in sorted_aa]
                
                # Create the bar plot using Matplotlib
                fig, ax = plt.subplots(figsize=(10, 6)) # Make the figure a bit wider
                
                # Create horizontal bar chart
                ax.barh(labels, percentages, color='skyblue') 
                
                # Add labels and title
                ax.set_xlabel("Percentage (%)")
                ax.set_title("Amino Acid Composition of Toxin")
                ax.invert_yaxis() # Put highest percentage at the top
                
                # Add percentage values on the bars
                for index, value in enumerate(percentages):
                    ax.text(value + 0.5, index, f'{value:.1f}%', va='center') # Adjust x for label position
                
                # Show it in Streamlit, ensuring layout is tight
                st.pyplot(fig, use_container_width=True) # use_container_width makes it responsive
                
            else:
                st.error("Failed to fetch sequence details.")
else:
    st.info("👈 Start by searching for a snake in the sidebar.")