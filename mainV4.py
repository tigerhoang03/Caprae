import streamlit as st
import pandas as pd
import numpy as np

# Set the page configuration
st.set_page_config(layout="wide")

# Cache data loading
@st.cache_data
def loadData():
    df = pd.read_excel("caprae.xlsx")
    return df.apply(lambda col: col.astype(str) if col.dtype == 'object' else col)

df = loadData()

# Extract the annotation row for detailed mutation information
annotationRow = df[df['Unnamed: 0'] == 'annotation']

# Show the full data in an expander section
with st.expander("Full Data:"):
    st.write(df)

# Extract root sequence
root = df.iloc[0, 2:]
rootSeq = pd.DataFrame(root).T
st.write("Root Sequence", rootSeq)

# Extract all isolates options and group options
isolateOptions = df["Unnamed: 0"].unique()
isolateOptions = [option for option in isolateOptions if option not in ['root', 'MQ', 'annotation']]
isolateOptions = pd.DataFrame(isolateOptions, columns=['Isolate Options'])

# Extract group options
groupOptions = ['All'] + list(df['Group'].dropna().unique()) if 'Group' in df.columns else ['All']
groupOptions = [group for group in groupOptions if group != 'nan']

# Initialize session state for storing selections if it doesn't exist
if 'selected_isolates' not in st.session_state:
    st.session_state.selected_isolates = []

# Group selection outside the form (will auto-update)
st.sidebar.subheader("Filter isolates by group")
selectedGroup = st.sidebar.selectbox("Choose Group", groupOptions)

# Filter isolates by selected group
if selectedGroup != 'All' and 'Group' in df.columns:
    filtered_options = df[df['Group'] == selectedGroup]["Unnamed: 0"].unique()
    filtered_options = [option for option in filtered_options if option not in ['root', 'MQ', 'annotation']]
else:
    filtered_options = isolateOptions['Isolate Options'].tolist()

# Ensure currently selected isolates are always in the options
all_options = list(set(filtered_options + st.session_state.selected_isolates))
all_options.sort()  # Keep the list ordered

# Create a form for just the isolate selections
with st.sidebar.form("selection_form"):
    st.subheader("Select individual isolates")
    # Store selections temporarily
    temp_selections = st.multiselect(
        "Choose Isolate(s) (type to search)", 
        options=all_options,
        default=st.session_state.selected_isolates,
        help="Select multiple isolates and click 'Apply Selections' when done"
    )
    
    # Add a submit button
    submitted = st.form_submit_button("Apply Selections")
    if submitted:
        st.session_state.selected_isolates = temp_selections

# Select threshold for frequency (outside the form)
threshold = st.sidebar.slider("Threshold for displaying mutations", 0.0, 1.0, 0.5)

# Use the confirmed selections for analysis
selected_isolates = st.session_state.selected_isolates

# Cache mutation summary computation
@st.cache_data
def getMutationSummary(filteredDf, root, selectedIsolates, threshold, annotationRow):
    mutationSummary = []
    for col in root.index:
        rootBase = root[col]
        
        selectedBases = filteredDf[col].value_counts()
        baseSummary = []
        otherBasesSummary = []
        
        # Calculate the frequency of each base and compare to the threshold
        for base, count in selectedBases.items():
            if base not in ['A', 'C', 'T', 'G']:
                continue
            frequency = count / len(selectedIsolates)
            if frequency >= threshold and base != rootBase:  # Append to the summary if base is a mutation from root
                baseSummary.append((base, frequency, count))
            else:  # Track bases that didn't meet the threshold or did not change from root
                otherBasesSummary.append((base, frequency, count))
        
        # Extract annotation details if available
        if col in annotationRow.columns:
            annotationDetails = annotationRow[col].values[0]
            mutation, gene, substitution, locus = annotationDetails.split(',') if ',' in annotationDetails else (annotationDetails, "Not annotated", "Not annotated", "Not annotated")
        else:
            mutation = gene = substitution = locus = "Not annotated"
        
        if len(baseSummary) > 0:
            mutantBase = ", ".join([base for base, _, _ in baseSummary])
            otherBasesFreq = ", ".join([f"{base}: {freq:.2f} (Count: {count})" for base, freq, count in otherBasesSummary])
            
            mutationSummary.append({
                'Location': col,
                'Root Base': rootBase,
                'Mutant Base': mutantBase,
                'Frequency': ", ".join([f"{freq:.2f}" for _, freq, _ in baseSummary]),
                'Count': ", ".join([str(count) for _, _, count in baseSummary]),
                'Other Bases Below Threshold (Frequency and Count)': otherBasesFreq,
                'Mutation': mutation.strip(),
                'Gene': gene.strip(),
                'Substitution': substitution.strip(),
                'Locus': locus.strip()
            })
    return pd.DataFrame(mutationSummary)

# Filter the selected isolates data
if selected_isolates:
    filteredDf = df[df["Unnamed: 0"].isin(selected_isolates)]

    # Add a search bar to filter MTBC0 positions
    searchPosition = st.text_input("Search by MTBC0 position:")
    if searchPosition:
        filteredColumns = ['Unnamed: 0'] + [col for col in filteredDf.columns if searchPosition in col]
        if len(filteredColumns) > 1:
            st.write(f"Highlighted Data by MTBC0 Position ({searchPosition}):", filteredDf[filteredColumns])
        else:
            st.write(f"No matching MTBC0 positions found for: {searchPosition}")
    st.write("Filtered Isolates Data:", filteredDf)

    # Get mutation summary
    mutationSummaryDf = getMutationSummary(filteredDf, root, selected_isolates, threshold, annotationRow)
    st.write(f"Mutation Summary: ({len(mutationSummaryDf)} rows displayed)", mutationSummaryDf)
else:
    st.write("No isolates selected.")