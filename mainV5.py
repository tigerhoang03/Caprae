import streamlit as st
import pandas as pd
import numpy as np

# Set the page configuration
st.set_page_config(layout="wide")

# Cache data loading
@st.cache_data(ttl=3600)
def loadData():
    df = pd.read_excel("caprae.xlsx")
    # Check if 'Group' column exists, if not, create it with default group
    if 'Group' not in df.columns:
        df['Group'] = 'All Isolates'  # Default group name when no groups exist
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

# Check if there are multiple groups before showing group selection
unique_groups = list(df['Group'].dropna().unique())
unique_groups = [group for group in unique_groups if group != 'nan' and group != 'All Isolates' and pd.notna(group)]
has_multiple_groups = len(unique_groups) > 0

# Initialize session state for storing selections if it doesn't exist
if 'selected_isolates' not in st.session_state:
    st.session_state.selected_isolates = []
if 'selection_mode' not in st.session_state:
    st.session_state.selection_mode = 'Individual Selection'
if 'previous_mode' not in st.session_state:
    st.session_state.previous_mode = 'Individual Selection'

with st.sidebar:
    # Sidebar selection mode
    st.subheader("Selection Mode")
    selection_mode = st.radio(
        "Choose how to select isolates:",
        options=['Individual Selection', 'Group Selection'],
        key='selection_mode'
    )

    # Clear selections when switching modes
    if st.session_state.previous_mode != selection_mode:
        st.session_state.selected_isolates = []
        st.session_state.previous_mode = selection_mode

    if has_multiple_groups:
        groupOptions = ['All'] + unique_groups
    else:
        groupOptions = ['All Isolates']
        
    if selection_mode == 'Group Selection':
        # Group selection interface
        st.subheader("Select Groups")
        selected_groups = st.multiselect(
            "Choose Group(s)",
            options=groupOptions if has_multiple_groups else ['All Isolates'],
            default=[groupOptions[0]]
        )
        
        # Update selected isolates based on group selection
        if selected_groups:
            if 'All' in selected_groups:
                filtered_options = isolateOptions['Isolate Options'].tolist()
            else:
                filtered_options = []
                for group in selected_groups:
                    group_isolates = df[df['Group'] == group]["Unnamed: 0"].unique()
                    filtered_options.extend([opt for opt in group_isolates 
                                          if opt not in ['root', 'MQ', 'annotation']])
            st.session_state.selected_isolates = filtered_options
            
    else:  # Individual Selection
        if has_multiple_groups:
            # Group filter for individual selection
            st.sidebar.subheader("Filter isolates by group")
            filter_group = st.sidebar.selectbox("Filter by Group", groupOptions)
            
            # Filter isolates by selected group
            if filter_group != 'All':
                filtered_options = df[df['Group'] == filter_group]["Unnamed: 0"].unique()
                filtered_options = [option for option in filtered_options 
                                  if option not in ['root', 'MQ', 'annotation']]
            else:
                filtered_options = isolateOptions['Isolate Options'].tolist()
        else:
            filtered_options = isolateOptions['Isolate Options'].tolist()

        # Ensure currently selected isolates are always in the options
        all_options = list(set(filtered_options + st.session_state.selected_isolates))
        all_options.sort()

        # Create a form for individual isolate selections
        with st.sidebar.form("selection_form"):
            st.subheader("Select individual isolates")
            
            help_text = ("Select multiple isolates and click 'Apply Selections' when done" 
                        if has_multiple_groups 
                        else "Select isolates and click 'Apply Selections' when done")
            
            # Store selections temporarily
            temp_selections = st.multiselect(
                "Choose Isolate(s) (type to search)", 
                options=all_options,
                default=st.session_state.selected_isolates,
                help=help_text
            )
            
            # Add a submit button
            submitted = st.form_submit_button("Apply Selections")
            if submitted:
                st.session_state.selected_isolates = temp_selections

    # Select threshold for frequency (outside the form)
    threshold = st.slider("Threshold for displaying mutations", 0.0, 1.0, 0.5)

# Use the confirmed selections for analysis
selected_isolates = st.session_state.selected_isolates

def getSNPTypeStats(mutationSummaryDf):
    # Fill NaN with "Not annotated" before counting
    counts = mutationSummaryDf['Substitution'].fillna("Not annotated").value_counts()
    total = len(mutationSummaryDf)
    
    # Calculate percentages
    percentages = (counts/total * 100).round(1).astype(str) + '%'
    
    # Add total row
    counts['Total'] = total
    percentages['Total'] = '100%'
    
    return {
        'Count': counts,
        'Percentage': percentages
    }

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
            if ',' in annotationDetails:
                mutation, gene, locus, substitution = annotationDetails.split(',')
            else:
                mutation = annotationDetails
                gene = locus = substitution = "Not annotated"
        else:
            mutation = gene = locus = substitution = "Not annotated"
        
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
                'Locus': locus.strip(),
                'Substitution': substitution.strip()
            })
    return pd.DataFrame(mutationSummary)

# Function to identify mutation type
def get_mutation_type(mutation_name):
    """Identify the type of mutation based on its name."""
    mutation_name = mutation_name.lower()
    if 'non' in mutation_name or 'missense' in mutation_name:
        return 'nonsynonymous'
    if 'syn' in mutation_name or 'silent' in mutation_name:
        return 'synonymous'
    return 'other'

# Display results based on selections
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
    
    # Add SNP distribution statistics in an expander
    with st.expander("View SNP Type Distribution Across Groups"):
        st.subheader("Distribution of SNP Types in M. caprae Lineages")
        
        # Calculate stats for each group
        group_stats = {}
        all_substitution_types = set()
        
        # First pass: collect all unique substitution types across all groups
        for group in unique_groups:
            group_isolates = df[df['Group'] == group]["Unnamed: 0"].unique()
            group_isolates = [x for x in group_isolates if x not in ['root', 'MQ', 'annotation']]
            
            if group_isolates:
                group_df = df[df["Unnamed: 0"].isin(group_isolates)]
                # Now using threshold instead of 0.0
                group_mutations = getMutationSummary(group_df, root, group_isolates, threshold, annotationRow)
                stats = getSNPTypeStats(group_mutations)
                group_stats[group] = stats
                all_substitution_types.update(stats['Count'].index)
        
        # Create the summary table with all discovered substitution types
        if group_stats:
            # Organize types: regular types first, then Not annotated, then Total
            regular_types = sorted(list(all_substitution_types - {'Not annotated', 'Total'}))
            ordered_types = regular_types + ['Not annotated', 'Total']
            
            # Initialize DataFrame with substitution types as index
            columns = []
            data = []
            
            # Add count columns
            for group in unique_groups:
                col_name = f"{group} Count"
                columns.append(col_name)
                counts = group_stats[group]['Count']
                data.append([counts.get(sub_type, 0) for sub_type in ordered_types])
            
            # Add percentage columns
            for group in unique_groups:
                col_name = f"{group} %"
                columns.append(col_name)
                percentages = group_stats[group]['Percentage']
                data.append([percentages.get(sub_type, '0.0%') for sub_type in ordered_types])
            
            # Create DataFrame
            stats_df = pd.DataFrame(data, columns=ordered_types, index=columns).T
            st.write(stats_df)
            
            # Calculate and display ratio table
            st.subheader("Non-synonymous to Synonymous Mutation Ratios")
            
            # Get the actual mutation types from the data (excluding 'Not annotated' and 'Total')
            mutation_types = [t for t in stats_df.index if t not in ['Not annotated', 'Total']]
            
            # Group mutation types
            nonsynonymous_types = [t for t in mutation_types if get_mutation_type(t) == 'nonsynonymous']
            synonymous_types = [t for t in mutation_types if get_mutation_type(t) == 'synonymous']
            
            ratio_data = []
            for group in unique_groups:
                # Sum counts for each category
                nonsynonymous_count = sum(float(stats_df.loc[t, f"{group} Count"]) 
                                        for t in nonsynonymous_types)
                synonymous_count = sum(float(stats_df.loc[t, f"{group} Count"]) 
                                     for t in synonymous_types)
                
                ratio = round(nonsynonymous_count / synonymous_count, 2) if synonymous_count > 0 else float('inf')
                
                ratio_data.append({
                    'Lineage': group,
                    'Non-synonymous': int(nonsynonymous_count),
                    'Synonymous': int(synonymous_count),
                    'Ratio (dN/dS)': ratio
                })
                
            ratio_df = pd.DataFrame(ratio_data)
            st.write(ratio_df)
        else:
            st.write("No group statistics available.")
else:
    st.write("No isolates selected.")