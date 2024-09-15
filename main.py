import streamlit as st
import pandas as pd

# Load the data
df = pd.read_excel("caprae.xlsx")
root_sequence = df.iloc[0, 2:]  # Extract root sequence from the first row and from the third column onward
st.write("Root Sequence:")
st.dataframe(pd.DataFrame([root_sequence]))  # Display the root sequence

# Display full data for reference
st.write("Full Data: ")
st.dataframe(df)

# Extract isolate options
isolate_options = df['Unnamed: 0'].unique()
isolate_options = pd.DataFrame(isolate_options[1:262], columns=['Isolate Options'])  # Adjust slicing dynamically if needed

# Sidebar for selecting isolates and groups
st.sidebar.title("Selection Options")

# Section 1: Individual isolate selection
st.sidebar.subheader("Select Individual Isolates")
selected_isolates = st.sidebar.multiselect("Choose Isolate(s)", isolate_options['Isolate Options'])

# Display the selected isolate names
if selected_isolates:
    st.write(f"Selected Individual Isolates: {selected_isolates}")
    
    # Filter and display the sequences for the selected isolates
    filtered_df_individual = df[df['Unnamed: 0'].isin(selected_isolates)].iloc[:, 2:]  # Sequence columns only
    st.write(f"Filtered Data for Selected Isolates:")
    st.dataframe(filtered_df_individual)
else:
    st.write("Please select isolates from the sidebar.")

# Section 2: Group selection
st.sidebar.subheader("Select Group")
group_options = df['Group'].dropna().unique()
selected_group = st.sidebar.selectbox("Choose Group", group_options)

# Display the selected group's sequences
if selected_group:
    st.write(f"Selected Group: {selected_group}")
    
    # Filter the DataFrame for the selected group
    filtered_df_group = df[df['Group'] == selected_group].iloc[:, 2:]  # Sequence columns only
    st.write(f"Filtered Data for Group {selected_group}:")
    st.dataframe(filtered_df_group)
    
    # Calculate frequency of mismatches with the root sequence
    mismatches = (filtered_df_group != root_sequence).sum(axis=0)  # Count mismatches per column
    total_isolates_in_group = filtered_df_group.shape[0]  # Number of isolates in the group
    
    # Calculate fraction of mismatches (mismatches/total isolates in the group)
    mismatch_fraction = mismatches / total_isolates_in_group
    
    # Display the mismatch frequency for the group
    st.write("Mismatch Frequency for Selected Group:")
    mismatch_df = pd.DataFrame({
        'Marker': root_sequence.index,
        'Mismatches': mismatches,
        'Total Isolates': total_isolates_in_group,
        'Mismatch Fraction': mismatch_fraction
    })
    st.dataframe(mismatch_df)
    
    # Now, calculate for the entire dataset (all isolates)
    all_isolates_df = df.iloc[1:, 2:]  # Exclude the root and take only sequence columns
    total_isolates = all_isolates_df.shape[0]
    mismatches_all = (all_isolates_df != root_sequence).sum(axis=0)
    mismatch_fraction_all = mismatches_all / total_isolates
    
    # Display the mismatch frequency for all isolates
    st.write("Mismatch Frequency for All Isolates:")
    mismatch_all_df = pd.DataFrame({
        'Marker': root_sequence.index,
        'Mismatches': mismatches_all,
        'Total Isolates': total_isolates,
        'Mismatch Fraction': mismatch_fraction_all
    })
    st.dataframe(mismatch_all_df)
else:
    st.write("Please select a group from the sidebar.")
