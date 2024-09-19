import streamlit as st
import pandas as pd

# Load the data
df = pd.read_excel("caprae.xlsx")
root_sequence = df.iloc[0, 2:]  # Extract root sequence from the first row and from the third column onward

# Collapsible section for Root Sequence
with st.expander("Root Sequence:"):
    st.dataframe(pd.DataFrame([root_sequence]))  # Display the root sequence

# Collapsible section for Full Data
with st.expander("Full Data:"):
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
    # st.write(f"Selected Individual Isolates: {selected_isolates}")
    
    # Filter and display the sequences for the selected isolates
    filtered_df_individual = df[df['Unnamed: 0'].isin(selected_isolates)].iloc[:, 2:]  # Sequence columns only
    
    # Collapsible section for selected isolates' data
    with st.expander(f"Filtered Data for Selected Isolates ({len(selected_isolates)} isolates):"):
        st.dataframe(filtered_df_individual)

    # Calculate frequency of mismatches with the root sequence for selected isolates
    mismatches_individual = (filtered_df_individual != root_sequence).sum(axis=0)  # Count mismatches per column
    total_isolates_individual = filtered_df_individual.shape[0]  # Number of selected isolates
    
    # Calculate fraction of mismatches (mismatches/total selected isolates)
    mismatch_fraction_individual = mismatches_individual / total_isolates_individual
    
    # Collapsible section for mismatch frequency for the selected individual isolates
    with st.expander(f"Mismatch Frequency for Selected Individual Isolates ({total_isolates_individual} total):"):
        mismatch_individual_df = pd.DataFrame({
            'Mismatches': mismatches_individual,
            # 'Total Isolates': total_isolates_individual,
            'Mismatch Fraction': mismatch_fraction_individual
        })
        st.dataframe(mismatch_individual_df)
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
    
    # Collapsible section for group data
    with st.expander(f"Filtered Data for Group {selected_group}:"):
        st.dataframe(filtered_df_group)
    
    # Calculate frequency of mismatches with the root sequence for the selected group
    mismatches_group = (filtered_df_group != root_sequence).sum(axis=0)  # Count mismatches per column
    total_isolates_in_group = filtered_df_group.shape[0]  # Number of isolates in the group
    
    # Calculate fraction of mismatches (mismatches/total isolates in the group)
    mismatch_fraction_group = mismatches_group / total_isolates_in_group
    
    # Collapsible section for mismatch frequency for the group
    with st.expander(f"Mismatch Frequency for Selected Group ({total_isolates_in_group} total):"):
        mismatch_group_df = pd.DataFrame({
            'Mismatches': mismatches_group,
            # 'Total Isolates': total_isolates_in_group,
            'Mismatch Fraction': mismatch_fraction_group
        })
        st.dataframe(mismatch_group_df)
    
    # Now, calculate for the entire dataset (all isolates)
    all_isolates_df = df.iloc[1:, 2:]  # Exclude the root and take only sequence columns
    total_isolates = all_isolates_df.shape[0]
    mismatches_all = (all_isolates_df != root_sequence).sum(axis=0)  # Count mismatches for all isolates
    mismatch_fraction_all = mismatches_all / total_isolates
    
    # Collapsible section for mismatch frequency for all isolates
    with st.expander(f"Mismatch Frequency for All {total_isolates} Isolates:"):
        mismatch_all_df = pd.DataFrame({
            'Mismatches': mismatches_all,
            # 'Total Isolates': total_isolates,
            'Mismatch Fraction': mismatch_fraction_all
        })
        st.dataframe(mismatch_all_df)
else:
    st.write("Please select a group from the sidebar.")
