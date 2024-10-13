import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(layout="wide")

df = pd.read_excel("caprae.xlsx")
df = df.apply(lambda col: col.astype(str) if col.dtype == 'object' else col)

with st.expander("Full Data:"):
    st.write(df)

root = df.iloc[0, 2:]
rootSeq = pd.DataFrame(root).T 
st.write("Root Sequence", rootSeq) #this converts the sequence to a numpy array

#extract all isolates
isolateOptions = df["Unnamed: 0"].unique()
isolateOptions = pd.DataFrame(isolateOptions[1:262], columns=['Isolate Options'])  # Adjust slicing dynamically if needed, this removes 'root' as in this column root is an element
# st.write("Isolate options", isolateOptions)

groupOptions = df["Group"].unique()
groupOptions = pd.DataFrame(groupOptions[1:262], columns=['Isolate Groups'])
# st.write("Group options", groupOptions)

#"""--------------------------Sidebar selection -------------------------------"""
st.sidebar.title("Selection")
st.sidebar.subheader("Select individual isolates")

#"isolate options is the column name of the isolateOptions Dataframe"
selectedIsolates = st.sidebar.multiselect("Choose Isolate(s)", isolateOptions['Isolate Options']) 

if selectedIsolates:
    #basically checks the first column of the originial df against the selectedIsolates selected by user and it makes a dataframe from the filtered rows selecting all rows and all columns after 2
    #.isin(selectedIsolates): checks whether the values in the 'Unnamed: 0' column are in the selectedIsolates list. 
    # It returns a boolean Series (True or False for each row) which rows correspond to the isolates selected by the user.
    #df['Unnamed: 0'].isin(selectedIsolates) returns the true or fall sequence (vertical sequence remember)
    #the second df is what returns the actual sequence and is why we use .iloc to get the sequences
    isolateSeq = df[df['Unnamed: 0'].isin(selectedIsolates)]
    markerCol = isolateSeq.columns[2:]
    
    # Adjusting column widths to use more space
    col1, col2 = st.columns([5, 5])  # Increase the width ratio to use more space horizontally

    # SEARCH FUNCTIONALITY FOR `isolateSeq` DATAFRAME
    with col1:
        search_isolate = st.text_input("Search Selected Isolates (enter column name or part of name):")
        if search_isolate:
            # Filter columns of isolateSeq DataFrame based on search term
            matching_isolate_columns = [col for col in isolateSeq.columns if search_isolate.lower() in col.lower()]
            if matching_isolate_columns:
                st.subheader(f"Filtered Data for Selected Isolates ({len(isolateSeq)} isolates) - Matching '{search_isolate}'")
                st.dataframe(isolateSeq[matching_isolate_columns])
            else:
                st.write(f"No columns matching '{search_isolate}' found.")
        else:
            # Display the full isolateSeq DataFrame if no search is applied
            st.subheader(f"Filtered Data for Selected Isolates ({len(isolateSeq)} isolates)")
            st.dataframe(isolateSeq)

    # st.write(markerCol)
    differences = pd.DataFrame()
    
    for marker in markerCol:
        # Get the values for this marker across the selected isolates
        isolateValues = isolateSeq[marker] #this is querying every marker and getting the columns for every iteration
        
        # If there's more than one unique value, it means there are differences
        if len(isolateValues.unique()) > 1:
            differences[marker] = isolateValues.values #appends the columns with differences between all selected isolates
    
    #filters out the columns of the  (marker sites) where there are actual differences
    # rootFiltered = df[differences.columns]
    
    dfFiltered = df.loc[263, differences.columns] #selcts rows and columns by label, 1 is considered a row label
    differences.loc["Annotations"] = dfFiltered.values #adds the values of all the different columns to the differences "annotations" row
    
    # Add an "Isolate" column to label the rows (selected isolates and annotation)
    differences.insert(0, 'Isolate', list(isolateSeq['Unnamed: 0']) + ['Annotation'])

    # SEARCH FUNCTIONALITY FOR `differences` DATAFRAME
    with col2:
        search_diff = st.text_input("Search Differences table (enter column name or part of name):")
        if search_diff:
            # Filter columns of differences DataFrame based on search term
            matching_diff_columns = [col for col in differences.columns if search_diff.lower() in col.lower()]
            if matching_diff_columns:
                st.subheader(f"Differences in selected Isolates (Including Annotations) - Matching '{search_diff}'")
                st.dataframe(differences[["Isolate"] + matching_diff_columns])  # Include Isolate column and matching columns
            else:
                st.write(f"No columns matching '{search_diff}' found.")
        else:
            # Display the full differences DataFrame if no search is applied
            st.subheader("Differences in selected Isolates (Including Annotations)")
            st.dataframe(differences)

### Calculate Frequency of Nucleotide Bases ###
    # We exclude the "Annotation" row for the frequency calculation
    freq_dict = {}  # Dictionary to store frequency data for all markers

    for marker in differences.columns[1:]:  # Skip the 'Isolate' column
        # Get the nucleotide counts (excluding the "Annotation" row)
        base_counts = differences.loc[differences['Isolate'] != 'Annotation', marker].value_counts(normalize=True)
        base_counts = base_counts * 100  # Convert to percentage
        # Store the base counts in a dictionary
        freq_dict[marker] = base_counts.reindex(['A', 'T', 'C', 'G'], fill_value=0)  # Fill missing bases with 0
    
    # Use pd.concat to join all the columns at once
    freq_df = pd.concat(freq_dict, axis=1)

    # Transpose the frequency DataFrame for better readability (markers as rows, bases as columns)
    freq_df = freq_df.T

    # Display frequency and differences using tabs
    tab1, tab2 = st.tabs(["Intra-Base Frequency", "Differences Table"])

    with tab1:
        st.subheader("Nucleotide Base Frequency (%)")
        st.write("This excludes non-base pair in our calculation")
        st.dataframe(freq_df)

    with tab2:
        st.subheader("Differences in selected Isolates and Location of mutation (Including Annotations)")
        st.dataframe(differences)
