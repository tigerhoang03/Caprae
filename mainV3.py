import streamlit as st
import pandas as pd
import numpy as np

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
    
    with st.expander(f"Filtered Data for Selected Isolates ({len(isolateSeq)} isolates):"):
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
    # Display the filtered markers with differences and the annotation
    with st.expander(f"Differences in selected Isolates and Location of mutation (Including Annotations)"):
        st.dataframe(differences)
        
        
### Calculate Frequency of Nucleotide Bases ###
    # We exclude the "Annotation" row for the frequency calculation
    freq_df = pd.DataFrame()

    for marker in differences.columns[1:]:  # Skip the 'Isolate' column
        # Get the nucleotide counts (excluding the "Annotation" row)
        base_counts = differences.loc[differences['Isolate'] != 'Annotation', marker].value_counts(normalize=True)
        
        # Normalize to get the frequency (i.e., # of bases / total bases)
        base_counts = base_counts * 100  # Convert to percentage

        # Include only A, T, C, G in the output
        freq_df[marker] = base_counts.reindex(['A', 'T', 'C', 'G'], fill_value=0)  # Fill missing bases with 0
    
    # Transpose the frequency DataFrame for better readability (markers as rows, bases as columns)
    freq_df = freq_df.T

    # Display the frequency table in a new dynamic window
    with st.expander("Base Frequency Across Selected Isolates"):
        st.write("Nucleotide Base Frequency (%) for Markers with Variations (excluding non-base paris)")
        st.dataframe(freq_df)
 

