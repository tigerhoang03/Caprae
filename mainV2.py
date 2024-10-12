import streamlit as st
import pandas as pd


def frequency(df):
    for i in range(len(df)):
        st.write(pd.DataFrame(df.iloc[i]).T)
        st.write(rootSeq)
        st.write(df.iloc[:, 2:])
        rootSeq = rootSeq.squeeze()
        freqDiff = (df.iloc[:, 2:] != rootSeq).sum(axis=0) 
    #sums column wise differences for the selected isolates
    return freqDiff
    
# Load and display data
df = pd.read_excel("caprae.xlsx")
with st.expander("Data"):
    st.write(df)

# Extract the root sequence
rootSeq = pd.DataFrame(df.iloc[0, 2:]).T 
st.write("Root Sequence:", rootSeq)

# Slices unnamed column for unique isolates
isolateOpt = df["Unnamed: 0"].unique()
# Filter out non-isolates and set column to "Isolate Options"
isolateOpt = pd.DataFrame(isolateOpt[1:262], columns=["Isolate Options"])

# Sidebar menu for selecting isolates
st.sidebar.title("Isolate Selection")
selected_isolates = st.sidebar.multiselect("Select Isolates to Compare", isolateOpt["Isolate Options"])

# Filter the dataframe based on the selected isolates
if selected_isolates:
    #checks if each "unnamed: 0" values are in the selected_isolates list
    selected_isolate_df = df[df["Unnamed: 0"].isin(selected_isolates)] 
    
    # Display selected isolates
    st.write("Selected Isolates", selected_isolate_df)
    
    freqDiff = frequency(selected_isolate_df)
    
    
    
else:
    st.write("No isolates selected. Please select isolates to compare.")
