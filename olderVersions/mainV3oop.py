import streamlit as st
import pandas as pd
import numpy as np

class NucleotideAnalysisApp:
    def __init__(self, df):
        self.df = df.apply(lambda col: col.astype(str) if col.dtype == 'object' else col)
        self.differences = pd.DataFrame()
        self.isolateSeq = pd.DataFrame()
        self.freq_df = pd.DataFrame()
    
    def display_full_data(self):
        """Displays the full data"""
        with st.expander("Full Data:"):
            st.write(self.df)
    
    def display_root_sequence(self):
        """Displays the root sequence (first row)"""
        root = self.df.iloc[0, 2:]
        rootSeq = pd.DataFrame(root).T
        st.write("Root Sequence", rootSeq) # this converts the sequence to a numpy array

    def sidebar_selection(self):
        """Handles sidebar selection for isolates"""
        isolateOptions = self.df["Unnamed: 0"].unique()
        isolateOptions = pd.DataFrame(isolateOptions[1:262], columns=['Isolate Options'])
        
        st.sidebar.title("Selection")
        st.sidebar.subheader("Select individual isolates")
        selectedIsolates = st.sidebar.multiselect("Choose Isolate(s)", isolateOptions['Isolate Options'])
        return selectedIsolates
    
    def filter_isolates(self, selectedIsolates):
        """Filters isolates based on user selection"""
        if selectedIsolates:
            self.isolateSeq = self.df[self.df['Unnamed: 0'].isin(selectedIsolates)]
            markerCol = self.isolateSeq.columns[2:]
            
            with st.expander(f"Filtered Sequence for Selected Isolates ({len(self.isolateSeq)} isolates):"):
                st.dataframe(self.isolateSeq)

            # Find markers with differences across isolates
            for marker in markerCol:
                isolateValues = self.isolateSeq[marker]
                if len(isolateValues.unique()) > 1:
                    self.differences[marker] = isolateValues.values

    def add_annotations(self):
        """Adds the annotation row to the differences DataFrame"""
        dfFiltered = self.df.loc[263, self.differences.columns]  # Selects annotation row by label
        self.differences.loc["Annotations"] = dfFiltered.values  # Adds annotation row
        self.differences.insert(0, 'Isolate', list(self.isolateSeq['Unnamed: 0']) + ['Annotation'])

        # Display filtered markers with differences and annotations
        with st.expander("Differences in selected Isolates and Location of mutation (Including Annotations)"):
            st.dataframe(self.differences)

    def calculate_frequency(self):
        """Calculates the frequency of nucleotide bases"""
        self.freq_df = pd.DataFrame()

        for marker in self.differences.columns[1:]:  # Skip the 'Isolate' column
            base_counts = self.differences.loc[self.differences['Isolate'] != 'Annotation', marker].value_counts(normalize=True)
            base_counts = base_counts * 100  # Convert to percentage
            self.freq_df[marker] = base_counts.reindex(['A', 'T', 'C', 'G'], fill_value=0)  # Include only A, T, C, G
        
        self.freq_df = self.freq_df.T  # Transpose for better readability
        with st.expander("Base Frequency Across Selected Isolates"):
            st.write("Nucleotide Base Frequency (%) for Markers with Variations (excluding non-base pairs)")
            st.dataframe(self.freq_df)

def main():
    # Load the Excel file
    df = pd.read_excel("caprae.xlsx")
    
    # Instantiate the NucleotideAnalysisApp class
    app = NucleotideAnalysisApp(df)

    # Display the full data and root sequence
    app.display_full_data()
    app.display_root_sequence()

    # Sidebar selection of isolates
    selectedIsolates = app.sidebar_selection()

    # Filter isolates and calculate differences
    if selectedIsolates:
        app.filter_isolates(selectedIsolates)
        app.add_annotations()
        app.calculate_frequency()

if __name__ == "__main__":
    main()
