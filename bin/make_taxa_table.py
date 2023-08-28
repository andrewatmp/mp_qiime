#!/usr/bin/env python3

import plotly.figure_factory as ff
import pandas as pd
import argparse
import plotly.graph_objects as go

def extract_genus_species(taxonomic_string):
    elements = taxonomic_string.split(';')
    genus = species = None
    for element in elements:
        if element.startswith("g__"):
            genus = element[3:]
        elif element.startswith("s__"):
            species = element[3:]
    return f"{genus} {species}" if genus and species else None


def make_taxa_table(input_file, output_file):

    # Load csv file into dataframe
    df = pd.read_csv(input_file)

    # Transpose dataframe
    df_transposed = df.transpose()

    # Reset the index for easier manipulation
    df_transposed.reset_index(inplace=True)

    # Rename the columns for clarity
    df_transposed.columns = ['Taxonomic_Classification', 'Read_Count']

    # Drop the first row as it is a placeholder ("subseq")
    df_transposed = df_transposed.drop(0)

    # Extract genus and species
    df_transposed['Genus_Species'] = df_transposed['Taxonomic_Classification'].apply(extract_genus_species)

    # Drop rows where Genus_Species could not be extracted (i.e., is None)
    df_transposed.dropna(subset=['Genus_Species'], inplace=True)

    # Convert the 'Read_Count' column to numeric for calculations
    df_transposed['Read_Count'] = pd.to_numeric(df_transposed['Read_Count'], errors='coerce')

    # Calculate the total read count
    total_read_count = df_transposed['Read_Count'].sum()

    # Calculate the relative abundance
    df_transposed['Relative_Abundance'] = df_transposed['Read_Count'] / total_read_count * 100

    # Create a Plotly table
    fig = go.Figure(data=[go.Table(
    columnwidth = [2,1],
    header=dict(values=['Species Name', 'Relative Abundance (%)'],
                line_color='darkslategray',
                fill_color='lightskyblue',
                align='center',
                font=dict(color='black', size=14)),
    cells=dict(values=[df_transposed['Genus_Species'], df_transposed['Relative_Abundance'].round(2)],
                line_color='darkslategray',
                fill_color='white',
                align='center',
                font=dict(color='darkslategray', size=12)))
    ])

    # Update layout and title
    fig.update_layout(title='Taxonomic Classification: Relative Abundance')
    fig.update_layout(width=800)
    fig.write_html(output_file)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create AMR table.')
    parser.add_argument('input_file', type=str, help='Input file name')
    parser.add_argument('output_file', type=str, help='Output file name')
    args = parser.parse_args()

    make_taxa_table(args.input_file, args.output_file)


