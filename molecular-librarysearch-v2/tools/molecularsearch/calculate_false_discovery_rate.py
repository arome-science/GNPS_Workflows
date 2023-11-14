"""
Based on Kerstin Scheubert et al, "Significance estimation for large scale metabolomics annotations by spectral matching"
https://www.nature.com/articles/s41467-017-01318-5
"""

import argparse
import pandas as pd


def calculate_false_discovery_rate(input_filename: str, output_filename: str):
    annotations = pd.read_csv(input_filename, header=0, sep='\t')
    scans = annotations['#Scan#'].unique()
    for scan in scans:
        selection = annotations['#Scan#'] == scan
        library_names = annotations.loc[selection, 'LibraryName']
        false_discovery_rates = pd.Series(index=library_names.index)
        for library in library_names.unique():
            if library.lower().startswith('decoy'):
                decoy_annotations = library_names[library_names == library]
                target_and_decoy_selection = library_names.apply(lambda x: x.endswith(library[6:]))
                target_and_decoy_annotations = library_names[target_and_decoy_selection]
                if len(target_and_decoy_annotations) > 0:
                    false_discovery_rate = len(decoy_annotations) / len(target_and_decoy_annotations)
                    false_discovery_rates[target_and_decoy_selection] = false_discovery_rate

        annotations.loc[selection, 'FalseDiscoveryRate'] = false_discovery_rates

    annotations.to_csv(output_filename, index=False, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add False Discovery Rate values')
    parser.add_argument('--input', help='TSV file with annotations', required=True)
    parser.add_argument('--output', help='TSV file with FDR', required=True)
    args = parser.parse_args()
    calculate_false_discovery_rate(args.input, args.output)
