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
        num_target_libraries = 0
        num_decoy_libraries = 0
        for library in library_names:
            if library.lower().startswith('decoy'):
                num_decoy_libraries += 1
            else:
                num_target_libraries += 1

        false_discovery_rate = num_decoy_libraries / (num_target_libraries + num_decoy_libraries)\
            if num_target_libraries > 0 or num_decoy_libraries > 0 else None
        annotations.loc[selection, 'FalseDiscoveryRate'] = false_discovery_rate

    annotations.to_csv(output_filename, index=False, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add False Discovery Rate values')
    parser.add_argument('--input', help='TSV file with annotations', required=True)
    parser.add_argument('--output', help='TSV file with FDR', required=True)
    args = parser.parse_args()
    calculate_false_discovery_rate(args.input, args.output)
