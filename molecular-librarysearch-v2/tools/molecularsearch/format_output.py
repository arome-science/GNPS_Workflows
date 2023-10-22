import argparse
import pandas as pd

from os.path import basename
from spectrum.spectrum import read_mgf
from spectrum.utils import expand_directories
from typing import List


def filter_by_library(data: pd.DataFrame) -> pd.DataFrame:
    common_columns = ['#Scan#', 'Charge', 'SpecMZ', 'SpectrumFile', 'p-value']
    annotation_columns = ['Peptide', 'ExactMass', 'SharedPeaks', 'SpectrumID', 'Score', 'MassDiff', 'Protein',
                          'mzErrorPPM', 'FalseDiscoveryRate', 'IonMode', 'Instrument', 'Adduct', 'INCHI']

    filtered_rows = []
    for _, annotation_group in data.groupby(by=['#Scan#']):
        row = annotation_group.iloc[0][common_columns]
        unique_libraries = annotation_group['LibraryName'].unique()
        for library in unique_libraries:
            library_group = annotation_group[annotation_group['LibraryName'] == library] if isinstance(library, str) else annotation_group
            sorted_group = library_group.sort_values(by='MQScore', ascending=False)
            for column in annotation_columns:
                # library_name = rename_nist_library(library)
                row[f'{library}: {column}'] = sorted_group.iloc[0][column]
        filtered_rows.append(row)

    # filtered_data = pd.concat(filtered_rows, ignore_index=True)
    filtered_data = pd.DataFrame(data=filtered_rows)
    return filtered_data


def add_library_info(annotations: pd.DataFrame, library_filename: str):
    library = read_mgf(library_filename, id_field='SPECTRUMID')
    library_name = basename(library_filename)
    for index, annotation in annotations[annotations['LibraryName'] == library_name].iterrows():
        spectrum_id = annotation['SpectrumID']
        if not spectrum_id:
            continue

        spectrum = library.get(spectrum_id)
        if not spectrum:
            continue

        properties = spectrum.properties
        annotations.loc[index, 'IonMode'] = properties.get('IONMODE', properties.get('ION_MODE'))
        annotations.loc[index, 'Instrument'] = properties.get('SOURCE_INSTRUMENT', properties.get('INSTRUMENT'))
        annotations.loc[index, 'Adduct'] = properties.get('PRECURSOR_TYPE')
        annotations.loc[index, 'INCHI'] = properties.get('INCHIKEY')
        annotations.loc[index, 'ExactMass'] = properties.get('EXACTMASS')


def format_output(annotations_filename: str, library_filenames: List[str], output_filename: str):
    annotations = pd.read_csv(annotations_filename, header=0, sep='\t')
    annotations.rename(mapper={
        'LibrarySpectrumID': 'SpectrumID',
        'CompoundName': 'Peptide',
        'ParentMassDiff': 'MassDiff',
        'LibSearchSharedPeaks': 'SharedPeaks',
        'MQScore': 'Score'
    }, axis=1, inplace=True)

    for library_filename in library_filenames:
        add_library_info(annotations, library_filename)

    annotations = filter_by_library(annotations)

    annotations.to_csv(output_filename, index=False, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate annotations in the GNPS format')
    parser.add_argument('--annotations', help='TSV file with MSPepSearch annotations', required=True)
    parser.add_argument('--libraries', help='MSP files with library spectra', nargs='+', required=True)
    parser.add_argument('--output', help='Output filename', required=True)
    args = parser.parse_args()
    args.libraries = expand_directories(args.libraries)
    format_output(args.annotations, args.libraries, args.output)
