import argparse
import pandas as pd

from os.path import basename
from spectrum.spectrum import read_mgf
from spectrum.utils import expand_directories
from tqdm import tqdm
from typing import List


def filter_by_library(data: pd.DataFrame) -> pd.DataFrame:
    common_columns = ['#Scan#', 'Charge', 'SpecMZ', 'SpectrumFile', 'p-value']
    annotation_columns = ['Peptide', 'Mass', 'SharedPeaks', 'Id', 'Score', 'MassDiff', 'RTdiff', 'RTmatch', 'Protein',
                          'mzErrorPPM', 'FalseDiscoveryRate', 'IonMode', 'Instrument', 'Prec.Type', 'InChIKey']

    filtered_rows = []
    for index, annotation_group in tqdm(data.groupby(by=['#Scan#']), total=len(data['#Scan#'].unique())):
        row = annotation_group.iloc[0][common_columns]
        unique_libraries = annotation_group['LibraryName'].unique()
        for library in unique_libraries:
            library_group = annotation_group[annotation_group['LibraryName'] == library] if isinstance(library, str) else annotation_group
            sorted_group = library_group.sort_values(by='Score', ascending=False)
            for column in annotation_columns:
                # library_name = rename_nist_library(library)
                row[f'{library}: {column}'] = sorted_group.iloc[0].get(column)
        filtered_rows.append(row)

    # filtered_data = pd.concat(filtered_rows, ignore_index=True)
    filtered_data = pd.DataFrame(data=filtered_rows)
    return filtered_data


# def add_library_info(annotations: pd.DataFrame, library_filename: str):
#     library = read_mgf(library_filename, id_field='SPECTRUMID')
#     library_name = basename(library_filename)
#     for index, annotation in annotations[annotations['LibraryName'] == library_name].iterrows():
#         spectrum_id = annotation['Id']
#         if not spectrum_id or pd.isna(spectrum_id):
#             raise ValueError('Cannot find library spectrum ID')
#
#         spectrum = library.get(spectrum_id)
#         if not spectrum:
#             continue
#
#         properties = spectrum.properties
#         annotations.loc[index, 'IonMode'] = properties.get('IONMODE', properties.get('ION_MODE'))
#         annotations.loc[index, 'Instrument'] = properties.get('SOURCE_INSTRUMENT', properties.get('INSTRUMENT'))
#         annotations.loc[index, 'Prec.Type'] = properties.get('PRECURSOR_TYPE')
#         annotations.loc[index, 'InChIKey'] = properties.get('INCHIKEY')
#         annotations.loc[index, 'Mass'] = properties.get('EXACTMASS')


def format_output(annotations_filename: str, library_filenames: List[str], output_filename: str):
    annotations = pd.read_csv(annotations_filename, header=0, sep='\t')
    annotations.rename(mapper={
        'LibrarySpectrumID': 'Id',
        'CompoundName': 'Peptide',
        'ParentMassDiff': 'MassDiff',
        'LibSearchSharedPeaks': 'SharedPeaks',
        'MQScore': 'Score',
        'ExactMass': 'Mass'
    }, axis=1, inplace=True)

    # for library_filename in library_filenames:
    #     add_library_info(annotations, library_filename)

    annotations = filter_by_library(annotations)

    annotations['MQScore'] = annotations[[c for c in annotations.columns if c.endswith(': Score')]].apply(
        lambda x: ', '.join(x.astype(str)), axis=1)
    annotations['SpectrumID'] = annotations[[c for c in annotations.columns if c.endswith(': Id')]].apply(
        lambda x: ', '.join(x.astype(str)), axis=1)
    annotations['Compound_Name'] = annotations[[c for c in annotations.columns if c.endswith(': Peptide')]].apply(
        lambda x: ', '.join(x.astype(str)), axis=1)
    annotations['Adduct'] = annotations[[c for c in annotations.columns if c.endswith(': Prec.Type')]].apply(
        lambda x: ', '.join(x.astype(str)), axis=1)
    annotations['INCHI'] = annotations[[c for c in annotations.columns if c.endswith(': InChIKey')]].apply(
        lambda x: ', '.join(x.astype(str)), axis=1)
    annotations['MZErrorPPM'] = annotations[[c for c in annotations.columns if c.endswith(': mzErrorPPM')]].apply(
        lambda x: ', '.join(x.astype(str)), axis=1)
    annotations['RTdiff'] = annotations[[c for c in annotations.columns if c.endswith(': RTdiff')]].apply(
        lambda x: ', '.join(x.astype(str)), axis=1)

    annotations['Smiles'] = None
    annotations['Ion_Source'] = None
    annotations['Instrument'] = None
    annotations['Compound_Source'] = None
    annotations['PI'] = None
    annotations['Data_Collector'] = None
    annotations['IonMode'] = None
    annotations['MassDiff'] = None
    annotations['SharedPeaks'] = None
    annotations['tags'] = None
    annotations['Library_Class'] = None

    annotations.to_csv(output_filename, index=False, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate annotations in the GNPS format')
    parser.add_argument('--annotations', help='TSV file with MSPepSearch annotations', required=True)
    parser.add_argument('--libraries', help='MSP files with library spectra', nargs='+', required=True)
    parser.add_argument('--output', help='Output filename', required=True)
    args = parser.parse_args()
    args.libraries = expand_directories(args.libraries)
    format_output(args.annotations, args.libraries, args.output)
