import argparse
import pandas as pd

from os.path import basename
from spectrum.spectrum import read_mgf
from spectrum.utils import expand_directories
from tqdm import tqdm
from typing import List


def filter_by_library(data: pd.DataFrame, retention_time_tolerance: float = 0.5, order_by: str = 'Score') -> pd.DataFrame:
    common_columns = ['#Scan#', 'Charge', 'SpecMZ', 'SpectrumFile', 'p-value']
    annotation_columns = ['Peptide', 'Mass', 'SharedPeaks', 'Id', 'Score', 'Clean Entropy', 'RetIndex', 'MassDiff', 'RTdiff', 'RTmatch', 'Protein',
                          'mzErrorPPM', 'FalseDiscoveryRate', 'IonMode', 'Instrument', 'Prec.Type', 'InChIKey',
                          'RetIndexDiff']

    filtered_rows = []
    for index, annotation_group in tqdm(data.groupby(by=['#Scan#']), total=len(data['#Scan#'].unique())):
        row = annotation_group.iloc[0][common_columns]
        unique_libraries = annotation_group['LibraryName'].unique()
        for library in unique_libraries:
            library_group = annotation_group[annotation_group['LibraryName'] == library] if isinstance(library, str) else annotation_group
            library_group['in_ret_time_tolerance'] = library_group['RTdiff'].apply(lambda x: x < retention_time_tolerance) \
                if 'RTdiff' in library_group.columns and library.startswith('LEVEL1') else None
            sorted_group = library_group.sort_values(by=['in_ret_time_tolerance', order_by], ascending=False)
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


def merge_columns(data: pd.DataFrame) -> pd.Series:
    return data.apply(lambda x: ', '.join(x.astype(str)) if len(x.dropna()) > 0 else ' ', axis=1)


def format_output(annotations_filename: str, library_filenames: List[str], output_filename: str,
                  retention_time_tolerance: float = 0.5, order_by: str = 'Score'):
    annotations = pd.read_csv(annotations_filename, header=0, sep='\t')
    annotations.rename(mapper={
        'LibrarySpectrumID': 'Id',
        'CompoundName': 'Peptide',
        'ParentMassDiff': 'MassDiff',
        'LibSearchSharedPeaks': 'SharedPeaks',
        'MQScore': 'Score',
        'CleanEntropy': 'Clean Entropy',
        'ExactMass': 'Mass'
    }, axis=1, inplace=True)

    # for library_filename in library_filenames:
    #     add_library_info(annotations, library_filename)

    annotations = filter_by_library(annotations, retention_time_tolerance, order_by)

    annotations['MQScore'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': Score')]])
    annotations['Entropy'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': Clean Entropy')]])
    annotations['SpectrumID'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': Id')]])
    annotations['Compound_Name'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': Peptide')]])
    annotations['Adduct'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': Prec.Type')]])
    annotations['INCHI'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': InChIKey')]])
    annotations['MZErrorPPM'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': mzErrorPPM')]])
    annotations['RTdiff'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': RTdiff')]])
    annotations['RetIndex'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': RetIndex')]])
    annotations['RetIndexDiff'] = merge_columns(annotations[[c for c in annotations.columns if c.endswith(': RetIndexDiff')]])

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


def map_order_by(order_by: str):
    if order_by == 'MQScore':
        return 'Score'
    if order_by == 'CleanEntropy':
        return 'Clean Entropy'
    return order_by


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate annotations in the GNPS format')
    parser.add_argument('--annotations', help='TSV file with MSPepSearch annotations', required=True)
    parser.add_argument('--libraries', help='MSP files with library spectra', nargs='+', required=True)
    parser.add_argument('--output', help='Output filename', required=True)
    parser.add_argument('--retention-time-tolerance', help='Retention time tolerance (min)', type=float,
                        default=0.5)
    parser.add_argument('--order-by', help='Order by', default='MQScore')
    args = parser.parse_args()
    args.libraries = expand_directories(args.libraries)
    format_output(args.annotations, args.libraries, args.output, args.retention_time_tolerance,
                  map_order_by(args.order_by))
