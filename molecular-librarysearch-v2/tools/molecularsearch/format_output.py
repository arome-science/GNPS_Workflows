import argparse
import pandas as pd

from spectrum.spectrum import read_mgf
from typing import List


def add_library_info(annotations: pd.DataFrame, library_filename: str):
    library = read_mgf(library_filename)
    for index, annotation in annotations[annotations['LibraryName'] == library_filename].iterrows():
        spectrum_id = annotation['SpectrumID']
        if not spectrum_id:
            continue

        spectrum = library.get(spectrum_id)
        if not spectrum:
            continue

        properties = spectrum.properties
        annotation.loc[index, 'IonMode'] = properties.get('IONMODE', properties.get('ION_MODE'))
        annotation.loc[index, 'Instrument'] = properties.get('SOURCE_INSTRUMENT', properties.get('INSTRUMENT'))
        annotations.loc[index, 'Adduct'] = properties.get('PRECURSOR_TYPE')
        annotation.loc[index, 'INCHI'] = properties.get('INCHIKEY')
        annotation.loc[index, 'ExactMass'] = properties.get('EXACTMASS')


def format_output(annotations_filename: str, spectra_filename: str, library_filenames: List[str], output_filename: str):
    annotations = pd.read_csv(annotations_filename, header=0)
    annotations.rename(mapper={
        'LibrarySpectrumID': 'SpectrumID',
        'Compound_Name': 'CompoundName',
        'ParentMassDiff': 'MassDiff',
        'LibSearchSharedPeaks': 'SharedPeaks'
    }, axis=1, inplace=True)

    for library_filename in library_filenames:
        add_library_info(annotations, library_filename)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate annotations in the GNPS format')
    parser.add_argument('--annotations', help='TSV file with MSPepSearch annotations', required=True)
    parser.add_argument('--spectra', help='MSP file with spectra', required=True)
    parser.add_argument('--libraries', help='MSP files with library spectra', nargs='+', required=True)
    parser.add_argument('--output', help='Output filename')
    args = parser.parse_args()
    format_output(args.annotations, args.spectra, args.libraries, args.output)
