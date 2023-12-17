import argparse
import pandas as pd

from os import listdir
from os.path import basename, isdir, isfile, join, splitext
from spectrum.spectrum import read_mgf, read_msp
from tqdm import tqdm
from typing import List


def add_retention_time_differences(annotation_file: str, spectrum_file: str, library_files: List[str], output_file: str,
                                   tolerance: float = 0.5, retention_time_matches_only: bool = False):
    property_mapping = {
        'IONMODE': 'IonMode',
        'ION_MODE': 'IonMode',
        'SOURCE_INSTRUMENT': 'Instrument',
        'INSTRUMENT': 'Instrument',
        'PRECURSOR_TYPE': 'Prec.Type',
        'INCHIKEY': 'InChIKey',
        'EXACTMASS': 'Mass',
        'EXACT_MASS': 'Mass'
    }

    try:
        annotations = pd.read_csv(annotation_file, header=0, sep='\t')
    except pd.errors.EmptyDataError:
        return
    spectra = read_mgf(spectrum_file)
    # spectra_dict = {s.properties['name']: s for s in spectra}

    libraries = {}
    for filename in library_files:
        library_name = basename(filename)
        libraries[library_name] = read_mgf(filename, id_field='SPECTRUMID')

    annotations_with_rt = annotations.copy()
    drop_indices = []
    for index, row in tqdm(annotations.iterrows(), total=len(annotations)):
        query_id = str(row['#Scan#'])
        # query_spectrum = spectra.get(query_id)
        query_spectrum = spectra[query_id] if query_id in spectra else None
        if query_spectrum is None:
            continue

        # query_ret_time = query_spectrum.properties.get('RTINSECONDS')
        # if query_ret_time is not None:
        #     query_ret_time = float(query_ret_time)
        query_ret_time = float(
            query_spectrum.properties['RTINSECONDS']) if 'RTINSECONDS' in query_spectrum.properties else None

        library_name = row['LibraryName']
        # library = libraries.get(library_name)
        library = libraries[library_name] if library_name in libraries else None
        if library is None:
            continue

        library_id = str(row['LibrarySpectrumID'])
        # library_spectrum = library.get(library_id)
        library_spectrum = library[library_id] if library_id in library else None
        if library_spectrum is None:
            continue

        # if query_spectrum.get_polarity_sign() != library_spectrum.get_polarity_sign():
        #     drop_indices.append(index)
        #     continue

        # library_ret_time = library_spectrum.properties.get('RTINSECONDS')
        # if library_ret_time is not None:
        #     library_ret_time = float(library_ret_time)
        library_ret_time = float(
            library_spectrum.properties['RTINSECONDS']) if 'RTINSECONDS' in library_spectrum.properties else None

        if query_ret_time is not None and library_ret_time is not None:
            ret_time_difference = abs(library_ret_time - query_ret_time) / 60.0
            annotations_with_rt.loc[index, 'RTdiff'] = ret_time_difference
            annotations_with_rt.loc[index, 'RTmatch'] = ret_time_difference < tolerance
            if retention_time_matches_only and ret_time_difference > tolerance:
                drop_indices.append(index)
        elif retention_time_matches_only:
            drop_indices.append(index)

        for key in library_spectrum.properties:
            if key in property_mapping:
                annotations_with_rt.loc[index, property_mapping[key]] = library_spectrum.properties[key]

    if len(drop_indices) > 0:
        annotations_with_rt.drop(labels=drop_indices, axis='index', inplace=True)

    annotations_with_rt.to_csv(output_file, index=False, sep='\t')


def expand_directory(folder_name: str) -> List[str]:
    filenames = [f for f in listdir(folder_name) if f.lower().endswith('.msp') or f.lower().endswith('.mgf')]
    paths = [join(folder_name, f) for f in filenames]
    return [p for p in paths if isfile(p)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate HTML with spectral plots')
    parser.add_argument('--annotations', help='TSV file with MSPepSearch annotations', required=True)
    parser.add_argument('--spectra', help='MSP file with spectra', required=True)
    parser.add_argument('--libraries', help='MSP files with library spectra', nargs='+', required=True)
    parser.add_argument('--output', help='Output filename')
    parser.add_argument('--tolerance', help='Retention time tolerance', type=float, default=0.5)
    parser.add_argument('--retention-time-matches-only',
                        help='If true, only matches with the retention time error below the threshold are kept',
                        type=lambda x: x.lower() in ('1', 'yes', 'true'),
                        default=False)
    args = parser.parse_args()

    if isdir(args.spectra):
        files = expand_directory(args.spectra)
        assert len(files) > 0, f'Cannot find spectrum files in folder {args.spectra}'
        args.spectra = files[0]

    libraries = []
    for lib in args.libraries:
        if isdir(lib):
            libraries += expand_directory(lib)
        else:
            libraries.append(lib)
    args.libraries = libraries

    add_retention_time_differences(args.annotations, args.spectra, args.libraries, args.output, args.tolerance,
                                   args.retention_time_matches_only)
