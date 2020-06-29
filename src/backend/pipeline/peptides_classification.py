import os
import logging
import re

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

def get_file_name(path):
    return os.path.splitext(path)[0].split('/')[-1]


def load_db(db_path):
    from DatabaseRecord import Record

    # a dictionary mapping each sequence to a list of its metadata from the db (ASAP format)
    sequence2db_record = {}

    db_name = get_file_name(db_path)
    with open(db_path) as f:
        i = -1
        for header in f:
            i += 1
            if not header:
                continue

            if not header.startswith('>'):
                raise ValueError(f'Line number {i+1} is illegal (header lines should start with ">" character, '
                                 f'not with "{header[0]}").')


            if header.startswith('>sp') or header.startswith('>tr'):
                # (default) human db record (header + sequence)
                logger.info(f'Human db record was detected:\n{header}')
                logger.info(f'Skipping records until ASAP records will be found...')

                empty_lines = 0
                while True:
                    # skip rows until next header is reached
                    header = f.readline()
                    i += 1
                    if header.startswith('>') and not header.startswith('>sp') and not header.startswith('>tr'):
                        logger.info(f'Finished skipping human db records. Next header is:\n{header}')
                        break
                    if not header.strip():
                        empty_lines += 1
                    if empty_lines > 2:
                        # to avoid infinite loop..
                        break

            header_tokens = header.lstrip('>').rstrip().split('|')
            sequence = f.readline().rstrip()
            i += 1

            if sequence in sequence2db_record:
                # TODO: investigate why it happens...?
                previous_header = sequence2db_record[sequence]
                previous_header_tokens = [previous_header.chain_type, previous_header.iso_type,
                                          previous_header.cdr3, previous_header.IGHV,
                                          previous_header.IGHD, previous_header.IGHJ,
                                          previous_header.counts, ','.join(str(x) for x in previous_header.reps_intensities)]
                logger.warning(f'The following sequence appears more than once in the uploaded database ({db_name}):\n{sequence}')
                logger.warning(f'Previous header: {previous_header_tokens}')
                logger.warning(f'Current header:  {header_tokens}')  # there are 2 spaces in purpose (aligned better).
                logger.warning('Skipping current record...\n')

                continue

            sequence2db_record[sequence] = Record(db_name, sequence, *header_tokens)

    return sequence2db_record


def create_ambiguous_pattern(protein):
    return re.sub('I|L', '[IL]', protein)


def map_peptide_to_dbs(peptide, db_name2db, minimum_overlap=1):
    """
    :param peptide: a string to blast against all dbs
    :param db_name2db: a dictionary mapping db_name to the db it self (db is a dictionary that maps sequence to db_record)
    :param minimum_overlap: an integer representing the minimum overlap needed for a peptide to be considered as cdr3 informative
    :return: 0, None -> non_informative;
             1, relevant_record -> informative;
             2, relevant_record -> cdr3_informative;
    """

    relevant_record = ''
    cdr3s = set()

    for db_name in db_name2db:
        sequence2db_record = db_name2db[db_name]
        for sequence in sequence2db_record:
            db_record = sequence2db_record[sequence]
            pattern = create_ambiguous_pattern(peptide)
            match = re.search(pattern, db_record.seq)
            if match:
                relevant_record = db_record
                match_start, match_stop = match.span()
                cdr3_start, cdr3_stop = re.search(relevant_record.cdr3, relevant_record.seq).span()
                cdr3s.add(db_record.cdr3)
                if len(cdr3s) > 1:
                    return 0, None  # non-informative

    if not relevant_record:
        # no matching record was found in all dbs!
        return 0, None  # non-informative

    if cdr3_stop < match_start + minimum_overlap or cdr3_start + minimum_overlap > match_stop:
        return 1, relevant_record  # peptide does not overlap with cdr3 with at least $minimum_overlap positions

    # peptide overlaps with cdr3
    return 2, relevant_record


def classify_peptides(db_paths, filtered_peptides_info, non_informative_path, informative_path, cdr3_informative_path):

    db_name2db = {}
    for db_path in db_paths:
        db_name = get_file_name(db_path)
        logger.info(f'Loading {db_name}...')
        db_name2db[db_name] = load_db(db_path)

    non_informative_f = open(non_informative_path, 'w')
    non_informative_f.write(f'Peptide\n')

    informative_f = open(informative_path, 'w')
    informative_f.write(f'peptide,mappedSequence,cdr3,IGHV,IGHD,IGHJ,isotype,chainType,counts\n')

    cdr3_informative_f = open(cdr3_informative_path, 'w')
    cdr3_informative_f.write(f'peptide,mappedSequence,cdr3,IGHV,IGHD,IGHJ,isotype,chainType,counts\n')

    for i, peptide_info in enumerate(filtered_peptides_info):

        if i % 50 == 0:
            logger.info(f'Peptide {i} is being queried')

        peptide = peptide_info[0]
        informativity_level, record = map_peptide_to_dbs(peptide, db_name2db)
        if informativity_level == 0:
            non_informative_f.write(f'{peptide}\n')
            continue

        result = f'{peptide},{record.seq},{record.cdr3},{record.IGHV},{record.IGHD},' \
                 f'{record.IGHJ},{record.iso_type},{record.chain_type},{record.counts}\n'
        if informativity_level == 1:
            informative_f.write(result)
        else:
            # informativity_level == 2
            cdr3_informative_f.write(result)

    non_informative_f.close()
    informative_f.close()
    cdr3_informative_f.close()




# from peptides_filtration import filter_peptides
# filtered_peptides_info = filter_peptides('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/el_peptides.txt',
#                                 '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/ft_peptides.txt',
#                                 '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/filtered_peptides.txt', min_fold=5)
# classify_peptides(['/Users/Oren/Dropbox/Projects/PASA/exampleDB.txt'],
#                   filtered_peptides_info,
#                   '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/non_informative.csv',
#                   '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/informative.csv',
#                   '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/cdr3_informative.csv')


