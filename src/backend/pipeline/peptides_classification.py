import os
import logging
import re
import CONSTANTS as CONSTS

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

def get_file_name(path):
    return os.path.splitext(path)[0].split('/')[-1]


def get_next_relevant_header(f, header, i):
    logger.info(f'Human db record was detected:\n{header}')
    logger.info(f'Skipping records until ASAP records will be found...')
    while header:
        # skip rows until next non-human header is reached
        if header.startswith('>') and not header.startswith('>sp') and not header.startswith('>tr'):
            logger.info(f'Finished skipping human db records. Next header is:\n{header}')
            break
        header = f.readline()
        i += 1
    return header, i


def get_record_relative_intensity(db_record, total_intensity=1):
    avg_intensity = sum(db_record.reps_intensities) / len(db_record.reps_intensities)
    return avg_intensity/total_intensity


def load_db(db_path):
    from DatabaseRecord import Record

    # a dictionary mapping each sequence to a list of its metadata from the db (ASAP format)
    sequence2db_record = {}
    total_intensity_in_db = 0

    db_name = get_file_name(db_path)
    with open(db_path) as f:
        header = f.readline()
        i = 0
        while header:
            header = header.rstrip()
            if not header.startswith('>'):
                raise ValueError(f'''{get_file_name(db_path)} is illegal. Line {i+1} contains a badly formatted header: <br>
                                 {header.rstrip()}
                                 <br>(FASTA header should start with ">" character, not with "{header[0]}").''')

            if header.startswith('>sp') or header.startswith('>tr'):
                # (default) human db record (header + sequence)
                header, i = get_next_relevant_header(f, header, i)
                if not header:
                    assert sequence2db_record, f'No records (other than human db records) were detected in {get_file_name(db_path)}'
                    break

            header_tokens = header.lstrip('>').rstrip().split('|')
            sequence = f.readline().rstrip()
            header = f.readline()
            while header and not header.startswith('>'):
                sequence += header.rstrip()
                header = f.readline()
                i += 1

            # this block must be only after next header was read!
            if header_tokens[0] != 'IGH':
                logger.debug(f"Skipping record with non IGH chain\n{'>'+'|'.join(header_tokens)}")
                continue

            if sequence in sequence2db_record:
                # TODO: investigate in ASAP why it happens...? maybe MiXCR?
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

            try:
                sequence2db_record[sequence] = Record(db_name, sequence, *header_tokens)
                total_intensity_in_db += get_record_relative_intensity(sequence2db_record[sequence])
                pass
            except:
                raise ValueError(f'''{get_file_name(db_path)} is illegal. Failed to parse line {i+1} since it contains a badly formatted header: <br> 
                                    {'>'+'|'.join(header_tokens)} 
                                    <br>For further details see "?" button next to the corresponding input field.''')

    return sequence2db_record, total_intensity_in_db


def create_ambiguous_pattern(protein):
    return re.sub('I|L', '[IL]', protein)


def map_peptide_to_dbs(peptide, peptide_intensity_in_elution, db_name2db, db_name2total_intensity,
                       db_name2db_intensity_file_handler, multiple_clones_peptides_file_handler,
                       no_clones_peptides_file_handler, minimum_overlap=1, normalization_factor=1):
    """
    :param peptide: a string to blast against all dbs
    :param db_name2db: a dictionary mapping db_name to the db it self (db is a dictionary that maps sequence to db_record)
    :param minimum_overlap: an integer representing the minimum overlap needed for a peptide to be considered as cdr3 informative
    :return: 0, None, None -> non_informative;
             1, relevant_records, db_name2record_intensity -> informative;
             2, relevant_records, db_name2record_intensity -> cdr3_informative;
    """

    db_name2record_intensity = {}
    db_name2relevant_records = {}
    cdr3s = set()

    for db_name in db_name2db:
        sequence2db_record = db_name2db[db_name]
        for sequence in sequence2db_record:
            db_record = sequence2db_record[sequence]
            pattern = create_ambiguous_pattern(peptide)
            match = re.search(pattern, db_record.seq)
            if match:
                db_name2relevant_records[db_name] = db_name2relevant_records.get(db_name, []) + [db_record]
                record_relative_intensity = get_record_relative_intensity(db_record,
                                                                          db_name2total_intensity[db_name])
                db_name2record_intensity[db_name] = db_name2record_intensity.get(db_name, 0) + record_relative_intensity
                match_start, match_stop = match.span()
                cdr3_start, cdr3_stop = re.search(db_record.cdr3, db_record.seq).span()
                cdr3s.add(db_record.cdr3)
                if len(cdr3s) > 1:
                    if multiple_clones_peptides_file_handler:
                        # for debugging purposes
                        multiple_clones_peptides_file_handler.write(f'{peptide},{",".join(cdr3s)}\n')
                    return 0, None, None  # non-informative

    if not db_name2relevant_records:
        # no matching record was found in all dbs!
        if no_clones_peptides_file_handler:
            no_clones_peptides_file_handler.write(f'{peptide}\n')
        return 0, None, None  # non-informative

    if cdr3_stop < match_start + minimum_overlap or cdr3_start + minimum_overlap > match_stop:
        return 1, db_name2relevant_records, db_name2record_intensity  # peptide does not overlap with cdr3 with at least $minimum_overlap positions

    # peptide overlaps with cdr3
    for db_name in db_name2db_intensity_file_handler:
        if db_name not in db_name2relevant_records:
            continue
        relevant_records = db_name2relevant_records[db_name]
        total_db_intensity = db_name2total_intensity[db_name]
        total_peptide_intensity_in_db = sum(get_record_relative_intensity(record, total_db_intensity) for record in relevant_records)
        f = db_name2db_intensity_file_handler[db_name]
        f.write(f'{peptide},{normalization_factor * peptide_intensity_in_elution},'
                f'{normalization_factor * total_peptide_intensity_in_db}\n')

    return 2, db_name2relevant_records, db_name2record_intensity


def classify_peptides(db_paths, filtered_peptides_info, raw_output_dir):

    db_name2db = {}
    db_name2total_intensity = {}
    db_name2db_intensity_file_handler = {}
    for db_path in db_paths:
        db_name = get_file_name(db_path)
        logger.info(f'_'*80)
        logger.info(f'Loading {db_name}...')
        db_name2db[db_name], db_name2total_intensity[db_name] = load_db(db_path)
        db_name2db_intensity_file_handler[db_name] = open(f'{raw_output_dir}/proteomic_intensities_with_respect_to_{db_name}.csv', 'w')
        db_name2db_intensity_file_handler[db_name].write('peptide,peptide_intensity_in_elution,total_peptide_intensity_in_db\n')

    non_informative_path = f'{raw_output_dir}/{CONSTS.NON_INFORMATIVE_CSV_NAME}'
    non_informative_f = open(non_informative_path, 'w')
    non_informative_f.write(f'Peptide\n')

    informative_path = f'{raw_output_dir}/{CONSTS.INFORMATIVE_CSV_NAME}'
    informative_f = open(informative_path, 'w')
    informative_f.write(f'peptide,mappedSequence,cdr3,IGHV,IGHD,IGHJ,isotype,chainType,counts\n')

    cdr3_informative_path = f'{raw_output_dir}/{CONSTS.CDR3_INFORMATIVE_CSV_NAME}'
    cdr3_informative_f = open(cdr3_informative_path, 'w')
    cdr3_informative_f.write(f'peptide,mappedSequence,cdr3,IGHV,IGHD,IGHJ,isotype,chainType,counts\n')

    if not os.path.exists('/bioseq'):
        # debugging in local run
        multiple_clones_peptides_file_handler = open(f'{raw_output_dir}/multiple_clones.csv', 'w')
        no_clones_peptides_file_handler = open(f'{raw_output_dir}/no_clone.txt', 'w')
    else:
        multiple_clones_peptides_file_handler = no_clones_peptides_file_handler = None

    for i, peptide_info in enumerate(filtered_peptides_info):

        if i % 10 == 0:
            logger.info(f'Peptide {i} is being queried')

        peptide = peptide_info[0]
        peptide_intensity_in_elution = peptide_info[3]
        informativity_level, db_name2relevant_records, db_name2record_intensity = map_peptide_to_dbs(peptide, peptide_intensity_in_elution, db_name2db, db_name2total_intensity, db_name2db_intensity_file_handler, multiple_clones_peptides_file_handler, no_clones_peptides_file_handler)
        if informativity_level == 0:
            non_informative_f.write(f'{peptide}\n')
            continue

        for records in db_name2relevant_records.values():
            for record in records:
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
    if not os.path.exists('/bioseq'):
        # debugging in local run
        multiple_clones_peptides_file_handler.close()
        no_clones_peptides_file_handler.close()


# from peptides_filtration import filter_peptides
# filtered_peptides_info = filter_peptides('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/el_peptides.txt',
#                                 '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/ft_peptides.txt',
#                                 '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/', min_fold=5)
# classify_peptides(["/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/11.fatsa",
#                    "/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/6.fatsa",
#                    "/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/7.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/8.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/9.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/10.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/1.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/2.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/3.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/4.fatsa","/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs/5.fatsa"],
#                   filtered_peptides_info,
#                   '/Users/Oren/Dropbox/Projects/PASA/linux_outputs')
# x=load_db('/Users/Oren/Dropbox/Projects/PASA/data/asap_format_dbs/159214090862050017940890971071.fasta')

