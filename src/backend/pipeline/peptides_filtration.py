import logging
logger = logging.getLogger()

# def fail(error_msg, error_file_path='error.txt'):
#     # TODO: implement
#     with open(error_file_path, 'w') as error_f:
#         error_f.write(error_msg + '\n')
#     raise ValueError(error_msg)


def get_intensities_columns(peptides_path):
    total_intensity = reverese = contaminant = None
    with open(peptides_path) as f:
        splitted_header = f.readline().split('\t')

    for i, field in enumerate(splitted_header):
        if 'Intensity' in field and not total_intensity:
            total_intensity = i
        elif field == 'Reverse':
            reverese = i
        elif field == 'Potential contaminant':
            contaminant = i

    # detect replication intensities columns
    i = total_intensity + 1
    replications_intensities = []
    while 'Intensity' in splitted_header[i]:
        replications_intensities.append(i)
        i += 1

    num_of_replications = len(replications_intensities)
    if num_of_replications < 2:
        raise ValueError(f'Number of replications detected in {peptides_path} is {num_of_replications}. At least 2 are required.')

    if not total_intensity:
        raise ValueError('Intensity column could be find')

    return replications_intensities, reverese, contaminant


def get_peptide_to_avg_intensity(mq_output_path):
    replications_intensities, reverese_index, contaminant_index = get_intensities_columns(mq_output_path)

    peptides_to_avg_intensity = {}

    with open(mq_output_path) as f:
        f.readline()  # skip header
        for line in f:
            splitted_line = line.rstrip().split('\t')

            if splitted_line[reverese_index] == '+':
                # decoy
                continue

            if splitted_line[contaminant_index] == '+':
                # possible contaminant
                continue

            peptide = splitted_line[0]

            # sanity check
            if peptide in peptides_to_avg_intensity:
                raise ValueError(f'{peptide} exists more than once in {mq_output_path}')

            intensities = [float(splitted_line[i]) for i in replications_intensities]

            peptide_existence = [intensity > 0 for intensity in intensities]
            if sum(peptide_existence) < 2:
                # peptide was found in less than 2 reps
                continue

            peptides_to_avg_intensity[peptide] = sum(intensities)/len(peptide_existence)

    return peptides_to_avg_intensity


def get_peptide_to_relative_avg_intensity(peptides_to_intensity):
    # TODO: should the total intensity be divided by the sum of filtered intensities or by the sum of all intensities?
    #       if the latter, the sum will not be 1 because we discard some records that contribute to the sum
    sum_of_intensities = sum(peptides_to_intensity.values())
    return {peptide: peptides_to_intensity[peptide]/sum_of_intensities for peptide in peptides_to_intensity}


def filter_peptides(el_path, ft_path, output_dir, min_fold=5):

    logger.info('Parsing Elution peptides file...')
    el_peptides_to_intensity = get_peptide_to_avg_intensity(el_path)
    el_peptides_to_relative_intensity = get_peptide_to_relative_avg_intensity(el_peptides_to_intensity)

    logger.info('Parsing FT peptides file...')
    ft_peptides_to_intensity = get_peptide_to_avg_intensity(ft_path)
    ft_peptides_to_relative_intensity = get_peptide_to_relative_avg_intensity(ft_peptides_to_intensity)

    filtered_peptides_info = []
    junk_peptides_info = []
    for peptide in el_peptides_to_intensity:
        el_relative_intensity = el_peptides_to_intensity[peptide]
        ft_relative_intensity = ft_peptides_to_intensity.get(peptide, 0)

        # notice that ft_relative_intensity might be zero so be careful with divisions here...
        if el_relative_intensity > min_fold * ft_relative_intensity:
            filtered_peptides_info.append([peptide,
                                           el_peptides_to_intensity[peptide],
                                           ft_peptides_to_intensity.get(peptide, 0),
                                           el_peptides_to_relative_intensity[peptide],
                                           ft_peptides_to_relative_intensity.get(peptide, 0)])
        else:
            junk_peptides_info.append([peptide,
                                       str(el_peptides_to_intensity[peptide]),
                                       str(ft_peptides_to_intensity.get(peptide, 0)),
                                       str(el_peptides_to_relative_intensity[peptide]),
                                       str(ft_peptides_to_relative_intensity.get(peptide, 0))])

    filtered_peptides_path = f'{output_dir}/filtered_peptides.csv'
    with open(filtered_peptides_path, 'w') as f:
        f.write('Peptide_name,Average_intensity_elution,Average_intensity_flowthrough,Average frequency_elution,Average frequency_flowthrough\n')
        for element in filtered_peptides_info:
            f.write(','.join([str(x) for x in element]) + '\n')

    junk_peptides_path = f'{output_dir}/non_enriched_elution_peptides.csv'
    with open(junk_peptides_path, 'w') as f:
        f.write('Peptide_name,Average_intensity_elution,Average_intensity_flowthrough,Average frequency_elution,Average frequency_flowthrough\n')
        for element in junk_peptides_info:
            f.write(','.join([str(x) for x in element]) + '\n')

    return filtered_peptides_info


# filter_peptides('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/el_peptides.txt',
#        '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/ft_peptides.txt',
#        '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/filtered_peptides.txt')