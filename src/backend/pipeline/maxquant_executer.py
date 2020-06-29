import os
import shutil
import subprocess
import CONSTANTS as CONSTS
import logging
logger = logging.getLogger()


def prepare_dbs(wd, raw_dbs_path):
    maxquant_dbs_path = f'{wd}/dbs'
    logger.info(f'Copying {raw_dbs_path} TO {maxquant_dbs_path}')
    try:
        shutil.copytree(raw_dbs_path, maxquant_dbs_path)
    except FileExistsError:
        pass

    for db in os.listdir(maxquant_dbs_path):
        db_path = f'{maxquant_dbs_path}/{db}'
        logger.info(f'Preparing {db_path} for MaxQuant analysis...')
        if not db_path.endswith('fasta'):
            logger.info(f'Removing this file from MaxQuant DBs path (not in a fasta format):\n{db_path}')
            shutil.rmtree(db_path)
            return

        with open(db_path) as f:
            cleaned_headers = f.readline()  # header
            for line in f:
                if not line.startswith('>'):
                    # accumulate sequence rows until a header was detected
                    cleaned_headers += line.rstrip()
                else:
                    # a header was found. Don't forget new line for the previous sequence...
                    cleaned_headers += '\n' + line
                    # cleaned_headers += '\n' + line.split('|')[0].rstrip() + '\n'
            cleaned_headers += '\n'

        # re-writing db
        with open(db_path, 'w') as f:
            f.write(cleaned_headers)

    return maxquant_dbs_path


def get_fasta_files_tag(db_paths):
    tag = '   <fastaFiles>\n'
    for db_path in db_paths:
        tag += f'''      <FastaFileInfo>
         <fastaFilePath>{db_path}</fastaFilePath>
         <identifierParseRule>>(.*)</identifierParseRule>
         <descriptionParseRule>>(.*)</descriptionParseRule>
         <taxonomyParseRule></taxonomyParseRule>
         <variationParseRule></variationParseRule>
         <modificationParseRule></modificationParseRule>
         <taxonomyId></taxonomyId>
      </FastaFileInfo>\n'''
    tag += '   </fastaFiles>\n'
    return tag


def get_raw_files_tag(raw_paths):
    tag = '   <filePaths>\n'
    for raw_path in raw_paths:
        tag += f'''      <string>{raw_path}</string>\n'''
    tag += '   </filePaths>\n'
    return tag


def get_enzymes_tag(enzymes):
    tag = '         <enzymes>\n'
    for enzyme in enzymes:
        tag += f'''            <string>{enzyme}</string>\n'''
    tag += '         </enzymes>\n'
    return tag


def prepare_mqpar_file(db_paths, raw_paths, enzymes, mqpar_path, num_of_threads, mqpar_template_path='/bioseq/pasa/pipeline/mqpar_template.xml'):
    logger.info(f'Preparing mqpar.xml file for MaxQuant run. File will be saved at:\n{mqpar_path}')
    with open(mqpar_template_path) as f:
        mqpar_content = ''
        for line in f:
            # search for lines/tags that need to be manipulated
            if '<fastaFiles></fastaFiles>' in line:
                mqpar_content += get_fasta_files_tag(db_paths)
            elif '<filePaths></filePaths>' in line:
                mqpar_content += get_raw_files_tag(raw_paths)
            elif '<enzymes></enzymes>' in line:
                mqpar_content += get_enzymes_tag(enzymes)
            else:
                # "regular" line, copy as is
                mqpar_content += line
        mqpar_content = mqpar_content.replace('<numThreads></numThreads>', f'<numThreads>{num_of_threads}</numThreads>')

    with open(mqpar_path, 'w') as f:
        f.write(mqpar_content)


def run_maxquant(pasa_dbs_path, enzymes, data_type, job_id,
                 maxquant_current_analysis_path, num_of_threads=20):

    logger.info(f'Running MaxQuant for {data_type}')

    raw_paths = [f'{maxquant_current_analysis_path}/{raw_name}' for raw_name in sorted(os.listdir(maxquant_current_analysis_path))]
    logger.info(f'Raw files detected are:\n' + '\n'.join(raw_paths))
    assert len(raw_paths)

    maxquant_dbs_path = prepare_dbs(maxquant_current_analysis_path, pasa_dbs_path)
    mqpar_path = f'{maxquant_current_analysis_path}/mqpar.xml'
    db_paths = [f'{maxquant_dbs_path}/{db_name}' for db_name in os.listdir(maxquant_dbs_path)]
    logger.info(f'DB files detected are:\n' + '\n'.join(db_paths))
    assert len(db_paths)

    prepare_mqpar_file(db_paths, raw_paths, enzymes, mqpar_path, num_of_threads)

    cmds_file = f'{maxquant_current_analysis_path}/maxquant.cmds'
    with open(cmds_file, 'w') as f:
        f.write(f'module load mono/mono-5.18.0.225 python/python-3.6.7')
        f.write('!@#')
        f.write(f'mono /bioseq/pasa/MaxQuant.1.6.14/bin/MaxQuantCmd.exe {mqpar_path}')
        f.write('\t')
        f.write(f'mq_{data_type}_{job_id}')

    cmd = f'{CONSTS.Q_SUBMITTER_SCRIPT} {cmds_file} {maxquant_current_analysis_path} -q pupkotmpr --cpu {num_of_threads}'
    logger.info(f'Starting MaxQuant. Calling:\n{cmd}')
    subprocess.run(cmd, shell=True)

    logger.info(f'MaxQuant analysis for {data_type} data was fetched successfully.')



