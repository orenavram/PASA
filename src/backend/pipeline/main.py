import sys
import logging
import os
import shutil
from time import time
# initialize logger
while logging.root.handlers:
    logging.root.removeHandler(logging.root.handlers.pop())
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%d/%m/%Y %H:%M:%S')
logger = logging.getLogger()

logger.info(f'Starting {sys.argv[0]}. Executed command is:\n{" ".join(sys.argv)}')
sys.path.append('/bioseq/bioSequence_scripts_and_constants/')
sys.path.append('/bioseq/pasa/pipeline/')

import CONSTANTS as CONSTS
# from email_sender import send_email
from html_editor import edit_failure_html, edit_success_html
from maxquant_executer import run_maxquant
from auxiliaries import wait_for_maxquant, measure_time, edit_progress
from peptides_classification import classify_peptides
from peptides_filtration import filter_peptides
from results_ploting import plot_results


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser()
    # parser.add_argument('DBs_dir', type=lambda x: x.rstrip('/'),
    #                     help='A path to a folder in which the DB files are stored')
    parser.add_argument('working_dir', type=lambda x: x.rstrip('/'),
                        help='An output directory in which ALL results (i.e., intermediate and final) will be written to')
    parser.add_argument('--min-fold', type=float, default=5,
                        help='A minimal fold threshold for passing filtration. If a peptide appears n times in the FT, it should appear at least n x $min_fold times in the elution in order to be kept for downstream analysis')

    # MaxQuant running parameters. Either none or all have to be set.
    parser.add_argument('-mq', '--maxquant', action='store_true')
    parser.add_argument('--enzymes', default='Trypsin', type=lambda x: x.split(','),
                        help='Restriction enzyme names (separated by comma if more than one. no spaces!) for MaxQuant Analysis. '
                             'Possible options are any of the following: '
                             'ArgC,AspC,AspN,Chymotrypsin,Chymotrypsin+,D.P,GluC,GluN,LysC,LysC/P,LysN,Trypsin,Trypsin/P')

    parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')
    return parser.parse_args()


def notify_user(status, email, run_number, start, end):
    logger.info(f'Sending a notification email to {email}')
    results_location = f'{CONSTS.WEBSERVER_RESULTS_URL}/{run_number}/{CONSTS.RESULT_WEBPAGE_NAME}'
    msg = f'Hello,\n\nWe just wanted to let you know that your {CONSTS.WEBSERVER_NAME.upper()} analysis {status}!'
    if status == 'is ready':
        msg += f' (Took {measure_time(int(end - start))})\n\nResults can be found at: {results_location}.\n\n' \
            f'Please note that the results will be kept on the server for the next three months.\n\n' \
            f'Kind regards,\nPASA team'
    else:
        msg += f'. For further information please visit: {results_location}'
    logger.info(msg)
    try:
        send_email('mxout.tau.ac.il', 'TAU BioSequence <bioSequence@tauex.tau.ac.il>', email,
                   subject=f'{CONSTS.WEBSERVER_NAME.upper()} run number {run_number} {status}.', content=msg)
    except:
        logger.info(f'\nFailed sending notification to {email}\n')


def main(wd, min_fold, maxquant_analysis_is_needed, enzymes):

    local_run = not os.path.exists('/boiseq')
    start_time = time()
    maxquant_analysis_dir = f'{wd}/{CONSTS.MAXQUANT_DIR_NAME}'
    el_raw_path = f'{maxquant_analysis_dir}/el'
    ft_raw_path = f'{maxquant_analysis_dir}/ft'
    dbs_dir = f'{wd}/dbs'

    html_path = f'{wd}/{CONSTS.RESULT_WEBPAGE_NAME}'
    edit_progress(html_path, progress=4)

    output_dir = os.path.join(wd, CONSTS.OUTPUT_DIR_NAME)
    os.makedirs(output_dir, exist_ok=True)

    run_number = output_dir.split('/')[-2]  # e.g., /bioseq/data/results/pasa/1234/outputs  ->  1234

    textual_output_path = f'{output_dir}/text'
    os.makedirs(textual_output_path, exist_ok=True)
    elution_peptides_path = f'{wd}/el_peptides.txt'  # A path to a MaxQuant output of ELUTION
    flowthrough_peptides_path = f'{wd}/ft_peptides.txt'  # A path to a MaxQuant output of FLOWTHROUGH

    unpack_zipped_data(f'{wd}/db.zip', dbs_dir)
    logger.info(f'After unpacking, DBs dir contains: {os.listdir(dbs_dir)}')

    edit_progress(html_path, progress=7)

    logger.info('Searching for MaxQuant analysis results...')
    if maxquant_analysis_is_needed and \
            not os.path.exists(elution_peptides_path) or \
            not os.path.exists(flowthrough_peptides_path):
        '''A path (to a folder) in which MaxQuant analysis results will be stored.
        In addition, this folder should contain two folders "el" and "ft".
        In each of which, there should be raw files for the elution experiment and the flow-through, respectively.
        If not set, MaxQuant analysis will be skipped.'''
        os.makedirs(maxquant_analysis_dir, exist_ok=True)
        logger.info('No results were detected. Preparing for MaxQuant analysis...')

        # elution analysis
        if os.path.exists(f'{wd}/el.zip'):
            # "real run"
            unpack_zipped_data(f'{wd}/el.zip', el_raw_path)
        # else, it's example run. The data is already there.
        logger.info(f'el_raw_path content is: {os.listdir(el_raw_path)}')
        run_maxquant(dbs_dir, enzymes, 'elution', run_number, el_raw_path)
        edit_progress(html_path, progress=15)

        # flowthrough analysis
        if os.path.exists(f'{wd}/ft.zip'):
            # "real run"
            unpack_zipped_data(f'{wd}/ft.zip', ft_raw_path)
            edit_progress(html_path, progress=20)
        # else, it's example run. The data is already there.
        logger.info(f'ft_raw_path content is: {os.listdir(ft_raw_path)}')
        run_maxquant(dbs_dir, enzymes, 'flowthrough', run_number, ft_raw_path)

        time_passed = wait_for_maxquant(el_raw_path, html_path=html_path)
        time_passed = wait_for_maxquant(ft_raw_path, time_passed=time_passed, html_path=html_path)
        logger.info(f'MaxQuant is done. Total running time was {measure_time(time_passed)}.')
        edit_progress(html_path, progress=75)

        logger.info(f'Copying results to {output_dir}')
        shutil.copy(f'{el_raw_path}/combined/txt/peptides.txt', f'{wd}/el_peptides.txt')
        shutil.copy(f'{ft_raw_path}/combined/txt/peptides.txt', f'{wd}/ft_peptides.txt')

        # add peptides list also to output files
        shutil.copy(f'{el_raw_path}/combined/txt/peptides.txt', f'{output_dir}/text/el_peptides.txt')
        shutil.copy(f'{ft_raw_path}/combined/txt/peptides.txt', f'{output_dir}/text/ft_peptides.txt')
        logger.info('MaxQuant analysis results are ready!')
    else:
        logger.info('_'*100)
        logger.info('MaxQuant analysis results were detected! Skipping...')
        logger.info('_'*100)

    edit_progress(html_path, progress=80)

    # filter maxquant resulted peptides
    logger.info('Filtering MaxQuant results...')
    filtered_peptides_info = filter_peptides(elution_peptides_path, flowthrough_peptides_path, textual_output_path, min_fold)
    edit_progress(html_path, progress=85)

    # classify filtered peptides
    db_paths = [f'{dbs_dir}/{db_name}' for db_name in os.listdir(dbs_dir) if not db_name.startswith('.')]  # ignoe system files
    if not os.path.exists(f'{textual_output_path}/{CONSTS.CDR3_INFORMATIVE_CSV_NAME}'):
        logger.info('No classifications were detected. Classifying filtered results (non-info/info/cdr3 info)...')
        classify_peptides(db_paths, filtered_peptides_info, textual_output_path)
    else:
        logger.info('_'*100)
        logger.info('Classification analysis results were detected! Skipping...')
        logger.info('_'*100)
    edit_progress(html_path, progress=90)

    # plot results
    plots_output_dir = f'{output_dir}/figures'
    if not os.path.exists(plots_output_dir) or local_run:
        os.makedirs(plots_output_dir, exist_ok=True)
        logger.info('Plotting cdr3-informative peptides...')
        plot_results(textual_output_path, plots_output_dir)
    else:
        logger.info('Plots folder is ready. Skipping...')
    edit_progress(html_path, progress=95)

    if not os.path.exists(f'{output_dir}.zip'):
        logger.info('Zipping results...')
        shutil.make_archive(output_dir,
                            'zip',
                            wd,
                            CONSTS.OUTPUT_DIR_NAME)

    edit_progress(html_path, progress=100, active=False)

    # edit html
    logger.info('Editing success html...')
    edit_success_html(html_path, run_number, output_dir, dbs_dir)

    if os.path.exists(f'{wd}/user_email.txt'):
        with open(f'{wd}/user_email.txt') as f:
            email = f.read().rstrip()
        notify_user('is ready', email, run_number, start_time, time())

    logger.info('_'*100)
    logger.info(f'PASA is done!')


def unpack_zipped_data(zipped_folder, target):
    logger.info(zipped_folder)
    print(zipped_folder)
    shutil.unpack_archive(zipped_folder, extract_dir=target)
    for file in os.listdir(target):
        file_path = f'{target}/{file}'
        try:
            shutil.unpack_archive(file_path, extract_dir=target)
            logger.info(f'Removing {file_path}')
            print(f'Removing {file_path}')
            os.remove(file_path)  # remove "leftovers"
        except:
            pass

    # make sure there are no irrelevant files after dbs extraction
    for file in os.listdir(target):
        file_path = f'{target}/{file}'
        if file.startswith(('_', '.')):
            logger.info(f'Removing {file_path}')
            print(f'Removing {file_path}')
            try:
                os.remove(file_path)  # ignore system files, such as __MACOSX
            except:
                pass
            try:
                shutil.rmtree(file_path)  # ignore system files, such as __MACOSX
            except:
                pass


if __name__ == '__main__':
    try:
        args = parse_arguments()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)

        logger.info(f'Starting a new job of {CONSTS.WEBSERVER_NAME.upper()}...')

        logger.info(args.maxquant)
        main(args.working_dir, args.min_fold,
             args.maxquant, args.enzymes)

    except Exception as e:
        html_path = f'{args.working_dir}/{CONSTS.RESULT_WEBPAGE_NAME}'
        run_number = args.working_dir.split('/')[-1]  # e.g., /bioseq/data/results/pasa/1234  ->  1234
        edit_failure_html(html_path, run_number, e.args[0])
