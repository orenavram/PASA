import os
import subprocess
import logging
logger = logging.getLogger()
from time import sleep

def measure_time(total):
    hours = total // 3600
    minutes = (total% 3600) // 60
    seconds = total % 60
    if hours != 0:
        return f'{hours}:{minutes:02}:{seconds:02} hours'
    elif minutes != 0:
        return f'{minutes}:{seconds:02} minutes'
    else:
        return f'{seconds} seconds'


def wait_for_maxquant(maxquant_wd, time_passed=0, sleeping_time=5, html_path='',
                      total_bar_for_maxquant=50, initial_shift=20):
    maxquant_running_times_file = f'{maxquant_wd}/combined/proc/#runningTimes.txt'
    maxquant_result_file = f'{maxquant_wd}/combined/txt/peptides.txt'
    maxquant_done_file = f'{maxquant_wd}/combined/proc/Finish_writing_tables 11.finished.txt'
    logger.info(f'Waiting for results at: {maxquant_result_file}')
    logger.info(f'Progress can be viewed at: {maxquant_running_times_file}')
    while True:
        sleep(sleeping_time)
        if html_path:
            try:
                # b'  54 maxquant_analysis/el/combined/proc/#runningTimes.txt\n  54 maxquant_analysis/ft/combined/proc/#runningTimes.txt\n 108 total\n'
                total_steps_done = sum(int(x) for x in subprocess.check_output(f'wc -l {maxquant_wd}/../*/combined/proc/#runningTimes.txt',
                                                                               shell=True).decode().rstrip().split()[:3][::2])
                # "#runningTimes.txt" file contains 54 rows when it's done. 2 files contains 108...
                edit_progress(html_path, progress=total_bar_for_maxquant * total_steps_done / 108 + initial_shift)
            except:
                logger.info('Could not find any #runningTimes.txt file to estimate progress...')
                pass
        time_passed += sleeping_time
        logger.info(f'\t{measure_time(time_passed)} have passed since started waiting.')
        if len([file for file in os.listdir(maxquant_wd) if file.endswith('ER')]) > 0:
            # ER file is ready, namely, job is done. Make sure it was successfully...
            if os.path.exists(maxquant_done_file) and os.path.exists(maxquant_result_file):
                break
            raise Exception('MaxQuant analysis failed.')

    return time_passed


def edit_progress(html_path, progress=None, active=True):
    result = ''
    if not os.path.exists('/bioseq'):
        # local run
        return

    with open(html_path) as f:
        for line in f:
            if 'progress-bar' in line:
                if progress:
                    # <div class="progress-bar ... style="width:0%">
                    line = line.split('style')[0]
                    line += f'style="width:{progress}%">\n'
                if not active:
                    # <div class="progress-bar progress-bar-striped bg-success active" ...
                    line = line.replace(' active', '')
            result += line

    with open(html_path, 'w') as f:
        f.write(result)

    # TODO: uncomment in production
    # sleep(3)
