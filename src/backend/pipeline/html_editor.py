import os
import subprocess
import CONSTANTS as CONSTS
import logging

from auxiliaries import edit_progress

logger = logging.getLogger()

def add_closing_html_tags(html_path, CONSTS, run_number):
    with open(html_path, 'a') as f:
        f.write(f'''
                <hr>
                <h4 class=footer>
                    <p align=\'center\'>Questions and comments are welcome! Please
                        <span class="admin_link"> 
                        <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number%20{run_number}">contact us</a> 
                        </span>
                    </p>
                </h4>
                <div id="bottom_links" align="center"><span class="bottom_link">
                <a href="{CONSTS.WEBSERVER_URL}" target="_blank">Home</a> 
                 | 
                <a href="{CONSTS.WEBSERVER_URL}/overview.html" target="_blank">Overview</a> 
            </span>
        </div>
        <br><br><br>
    </body>
</html>
                ''')
        f.flush()

    # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
    from time import sleep
    sleep(max(20, 2 * CONSTS.RELOAD_INTERVAL))
    with open(html_path) as f:
        html_content = f.read()
    html_content = html_content.replace(CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')
    with open(html_path, 'w') as f:
        f.write(html_content)


def get_str_of_table_row_with_result_links(output_dir, file_name, str_to_show_on_html, final_output_dir_name=CONSTS.OUTPUT_DIR_NAME):
    result = f'''\t<tr>
                            <td>'''
    if os.path.exists(f'{output_dir}/{file_name}'):
        result += f'''
                            <a href="{final_output_dir_name}/figures/{file_name}" target="_blank">
                                {str_to_show_on_html}
                            </a>'''
    else: # disable link
        logger.info(f'{output_dir}/{file_name} is not available.')
        result += f'''
                                {str_to_show_on_html} - result is not available'''

    result += '''
                            </td>
                        </tr>\n\t\t\t\t\t'''

    return result


def edit_success_html(html_path, run_number, output_dir, dbs_dir):
    plots_output_dir = f'{output_dir}/figures'

    genetics_vs_proteomics_raw = [f'proteomic_intensities_with_respect_to_{os.path.splitext(file)[0]}.png'
                                  for file in sorted(os.listdir(dbs_dir))]

    result_figures_file_names = ['isotype_distribution.png',
                                 'cdr3_length.png',
                                 'v_usage.png',
                                 'd_usage.png',
                                 'j_usage.png',
                                 'vd_usage.png',
                                 'vj_usage.png',
                                 'dj_usage.png',
                                 'vdj_usage.png',
                                 *genetics_vs_proteomics_raw]

    genetics_vs_proteomics_png = []
    for file in genetics_vs_proteomics_raw:
        file = os.path.splitext(file)[0]
        file = file.replace('_', ' ')
        file = file[0].upper() + file[1:]
        genetics_vs_proteomics_png.append(f'{file} (scatter plot)')
    strs_to_show_on_html = ['Isotype distribution (pie chart)',
                            'CDR3 length distribution (bar plot)',
                            'V usage distribution (bar plot)',
                            'D usage distribution (bar plot)',
                            'J usage distribution (bar plot)',
                            'VD usage distribution (bar plot)',
                            'VJ usage distribution (bar plot)',
                            'DJ usage distribution (bar plot)',
                            'VDJ usage distribution (bar plot)',
                            *genetics_vs_proteomics_png]

    html_text = ''

    try: # Look for the initial file (generated by the cgi) so we can read and parse it.
        with open(html_path) as f:
            # html_text = f.read()
            for line in f:
                html_text += line
                if '<!--result-->' in line:
                    break
        html_text = html_text.replace('RUNNING', 'FINISHED').replace(f'{CONSTS.WEBSERVER_NAME.upper()} is now processing your request. This page will be automatically updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, the output will appear below. ', '')
    except FileNotFoundError as e:
        e.strerror = f"Couldn't find html file at: {html_path}"
        logger.warning(e.strerror)
        raise e

    html_text += f'''\t\t\t<div class="container" style="{CONSTS.CONTAINER_STYLE}">\n
                <h2>RESULTS:</h2>
                <table class="table">
                    <thead>
                        <tr><th><h3>Quick access to graphical results:</h3></th></tr>
                    </thead>
                    <tbody>
                        <h3><b><a href='{CONSTS.OUTPUT_DIR_NAME}.zip' target='_blank'>
                            Download results (both textual & visual)
                        </a></b></h3>\n\t\t\t\t\t'''

    for i in range(len(result_figures_file_names)):
        html_text += get_str_of_table_row_with_result_links(plots_output_dir, result_figures_file_names[i], strs_to_show_on_html[i])

    html_text += f'''</tbody>
                </table>'''
    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()

    add_closing_html_tags(html_path, CONSTS, run_number)


def edit_failure_html(html_path, run_number, msg):
    edit_progress(html_path, active=False)
    html_text = ''
    try:
        with open(html_path) as f:
            html_text = f.read()
        # The initial file exists (generate by the cgi) so we can read and parse it.
        html_text = html_text.replace('RUNNING', 'FAILED').replace(f'{CONSTS.WEBSERVER_NAME.upper()} is now processing your request. This page will be automatically updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, several links to the output files will appear below. ', '')
    except FileNotFoundError:
        logger.warning(f"Couldn't find html prefix at: {html_path}")

    html_text += f'''<div class="container" style="{CONSTS.CONTAINER_STYLE}">
            <br><br><br>
            <h2 align="justify">
            {CONSTS.WEBSERVER_NAME.upper()} failed due to the following reason:<br>
            <font color="red">{msg}</font><br><br>
            </h2>
            <h4>
            Please try to re-run your job after careful verification of your input or <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number:%20{run_number}">contact us</a> for further information
            </h4>
            <br><br>'''

    with open(html_path, 'w') as f:
        f.write(html_text)
        f.flush()

    logger.warning(f'{CONSTS.WEBSERVER_NAME.upper()} failed due to the following reason: {msg}')

    add_closing_html_tags(html_path, CONSTS, run_number)

