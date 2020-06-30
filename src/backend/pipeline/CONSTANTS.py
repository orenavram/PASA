#!/data/shared/python/anaconda3-5.1.0/bin/python3.6

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
# ADMIN_EMAIL = 'Yariv Wine <yariv.wine@gmail.com>'
SMTP_SERVER = 'mxout.tau.ac.il'
OWNER_EMAIL = 'orenavram@gmail.com'

# general variables
SERVERS_RESULTS_DIR = '/bioseq/data/results'
RELOAD_INTERVAL = 5
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

WEBSERVER_NAME = 'pasa'
WEBSERVER_URL = f'https://{WEBSERVER_NAME}.tau.ac.il'
WEBSERVER_PROJECT_PATH = f'/bioseq/{WEBSERVER_NAME}'
WEBSERVER_TITLE = 'Proteomics Analysis of Serum Antibodies'

WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, WEBSERVER_NAME)
WEBSERVER_HTML_DIR = f'/data/www/html/{WEBSERVER_NAME}'

WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')

Q_SUBMITTER_SCRIPT = '/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py'
MAIN_SCRIPT = os.path.join(f'{WEBSERVER_PROJECT_PATH}/pipeline/main.py')

OUTPUT_DIR_NAME = 'outputs'
RESULT_WEBPAGE_NAME = 'results.html'

#path to example data
EXAMPLE_DATA_PATH = f'{WEBSERVER_PROJECT_PATH}/example_data'
EXAMPLE_DB_NAME = 'db.zip'
EXAMPLE_ELUTION_RAW = 'el.zip'
EXAMPLE_FLOWTHROUGH_RAW = 'ft.zip'
ELUTION_FILE_NAMES = ['elution_data_1.raw', 'elution_data_2.raw', 'elution_data_3.raw']
FLOWTHROUGH_FILE_NAMES = ['flowthrough_data_1.raw', 'flowthrough_data_2.raw', 'flowthrough_data_3.raw']

# external programs
MAXQUANT_EXECUTABLE_PATH = f'/bioseq/{WEBSERVER_NAME}/MaxQuant.1.6.14/bin/MaxQuantCmd.exe'

# html styling
WEBSERVER_JUMBOTRON = f'&nbsp;&nbsp;&nbsp;&nbsp;<span id="server-title">{WEBSERVER_NAME}</span>&nbsp;&nbsp;&nbsp;&nbsp;<span id="sub-title">{WEBSERVER_TITLE}</span>'

CONTAINER_WIDTH = 'width: 850px'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
FONT_STYLE = 'font-family: Georgia'

CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}; {FONT_STYLE}'

PROGRESS_BAR_TAG = '''<div class="progress">
                <div class="progress-bar progress-bar-striped active" role="progressbar" aria-valuemin="0" aria-valuemax="100" style="width:1%">
                <!--MaxQuant Analysis ; Peptides Filtration ; Peptides Classification ; Plots Generation-->
                </div>
            </div>\t\t'''

# miscellaneous
MAXQUANT_DIR_NAME = 'maxquant_analysis'
