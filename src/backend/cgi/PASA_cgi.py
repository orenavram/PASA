#!/powerapps/share/centos7/Python-3.6.7/bin/python
import os
import shutil
import sys
import cgi
import cgitb
import subprocess
from time import time, ctime
from random import randint

if os.path.exists('/bioseq/pasa'): # remote run
    sys.path.insert(0, '/bioseq/pasa/pipeline')
    sys.path.insert(1, '/bioseq/bioSequence_scripts_and_constants')

import CONSTANTS as CONSTS  # from /bioseq/pasa/pipeline
from auxiliaries import edit_progress  # from /bioseq/pasa/pipeline
from email_sender import send_email  # from /bioseq/bioSequence_scripts_and_constants


def write_to_debug_file(cgi_debug_path_f, content):
    cgi_debug_path_f.write(f'{ctime()}: {content}\n')


def write_html_prefix(output_path, run_number):
    with open(output_path, 'w') as output_path_f:
        # html prefix. When creating a new one, copy it from the index page and replace the relevant values (i.e., web server name, etc...)
        output_path_f.write(f'''
<!DOCTYPE html>
<html lang="en">
    <head>
        <title>{CONSTS.WEBSERVER_NAME.upper()} Job #{run_number}</title>
        <meta http-equiv="cache-control" content="no-cache, must-revalidate, post-check=0, pre-check=0" />
        <meta http-equiv="cache-control" content="max-age=0" />
        <meta http-equiv="expires" content="0" />
        <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
        <meta http-equiv="pragma" content="no-cache" />
        {CONSTS.RELOAD_TAGS}
            
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
        <link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">
        <link rel="stylesheet" href="{CONSTS.WEBSERVER_URL}/css/general.css">
        <link rel="stylesheet" href="{CONSTS.WEBSERVER_URL}/css/nav.css">
        <link rel="shortcut icon" type="image/x-icon" href="{CONSTS.WEBSERVER_URL}/pics/logo.png" />

        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
        <script src="https://code.jquery.com/jquery-1.10.2.js"></script>
    </head>
    <body>
        <nav role="navigation" class="navbar navbar-inverse navbar-fixed-top" id="nav">
            <div class="jumbotron" id="jumbo">
                <div class="container">
                    <div class="row" id="title-row">
                        <div class="col-md-1">
                        </div>
                        <div class="col-md-1">
                            <img src="{CONSTS.WEBSERVER_URL}/pics/logo.png" id="antibody_image" class="img-rounded">
                        </div>
                        <div class="col-md-10">
                            <span id="pasa-title">PASA</span>&nbsp;&nbsp;&nbsp;&nbsp;<span id="sub-title"><b>Proteomics Analysis of Serum Antibodies </b></span><br>
                        </div>
                    </div>
                </div>
            </div>
        </nav>
        <div id="behind-nav-bar-results">
        </div>
        <br><br>
        <div class="container" style="font-size: 17px; {CONSTS.CONTAINER_STYLE}"  align="justify"> 
            <H1 align=center>Job status: <FONT color='red'>QUEUED</FONT></h1>
            <br>
            {CONSTS.PROGRESS_BAR_TAG}
            {CONSTS.WEBSERVER_NAME.upper()} is now processing your request. This page will be automatically updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, the output will appear below. A link to this page was sent to your email in case you wish to view these results at a later time without recalculating them. Please note that the results will be kept in the server for 3 months.
            <br><br>
        </div>\n\t\t''')
        output_path_f.flush()


def append_running_parameters_to_html(output_html_path, elution_file_name, flowthrough_file_name,
                                      el_peptides_file_name, ft_peptides_file_name, database_name,
                                      digestion_enzymes, enrichment_threshold, job_title):

    with open(output_html_path, 'a') as output_path_f:
        output_path_f.write(f'<div class="container" style="{CONSTS.CONTAINER_STYLE}">\n')

        # optional param row
        if job_title != '':
            output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
            output_path_f.write(f'<b>Job title: </b>{job_title}\n')
            output_path_f.write('</div></div>\n')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
        output_path_f.write(f'<b>BCR-Seq dataset: </b>{database_name}\n')
        output_path_f.write('</div></div>\n')

        if el_peptides_file_name != None:
            output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
            output_path_f.write(f'<b>Elution peptides list: </b>{el_peptides_file_name}\n')
            output_path_f.write('</div></div>\n')
            output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
            output_path_f.write(f'<b>Flow-through peptides list: </b>{ft_peptides_file_name}\n')
            output_path_f.write('</div></div>\n')
        else:
            output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
            output_path_f.write(f'<b>Elution raw dataset: </b>{elution_file_name}\n')
            output_path_f.write('</div></div>\n')

            output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
            output_path_f.write(f'<b>Flow-through raw dataset: </b>{flowthrough_file_name}\n')
            output_path_f.write('</div></div>\n')

            output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
            output_path_f.write(f'<b>Digestion enzyme: </b>{digestion_enzymes}\n')
            output_path_f.write('</div></div>\n')

        # mandatory param row
        output_path_f.write('<div class="row" style="font-size: 20px;"><div class="col-md-12">\n')
        output_path_f.write(f'<b>Min enrichment ratio: </b>{enrichment_threshold}\n')
        output_path_f.write('</div></div>\n')

        output_path_f.write('</div><br>\n')
        output_path_f.write('\n\n<!--result-->\n\n\t\t')

        output_path_f.flush()


def write_cmds_file(cmds_file, parameters, run_number):
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
    new_line_delimiter = '!@#'

    required_modules = ' '.join(['mono/mono-5.18.0.225', 'python/python-3.6.7'])
    with open(cmds_file, 'w') as f:
        f.write(f'module load {required_modules};')
        f.write(new_line_delimiter)
        f.write(f'python {CONSTS.MAIN_SCRIPT} {parameters}\tpasa_{run_number}')
        f.write('\n')


def save_file_to_disk(cgi_debug_path_f, form, wd, file_type):
    write_to_debug_file(cgi_debug_path_f, f'{file_type} in form: {file_type in form}')
    if form[file_type].filename in os.listdir(wd):
        write_to_debug_file(cgi_debug_path_f, f'changing file name of :\n{form[file_type].filename}')
        # form[file_type].filename = get_alternative_name(file_type, form, wd)

    file_name = form[file_type].filename
    write_to_debug_file(cgi_debug_path_f, f'file name is:\n{file_name}')
    # data = form[file_type].value
    write_to_debug_file(cgi_debug_path_f, f'{file_name} first 100 chars are: {form[file_type].value[:100]}\n')
    data_path = os.path.join(f'{wd}/{file_type}{os.path.splitext(file_name)[-1]}')
    with open(data_path, 'wb') as data_f:
        data_f.write(form[file_type].value)
    write_to_debug_file(cgi_debug_path_f, f'Uploaded data was saved to {data_path} successfully\n')
    return data_path, file_name


def get_alternative_name(file_type, form, wd):
    new_name = form[file_type].filename
    i = 1
    while new_name in os.listdir(wd):
        prefix, suffix = os.path.splitext(form[file_type].filename)
        new_name = f'{prefix}_{i}{suffix}'
        i += 1
    return new_name


def run_cgi():

    # prints detailed error report on BROWSER when backend crashes
    # This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
    cgitb.enable()

    # print_hello_world() # for debugging
    form = cgi.FieldStorage()  # extract POSTed object

    # adding 20 random digits to prevent users guess a number and see data that are not their's
    run_number = str(round(time())) + str(randint(10 ** 19, 10 ** 20 - 1))
    # run_number = 'test'
    output_url = f'{CONSTS.WEBSERVER_RESULTS_URL}/{run_number}/{CONSTS.RESULT_WEBPAGE_NAME}'

    if form['example_page'].value == 'yes':
        run_number = 'example'  # str(round(time())) + str(randint(1000,9999)) # adding 4 random figures to prevent users see data that are not their's
        output_url = os.path.join(f'{CONSTS.WEBSERVER_URL}/{run_number}/results.html')

    wd = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, run_number)
    os.makedirs(wd, exist_ok=True)

    output_html_path = os.path.join(wd, CONSTS.RESULT_WEBPAGE_NAME)
    cgi_debug_path = os.path.join(wd, 'cgi_debug.txt')

    page_is_ready = os.path.exists(output_html_path)
    if not page_is_ready:
        write_html_prefix(output_html_path, run_number)  # html's prefix must be written BEFORE redirecting...

    print(f'Location: {output_url}')  # Redirects to the results url. MUST appear before any other print.
    print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
    sys.stdout.flush()  # must be flushed immediately!!!

    try:
        # first, make sure the submission is from a "real user"
        exit_if_bot(form, wd)

        # a file handler for cgi log
        cgi_debug_path_f = open(cgi_debug_path, 'w')

        write_to_debug_file(cgi_debug_path_f, f'{"#" * 100}\n{ctime()}: A new CGI request has been received!\n')

        # extract form's values:
        email = ''
        if 'email' in form and form['email'].value != '':
            email = form['email'].value.strip()

        # Send me a notification email every time there's a new request
        send_email(smtp_server=CONSTS.SMTP_SERVER,
                   sender=CONSTS.ADMIN_EMAIL,
                   receiver=f'{CONSTS.OWNER_EMAIL}',
                   subject=f'{CONSTS.WEBSERVER_NAME.upper()} - A new job has been submitted: {run_number}',
                   content=f'{os.path.join(CONSTS.WEBSERVER_URL, "results", run_number, "cgi_debug.txt")}\n'
                   f'{os.path.join(CONSTS.WEBSERVER_URL, "results", run_number, CONSTS.RESULT_WEBPAGE_NAME)}\n'
                   f'{email if email else "NO EMAIL"}')

        peek_form(cgi_debug_path_f, form)

        job_title = ''
        if 'job_title' in form and form['job_title'].value != '':
            # 'job_title' is NOT in form when the text field is disabled...
            job_title = form['job_title'].value.strip()

        elution_file_name = None
        flowthrough_file_name = None
        el_peptides_file_name = None
        ft_peptides_file_name = None
        digestion_enzyme = None
        maxquant_analysis_is_needed = True
        if form['example_page'].value == 'yes':  # example data
            write_to_debug_file(cgi_debug_path_f, f'Linking example data FROM {CONSTS.EXAMPLE_DATA_PATH} TO {wd}\n')
            copy_example_data(wd, cgi_debug_path_f)
            database_file_name = CONSTS.EXAMPLE_DB_NAME
            elution_file_name = 'Elution example data'  # CONSTS.ELUTION_FILE_NAMES
            flowthrough_file_name = 'Flow-through example data'  # CONSTS.FLOWTHROUGH_FILE_NAMES
            digestion_enzyme = 'Trypsin'  # TODO: check what is recieved when sending several values
            enrichment_threshold = 5
        else:
            write_to_debug_file(cgi_debug_path_f, f'\n{"#"*80}\nuploading data\n')
            database_file_path, database_file_name = save_file_to_disk(cgi_debug_path_f, form, wd, 'db')
            enrichment_threshold = form['enrichment_threshold'].value.strip()
            if 'el' in form and form['el'].filename != '':
                elution_file_path, elution_file_name = save_file_to_disk(cgi_debug_path_f, form, wd, 'el')
                flowthrough_file_path, flowthrough_file_name = save_file_to_disk(cgi_debug_path_f, form, wd, 'ft')
                digestion_enzyme = form['enzyme'].value.strip()  # TODO: check what is recieved when sending several values
            else:
                maxquant_analysis_is_needed = False
                _, el_peptides_file_name = save_file_to_disk(cgi_debug_path_f, form, wd, 'el_peptides')
                _, ft_peptides_file_name = save_file_to_disk(cgi_debug_path_f, form, wd, 'ft_peptides')


        write_to_debug_file(cgi_debug_path_f, f'ls of {wd} yields:\n{os.listdir(wd)}\n')

        write_to_debug_file(cgi_debug_path_f, f'{ctime()}: Writing running parameters to html...\n')

        if not page_is_ready:
            append_running_parameters_to_html(output_html_path, elution_file_name, flowthrough_file_name,
                                              el_peptides_file_name, ft_peptides_file_name, database_file_name,
                                              digestion_enzyme, enrichment_threshold, job_title)

        write_to_debug_file(cgi_debug_path_f, f'{ctime()}: Running parameters were written to html successfully.\n')

        parameters = f'{wd} --min-fold {enrichment_threshold}'
        if maxquant_analysis_is_needed:
            parameters += f' --enzymes {form["enzyme"].value} -mq'

        cmds_file = os.path.join(wd, 'qsub.cmds')
        write_cmds_file(cmds_file, parameters, run_number)

        job_id_file = os.path.join(wd, 'job_id.txt')

        submission_cmd = f'{CONSTS.Q_SUBMITTER_SCRIPT} {cmds_file} {wd} -q pupkowebr --verbose > {job_id_file}'

        if not page_is_ready:
            write_to_debug_file(cgi_debug_path_f, f'\nSUBMITTING JOB TO QUEUE:\n{submission_cmd}\n')
            subprocess.call(submission_cmd, shell=True)
        else:
            write_to_debug_file(cgi_debug_path_f, f'\nPage already exists! no need to run the analysis again\n')

        user_email_file = os.path.join(wd, CONSTS.EMAIL_FILE_NAME)
        if email != '':
            with open(user_email_file, 'w') as email_f:
                email_f.write(f'{email}\n')

            try:
                # Send the user a notification email for their submission
                notify_user(job_title, database_file_name,
                            elution_file_name, flowthrough_file_name,
                            el_peptides_file_name, ft_peptides_file_name,
                            digestion_enzyme, enrichment_threshold, run_number, email)
            except:
                write_to_debug_file(cgi_debug_path_f, f'\nFailed sending notification to {email}\n')

        else:
            try:
                os.remove(user_email_file)      #for example mode
            except OSError:
                pass

        job_title_file = os.path.join(wd, 'job_title.txt')
        if job_title != '':
            with open(job_title_file, 'w') as job_f:
                job_f.write(f'{job_title}\n')
        else:
            try:
                os.remove(job_title_file)      #for example mode
            except OSError:
                pass

        write_to_debug_file(cgi_debug_path_f, f'\n\n{"#"*50}\nCGI finished running!\n{"#"*50}\n')

        cgi_debug_path_f.close()

    except Exception as e:
        edit_progress(output_html_path, active=False)

        msg = 'CGI crashed before the job was submitted :('
        with open(output_html_path) as f:
            html_content = f.read()
        html_content = html_content.replace('QUEUED', 'FAILED').replace('RUNNING', 'FAILED')
        html_content += f'<br><br><br><center><h2><font color="red">{msg}</font><br><br>Please try to re-run your job or <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number%20{run_number}">contact us</a> for further information</h2></center><br><br>\n</body>\n</html>\n'
        with open(output_html_path, 'w') as f:
            f.write(html_content)

        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        with open(cgi_debug_path, 'w') as cgi_debug_path_f:
            write_to_debug_file(cgi_debug_path_f, f'\n{"$"*100}\n\n{msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\ne.args[0]: {e.args[0]}\n\n{"$"*100}')

        # Send me a notification email every time there's a failure
        try:
            email = form['email'].value.strip() if form['email'].value.strip() else 'NO_EMAIL'
        except:
            email = 'NO_EMAIL'
        send_email(smtp_server=CONSTS.SMTP_SERVER,
                   sender=CONSTS.ADMIN_EMAIL,
                   receiver=f'{CONSTS.OWNER_EMAIL}',
                   subject=f'{CONSTS.WEBSERVER_NAME.upper()} job {run_number} by {email} has been failed!',
                   content=f"{email}\n\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, CONSTS.RESULT_WEBPAGE_NAME)}\n"
                            f"\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, 'cgi_debug.txt')}")

        # logger.info(f'Waiting {2*CONSTS.RELOAD_INTERVAL} seconds to remove html refreshing headers...')
        # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
        from time import sleep

        sleep(2 * CONSTS.RELOAD_INTERVAL)
        with open(output_html_path) as f:
            html_content = f.read()
        html_content = html_content.replace(CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')
        with open(output_html_path, 'w') as f:
            f.write(html_content)

    # logging submission
    with open(CONSTS.SUBMISSIONS_LOG, 'a') as f:
        f.write(f'{email}\t{run_number}\t{ctime()}\n')

    with open(cgi_debug_path, 'a') as f:  # for cgi debugging
        f.write(f'{ctime()}: Submission was documented in \n')


def notify_user(job_title, database_file_name,
                elution_file_name, flowthrough_file_name,
                el_peptides_file_name, ft_peptides_file_name,
                digestion_enzyme, enrichment_threshold, run_number, email):
    job_name = f'Job title: {job_title}\n' if job_title else ''
    notification_content = f'Your submission details are:\n\n{job_name}BCR-Seq dataset: {database_file_name}\n'

    if el_peptides_file_name != None:
        notification_content += f'Elution peptides list: {el_peptides_file_name}\nFlow-through peptides list: {ft_peptides_file_name}'
    else:
        notification_content += f'Elution raw dataset: {elution_file_name}\nFlow-through raw dataset: {flowthrough_file_name}\n'

    notification_content += f'Digestion enzyme: {digestion_enzyme}\nMin enrichment ratio: {enrichment_threshold}\n\n'

    notification_content += f'Once the analysis will be ready, we will let you know! Meanwhile, you can track the ' \
        f'progress of your job at:\n{CONSTS.WEBSERVER_URL}/results/{run_number}/{CONSTS.RESULT_WEBPAGE_NAME}\n\n'

    send_email(smtp_server=CONSTS.SMTP_SERVER,
               sender=CONSTS.ADMIN_EMAIL,
               receiver=f'{email}',
               subject=f'{CONSTS.WEBSERVER_NAME.upper()} - your job has been submitted! (Run number: {run_number})',
               content=notification_content)


def copy_example_data(target_dir, cgi_debug_path_f):
    maxquant_analysis_dir = f'{target_dir}/{CONSTS.MAXQUANT_DIR_NAME}'
    os.makedirs(maxquant_analysis_dir, exist_ok=True)
    db_path = f'{CONSTS.EXAMPLE_DATA_PATH}/{CONSTS.EXAMPLE_DB_NAME}'

    for data_type in ['el', 'ft']:
        raw_path = f'{maxquant_analysis_dir}/{data_type}'
        write_to_debug_file(cgi_debug_path_f, f'Creating {raw_path}...')
        os.makedirs(raw_path, exist_ok=True)
        write_to_debug_file(cgi_debug_path_f, f'Copying {db_path} to {raw_path}')
        shutil.copy(db_path, f'{target_dir}/{CONSTS.EXAMPLE_DB_NAME}')
        for file in os.listdir(f'{CONSTS.EXAMPLE_DATA_PATH}/{data_type}'):
            alias = f'{raw_path}/{file}'
            write_to_debug_file(cgi_debug_path_f, f'Soft linking {alias}')
            #TODO: copy example data as zip
            soft_link_cmd = f'ln -sf {CONSTS.EXAMPLE_DATA_PATH}/{data_type}/{file} {alias}\n'
            write_to_debug_file(cgi_debug_path_f, f'{soft_link_cmd}')
            subprocess.run(soft_link_cmd, shell=True)


def exit_if_bot(form, wd):
    # email field should ALWAYS exist in the form (even if it's empty!!).
    # If it's not there, someone sent a request not via the website so they should be blocked.
    # confirm_email is a hidden field that only spamming bots might fill in...
    if 'email' not in form or ('confirm_email' in form and form['confirm_email'].value != ''):
        shutil.rmtree(wd)
        exit()


def peek_form(cgi_debug_path_f, form):
    # for debugging
    sorted_form_keys = sorted(form.keys())
    write_to_debug_file(cgi_debug_path_f,
                        f'These are the keys that the CGI received:\n{"; ".join(sorted_form_keys)}\n\n')
    write_to_debug_file(cgi_debug_path_f, 'Form values are:\n')
    for key in sorted_form_keys:
        if not key.startswith(('el', 'ft', 'db')):
            if 'peptides' in key:
                write_to_debug_file(cgi_debug_path_f, f'100 first characters of {key} = {form[key].value[:100]}\n')
            else:
                write_to_debug_file(cgi_debug_path_f, f'{key} = {form[key]}\n')
    # for key in sorted_form_keys:
    #     if 'data' in key:
    #         write_to_debug_file(cgi_debug_path_f, f'100 first characters of {key} = {form[key].value[:100]}\n')
    cgi_debug_path_f.flush()


if __name__ == '__main__':
    run_cgi()