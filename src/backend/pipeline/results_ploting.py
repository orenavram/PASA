import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import CONSTANTS as CONSTS
import logging
logger = logging.getLogger()
from adjustText import adjust_text

def upload_cdr3_data(data_path, token_numbers, delimiter=',', contains_header=True, manipulate=lambda x: x):
    frequency_counter = {}

    with open(data_path) as f:
        if contains_header:
            # skip header
            f.readline()

        for line in f:
            tokens = line.strip().split(delimiter)
            token = ' '.join([manipulate(tokens[token_number]) for token_number in token_numbers])

            if not token:
                token = 'NA'

            token = token.replace('IGH', '')  # show combinations as "V5 D3 J4" instead of "IGHV5 IGHD3 IGHJ4"
            frequency_counter[token] = frequency_counter.get(token, 0) + 1

    return frequency_counter



def generate_pie_chart(data_path, out_path, token_number):

    isotype_frequency_counter = upload_cdr3_data(data_path, token_number)
    if not isotype_frequency_counter:
        logger.info(f'An empty file was detected. Nothing to plot for {data_path}')
        return

    isotypes = sorted(isotype_frequency_counter, key=isotype_frequency_counter.get, reverse=True)
    portions = [isotype_frequency_counter[isotype] for isotype in isotypes]

    portions_percents = [100*portions[i]/sum(portions) for i in range(len(portions)) if portions[i] != 0]
    isotypes = [isotypes[i] for i in range(len(isotypes)) if portions[i] != 0]
    labels = [f'{isotypes[i]} ({portions_percents[i]:.3f}%)' for i in range(len(portions_percents))]

    patches, texts = plt.pie(portions_percents, colors=sns.color_palette('colorblind'), counterclock=False)

    plt.legend(patches, labels, loc="best")
    # plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.title('Isotype Distribution')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_barplot(data_path, out_path, token_numbers, manipulate=lambda x: x, sorting_key=lambda x: x, x_label='', y_label='% of total\n', as_proportions=True, rotation=0):

    logger.info(f'Ploting proteomics_vs_genetics_plot of {out_path}')
    frequency_counter = upload_cdr3_data(data_path, token_numbers, manipulate=manipulate)

    if not frequency_counter:
        logger.info(f'An empty file was detected. Nothing to plot for {data_path}')
        return

    x_values = sorted(frequency_counter, key=sorting_key)
    y_values = [frequency_counter[x] for x in x_values]

    if as_proportions:
        sum_y_values = sum(y_values)/100 # to get y's in percents
        y_values = [round(y/sum_y_values, 2) for y in y_values]

    barplot = sns.barplot(x_values, y_values, color='mediumseagreen')
    #print(x_label, len(x_values))
    # if len(x_values) < 20:
    #     fontsize = 10
    # elif len(x_values) < 60:
    #     fontsize = 5
    # else:
    #     fontsize = 2
    barplot.set_xticklabels(x_values, rotation=rotation)#, fontsize=fontsize)

    barplot.set_xlabel(x_label)
    barplot.set_ylabel(y_label)
    # barplot.set_title(title)

    ylim = [0, 1.2*max(y_values)]
    barplot.set_ylim(ylim)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def proteomics_vs_genetics_plot(data_path, figure_path, normalization_factor=1e6):

    df = pd.read_csv(data_path)
    # db_name = os.path.splitext(data_path.split('/')[-1])[0]
    if not len(df):
        logger.info(f'An empty file was detected. Nothing to plot for {data_path}')
        return

    df.iloc[:, 1] = df.iloc[:, 1]*normalization_factor
    df.iloc[:, 2] = df.iloc[:, 2]*normalization_factor
    muted_red = sns.color_palette('muted')[3]
    max_lim = df.iloc[:, 1:3].max().max() + 10
    plt.plot([0, max_lim], [0,  max_lim], linestyle=':', color='grey')
    plt.scatter(df.iloc[:, 1], df.iloc[:, 2], color=muted_red)
    plt.xlim([0, max_lim])
    plt.ylim([0, max_lim])
    plt.xlabel(f'\nBCR-Seq relative intensity')
    plt.ylabel(f'Proteomics relative intensity\n')
    texts = []
    for i in range(len(df)):
        texts.append(plt.text(df.iloc[i, 1], df.iloc[i, 2], df.iloc[i, 0],
                              horizontalalignment='left', color='black',
                              weight='semibold', fontsize=6))
    # adjust_text(texts, autoalign=True)
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()


def plot_results(textual_output_path, output_dir):

    cdr3_informative_path = f'{textual_output_path}/{CONSTS.CDR3_INFORMATIVE_CSV_NAME}'
    try:
        generate_pie_chart(cdr3_informative_path,
                           f'{output_dir}/isotype_distribution.png',
                           [6])  # isotype token. See peptides_classification.py.
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/cdr3_length.png',
                     [2],  # cdr3 token. See peptides_classification.py.
                     lambda x: str(len(x)), int, x_label='\nCDR3 length (AA level)')
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/v_usage.png',
                     [3], x_label='\nV family')  # IGHV token. See peptides_classification.py.
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/d_usage.png',
                     [4], x_label='\nD family')  # IGHD token. See peptides_classification.py.
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/j_usage.png',
                     [5], x_label='\nJ family')  # IGHJ token. See peptides_classification.py.
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/vd_usage.png',
                     [3, 4], rotation=90, x_label='\nVD family combination')
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/vj_usage.png',
                     [3, 5], rotation=90, x_label='\nVJ family combination')
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/dj_usage.png',
                     [4, 5], rotation=90, x_label='\nDJ family combination')
    except:
        pass

    try:
        plot_barplot(cdr3_informative_path, f'{output_dir}/vdj_usage.png',
                     [3, 4, 5], rotation=90, x_label='\nVDJ family combination')
    except:
        pass

    for raw_file_name in os.listdir(textual_output_path):
        if 'proteomic_intensities' not in raw_file_name:
            continue
        figure_name = f'{os.path.splitext(raw_file_name)[0]}.png'
        logger.info(f'Ploting {figure_name}')
        text_file_path = os.path.join(textual_output_path, raw_file_name)
        try:
            proteomics_vs_genetics_plot(text_file_path, f'{output_dir}/{figure_name}')
        except:
            logger.warning(f'Failed to plot {figure_name}')


# generate_pie_chart('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/outputs/text/informative.csv',
#                    '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/colorblind_pie_chart.png', [6])

# plot_barplot('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/outputs/text/informative.csv',
#              '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/cdr3_length.png',
#              [2],  # cdr3 token. See peptides_classification.py.
#              lambda x: str(len(x)))
#
#
# plot_barplot('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/informative.csv',
#              '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/v_subgroups.png',
#              [3])    # IGHV token. See peptides_classification.py.
#
# plot_barplot('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/informative.csv',
#              '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/vd_subgroups.png',
#              [3, 4], rotation=90)    # IGHV token. See peptides_classification.py.
#
# plot_barplot('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/informative.csv',
#              '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/vdj_subgroups.png',
#              [3, 4, 5], rotation=90)    # IGHV token. See peptides_classification.py.
#
#
# proteomics_vs_genetics_plot('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/informative.csv',
#                    '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/proteomics_vs_genetics.png')


