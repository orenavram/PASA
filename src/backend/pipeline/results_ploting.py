import matplotlib.pyplot as plt
import seaborn as sns


def upload_data(data_path, token_numbers, delimiter=',', contains_header=True, manipulate=lambda x: x):
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

            frequency_counter[token] = frequency_counter.get(token, 0) + 1

    return frequency_counter


def generate_pie_chart(data_path, out_path, token_number):

    isotype_frequency_counter = upload_data(data_path, token_number)

    isotypes = sorted(isotype_frequency_counter, key=isotype_frequency_counter.get, reverse=True)
    portions = [isotype_frequency_counter[isotype] for isotype in isotypes]

    portions_percents = [100*portions[i]/sum(portions) for i in range(len(portions)) if portions[i]!=0]
    isotypes = [isotypes[i] for i in range(len(isotypes)) if portions[i]!=0]
    labels = [f'{isotypes[i]} ({portions_percents[i]:.3f}%)' for i in range(len(portions_percents))]

    patches, texts = plt.pie(portions_percents, counterclock=False)
    plt.legend(patches, labels, loc="best")
    # plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    plt.title('Isotype Distribution')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_barplot(data_path, out_path, token_numbers, manipulate=lambda x: x, x_label='', y_label='% of total\n', as_proportions=True, rotation=0):

    frequency_counter = upload_data(data_path, token_numbers, manipulate=manipulate)

    x_values = sorted(frequency_counter)
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


def proteomics_vs_genetics_plot(data_path, out_path, token_number):
    frequency_counter = upload_data(data_path, token_number)
    pass


def plot_results(csv_file, output_dir):

    proteomics_vs_genetics_plot(csv_file, f'{output_dir}/proteomics_vs_genetics.png',
                                [2])  # cdr3 token. See peptides_classification.py.

    plot_barplot(csv_file, f'{output_dir}/cdr3_length.png',
                 [2],  # cdr3 token. See peptides_classification.py.
                 lambda x: str(len(x)))

    # TODO: separate plot for each chain type
    # plot_barplot(csv_file, f'{output_dir}/IGH_cdr3_length.png',
    #              [2],  # cdr3 token. See peptides_classification.py.
    #              lambda x: str(len(x)))
    #
    # plot_barplot(csv_file, f'{output_dir}/IGL_cdr3_length.png',
    #              [2],  # cdr3 token. See peptides_classification.py.
    #              lambda x: str(len(x)))
    #
    # plot_barplot(csv_file, f'{output_dir}/IGK_cdr3_length.png',
    #              [2],  # cdr3 token. See peptides_classification.py.
    #              lambda x: str(len(x)))

    plot_barplot(csv_file, f'{output_dir}/subgroups_v.png',
                 [3])  # IGHV token. See peptides_classification.py.

    plot_barplot(csv_file, f'{output_dir}/subgroups_d.png',
                 [4])  # IGHD token. See peptides_classification.py.

    plot_barplot(csv_file, f'{output_dir}/subgroups_j.png',
                 [5])  # IGHJ token. See peptides_classification.py.

    plot_barplot(csv_file, f'{output_dir}/subgroups_vd.png',
                 [3, 4], rotation=90)

    plot_barplot(csv_file, f'{output_dir}/subgroups_vj.png',
                 [3, 5], rotation=90)

    plot_barplot(csv_file, f'{output_dir}/subgroups_dj.png',
                 [4, 5], rotation=90)

    plot_barplot(csv_file, f'{output_dir}/subgroups_vdj.png',
                 [3, 4, 5], rotation=90)  # IGHV token. See peptides_classification.py.

    generate_pie_chart(csv_file, f'{output_dir}/pie_chart.png',
                       [6])  # isotype token. See peptides_classification.py.


# generate_pie_chart('/Users/Oren/Dropbox/Projects/PASA/linux_outputs/outputs/text/informative.csv',
                   # '/Users/Oren/Dropbox/Projects/PASA/linux_outputs/pie_chart.png')

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


