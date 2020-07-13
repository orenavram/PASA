import os
dbs_dir = '/Users/Oren/Dropbox/Projects/PASA/data/aya_dbs_trimmed_header'

for file in os.listdir(dbs_dir):
    if not file.endswith('fatsa'):
        continue
    with open(f'{dbs_dir}/{file}') as f:
        i = 1
        cleaned_headers = f.readline().split('|')[0].rstrip() + f'_{i}\n'
        for line in f:
            if line.startswith('>'):
                i += 1
                cleaned_headers += '\n' + line.split('|')[0].rstrip() + f'_{i}\n'

            else:
                cleaned_headers += line.rstrip()

        cleaned_headers += '\n'

    with open(f'{dbs_dir}/{file}', 'w') as f:
        f.write(cleaned_headers)

