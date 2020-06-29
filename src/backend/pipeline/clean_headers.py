import os
dbs_dir = '/Users/Oren/Dropbox/Projects/PASA/data/tmp'

for file in os.listdir(dbs_dir):
    if not file.endswith('fasta'):
        continue
    cleaned_headers = ''
    with open(f'{dbs_dir}/{file}') as f:
        cleaned_headers = f.readline()
        for line in f:
            if line.startswith('>'):
                cleaned_headers += '\n' + line  # .split('|')[0].rstrip() + '\n'
            else:
                cleaned_headers += line.rstrip()

        cleaned_headers += '\n'

    with open(f'{dbs_dir}/{file}', 'w') as f:
        f.write(cleaned_headers)

