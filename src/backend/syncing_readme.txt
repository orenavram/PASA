rsync -avzr --exclude *pycache* /Users/Oren/Dropbox/Projects/PASA/src/backend/pipeline/* bioseq@powerweb1.tau.ac.il:/bioseq/pasa/pipeline/
rsync -avzr --exclude *pycache* /Users/Oren/Dropbox/Projects/PASA/src/backend/cgi/* bioseq@powerweb1.tau.ac.il:/var/www/cgi-bin/pasa/
rsync -avzr --exclude *pycache* /Users/Oren/Dropbox/Projects/PASA/src/frontend/* bioseq@powerweb1.tau.ac.il:/var/www/html/pasa/


# find all runs with a specific email
grep <email> /bioseq/data/results/microbializer/*/user_email.txt | awk -F ':' '{print $1}' | xargs ls -lt
