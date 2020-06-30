import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()
legal_iso_types = ('M', 'G', 'A', 'A1', 'A2', 'E', 'D', 'unknown')
legal_chain_types = ('IGH', 'IGK', 'IGL')

class Record:


    def __init__(self, db_name, seq, chain_type, iso_type, cdr3, IGHV, IGHD, IGHJ, counts, reps_intensities):
        if chain_type not in legal_chain_types or iso_type not in legal_iso_types:
            raise ValueError

        self.db_num = db_name
        self.seq = seq
        self.chain_type = chain_type
        self.iso_type = iso_type
        self.cdr3 = cdr3
        self.IGHV = IGHV
        self.IGHD = IGHD
        self.IGHJ = IGHJ
        self.counts = counts
        self.reps_intensities = [int(x) for x in reps_intensities.split(',')]

    def __repr__(self):
        return '<' + str([self.db_num, self.seq, self.chain_type, self.iso_type, self.cdr3,
                    self.IGHV, self.IGHD, self.IGHJ, self.counts, self.reps_intensities])[1:-1] + '>'





