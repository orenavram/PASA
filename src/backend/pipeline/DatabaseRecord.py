class Record:

    def __init__(self, db_name, seq, chain_type, iso_type, cdr3, IGHV, IGHD, IGHJ, counts, reps_intensities):
        self.db_num = db_name
        self.seq = seq
        # self.ibd = ibd
        self.chain_type = chain_type
        self.iso_type = iso_type
        self.cdr3 = cdr3
        # self.cdr3_start = -1
        # self.cdr3_end = -1
        self.IGHV = IGHV
        self.IGHD = IGHD
        self.IGHJ = IGHJ
        self.counts = counts
        self.reps_intensities = [int(x) for x in reps_intensities.split(',')]
        # self.avg_counts = counts / len(self.reps_intensities) if reps_intensities != None else self.counts
        # self.enrichment = 0

    def __repr__(self):
        return '<' + str([self.db_num, self.seq, self.chain_type, self.iso_type, self.cdr3,
                    self.IGHV, self.IGHD, self.IGHJ, self.counts, self.reps_intensities])[1:-1] + '>'





