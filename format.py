
def main():
    import numpy as np
    import pandas as pd
    from difflib import SequenceMatcher
    from functools import partial

    df = pd.read_csv('/Users/patricknaylor/Desktop/kaggle/kaggle-NOVO/data/train.csv')
    corr = pd.read_csv('/Users/patricknaylor/Desktop/kaggle/kaggle-NOVO/data/train_updates_20220929.csv')
    corr_ids = list(corr['seq_id'])
    df = df[~df['seq_id'].isin(corr_ids)]

    wild_type = 'VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK'

    def similar(s, c, w, b, side): 
        if side == 'start':
            return SequenceMatcher(None, s[c][:b], w).ratio()
        elif side == 'end':
            return SequenceMatcher(None, s[c][-b:], w).ratio()
        else:
            raise ValueError('side must be `start` or `end`')
        
    def match_columns(wild, train):
        bounds = [1, 3, 5, 15, 30]
        for bound in bounds:

            start_wild = wild[:bound]
            end_wild = wild[-bound:]
            train[f'first{bound}'] = train.apply(partial(similar, c='protein_sequence', w=start_wild, b=bound, side='start'), axis=1)
            train[f'last{bound}'] = train.apply(partial(similar, c='protein_sequence', w=start_wild, b=bound, side='end'), axis=1)

        return train
	col = wild_type
	length = len(col)
	tenths = [int((length/10)*(x+1)) for x in range(10)]
	wild_list = []
	for idx, ten in enumerate(tenths):
	    col_tenth = col[idx:ten]
	    col_dict = {}
	    for letter in col_tenth:
	        if letter in col_dict:
	            col_dict[letter] += 1
	        else:
	            col_dict[letter] = 1
	    for letter in letters:
	        if letter not in col_dict:
	            col_dict[letter] = 0
	    for key, value in col_dict.items():
	        col_dict[key] = value/(length/10)
	    wild_list.append(col_dict)

    def tenths(train, letters):
        col = train['protein_sequence']
        length = len(col)
        tenths = [int((length/10)*(x+1)) for x in range(10)]
        tenth_list = []
        for idx, ten in enumerate(tenths):
            col_tenth = col[idx:ten]
            col_dict = {}
            for letter in col_tenth:
                if letter in col_dict:
                    col_dict[letter] += 1
                else:
                    col_dict[letter] = 1
            for letter in letters:
                if letter not in col_dict:
                    col_dict[letter] = 0
            for key, value in col_dict.items():
                col_dict[key] = value/(length/10)
            tenth_list.append(col_dict)
        labels = np.arange(0, 110, 10)
        for idx, (t, w) in enumerate(zip(tenth_list, wild_list)):
            for l in letters:
                train[f'{labels[idx]}-{labels[idx+1]}{l}_count'] = t[l] - w[l]
        return train

    df = df.apply(partial(tenths, letters=letters), axis=1)
    df = match_columns(wild_type, df)

    df['seq_len'] = df['protein_sequence'].str.len()
    df['seq_len'] = df['seq_len']/df['seq_len'].max()
    df['rel_stab'] = df['tm']/df['tm'].max()

if __name__ == '__main__':
    main()