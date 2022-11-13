def main()

    import numpy as np
    import pandas as pd
    from difflib import SequenceMatcher
    from functools import partial

    df = pd.read_csv('/Users/patricknaylor/Desktop/kaggle/kaggle-NOVO/data/train.csv')
    corr = pd.read_csv('/Users/patricknaylor/Desktop/kaggle/kaggle-NOVO/data/train_updates_20220929.csv')
    corr_ids = list(corr['seq_id'])
    df = df[~df['seq_id'].isin(corr_ids)]

    wild_type = 'VPVNPEPDATSVENVALKTGSGDSQSDPIKADLEVKGQSALPFDVDCWAILCKGAPNVLQRVNEKTKNSNRDRSGANKGPFKDPQKWGIKALPPKNPSWSAQDFKSPEEYAFASSLQGGTNAILAPVNLASQNSQGGVLNGFYSANKVAQFDPSKPQQTKGTWFQITKFTGAAGPYCKALGSNDKSVCDKNKNIAGDWGFDPAKWAYQYDEKNNKFNYVGK'
    print(df.head())
    letters = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def similar(s, c, w, b, side=None): 
        if side == 'start':
            return SequenceMatcher(None, s[c][:b], w).ratio()
        elif side == 'end':
            return SequenceMatcher(None, s[c][-b:], w).ratio()
        elif side == None:
            return SequenceMatcher(None, s[c][:], w).ratio()
        else:
            raise ValueError('side must be `start`, `end`, or None')
            
    def match_columns(wild, train):
        bounds = [1, 3, 5, 15, 30]
        for bound in bounds:
            start_wild = wild[:bound]
            end_wild = wild[-bound:]
            train[f'first{bound}'] = train.apply(partial(similar, c='protein_sequence', w=start_wild, b=bound, side='start'), axis=1)
            train[f'last{bound}'] = train.apply(partial(similar, c='protein_sequence', w=start_wild, b=bound, side='end'), axis=1)
        train['sim_ratio'] = train.apply(partial(similar, c='protein_sequence', w=start_wild, b=bound), axis=1)
        return train

    def tenth_strings(s, c, num, tenth):
        return s[c][num*s[tenth] : (num+1)*s[tenth]]

    def tenths_letters(train, letters):
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
            
        
            
            

    df['seq_len'] = df['protein_sequence'].str.len()
    df['seq_len_norm'] = df['seq_len']/df['seq_len'].max()

    df['len_tenth'] = df['seq_len']//10
    wild_tenth = len(wild_type)//10
    wild_tenths = []
    tenth_labels = []
    labels = np.arange(0, 10)
    for label in labels:
        df[f'{label}_tenth_str'] = df.apply(partial(tenth_strings, c='protein_sequence', num=label, tenth='len_tenth'), axis=1)
        tenth_labels.append(f'{label}_tenth_str')
        wild_tenths.append(wild_type[label*wild_tenth: (label+1)*wild_tenth])
        
    for lab, wild in zip(tenth_labels, wild_tenths):
        df[f'{lab}_match'] = df.apply(partial(similar, c=lab, w=wild, b=0, side=None), axis=1)

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

    df = df.apply(partial(tenths_letters, letters=letters), axis=1)

    df = match_columns(wild_type, df)

    df['rel_stab'] = df['tm']/df['tm'].max()

    rel_stab = df['rel_stab']

    removals = ['seq_id', 'protein_sequence', 'pH', 'data_source', 'tm', 'seq_len', 'len_tenth', '0_tenth_str', '1_tenth_str', '2_tenth_str', '3_tenth_str', '4_tenth_str', '5_tenth_str', '6_tenth_str', '7_tenth_str', '8_tenth_str', '9_tenth_str', 'rel_stab']
    df = df.drop(removals, axis=1)

    df.to_csv('/Users/patricknaylor/Desktop/kaggle/kaggle-NOVO/data/variables_1.csv')
    rel_stab.to_csv('/Users/patricknaylor/Desktop/kaggle/kaggle-NOVO/data/target.csv')



if __name__ == '__main__':
    main()