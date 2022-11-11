
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

    df = match_columns(wild_type, df)
    print(df)


if __name__ == '__main__':
    main()