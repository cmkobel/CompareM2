#!/usr/bin/env python

__author__ = "Carl M. Kobel"


import pandas as pd


def pretty_print_snakemake_summary(path): 
    summary = pd.read_csv(path, delimiter="\t") 
    summary = summary.dropna(subset=['rule'])
    #print(summary)
    
    # Hide stuff
    summary = summary.loc[~summary['rule'].str.endswith('_download')]
    summary = summary.loc[~summary['rule'].str.contains('^copy$')]
    summary = summary.loc[~summary['rule'].str.contains('^compute_snp_dists$')]
    
    summary = summary.groupby(['rule','status']).size().reset_index().rename(columns={0:'count'})
    #print(summary)
    
    summary = summary.pivot(index = "rule", columns = "status", values = "count").rename_axis(None, axis=1).reset_index()
    # print(summary)
    
    summary = summary.fillna(0)
    summary['complete'] = summary['ok'] / (summary['ok']+summary['missing'])
    summary['complete %'] = (summary['complete'] * 100).astype(str) + " %"
    
    
    completed = summary.loc[summary['complete'] == 1.0]
    non_completed = summary.loc[summary['complete'] != 1.0]
    
    print("Completed rules:")
    print(completed[['rule', 'complete %']])
    print()
    print("Non-completed rules:")
    print(non_completed[['rule', 'complete %']])
    
    
    #print(summary)

    




# For quick development.
if __name__ == "__main__":
    pretty_print_snakemake_summary(".status.tsv")
