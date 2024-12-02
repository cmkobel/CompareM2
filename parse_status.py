#!/usr/bin/env python

__author__ = "Carl M. Kobel"


import os
import pandas as pd


def present_status(path = ".status.tsv", delete_after_use = True): 
    
    print()
    summary = pd.read_csv(path, delimiter="\t") 
    summary = summary.dropna(subset=['rule'])
    #print(summary)
    
    # Hide irrelevant stuff
    summary = summary.loc[~summary["rule"].str.endswith("_download")]
    summary = summary.loc[~summary["rule"].str.contains("^copy$")]
    summary = summary.loc[~summary["rule"].str.contains("^compute_snp_dists$")]
    summary = summary.loc[~summary["rule"].str.contains("^report_env$")]
    summary = summary.loc[~summary["rule"].str.contains("^bakta_env$")]
    summary = summary.loc[~summary["rule"].str.contains("^annotate$")]
    
    # Count occurrences of statuses and pivot.
    summary = summary.groupby(["rule", "status"]).size().reset_index().rename(columns={0:"count"})
    summary = summary.pivot(index = "rule", columns = "status", values = "count").rename_axis(None, axis=1).reset_index()
    summary = summary.fillna(0)
    #print(summary)
    
    # Make sure that required columns are available after pivoting, which may not be the case when _all_ are missing or ok.
    try:
        summary["ok"]
        summary["missing"]
    except KeyError as e:
        # In case any key is missing, fill it with zeros. Could also have been a passive coalesce action with no ifs. But I guess python makes it beautiful even though it is a bit complex.
        for i in e.args:
            print(f"debug: filling key \"{i}\" with zeros.")
            summary[i] = 0
    
    # Calculate interpretative statistic
    summary["complete"] = summary["ok"] / (summary["ok"]+summary["missing"])
    
    # In case NaNs are introduced (division by zero) fill with zeros.
    summary = summary.fillna(0)
    #print(summary["complete"])
    summary["complete %"] = (summary["complete"] * 100).astype(int).astype(str) + " %"

    # Present
    print("Below is an overview of the completion of each rule.")
    print()
    print("  CompareM2 status")
    print("  ----------------")
    print(summary[["rule", "complete %"]].to_string(index = False))
    print("//")
    
    try:
        if delete_after_use:
            os.remove(path)
    except FileNotFoundError:
        pass # ignore
    

# For quick development.
if __name__ == "__main__":
    present_status()
