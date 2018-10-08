#!/usr/bin/env python

""" Collection of useful command line scripts """
import fire
import pandas as pd
import pandas_profiling

def describe(file):
    df = pd.read_csv(file, header=1, sep='\s', na_values=["-9"])
    # pandas_profiling.ProfileReport(df)
    print(df.describe().transpose())

if __name__ == '__main__':
  fire.Fire()
