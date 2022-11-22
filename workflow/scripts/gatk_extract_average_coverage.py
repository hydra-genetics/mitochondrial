

import pandas as pd
from math import floor

def extract_avg_cov(infile, mean_outfile, median_outfile):
    df = pd.read_csv(infile, skiprows=6, sep='\t', nrows=1)
    mean_df = df.get(['MEAN_COVERAGE']).apply(floor)
    mean_df.to_csv(mean_outfile, index=False, header=False)
    median_df = df.get(['MEDIAN_COVERAGE'])
    median_df.to_csv(median_outfile, index=False, header=False)
    

extract_avg_cov(snakemake.input.metrics, snakemake.output.mean, snakemake.output.median)

