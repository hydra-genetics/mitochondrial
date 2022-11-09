__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from pysam import AlignmentFile

min_version("6.8.0")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")


validate(config, schema="../schemas/config.schema.yaml")
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

samples = pd.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

### Read and validate units file

units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

### Set wildcard constraints


wildcard_constraints:
    sample="|".join(samples.index),
    type="N|T|R",


def get_pg_info(bam):

    bamfile = AlignmentFile(bam, 'rb')
    pg_header = bamfile.header['PG']

    for d in pg_header:
        if d['ID'] == 'bwa':
            version = d['VN']
            bwa_cmd_line = d['CL']

    return(version, bwa_cmd_line)


def compile_output_list(wildcards):

    files = {
        "mitochondrial/gatk_sort_sam": [
            "bam"
        ],
        "mitochondrial/gatk_collect_wgs_metrics": [
            "metrics.txt", "theoretical_sensitivity.txt", "mean_coverage.txt", "median_coverage.txt",
        ],
        "mitochondrial/gatk_mutect2": [
            "vcf"
        ],

    }

    output_files = [
        "%s/%s_%s_%s.%s" % (prefix, sample, unit_type, ref, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for ref in ['mt', 'mt_shifted']
        for suffix in files[prefix]
    ]

    
    return output_files
