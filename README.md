# CoDeS3D: Contextualising Developmental SNPs in Three Dimensions

## Installation

Clone the repository using the following command:

```
https://github.com/Genome3d/codes3d.git
cd codes3d/
```

## Setup
(The following guide has been tested on Ubuntu 14.04+.)

Install required apt packages:
```
sudo apt update
sudo apt install -y libz-dev libxslt1-dev libpq-dev python3-pip python3-dev python3-virtualenv bedtools
pip3 install --upgrade pip virtualenv setuptools
```


#### Create environment with conda
Ensure that you have [conda](https://docs.conda.io/en/latest/) installed. Then create a conda environment from the root directory. In this example, the environment is created in the envs directory and has all the required dependencies. (Size of environment ~1.8GB.)
```
conda env create --prefix ./codes3d_env --file environment.yaml
```

Then activate the environment with
```
conda activate codes3d_env/
```
To deactivate,
```
conda deactivate
```


#### Data installation
You will also need to install Hi-C libraries and data necessary to calculate eQTLs. You can find detailed instructions and helper scripts in the `download_data` directory's [README](download_data/README.md) file.

## Basic Usage

The CoDeS3D interface is heavily inspired by the Qiime interface (J Gregory Caporaso *et al*., Nature Methods, 2010; doi:10.1038/nmeth.f.303). Running the CoDeS3D script in the codes3d directory will drop the user into the CoDeS3D shell, in which all CoDeS3D scripts are accessible from anywhere in the system, e.g.

```
(/mnt/projects/codes3d/codes3d_env) /m/p/codes3d$ python codes3d/codes3d.py -h
usage: codes3d.py [-h] [-s SNP_INPUT [SNP_INPUT ...]] [-g GENE_INPUT [GENE_INPUT ...]] [--snps-within-gene SNPS_WITHIN_GENE [SNPS_WITHIN_GENE ...]] [-o OUTPUT_DIR] [--multi-test MULTI_TEST] [--pval-threshold PVAL_THRESHOLD]
                  [-f FDR_THRESHOLD] [--maf-threshold MAF_THRESHOLD] [-p NUM_PROCESSES] [--no-afc] [--afc-bootstrap AFC_BOOTSTRAP] [-n INCLUDE_CELL_LINES [INCLUDE_CELL_LINES ...]] [-x EXCLUDE_CELL_LINES [EXCLUDE_CELL_LINES ...]]
                  [--list-hic-libraries] [--match-tissues ...] [--list-tissue-tags] [-t TISSUES [TISSUES ...]] [--eqtl-project EQTL_PROJECT] [--list-eqtl-db] [-r RESTRICTION_ENZYMES [RESTRICTION_ENZYMES ...]] [--list-enzymes]
                  [--list-eqtl-tissues] [-c CONFIG] [--output-format OUTPUT_FORMAT] [--do-not-produce-summary] [--suppress-intermediate-files] [--non-spatial] [--gene-list GENE_LIST [GENE_LIST ...]] [--gtex-cis]

CoDeS3D maps gene regulatory landscape using chromatin
        interaction and eQTL data.

optional arguments:
  -h, --help            show this help message and exit
  -s SNP_INPUT [SNP_INPUT ...], --snp-input SNP_INPUT [SNP_INPUT ...]
                        The dbSNP IDs or loci of SNPs of interest in the format 'chr<x>:<locus>'. Use this flag to identify eGenes associated with the SNPs of interest.
  -g GENE_INPUT [GENE_INPUT ...], --gene-input GENE_INPUT [GENE_INPUT ...]
                        The symbols, Ensembl IDs or loci of genes interest in the format 'chr<x>:<start>-<end>'. Use this flag to identify eQTLs associated with the gene of interest.
  --snps-within-gene SNPS_WITHIN_GENE [SNPS_WITHIN_GENE ...]
                        A gene symbol, Ensembl ID or location in the format 'chr<x>:<start>-<end>'. Use this flag to identify eGenes associated with the SNPs located within the gene of interest.
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        The directory in which to output results.
  --multi-test MULTI_TEST
                        Options for BH multiple-testing: ['snp', 'tissue', 'multi']. 'snp': corrects for genes associated with a given SNP in a given tissue. 'tissue': corrects for all associations in a given tissue. 'multi': corrects
                        for all associations across all tissues tested.
  --pval-threshold PVAL_THRESHOLD
                        Maximum p value for mapping eQTLs. Default: 1.
  -f FDR_THRESHOLD, --fdr-threshold FDR_THRESHOLD
                        The FDR threshold to consider an eQTL statistically significant (default: 0.05).
  --maf-threshold MAF_THRESHOLD
                        Minimum MAF for variants to include. Default: 0.1.
  -p NUM_PROCESSES, --num-processes NUM_PROCESSES
                        Number of CPUs to use (default: half the number of CPUs).
  --no-afc              Do not calculate allelic fold change (aFC). If true, eQTL beta (normalised effect size) is calculated instead. (default: False).
  --afc-bootstrap AFC_BOOTSTRAP
                        Number of bootstrap for aFC calculation (default: 1000).
  -n INCLUDE_CELL_LINES [INCLUDE_CELL_LINES ...], --include-cell-lines INCLUDE_CELL_LINES [INCLUDE_CELL_LINES ...]
                        Space-separated list of cell lines to include (others will be ignored). NOTE: Mutually exclusive with EXCLUDE_CELL_LINES. --match-tissues takes precedence.
  -x EXCLUDE_CELL_LINES [EXCLUDE_CELL_LINES ...], --exclude-cell-lines EXCLUDE_CELL_LINES [EXCLUDE_CELL_LINES ...]
                        Space-separated list of cell lines to exclude (others will be included). NOTE: Mutually exclusive with INCLUDE_CELL_LINES. --match-tissues and -n take precedence.
  --list-hic-libraries  List available Hi-C libraries.
  --match-tissues ...   Try to match eQTL and Hi-C tissue types using space-separated tags. When using this, make sure that it is the last tag. Note that tags are combined with the AND logic. Prepend "-" to a tag if you want it to be
                        excluded. Use `--list-tissues-tags' for possible tags.
  --list-tissue-tags    List tags to be used with `--match-tissues'.
  -t TISSUES [TISSUES ...], --tissues TISSUES [TISSUES ...]
                        Space-separated list of eQTL tissues to query. Note that tissues are case-sensitive and must be from the same eQTL projects. Default is all tissues from the GTEx project. Use 'codes3d.py --list-eqtl-tissues'
                        for a list of installed tissues.
  --eqtl-project EQTL_PROJECT
                        The eQTL project to query. Default: GTEx. 'use 'codes3d.py --list-eqtl-db' to list available databases.
  --list-eqtl-db        List available eQTL projects to query.
  -r RESTRICTION_ENZYMES [RESTRICTION_ENZYMES ...], --restriction-enzymes RESTRICTION_ENZYMES [RESTRICTION_ENZYMES ...]
                        Space-separated list of restriction enzymes used in Hi-C data.' Use 'list-enzymes' to see available enzymes.
  --list-enzymes        List restriction enzymes used to prepare installed Hi-C libraries.
  --list-eqtl-tissues   List available eQTL tissues to query.
  -c CONFIG, --config CONFIG
                        The configuration file to use (default: docs/codes3d.conf).
  --output-format OUTPUT_FORMAT
                        Determines columns of output. 'short': snp|gencode_id|gene|tissue|adj_pval| [log2_aFC|log2_aFC_lower|log2_aFC_upper] or [beta|beta_se]| maf|interaction_type|hic_score. 'medium' adds:
                        eqtl_pval|snp_chr|snp_locus|ref|alt|gene_chr| gene_start|gene_end|distance. 'long' adds:cell_lines|cell_line_hic_scores| expression|max_expressed_tissue|max_expression| min_expressed_tissue|min_expression.
                        (default: long).
  --do-not-produce-summary
                        Do not produce final summary file, stop process after mapping eQTLs. Helpful if running batches (default: False).
  --suppress-intermediate-files
                        Do not produce intermediate files. These can be used to run the pipeline from an intermediate stage in the event of interruption (default: False).
  --non-spatial         Map non-spatial eQTLs.
  --gene-list GENE_LIST [GENE_LIST ...]
                        List of genes for non-spatial eQTL mapping.
  --gtex-cis            Retrieve spatially unconstrained cis-eQTLs as calculated in GTEx. To be used in combination with 'non-spatial'.
```

If you want to make CoDeS3D accessible from anywhere in your system, you can add the codes3d directory to your PATH by appending this line to the bottom of your `.bash_profile`/`.bashrc` file in your home directory:

```
export PATH=$PATH:/path/to/CoDeS3D # Replace "/path/to" with the path where the CoDeS3D program is stored
```

Alternatively, you can make a link to the program in a directory which is already in your path (e.g. `/usr/local/bin`) using the following commands:

```
cd /usr/local/bin
sudo ln -s /path/to/CoDeS3D CoDeS3D # Replace "/path/to" with the path where the CoDeS3D program is stored
```

For details on the expected inputs of any CoDeS3D script, simply run the script with the `-h` argument, as in the example above.
