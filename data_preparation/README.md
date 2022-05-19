# CoDeS3D: Contextualising Developmental SNPs in Three Dimensions

## Database setup
CoDeS3D requires a PostgreSql database to run. This is because it integrates different types and sources of data, and using a database improves efficiency.
- Install [PostgreSql](https://www.postgresql.org) if you haven't already.
- If you plan to use more than a few Hi-C or eQTL datasets, you may want to create a tablespace for your database. To do that,
    1. Create a directory in location with adequate disk space
    ```
    mkdir /path/to/tablespace
    chown -R postgres /path/to/tablespace
    chrgrp -R postgres /path/to/tablespace
    ```

    2. Start psql as postgres and create tablespace
    ```
    psql -U postgres -h localhost
    CREATE TABLESPACE tblspace_codes3d LOCATION '</path/to/tablespace>';
    ```

    3. Confirm that a tablespace is created
    ```
    \db+
    ```
- Create the `codes3d` user,
```
CREATE USER codes3d WITH PASSWORD 'your password' CREATEDB;
CREATE DATABASE codes3d OWNER codes3d;
```

- Then create the `codes3d_commons` database,
```
CREATE DATABASE codes3d_commons OWNER codes3d tablespace tblspace_codes3d;
\q
```


## Hi-C data pre-processing and installation
- First process your raw Hi-C data with the [Juicer pipeline](https://github.com/aidenlab/juicer). Ensure that you use the GRCh38 reference human genome. The final outputs of the Hi-C processing should have the following columns: name of the paired-ends reads, strand (`str`), chromosome (`chr`), position (`pos`), restriction site fragment (`frag`), and mapping quality score (`score`) for both read pairs. (Read pairs with mapping quality score < 30 will be removed during installation.)

- Create a sub-directory for the processed Hi-C libraries in `codes3d/lib/hic/<restriction enzyme>`. After executing the subsequent steps, your Hi-C library directory should have a structure like the following:
```
codes3d/lib/hic/DpnII/
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.DpnII.fragments.bed
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.DpnII.fragments.db
├── GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
├── dna.fragments.db
└── hic_data
    ├── A549-ACE2_0h_Ho2020
    │   ├── A549_ACE2_0h_Ho2020_GSM4955378_merged_nodups.db
    │   └── A549_ACE2_0h_Ho2020_GSM4955379_merged_nodups.db
    |
    └── TeloHAEC_Lalonde2019
        ├── TeloHAEC_GSM3593256_merged_nodups.db
        └── TeloHAEC_GSM3593257_merged_nodups.db
```
- Digest the human genome with the restriction enzyme that was used to create your Hi-C library. See `python data_preparation/digest_genome.py -h` for help.

- Load the Hi-C libraries into sqlite3 or PostgreSql databases with the appropriate helper scripts: `data_preparation/hic_db_sqlite.py` for sqlite3 or `data_preparation/hic_db_postgres.py` for PostgreSql.

- Update the `lib/meta_info/meta_hic.txt` with information about the Hi-C cells. This file should have the following columns: `library`, `name`, `cell_type`, `tissue`, `disease` `sex`, `ethnicity`, `karyotype`, `age`, `enzyme`, and `tags`.

- Finally, load the Hi-C meta-info to the `codes3d_commons` PostgreSql database with `data_preparation/init_hic_meta.py`.


### eQTL data
- CoDeS3D requires genotype and expression datasets that are prepared using GTEx's [genotype](https://github.com/broadinstitute/gtex-pipeline/tree/master/genotype) and [RNA-seq](https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq) pipelines. In addition, follow steps 1, 2, and 3 of the GTEx's [QTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl) to produce your normalised expression and combined covariates files.

- Place all the required file in a sub-directory in `codes3d/lib/eqtls/`. Your file structure should look like this:
```
codes3d/lib/eqtls/GEUVADIS/
├── GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.lookup_table.txt.gz
├── GEUVADIS.445_samples.rnaseqc_tpm.gct.gz
├── GEUVADIS.445_samples.rnaseqc_tpm_median.gct.gz
├── GEUVADIS_covariates
│   └── GEUVADIS.445_samples.covariates.txt
├── GEUVADIS_expression_matrices
│   ├── GEUVADIS.445_samples.normalized_expression.bed.gz
│   └── GEUVADIS.445_samples.normalized_expression.bed.gz.tbi
└── GEUVADIS_genotypes
    ├── GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bed
    ├── GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.bim
    ├── GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.fam
    ├── GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.log
    ├── GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.vcf.gz
    └── GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.vcf.gz.tbi
```

- Update the `codes3d/lib/meta_info/meta_eqtls.txt` file to include the `name`,	`tissue`,	`project`, and `tags` of your eQTL project.

- Load the updated meta-info to the `codes3d_commons` database using the `data_preparation/init_eqtl_meta.py` helper script.

- Create a database for your eQTL project. The name should have the syntax `eqtls_<eQTL project name>`
```
CREATE DATABASE eqtls_<eQTL project name> owner codes3d tablespace <tablespace>;
```

- Create a `variants` table in your new eQTL database using
```
data_preparation/init_gene_variant_db.py
```

- Create `variant_lookup_<restriction fragment>` tables for the restriction enzymes that used to generate your Hi-C libraries
```
data_preparation/init_fragment_db.py
```

- Your eQTL database should have the following tables:

```
eqtls_geuvadis=> \d
                 List of relations
 Schema |          Name          | Type  |  Owner  
--------+------------------------+-------+---------
 public | variant_lookup_hindiii | table | codes3d
 public | variant_lookup_mboi    | table | codes3d
 public | variants               | table | codes3d
(3 rows)
```

- Finally, update `docs/codes3d.conf` to point CoDeS3D scripts to your eQTL datasets e.g.
```
.
.
.
[GEUVADIS]
eqtl_dir=../lib/eqtls/GEUVADIS/
GENE_FP: %(eqtl_dir)sGEUVADIS.445_samples.rnaseqc_tpm_median.gct.gz
GENOTYPES_DIR: %(eqtl_dir)sGEUVADIS_genotypes
GENOTYPES_FP: %(GENOTYPES_DIR)s/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.nodup.vcf.gz
COVARIATES_DIR: %(eqtl_dir)sGEUVADIS_covariates
EXPRESSION_DIR: %(eqtl_dir)sGEUVADIS_expression_matrices
.
.
```


## Reference tables
CoDeS3D needs genes and SNP reference files to run.
The gene reference file should be in `.BED` format: `chrom`, `start`, `end`,
`gene name`, `gene id (Ensembl)`.
Similarly, you will SNP reference files (e.g. human_9606_b151_GRCh38p7). Lastly,
CoDeS3D uses the [RsMergeArch.bcp.gz](https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/database/organism_data/RsMergeArch.bcp.gz) reference to resolve SNPs with multiple rsIDs. Your reference files should have the following structure.
```
codes3d/lib/reference_files/
├── genes
│   ├── gene_reference.bed
│   ├── gene_reference.db
│   └── gene_reference_sorted.txt
└── snps
    ├── human_9606_b151_GRCh38p7
    │   ├── bed_chr_1.bed
    │   ├── bed_chr_10.bed
    │   ├── bed_chr_11.bed
    │   ├── bed_chr_12.bed
    │   ├── bed_chr_13.bed
    │   ├── bed_chr_14.bed
    │   ├── bed_chr_15.bed
    │   ├── bed_chr_16.bed
    │   ├── bed_chr_17.bed
    │   ├── bed_chr_18.bed
    │   ├── bed_chr_19.bed
    │   ├── bed_chr_2.bed
    │   ├── bed_chr_20.bed
    │   ├── bed_chr_21.bed
    │   ├── bed_chr_22.bed
    │   ├── bed_chr_3.bed
    │   ├── bed_chr_4.bed
    │   ├── bed_chr_5.bed
    │   ├── bed_chr_6.bed
    │   ├── bed_chr_7.bed
    │   ├── bed_chr_8.bed
    │   ├── bed_chr_9.bed
    │   ├── bed_chr_MT.bed
    │   ├── bed_chr_X.bed
    │   └── bed_chr_Y.bed
    ├── human_9606_b151_GRCh38p7_RsMergeArch.bcp.gz
    └── human_9606_b151_GRCh38p7_RsMergeArch.pairs.bcp.gz
```

- Create a genes table in the `codes3d_commons` database. Make sure the `.BED`
file has gene coordinates in GrCh38.
```
data_preparation/init_gene_variant_db.py
```

- Create `gene_lookup_<restriction enzyme>` lookup tables in the `codes3d_commons`
database.
```
data_preparation/init_fragment_db.py
```

Your `codes3d_commons` database should look like the following:
```
codes3d_commons=> \d
               List of relations
 Schema |        Name         | Type  |  Owner  
--------+---------------------+-------+---------
 public | gene_lookup_hindiii | table | codes3d
 public | gene_lookup_mboi    | table | codes3d
 public | genes               | table | codes3d
 public | meta_eqtls          | table | codes3d
 public | meta_hic            | table | codes3d
(5 rows)
```
