#!/usr/bin/env python
from itertools import repeat
import multiprocessing
import istarmap
import tqdm
import time
import tempfile
import subprocess
import statsmodels.formula.api as smf
import statsmodels.api as sm
import scikits.bootstrap as boot
import pysam
import pandas as pd
import os
import numpy as np
import math
import gzip
import copy
import argparse
from itertools import repeat
import warnings
warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser()
    # REQUIRED
    parser.add_argument("--vcf", required=True, help="Genotype VCF")
    parser.add_argument("--pheno", required=True, help="Phenotype file")
    parser.add_argument(
        "--qtl", required=True,
        help="File containing QTL to calculate allelic fold " +
        "change for. Should contain tab separated columns 'pid' with " +
        "phenotype (gene) IDs and 'sid' with SNP IDs. Optionally can " +
        "include the columns 'sid_chr' and 'sid_pos', which will " +
        "facilitate tabix retrieval of genotypes, greatly reducing runtime.")
    parser.add_argument(
        "--geno", required=False, default="GT",
        help="Which field in VCF to use as the genotype. " +
        "By default 'GT' = genotype. Setting to 'DS' will use dosage rounded " +
        "to the nearest integer (IE 1.75 = 2 = 1|1).")
    parser.add_argument("--chr", type=str,
                        help="Limit to a specific chromosome.")
    parser.add_argument("--log_xform", type=int, required=True,
                        help="The data has been log transformed (1/0). If so, please set --log_base.")
    parser.add_argument("--output", "--o", required=True, help="Output file")

    # OPTIONAL
    parser.add_argument("--cov", help="Covariates file")
    parser.add_argument(
        "--matrix_o", help="Output the raw data matrix used to calculate aFC " +
        "for each QTL into the specific folder.")
    parser.add_argument(
        "--boot", default=100, type=int,
        help="Number of bootstraps to perform for effect size confidence interval. " +
        "Can be set to 0 to skip confidence interval calculation, which will greatly reduce runtimes.")
    parser.add_argument("--ecap", default=math.log(100, 2),
                        type=float, help="Absolute aFC cap in log2.")
    parser.add_argument(
        "--log_base", default=2, type=int,
        help="Base of log applied to data. If other than 2, data will be converted to log2.")
    parser.add_argument(
        "--min_samps", default=2, type=int,
        help="Minimum number of samples with genotype data required to calculate effect size, default = 2.")
    parser.add_argument(
        "--min_alleles", default=1, type=int,
        help="Minimum observations of each allele in data to calculate aFC, default = 1.")
    parser.add_argument("--count_o", default=0, type=int,
                        help="Output the observed allele counts, default = 0.")
    return parser.parse_args()


def main(
        eqtl_df,
        vcf_,
        expression_dir_,
        covariates_dir_,
        eqtl_project,
        output_dir_,
        boot_=0,
        num_processes=2,
        log_xform_=1,
        log_base_=2):
    warnings.simplefilter('ignore')
    if not 'pid' in eqtl_df.columns:
        eqtl_df['pid'] = eqtl_df['gencode_id']
    global log_xform, log_base, ecap, boot_val, geno, chrom, matrix_o, min_samps, min_alleles, count_o
    log_xform = log_xform_
    log_base = log_base_
    boot_val = boot_
    geno = "GT"
    chrom = None
    matrix_o = None
    ecap = math.log(100, 2)
    min_samps = 2
    min_alleles = 1
    count_o = 0

    global vcf, covariates_dir, expression_dir, output_dir
    vcf = vcf_
    covariates_dir = covariates_dir_
    expression_dir = expression_dir_
    output_dir = output_dir_

    start_time = time.time()

    version = "0.3"
    print("")
    print("########################################################")
    print("		 Welcome to aFC v%s" % (version))
    print("  Authors: Pejman Mohammadi (pmohammadi@nygenome.org),\n	   Stephane Castel (scastel@nygenome.org)")
    print("########################################################")

    print("RUN SETTINGS")
    #print("     Genotype VCF: %s" % (vcf))
    #print("     Phenotype File: %s" % (pheno))
    # if cov != None:
    #    print("     Covariate File: %s" % (cov))
    #print("     QTL File: %s" % (qtl))
    print("     Genotype Field: %s" % (geno))
    print("     Log Transformed: %d" % (log_xform))
    if log_xform == 1:
        print("     Log Base: %d" % (log_base))
    print("     Boostraps: %d" % (boot_val))
    print("     Minimum number of samples with genotype data: %d" %
          (min_samps))
    print("     Minimum number of allele observations: %d" % (min_alleles))
    if chrom != None:
        print("     Chromosome: %s" % (chrom))

    print("")

    tissues = eqtl_df.tissue.drop_duplicates().tolist()
    tissue_dfs = [eqtl_df[eqtl_df.tissue == tissue] for tissue in tissues]
    '''
    for i in range(len(tissue_dfs)):
        calculate_tissue_aFC(tissues[i], tissue_dfs[i], eqtl_project)
    '''
    with multiprocessing.Pool(num_processes) as pool:
        desc = '  * Calculating eQTL effect size (aFC)'
        bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'
        for _ in tqdm.tqdm(
                pool.istarmap(calculate_tissue_aFC,
                              zip(tissues,
                                  tissue_dfs,
                                  repeat(eqtl_project))),
                total=len(tissue_dfs),
                desc=desc, unit='tissues', ncols=80,
                bar_format=bar_format):
            pass
    df = []
    for fp in os.listdir(output_dir):
        if fp.endswith('aFC.txt'):
            df.append(pd.read_csv(os.path.join(output_dir, fp), sep='\t'))
            os.remove(os.path.join(output_dir, fp))
    del tissue_dfs
    df = pd.concat(df)
    '''
    if 'adj_pval' in df.columns:
        df = df[['sid', 'gencode_id', 'tissue', 'pval', 'adj_pval', 'log2_aFC',
                 'log2_aFC_lower', 'log2_aFC_upper', 'maf', 'sid_chr',
                 'sid_pos']]
        #df['pid'] = df['pid'].str.split('.', expand=True)[0]
        df.to_csv(os.path.join(output_dir, 'significant_eqtls.txt'),
                  sep='\t', index=False)
    '''
    duration = (time.time() - start_time) / 60
    print("  Time elasped: {:.2f} mins".format(duration))
    print("  Done.")
    print("")
    return df


def calculate_tissue_aFC(tissue, qtl, eqtl_project):

    pheno = ''
    cov = ''
    if eqtl_project.lower() == 'gtex': # TODO: rename files for consistency
        pheno = os.path.join(expression_dir, tissue +
                             '.v8.normalized_expression.bed.gz')
        cov = os.path.join(covariates_dir, tissue + '.v8.covariates.txt')
    else:
        pheno = os.path.join(expression_dir, tissue +
                             '.normalized_expression.bed.gz')
        cov = os.path.join(covariates_dir, tissue + '.covariates.txt')
    output = os.path.join(output_dir, tissue + '.aFC.txt')

    # get sample - column map from VCF
    #print("1. Loading VCF...")
    vcf_map = sample_column_map(vcf)
    tabix_vcf = pysam.Tabixfile(vcf, "r")

    global df_cov
    if cov != None:
     #   print("1b. Loading covariates...")
        df_cov = pd.read_csv(cov, sep="\t", index_col=False)
        if "ID" in df_cov.columns:
            cov_id_col = "ID"
        elif "id" in df_cov.columns:
            cov_id_col = "id"
        else:
            print("Could not find covariate ID column in covariates column. Please ensure that it is either labeled 'ID' or 'id'")
            quit()
    else:
        df_cov = pd.DataFrame(columns=['ID'])
        cov_id_col = "ID"
    df_cov[cov_id_col] = df_cov[cov_id_col].str.split('_').str[0]
    # 2 get sample - column map from phenotype file
    #print("2. Loading phenotype data...")
    pheno_map = sample_column_map(pheno, line_key="#", start_col=4)
    tabix_pheno = pysam.Tabixfile(pheno, "r")

    # 3 load fastQTL results
    #print("3. Loading fastQTL results...")
    df_qtl = qtl
    #df_qtl = pd.read_csv(args.qtl, sep="\t", index_col=False)

    #print("4. Retrieving eSNP positions...")
    set_esnp = set(df_qtl['sid'].tolist())
    dict_esnp = {}

    if "sid_chr" in df_qtl.columns and "sid_pos" in df_qtl.columns:
        # eSNP positions are specified in file
        for index, row in df_qtl.iterrows():
            if chrom == None or str(row['sid_chr']) == chrom:
                dict_esnp[row['sid']] = [row['sid_chr'], int(row['sid_pos'])]
    else:
        # retrieve SNP positions from the VCF (since these are not included in the fastQTL output)
        #print("     unpacking VCF...")
        # retrieve the SNP positions from the VCF
        tfile = tempfile.NamedTemporaryFile(delete=False)
        vcf_in = tfile.name
        tfile.close()
        if chrom != None:
            # retrieve only genotypes from the desired chromosome
            error = subprocess.call(
                "tabix "+vcf+" "+chrom+": | cut -f 1-3 > "+vcf_in, shell=True)
            if error != 0:
                print(
                    "     ERROR loading retrieving genotype data. Ensure tabix index exists and is current.")
                quit()
        else:
            error = subprocess.call(
                "gunzip -c "+vcf+" | cut -f 1-3 > "+vcf_in, shell=True)
            if error != 0:
                print("     ERROR loading retrieving genotype data.")
                quit()

        stream_in = open(vcf_in, "r")

        current_chr = ""

        for line in stream_in:
            if line[0:1] != "#":
                # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
                columns = line.rstrip().split("\t")

                if columns[0] != current_chr:
                    print("     chr: %s" % (columns[0]))
                    current_chr = columns[0]

                if columns[2] in set_esnp:
                    dict_esnp[columns[2]] = [columns[0], int(columns[1])]

        stream_in.close()

        os.remove(vcf_in)

    # determine how many total eQTL there are
    total_eqtl = 0
    for esnp in df_qtl['sid'].tolist():
        if esnp in dict_esnp:
            total_eqtl += 1

    # 5 retrieve phenotype positions
    #print("5. Retrieving ePhenotype positions...")
    set_epheno = set(df_qtl['pid'].tolist())

    stream_in = gzip.open(pheno, "r")

    dict_ephenotype = {}
    for line in stream_in:
        if isinstance(line, bytes) and not isinstance(line, str):
            line = line.decode()
        if line[0:1] != "#":
            columns = line.rstrip().split("\t")
            # Chr    start   end     ID
            if columns[3] in set_epheno:
                dict_ephenotype[columns[3]] = [columns[0], int(columns[1])]

    stream_in.close()

    # 6 calculate effect sizes
    #print("6. Calculating eQTL effect sizes...")
    stream_vcf = open

    completed = 0
    t = time.time()

    stream_out = open(output, "w")
    headers = ['log2_aFC', 'log2_aFC_lower', 'log2_aFC_upper']
    if count_o == 1:
        headers += ['ref_allele_count', 'alt_allele_count']
    stream_out.write("\t".join(df_qtl.columns.tolist()+headers)+"\n")

    for index, row in df_qtl.iterrows():
        # now retrieve the genotypes for the snp
        # only for those individuals with phenotype data
        dict_geno = {}
        line_written = False
        allele_counts = None

        if row['sid'] in dict_esnp:
            if row['pid'] in dict_ephenotype:

                esnp = dict_esnp[row['sid']]

                records = tabix_vcf.fetch(
                    esnp[0], esnp[1]-1, esnp[1])

                snp_found = 0
                for record in records:
                    cols = record.rstrip().split("\t")
                    if cols[2] == row['sid']:
                        gt_index = cols[8].split(":").index(geno)
                        snp_found = 1
                        for sample in pheno_map.keys():
                            sample_col = cols[vcf_map[sample]]
                            dict_geno[sample] = sample_col.split(":")[gt_index]
                if snp_found == 0:
                    # print("	  WARNING: eSNP %s not found in VCF" %
                    #      (row['sid']))
                    continue

                # assume phenotype is within a megabase of SNP
                ephenotype = dict_ephenotype[row['pid']]

                records = tabix_pheno.fetch(
                    ephenotype[0], ephenotype[1]-1, ephenotype[1]+1)

                dict_pheno = {}

                for record in records:
                    cols = record.rstrip().split("\t")
                    if cols[3] == row['pid']:
                        for sample in dict_geno.keys():
                            if log_xform == 1 and log_base != 2:
                                # if data has been log transformed but is not in base 2 convert it
                                dict_pheno[sample] = float(
                                    cols[pheno_map[sample]]) * math.log(log_base, 2)
                            else:
                                dict_pheno[sample] = float(
                                    cols[pheno_map[sample]])

                # make a dataframe with all covariates and genotype classes
                list_rows = []
                allele_counts = [0, 0]

                for sample in dict_geno.keys():
                    if geno == "GT":
                        # only include samples w/ complete genotype data (no '.')
                        if "." not in dict_geno[sample]:
                            list_rows.append([dict_geno[sample].count(
                                "1"), dict_pheno[sample]] + return_cov(sample))
                            allele_counts[0] += dict_geno[sample].count("0")
                            allele_counts[1] += dict_geno[sample].count("1")
                    elif geno == "DS":
                        list_rows.append(
                            [round(float(dict_geno[sample])), dict_pheno[sample]] + return_cov(sample))

                # Changed to only run effect size calc when more than minimum # samps w/ GT data and minimum number of observations for each allele
                if len(list_rows) >= min_samps and min(allele_counts) >= min_alleles:
                    df_test = pd.DataFrame(list_rows, columns=[
                                           'geno', 'pheno']+["cov_"+x for x in df_cov[cov_id_col].tolist()])

                    if matrix_o != None:
                        df_test.to_csv(
                            matrix_o+"/"+row['pid']+":"+row['sid']+".txt", sep="\t", index=False)

                    # correct for covariates
                    df_test = correct_covariates(df_test)

                    esize = effect_size(df_test)
                    line_out = row.tolist()+esize[0:3]
                    if count_o == 1:
                        line_out += allele_counts

                    stream_out.write("\t".join(map(str, line_out))+"\n")
                    line_written = True
                else:
                    if len(list_rows) == 0:
                        print("	  WARNING: no individual with genotype data for eQTL %s - %s" %
                              (row['pid'], row['sid']))
                    elif len(list_rows) < min_samps:
                        print("	  WARNING: only %d individual(s) with genotype data for eQTL %s - %s" %
                              (len(list_rows), row['pid'], row['sid']))
                    elif min(allele_counts) < min_alleles:
                        print("	  WARNING: only %d observations of minor allele for eQTL %s - %s" %
                              (min(allele_counts), row['pid'], row['sid']))

            completed += 1
            '''
            if completed % 100 == 0:
                print("     COMPLETED %d of %d = %f in %d seconds" % (
                    completed, total_eqtl,
                    float(completed)/float(total_eqtl),
                    time.time() - t))
                t = time.time()
            '''
        else:
            if row['pid'] != "nan" and chrom == None:
                print(
                    "	  WARNING: positional information not found for eSNP %s" % (row['sid']))

        # ensure that every eQTL has a line written so number of output lines = number of input eQTLs
        if ("sid_chr" in df_qtl.columns and (str(row['sid_chr']) == chrom or chrom == None)):
            if line_written == False:
                line_out = row.tolist()+[float('nan'),
                                         float('nan'), float('nan')]
                if count_o == 1:
                    if allele_counts != None:
                        line_out += allele_counts
                    else:
                        line_out += [float('nan'), float('nan')]
                stream_out.write("\t".join(map(str, line_out))+"\n")

    stream_out.close()


def return_cov(sample):
    global df_cov

    if sample in df_cov.columns:
        return(df_cov[sample].tolist())
    else:
        return([])


def sample_column_map(path, start_col=9, line_key="#CHR"):
    stream_in = gzip.open(path, "r")

    out_map = {}
    for line in stream_in:
        if isinstance(line, bytes) and not isinstance(line, str):
            line = line.decode()
        if line_key in line:
            line = line.rstrip().split("\t")
            for i in range(start_col, len(line)):
                out_map[line[i]] = i

            break

    stream_in.close()

    return(out_map)


def correct_covariates(df_test):
    global df_cov

    if len(df_cov.index) > 0:
        # correct for covariates
        # add genotype categorical covariates
        cov_homo_ref = [int(x == 0) for x in df_test['geno']]
        if sum(cov_homo_ref) > 0:
            df_test['cov_homo_ref'] = cov_homo_ref

        cov_homo_alt = [int(x == 2) for x in df_test['geno']]
        if sum(cov_homo_alt) > 0:
            df_test['cov_homo_alt'] = cov_homo_alt

        cov_ids = [x for x in df_test.columns if "cov_" in x]

        # convert categorical covariates to n-1 binary covariates
        new_cols = {}
        drop_cols = []

        for xcov in cov_ids:
            if df_test.dtypes[xcov] == object:
                values = list(set(df_test[xcov]))[1:]
                for xval in values:
                    xname = xcov+"_"+xval
                    new_cols[xname] = [int(x == xval) for x in df_test[xcov]]

                drop_cols.append(xcov)

        df_test.drop(drop_cols, axis=1, inplace=True)
        for xcov in new_cols.keys():
            df_test[xcov] = new_cols[xcov]
        cov_ids = [x for x in df_test.columns if "cov_" in x]

        # NOTE any variable that is a string will be treated as categorical - this is the same functionality as FASTQTL, so good
        # see: http://statsmodels.sourceforge.net/devel/example_formulas.html

        xformula = "pheno ~ "+"+".join(cov_ids)
        result = smf.ols(formula=xformula, data=df_test).fit()

        # use only significant (95% CI doesn't overlap 0) covariates to correct expression values
        # do not include intercept or genotypes in correction

        drop_covs = []
        for xcov in list(result.params.index):
            if xcov in df_test.columns:
                coefficient = result.params.loc[xcov]
                upper_ci = result.conf_int(0.05).loc[xcov][1]
                lower_ci = result.conf_int(0.05).loc[xcov][0]
                if (lower_ci <= 0 and upper_ci >= 0):
                    drop_covs.append(xcov)

        # drop insignificant covariates
        df_test.drop(drop_covs, axis=1, inplace=True)
        cov_ids = [x for x in df_test.columns if "cov_" in x]

        # redo regression without insignificant covs
        if len(cov_ids) > 0:
            xformula = "pheno ~ "+"+".join(cov_ids)
            result = smf.ols(formula=xformula, data=df_test).fit()

            df_test_corrected = copy.deepcopy(df_test)
            for xcov in list(result.params.index):
                coefficient = result.params.loc[xcov]
                if xcov == "Intercept" or xcov == "cov_homo_ref" or xcov == "cov_homo_alt":
                    df_test_corrected[xcov] = [0] * \
                        len(df_test_corrected.index)
                else:
                    df_test_corrected[xcov] = [
                        x * coefficient for x in df_test_corrected[xcov]]

            # add residual to dataframe
            df_test_corrected['pheno_cor'] = [
                row['pheno'] - sum(row[2:len(row)]) for index, row in df_test_corrected.iterrows()]

        else:
            # if none of the covariates are significant then just leave the values as is
            df_test_corrected = copy.deepcopy(df_test)
            df_test_corrected['pheno_cor'] = df_test_corrected['pheno']
    else:
        # covariates not provided
        df_test_corrected = copy.deepcopy(df_test)
        df_test_corrected['pheno_cor'] = df_test_corrected['pheno']

    return(df_test_corrected)


def effect_size(df_test):
    global args
    # calculate effect size
    esize = calculate_effect_size(df_test['geno'], df_test['pheno_cor'])

    # calculate 95% CI for effect size using BCa bootstrapping
    if boot_val > 0:
        try:
            ci = boot.ci((df_test['geno'].tolist(), df_test['pheno_cor'].tolist(
            )), statfunction=calculate_effect_size, alpha=0.05, n_samples=boot_val, method="bca")
        except (IndexError, ValueError):  # ValueError added for calculating CI on one sample
            ci = [float('nan'), float('nan')]
    else:
        ci = [float('nan'), float('nan')]

    return([esize, ci[0], ci[1]])


def calculate_effect_size(genos, phenos):
    global args

    # in cases where there is only a single genotype in the data return nan
    if len(set(genos)) == 1:
        return(float('nan'))

    if log_xform == 1:
        # M5 - for log2 transformed data

        # 1 need to prepare 4 estimates
        p_m = [
            np.mean(phenos[genos == 0]),
            np.mean(phenos[genos == 1]),
            np.mean(phenos[genos == 2])
        ]

        log2ratio_M2M0 = bound_basic(p_m[2] - p_m[0], -ecap, ecap)
        log2ratio_M1M2 = bound_basic(p_m[1] - p_m[2], -0.9999999, ecap)
        log2ratio_M1M0 = bound_basic(p_m[1] - p_m[0], -1, ecap)

        p_delta = [
            float('nan'),
            math.pow(2, log2ratio_M2M0),
            float(1) / (math.pow(2, log2ratio_M1M2+1) - 1),
            math.pow(2, log2ratio_M1M0+1) - 1,
            None,
        ]

        X = sm.add_constant(genos)
        result = sm.OLS(phenos, X).fit()
        result_coef = bound_basic(result.params[1]*2, -ecap, ecap)
        p_delta[4] = math.pow(2, result_coef)

        for x in range(1, 5):
            p_delta[x] = bound_basic(p_delta[x], math.pow(
                2, -ecap), math.pow(2, ecap))

        stdevs = {}

        # pick the estimate that minimizes residual variance
        for i in range(1, 5):
            #stdevs[i] = numpy.std([yi - calculate_expected_expr(p_delta[i], xi) for xi, yi in zip(genos, phenos)])
            stdevs[i] = np.std(
                phenos - np.log2((2 - genos) + (p_delta[i] * genos)))

        min_delta = min([x for x in stdevs.values() if math.isnan(x) == False])
        use_delta = 0
        for delta in range(1, 5):
            if stdevs[delta] == min_delta:
                use_delta = delta
                break
        p_delta[0] = float('nan')

        return(math.log(p_delta[use_delta], 2))
    else:
        # linear regression on untransformed data
        X = sm.add_constant(genos)
        result = sm.OLS(phenos, X).fit()
        # ensure intercept is positive
        b0 = bound_basic(result.params[0], np.finfo(float).eps, float('inf'))
        # calculate the effect size
        use_delta = (float(2 * result.params[1]) / float(b0)) + 1

        # bound delta between caps
        if use_delta < 0:
            if result.params[1] > 0:
                use_delta = math.pow(2, ecap)
            else:
                use_delta = math.pow(2, -ecap)

        # bound effect size between -ecap and ecap in log space
        use_delta_log = math.log(use_delta, 2)
        use_delta_log_bounded = bound_basic(
            use_delta_log, -ecap, ecap)

        return(use_delta_log_bounded)


def bound_basic(x, l, h):
    y = min([x, h])
    y = max([y, l])
    if math.isnan(x) == True:
        y = float('nan')
    return(y)


def calculate_expected_expr(delta, alt_alleles):
    if math.isnan(delta) == False:
        return(math.log((2 - alt_alleles) + (delta * alt_alleles), 2))
    else:
        return(float('nan'))


if __name__ == "__main__":
    global args

    args = parse_args()

    # main(args.qtl, args.vcf, args.phen, ags.cov, args.boot,
    #    args.log_xform, args.log_base, args.output)
