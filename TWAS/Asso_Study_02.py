#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import multiprocessing
import subprocess
import sys

from time import time

import numpy as np
import pandas as pd

from scipy.stats import chi2
from joblib import Parallel, delayed
from google.cloud import storage
import threading


def get_pval(z):
    return np.format_float_scientific(chi2.sf(z**2, 1), precision=15, exp_digits=0)


def get_V_cor(V_cov):
    V_cov = V_cov.copy()
    v = np.sqrt(np.diag(V_cov))
    outer_v = np.outer(v, v)
    V_cor = V_cov / outer_v
    V_cor[V_cov == 0] = 0
    return V_cor


def get_z_denom(V, w):
    return np.sqrt(np.linalg.multi_dot([w, V, w]))


def get_spred_zscore(V_cov, w, Z_gwas, snp_sd):
    Z_twas = snp_sd.dot(w * Z_gwas) / get_z_denom(V_cov, w)
    return Z_twas, get_pval(Z_twas)


def get_fusion_zscore(V_cov, w, Z_gwas, snp_sd=None):
    V_cor = get_V_cor(V_cov)
    Z_twas = np.vdot(Z_gwas, w) / get_z_denom(V_cor, w)
    return Z_twas, get_pval(Z_twas)


def get_burden_zscore(test_stat, get_zscore_args):
    if test_stat == "FUSION":
        return get_fusion_zscore(*get_zscore_args)
    if test_stat == "SPrediXcan":
        return get_spred_zscore(*get_zscore_args)


def write_to_file(stuff):
    storage_client = storage.Client("rome-repeat-genetics")
    bucket = storage_client.bucket("rome-repeat-genetics")
    blob = bucket.blob("output_logs.txt")
    with blob.open("w") as f:
        f.write(f"{stuff}\n")


def main():
    ###############################################################
    # time calculation
    start_time = time()

    ###############################################################
    # parse input arguments
    parser = argparse.ArgumentParser(description="Asso Study 02")
    parser.add_argument("--TIGAR_dir", type=str)
    parser.add_argument("--gene_anno", type=str, dest="annot_path")
    parser.add_argument("--chr", type=str, dest="chrm")
    parser.add_argument("--weight", type=str, dest="w_path")
    parser.add_argument("--Zscore", type=str, dest="z_path")
    parser.add_argument("--LD", type=str, dest="ld_path")
    parser.add_argument("--window", type=float)
    parser.add_argument("--weight_threshold", type=float, help="Weight threshold to include SNP in TWAS")
    parser.add_argument("--test_stat", type=str, help="specify 'FUSION', 'SPrediXcan', or 'both': Zscore test statistic to use")
    parser.add_argument("--thread", type=int)
    parser.add_argument("--out_dir", type=str)
    parser.add_argument("--out_twas_file", type=str)

    args = parser.parse_args()

    sys.path.append(args.TIGAR_dir)
    import TIGARutils as tg

    out_twas_path = args.out_dir + "/" + args.out_twas_file
    write_to_file(
        """********************************
        Input Arguments
        Gene annotation file specifying genes for TWAS: {annot_path}
        Chromosome: {chrm}
        cis-eQTL weight file: {w_path}
        GWAS summary statistics Z-score file: {z_path}
        Reference LD genotype covariance file: {ld_path}
        Gene training region SNP inclusion window: +-{window}
        SNP weight inclusion threshold: {weight_threshold}
        Test statistic to use: {test_stat_str}
        Number of threads: {thread}
        Output directory: {out_dir}
        Output TWAS results file: {out_path}
        ********************************""".format(
            **args.__dict__,
            test_stat_str="FUSION and SPrediXcan"
            if args.test_stat == "both"
            else args.test_stat,
            out_path=out_twas_path
        )
    )

    tg.print_args(args)


    write_to_file("Reading gene annotation file.")
    # Gene, TargetID, n_targets = tg.read_gene_annot_exp(args.annot_path, args.chrm)
    Gene, TargetID, n_targets = tg.read_gene_annot_exp(**args.__dict__)

    # read in headers for Weight and Zscore files; get the indices and dtypes for reading files into pandas
    write_to_file("Reading file headers.\n")
    weight_info = tg.weight_file_info(**args.__dict__)
    zscore_info = tg.zscore_file_info(**args.__dict__)

    # PREP OUTPUT - print output headers to files
    # print("Creating file: " + out_twas_path + "\n")
    # out_cols = ["CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName", "n_snps"]
    # if args.test_stat == "both":
    #     out_cols += ["FUSION_Z", "FUSION_PVAL", "SPred_Z", "SPred_PVAL"]
    # else:
    #     out_cols += ["Zscore", "PVALUE"]

    @tg.error_handler
    def thread_process(num):
        target = TargetID[num]
        write_to_file("num=" + str(num) + "\nTargetID=" + target)
        Gene_Info = Gene.iloc[[num]].reset_index(drop=True)

        # get start and end positions to tabix
        start = str(max(int(Gene_Info.GeneStart) - args.window, 0))
        end = str(int(Gene_Info.GeneEnd) + args.window)

        # check that both files have data for target
        write_to_file("reading query")
        tabix_query = tg.tabix_query_files(start, end, **args.__dict__)

        if not tabix_query:
            write_to_file(
                "No test SNPs with non-zero cis-eQTL weights and/or no test SNPs GWAS Zscore for TargetID: "
                + target
                + "."
            )
            return None

        # read in weight data for target, filtered by weight_threshold
        write_to_file("reading Weight")
        Weight = tg.read_tabix(start, end, target=target, semaphore_key="w", **weight_info)

        # read in Zscore data
        write_to_file("reading Zscore")
        Zscore = tg.read_tabix(start, end, semaphore_key="z", **zscore_info)

        # get flipped snpIDs
        Zscore["snpIDflip"] = tg.get_snpIDs(Zscore, flip=True)

        snp_overlap = np.intersect1d(Weight.snpID, Zscore[["snpID", "snpIDflip"]])

        if not snp_overlap.size:
            write_to_file(
                "No overlapping test SNPs that have magnitude of cis-eQTL weights greater than threshold value and with GWAS Zscore for TargetID: "
                + target
                + ".\n"
            )
            return None

        # filter out non-matching snpID rows
        write_to_file("step 1")
        Weight = Weight[Weight.snpID.isin(snp_overlap)]
        Zscore = Zscore[
            np.any(Zscore[["snpID", "snpIDflip"]].isin(snp_overlap), axis=1)
        ].reset_index(drop=True)

        # if not in Weight.snpIDs, assumed flipped; if flipped, flip Zscore sign
        flip = np.where(Zscore.snpID.isin(Weight.snpID.values), 1, -1)

        write_to_file("step 2")
        if not np.all(flip == 1):
            Zscore["snpID"] = np.where(flip == 1, Zscore.snpID, Zscore.snpIDflip)
            Zscore["Zscore"] = flip * Zscore["Zscore"]

        # drop unneeded columns
        Zscore = Zscore.drop(columns=["CHROM", "POS", "REF", "ALT", "snpIDflip"])

        # merge Zscore and Weight dataframes on snpIDs
        write_to_file("step 3")
        ZW = Weight.merge(
            Zscore[["snpID", "Zscore"]], left_on="snpID", right_on="snpID", how="inner"
        )

        snp_search_ids = ZW.snpID.values

        # Read in reference covariance matrix file by snpID
        write_to_file("step 4")
        MCOV = tg.get_ld_data(args.ld_path, snp_search_ids)

        if MCOV.empty:
            write_to_file(
                "No reference covariance information for target SNPs for TargetID: "
                + target
                + "\n"
            )
            return None

        # get the snp variance and covariance matrix
        write_to_file("step 5")
        snp_sd, V_cov = tg.get_ld_matrix(MCOV)

        ZW = ZW[ZW.snpID.isin(MCOV.snpID)]
        ZW = ZW.drop_duplicates(["snpID"], keep="first").reset_index(drop=True)
        n_snps = str(ZW.snpID.size)

        write_to_file("Running TWAS.\nN SNPs=" + n_snps)

        ### create output dataframe
        Result = Gene_Info.copy()
        Result["n_snps"] = n_snps

        ### calculate zscore(s), pvalue(s)
        write_to_file("step 6")
        get_zscore_args = [V_cov, ZW.ES.values, ZW.Zscore.values, snp_sd]

        if args.test_stat == "both":
            write_to_file("step 7")
            Result["FUSION_Z"], Result["FUSION_PVAL"] = get_fusion_zscore(*get_zscore_args)
            write_to_file("step 8")
            Result["SPred_Z"], Result["SPred_PVAL"] = get_spred_zscore(*get_zscore_args)

        else:
            Result["TWAS_Zscore"], Result["PVALUE"] = get_burden_zscore(
                args.test_stat, get_zscore_args
            )

        # write to file
        # Result.to_csv(out_twas_path, sep="\t", index=None, header=None, mode="a")

        write_to_file("Target TWAS completed.\n")
        return Result

    write_to_file("Starting TWAS for " + str(n_targets) + " target genes.\n")
    res = Parallel(n_jobs=args.thread)(delayed(thread_process)(num) for num in range(n_targets))

    # with multiprocessing.Pool(args.thread) as pool:
    #     # res = pool.imap(thread_process, range(n_targets))
    #     res = pool.map_async(thread_process, range(n_targets))
    #     res = res.get()
    #     pool.close()
    #     pool.join()

    pd.concat(res).to_csv(
        out_twas_path, sep="\t", index=None, header=True, mode="w"
    )
    # pool = multiprocessing.Pool(args.thread)
    # pool.imap(thread_process,[num for num in range(n_targets)])
    # pool.close()
    # pool.join()
    write_to_file("Done.")


    ###############################################################
    # time calculation
    elapsed_sec = time() - start_time
    elapsed_time = tg.format_elapsed_time(elapsed_sec)
    write_to_file("Computation time (DD:HH:MM:SS): " + elapsed_time)


###############################################################
# thread process
if __name__ == "__main__":
    main()
