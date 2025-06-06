#!/usr/bin/env python
import sys
import os
import argparse
from datetime import datetime
import csv
from pathlib import Path
import pandas as pd
import numpy as np
from simplesam import Reader
import gzip
import math

ebv_accs = ["NC_007605.1","NC_009334.1"]
other_accs = ["KC670213.1","XR_003525368.1","XR_003525370.1","KC670203.1"]

def get_gc_ratio(inputStr):
    compression_ratio = len(inputStr.replace("A","").replace("T",""))/len(inputStr)
    return compression_ratio

def get_kmer_ratio(inputStr, k):
    kmers = set()
    max_size = math.pow(4,k)
    if len(inputStr) < k:
        return 0
    for i in range(len(inputStr) - k):
        kmers.add(inputStr[i:i+k])
        if i%max_size == 0:
            if len(kmers) == max_size:
                break
    return len(kmers)/max_size

def load_map_info_from_sam(sam):
    sys.stderr.write("LOAD MAP INFO FROM: " + sam + "\n")
    map_details = []
    in_file = open(sam, 'r')
    in_sam = Reader(in_file)
    for x in in_sam:
        read_id, status, ref, pos, mapped_length = x.qname, x.flag, x.rname, x.pos, len(x)
        if status not in [0, 16]:
            continue
        mismatches, divergence =  x.tags["NM"], x.tags["de"]
        entry={"read_id": read_id, "ref": ref, "pos":pos, "mapped_length":int(mapped_length), "mismatches": int(mismatches), "identity": 1-(float(mismatches)/float(mapped_length)), "divergence":float(divergence)}
        entry["seq_length"] = len(x.seq)
        entry["gc_ratio"] = get_gc_ratio(x.seq)
        entry["5mer_ratio"] = get_kmer_ratio(x.seq, 5)
        for motif in ["A","C","G","T"]:
            entry[motif] = x.gapped('seq').count(motif)
        entry["mapped_prop"] = entry["mapped_length"] / entry["seq_length"]
        map_details.append(entry)
        if len(map_details) % 10000 == 0:
            sys.stderr.write("Processed " + str(len(map_details)) + "\n")
    df = pd.DataFrame(map_details)
    sys.stderr.write("Found " + str(df.shape)  + " entries\n")
    return df


def load_both_sam(host_sam, microbial_sam):
    host_df = load_map_info_from_sam(host_sam)
    host_df["sam"] = "host"
    microbial_df = load_map_info_from_sam(microbial_sam)
    microbial_df["sam"] = "microbial"
    sys.stderr.write("COMBINE DATAFRAMES\n")
    df = pd.concat([host_df, microbial_df], ignore_index=True)
    df.set_index("read_id")
    return df

def load_charon_output(path):
    sys.stderr.write("LOAD CHARON OUTPUT from " + path + "\n")
    df = pd.read_csv(path, sep="\t", index_col=False, header=None)
    entries = []
    for  i, row in df.iterrows():
        entry = {"status":row[0], "read_id":row[1], "classification": row[2], "length": row[3], "num_hashes":row[4],   "mean_quality": row[5], "confidence": row[6], "compression":row[7]}
        try:
            details = row[8].split(" ")
        except:
            print(row)
        for part in details:
            try:
                category, num_hits, prop_hits, prop_unique_hits, prob = part.split(":")
            except:
                continue
            entry[f"p_{category}"] = float(prob)
            entry[f"{category}_num_hits"] = int(num_hits)
            entry[f"{category}_prop"] = float(prop_hits)
            entry[f"{category}_unique_prop"] = float(prop_unique_hits)
        entries.append(entry)
    df =  pd.DataFrame(entries)
    df['classification'] = df['classification'].fillna("")
    for column in ["mean_quality", "length", "compression"]:
        m = df[column].mean()
        sd = df[column].std()
        df[f"{column}_num_stds"] = (df["column"]-m)/s
    return df

def generate_summary(df):
    sys.stderr.write("GENERATE SUMMARY\n")
    summary = {}

    #1. How many host, microbial, unclassified reads were there?
    g = df.groupby(["status","classification"]).count()
    if ("C","human") in g["read_id"].index:
        summary["num_host"] = g["read_id"]["C"]["human"]
    else:
        summary["num_host"] = 0
    if ("C","microbial") in g["read_id"].index:
        summary["num_microbial"] = g["read_id"]["C"]["microbial"]
    else:
        summary["num_microbial"] = 0
    if ("U","") in g["read_id"].index:
        summary["num_unclassified"] = g["read_id"]["U"][""]
    else:
        summary["num_unclassified"] = 0
    summary["total"] = summary["num_host"] + summary["num_microbial"] + summary["num_unclassified"]
    summary["classified"] = summary["num_host"] + summary["num_microbial"]

    #2. Scale these to proportions
    summary["prop_host"] = summary["num_host"]/summary["total"]
    summary["prop_microbial"] = summary["num_microbial"]/summary["total"]
    summary["prop_unclassified"] = summary["num_unclassified"]/summary["total"]

    #3. Of the host reads, what proportion map back to the host reference genome, or EBV (minimap2 T2T+EBV)?
    df_host = df[df["classification"] == "human"]
    host_total = df_host.shape[0]
    host_unmapped_df = df_host[df_host["unmapped"] == True]
    summary["num_host_unmapped"] = host_unmapped_df.shape[0]
    df_host = df_host[df_host["unmapped"] == False]
    host_ebv_df = df_host[df_host["ref"].isin(ebv_accs)]
    summary["num_host_map_ebv"] = host_ebv_df.shape[0]
    host_host_df = df_host[~df_host["ref"].isin(ebv_accs + other_accs)]
    summary["num_host_map_host"] = host_host_df.shape[0]
    if host_total > 0:
        summary["prop_host_unmapped"] = summary["num_host_unmapped"]/host_total
        summary["prop_host_map_ebv"] = summary["num_host_map_ebv"]/host_total
        summary["prop_host_map_host"] = summary["num_host_map_host"]/host_total
    else:
        summary["prop_host_unmapped"] = 0
        summary["prop_host_map_ebv"] = 0
        summary["prop_host_map_host"] = 0

    #4. For reads which classify as host and map to host, what is the (mean, median, max) length, quality, confidence, prop_unique_microbial, prop_unique_host, prop_microbial prop_host, num_hits_microbial, num_hits_host
    for column in ["length","mean_quality","confidence",'microbial_num_hits', 'microbial_prop','microbial_unique_prop', 'human_num_hits', 'human_prop', 'human_unique_prop', 'gc_ratio', '5mer_ratio', 'compression', 'mapped_prop']:
        summary[f"mean_{column}_host_map_host"] = host_host_df[column].mean()
        summary[f"median_{column}_host_map_host"] = host_host_df[column].median()
        summary[f"max_{column}_host_map_host"] = host_host_df[column].max()
        summary[f"min_{column}_host_map_host"] = host_host_df[column].min()

    #5. Of the microbial reads, what proportion map back to the host reference genome, or EBV?
    df_microbial = df[df["classification"] == "microbial"]
    microbial_total = df_microbial.shape[0]
    microbial_unmapped_df = df_microbial[df_microbial["unmapped"] == True]
    summary["num_microbial_unmapped"] = microbial_unmapped_df.shape[0]
    df_microbial = df_microbial[df_microbial["unmapped"] == False]
    microbial_ebv_df = df_microbial[df_microbial["ref"].isin(ebv_accs)]
    summary["num_microbial_map_ebv"] = microbial_ebv_df.shape[0]
    microbial_host_df = df_microbial[~df_microbial["ref"].isin(ebv_accs + other_accs)]
    summary["num_microbial_map_host"] = microbial_host_df.shape[0]
    microbial_host_strong_df = microbial_host_df[microbial_host_df["mapped_prop"]>0.8]
    summary["num_microbial_map_host_strong"] = microbial_host_strong_df.shape[0]

    if microbial_total > 0:
        summary["prop_microbial_unmapped"] = summary["num_microbial_unmapped"]/microbial_total
        summary["prop_microbial_map_ebv"] = summary["num_microbial_map_ebv"]/microbial_total
        summary["prop_microbial_map_host"] = summary["num_microbial_map_host"]/microbial_total
        summary["prop_microbial_map_host_strong"] = summary["num_microbial_map_host_strong"]/microbial_total
    else:
        summary["prop_microbial_unmapped"] = 0
        summary["prop_microbial_map_ebv"] = 0
        summary["prop_microbial_map_host"] = 0
        summary["prop_microbial_map_host_strong"] = 0

    #6. For reads which classify as microbial but map to host, what is the (mean, median, max) length, quality, confidence, prop_unique_microbial, prop_unique_host, prop_microbial prop_host, num_hits_microbial, num_hits_host
    for column in ["length","mean_quality","confidence",'microbial_num_hits', 'microbial_prop','microbial_unique_prop', 'human_num_hits', 'human_prop', 'human_unique_prop', 'gc_ratio', '5mer_ratio', 'compression', 'mapped_prop']:
        summary[f"mean_{column}_microbial_map_host"] = microbial_host_df[column].mean()
        summary[f"median_{column}_microbial_map_host"] = microbial_host_df[column].median()
        summary[f"max_{column}_microbial_map_host"] = microbial_host_df[column].max()
        summary[f"min_{column}_microbial_map_host"] = microbial_host_df[column].min()
    for column in ["length","mean_quality","confidence",'microbial_num_hits', 'microbial_prop','microbial_unique_prop', 'human_num_hits', 'human_prop', 'human_unique_prop', 'gc_ratio', '5mer_ratio', 'compression', 'mapped_prop']:
        summary[f"mean_{column}_microbial_map_host_strong"] = microbial_host_strong_df[column].mean()
        summary[f"median_{column}_microbial_map_host_strong"] = microbial_host_strong_df[column].median()
        summary[f"max_{column}_microbial_map_host_strong"] = microbial_host_strong_df[column].max()
        summary[f"min_{column}_microbial_map_host_strong"] = microbial_host_strong_df[column].min()

    #7. For reads which classify as microbial and do not classify as host or ebv , what is the (mean, median, max) length, quality, confidence, prop_unique_microbial, prop_unique_host, prop_microbial prop_host, num_hits_microbial, num_hits_host
    df_microbial_microbial = microbial_unmapped_df
    for column in ["length","mean_quality","confidence",'microbial_num_hits', 'microbial_prop','microbial_unique_prop', 'human_num_hits', 'human_prop', 'human_unique_prop', 'gc_ratio', '5mer_ratio', 'compression', 'mapped_prop']:
        summary[f"mean_{column}_microbial_map_microbial"] = df_microbial_microbial[column].mean()
        summary[f"median_{column}_microbial_map_microbial"] = df_microbial_microbial[column].median()
        summary[f"max_{column}_microbial_map_microbial"] = df_microbial_microbial[column].max()
        summary[f"min_{column}_microbial_map_microbial"] = df_microbial_microbial[column].min()

    return summary, microbial_host_df, host_unmapped_df

# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        dest="input",
        required=True,
        help="TSV output by charon",
    )
    parser.add_argument(
        "-p",
        dest="prefix",
        required=True,
        help="Prefix for output CSV files",
    )
    parser.add_argument(
        "--host_sam",
        dest="host_sam",
        required=True,
        help="SAM file of mapping results from host-extracted file",
    )
    parser.add_argument(
        "--microbial_sam",
        dest="microbial_sam",
        required=True,
        help="SAM file of mapping results from microbial-extracted file",
    )

    args = parser.parse_args()

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    full_file = Path(args.prefix + "_full.csv")
    if not full_file.is_file():
        mapped_df = load_both_sam(args.host_sam, args.microbial_sam)
        charon_df = load_charon_output(args.input)
        sys.stderr.write("COMBINE CHARON AND MAPPING DATA\n")
        charon_df.set_index("read_id")
        charon_df = charon_df.merge(mapped_df, how="left")
        charon_df["unmapped"] = charon_df["mapped_length"].isna()
        charon_df["file"] = args.input
        charon_df.to_csv(full_file, index=False)
    else:
        charon_df = pd.read_csv(full_file, index_col=None)
        charon_df['classification'] = charon_df['classification'].fillna("")

    summary,microbial_host_df, host_unmapped_df = generate_summary(charon_df)
    data_file = Path(args.prefix + "_microbial_data.csv")
    microbial_host_df.to_csv(data_file, index=False)
    data_file = Path(args.prefix + "_host_data.csv")
    host_unmapped_df.to_csv(data_file, index=False)

    # Save to CSV
    fieldnames = ["file"] + list(summary.keys())
    summary["file"] = args.input

    summary_file = Path(args.prefix + "_summary.csv")
    writer = None
    if summary_file.is_file():
        out_handle = open(summary_file, 'a', newline='')
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
    else:
        out_handle = open(summary_file, 'w', newline='')
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
        writer.writeheader()
    writer.writerow(summary)

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM END TIME: " + time + "\n")

    sys.exit(0)


if __name__ == "__main__":
    main()

