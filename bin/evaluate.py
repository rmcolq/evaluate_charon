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
from collections import defaultdict
import taxoniq
import copy

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
    processed_read_ids = set()
    sys.stderr.write("LOAD MAP INFO FROM: " + sam + "\n")
    if (os.path.getsize(sam) > 0):
        map_details = []
        in_file = open(sam, 'r')
        in_sam = Reader(in_file)
        for x in in_sam:
            read_id, status, ref, pos, ref_coords, mapped_length = x.qname, x.flag, x.rname, x.pos, x.coords, len(x)
            if status not in [0, 16]:
                continue
            if read_id in processed_read_ids:
                continue
            else:
                processed_read_ids.add(read_id)
            mismatches, divergence =  x.tags["NM"], x.tags["de"]
            entry={"read_id": read_id, "ref": ref, "pos":pos, "ref_start":ref_coords[0], "ref_end":ref_coords[-1], "mapped_length":int(mapped_length), "mismatches": int(mismatches), "identity": 1-(float(mismatches)/float(mapped_length)), "divergence":float(divergence)}
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
    else:
        columns = ["read_id","ref","pos","ref_start","ref_end","mapped_length","mismatches","identity","divergence","seq_length","gc_ratio","5mer_ratio","A","C","G","T","mapped_prop"]
        df = pd.DataFrame(columns=columns)
    sys.stderr.write("Found " + str(df.shape)  + " entries\n")
    return df

def load_blast_info(blast_results):
    processed_ids = set()
    sys.stderr.write("LOAD BLAST INFO FROM: " + blast_results + "\n")
    if (os.path.getsize(blast_results) > 0):
        default_ = {"read_id":None, "taxids":[], "names":[], "human_accs":[], "human":False, "pident":0, "top_hit":None}
        details = defaultdict(lambda:copy.deepcopy(default_))
        with open(blast_results, 'r') as f:
            for line in f:
                if len(line.strip())==0:
                    continue
                qseqid,sacc,sscinames,staxids,sstart,send,evalue,pident,length = line.strip().split()

                if f"{qseqid}_{staxids}" in processed_ids:
                    continue
                else:
                    processed_ids.add(f"{qseqid}_{staxids}")

                if details[qseqid]["top_hit"] is None:
                    details[qseqid]["top_hit"] = staxids
                if float(pident) < 90 and len(details[qseqid]["taxids"]) > 0:
                    continue
                if staxids not in details[qseqid]["taxids"]:
                    details[qseqid]["taxids"].append(staxids)
                    details[qseqid]["names"].append(sscinames)
                    details[qseqid]["read_id"] = qseqid
                    if staxids == "9606":
                        details[qseqid]["human"] = True
                        details[qseqid]["human_accs"].append(f"{sacc}:{sstart}-{send}")
                        details[qseqid]["pident"] = float(pident)
        for qseqid in details:
            details[qseqid]["taxids"] = ";".join(list(set(details[qseqid]["taxids"])))
            details[qseqid]["names"] = ";".join(list(set(details[qseqid]["names"])))
            details[qseqid]["human_accs"] = ";".join(list(set(details[qseqid]["human_accs"])))
        df = pd.DataFrame(details.values())
    else:
        columns = ["read_id", "taxids", "names", "human_accs", "human", "pident"]
        df = pd.DataFrame(columns=columns)
    df.set_index("read_id")
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
        df[f"{column}_num_stds"] = (df[column]-m)/sd
    return df

def load_tsv_output(path, classifier, df):
    sys.stderr.write("LOAD TSV OUTPUT from " + path + "\n")
    entries = []
    with open(path, newline='', sep="\t") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            entry = {"read_id":row["read_id"], classifier: row["classification"]}
            entries.append(entry)
    new_df =  pd.DataFrame(entries)
    new_df.set_index("read_id")
    df = df.merge(new_df, how="left")
    return df

def add_classified_counts_to_summary(df, summary, others=[]):
    #1. How many host, microbial, unclassified reads were there for charon?
    g = df.groupby(["status","classification"]).count()

    if ("C","human") in g["read_id"].index:
        summary["num_host_charon"] = g["read_id"]["C"]["human"]
    else:
        summary["num_host_charon"] = 0

    if ("C","microbial") in g["read_id"].index:
        summary["num_microbial_charon"] = g["read_id"]["C"]["microbial"]
    else:
        summary["num_microbial_charon"] = 0

    if ("U","") in g["read_id"].index:
        summary["num_unclassified_charon"] = g["read_id"]["U"][""]
    else:
        summary["num_unclassified_charon"] = 0

    summary["total"] = summary["num_host_charon"] + summary["num_microbial_charon"] + summary["num_unclassified_charon"]
    summary["classified"] = summary["num_host_charon"] + summary["num_microbial_charon"]

    #2. Scale these to proportions
    summary["prop_host_charon"] = summary["num_host_charon"]/summary["total"]
    summary["prop_microbial_charon"] = summary["num_microbial_charon"]/summary["total"]
    summary["prop_unclassified_charon"] = summary["num_unclassified_charon"]/summary["total"]

    # Process OTHER classifiers
    for classifier in others:
        #3. How many host, microbial, unclassified reads were there?
        g = df.groupby([classifier]).count()

        if ("human") in g["read_id"].index:
            summary[f"num_host_{classifier}"] = g["read_id"]["human"]
        else:
            summary[f"num_host_{classifier}"] = 0

        if ("microbial") in g["read_id"].index:
            summary[f"num_microbial_{classifier}"] = g["read_id"]["microbial"]
        else:
            summary[f"num_microbial_{classifier}"] = 0

        #4. Scale these to proportions
        assert summary["total"] == summary[f"num_host_{classifier}"] + summary[f"num_microbial_{classifier}"]
        summary[f"prop_host_{classifier}"] = summary[f"num_host_{classifier}"]/summary["total"]
        summary[f"prop_microbial_{classifier}"] = summary[f"num_microbial_{classifier}"]/summary["total"]

    return summary

def add_host_counts_to_summary(df, summary, classifier, prefix):
    df_host = df[df["classification"] == "human"]
    if classifier != "charon":
        df_host = df[df[classifier] == "human"]
    host_total = df_host.shape[0]

    #5. Of the host reads, what proportion map back to the host reference genome, or EBV (minimap2 T2T+EBV)?
    host_unmapped_df = df_host[df_host["unmapped"] == True]
    summary[f"num_host_unmapped_{classifier}"] = host_unmapped_df.shape[0]
    df_host = df_host[df_host["unmapped"] == False]

    host_ebv_df = df_host[df_host["ref"].isin(ebv_accs)]
    summary[f"num_host_map_ebv_{classifier}"] = host_ebv_df.shape[0]

    host_host_df = df_host[~df_host["ref"].isin(ebv_accs + other_accs)]
    summary[f"num_host_map_host_{classifier}"] = host_host_df.shape[0]

    if host_total > 0:
        summary[f"prop_host_unmapped_{classifier}"] = summary[f"num_host_unmapped_{classifier}"]/host_total
        summary[f"prop_host_map_ebv_{classifier}"] = summary[f"num_host_map_ebv_{classifier}"]/host_total
        summary[f"prop_host_map_host_{classifier}"] = summary[f"num_host_map_host_{classifier}"]/host_total
    else:
        summary[f"prop_host_unmapped_{classifier}"] = 0
        summary[f"prop_host_map_ebv_{classifier}"] = 0
        summary[f"prop_host_map_host_{classifier}"] = 0

    data_file = Path(f"{prefix}_{classifier}_host_data.csv")
    host_host_df.to_csv(data_file, index=False)
    return

def add_microbial_counts_to_summary(df, summary, classifier, prefix):
    df_microbial = df[df["classification"] == "human"]
    if classifier != "charon":
        df_microbial = df[df[classifier] == "microbial"]
    microbial_total = df_microbial.shape[0]

    #6. Of the microbial reads, what proportion map back to the host reference genome, or EBV?
    microbial_unmapped_df = df_microbial[df_microbial["unmapped"] == True]
    summary[f"num_microbial_unmapped_{classifier}"] = microbial_unmapped_df.shape[0]
    df_microbial = df_microbial[df_microbial["unmapped"] == False]

    microbial_ebv_df = df_microbial[df_microbial["ref"].isin(ebv_accs)]
    summary[f"num_microbial_map_ebv_{classifier}"] = microbial_ebv_df.shape[0]

    microbial_host_df = df_microbial[~df_microbial["ref"].isin(ebv_accs + other_accs)]
    summary[f"num_microbial_map_host_{classifier}"] = microbial_host_df.shape[0]

    microbial_host_verified_df = microbial_host_df[microbial_host_df["human"]==True]
    summary[f"num_microbial_map_host_verified_{classifier}"] = microbial_host_verified_df.shape[0]

    if microbial_total > 0:
        summary[f"prop_microbial_unmapped_{classifier}"] = summary[f"num_microbial_unmapped_{classifier}"]/microbial_total
        summary[f"prop_microbial_map_ebv_{classifier}"] = summary[f"num_microbial_map_ebv_{classifier}"]/microbial_total
        summary[f"prop_microbial_map_host_{classifier}"] = summary[f"num_microbial_map_host_{classifier}"]/microbial_total
        summary[f"prop_microbial_map_host_verified_{classifier}"] = summary[f"num_microbial_map_host_verified_{classifier}"]/microbial_total
    else:
        summary[f"prop_microbial_unmapped_{classifier}"] = 0
        summary[f"prop_microbial_map_ebv_{classifier}"] = 0
        summary[f"prop_microbial_map_host_{classifier}"] = 0
        summary[f"prop_microbial_map_host_verified_{classifier}"] = 0

    data_file = Path(f"{prefix}_{classifier}_microbial_data.csv")
    microbial_host_df.to_csv(data_file, index=False)
    return microbial_host_df

def check_related_taxa(microbial_host_df, classifier, prefix):
    #7. For reads which classify as microbial and minimap to host but do not have a blast human result, what taxa does blast return
    related_taxa = set()
    microbial_host_unverified_ids = microbial_host_df[microbial_host_df["human"]==False]["taxids"]
    for i in microbial_host_unverified_ids:
        related_taxa.update(i.split(";"))

    if len(related_taxa) > 0:
        taxa_file = Path(f"{prefix}_{classifier}_related_taxa.csv")
        with open(taxa_file, "w") as f:
            species_ids = []
            species_names = []
            for taxid in related_taxa:
                t = taxoniq.Taxon(int(taxid))
                for s in t.ranked_lineage:
                    if s.rank.name == "species":
                        if s.tax_id not in species_ids:
                            species_ids.append(s.tax_id)
                            species_names.append(s.scientific_name)
            species_ids = [str(id) for id in species_ids]
            f.write(f"{','.join(related_taxa)}\n")
            f.write(f"{','.join(species_ids)}\n")
            f.write(f"{','.join(species_names)}\n")
        sys.stderr.write(f"Found microbial taxa which are closely related to human for classifier {classifier}:\n{species_names}\n")

def check_human_accs(microbial_host_df, classifier, prefix):
    #8. For reads which classify as microbial and map to host and have a blast human result, what human accessions
    human_accs = set()
    microbial_host_verified_accs = microbial_host_df[microbial_host_df["human"]==True]["human_accs"]
    for i in microbial_host_verified_accs:
        human_accs.update(i.split(";"))

    if len(human_accs) > 0:
        accs_file = Path(f"{prefix}_{classifier}_human_accs.csv")
        with open(accs_file, "w") as f:
            f.write(",".join(human_accs))
        sys.stderr.write(f"Found human accessions which are classified as microbial for classifier {classifier}:\n{human_accs}\n")

def add_unclassified_to_summary(df, summary):
    #8. Collect basic stats for unclassified reads
    df_unclassified = df[df["status"] == "U"]
    unclassified_total = df_unclassified.shape[0]

    for column in ["length","mean_quality","confidence",'microbial_num_hits', 'microbial_prop','microbial_unique_prop', 'human_num_hits', 'human_prop', 'human_unique_prop']:
        summary[f"mean_{column}_unclassified_charon"] = df_unclassified[column].mean()
        summary[f"median_{column}_unclassified_charon"] = df_unclassified[column].median()
        summary[f"max_{column}_unclassified_charon"] = df_unclassified[column].max()
        summary[f"min_{column}_unclassified_charon"] = df_unclassified[column].min()

def generate_summary(df, prefix, others=[]):
    sys.stderr.write("GENERATE SUMMARY\n")
    summary = {}

    add_classified_counts_to_summary(df, summary, others=others)

    for classifier in ["charon"] + others:
        add_host_counts_to_summary(df, summary, classifier, prefix)

        microbial_host_df = add_microbial_counts_to_summary(df, summary, classifier, prefix)
        check_related_taxa(microbial_host_df, classifier, prefix)
        check_human_accs(microbial_host_df, classifier, prefix)

        if classifier == "charon":
            add_unclassified_to_summary(df, summary)

    return summary

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
        "-a",
        dest="additional",
        required=False,
        help="Additional TSV for another classifier",
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
    parser.add_argument(
        "--blast_result",
        dest="blast_result",
        required=True,
        help="TAB separated result from blastn showing top blast hits for microbial reads which map to T2T reference",
    )

    args = parser.parse_args()

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    full_file = Path(args.prefix + "_full.csv")
    if not full_file.is_file():
        mapped_df = load_both_sam(args.host_sam, args.microbial_sam)
        mapped_df.to_csv("mapped_df.csv")

        blast_df = load_blast_info(args.blast_result)
        blast_df.to_csv("blast_df.csv")

        combined_df = mapped_df.merge(blast_df, how="left")
        combined_df.to_csv("combined_df.csv")

        charon_df = load_charon_output(args.input)
        if args.additional:
            charon_df = load_tsv_output(args.additional, "deacon", charon_df)

        sys.stderr.write("COMBINE CHARON AND MAPPING DATA\n")
        charon_df.set_index("read_id")
        charon_df = charon_df.merge(combined_df, how="left")
        charon_df["unmapped"] = charon_df["mapped_length"].isna()

        charon_df["file"] = args.input
        charon_df.to_csv(full_file, index=False)
    else:
        charon_df = pd.read_csv(full_file, index_col=None)
        charon_df['classification'] = charon_df['classification'].fillna("")

    summary = generate_summary(charon_df)

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

