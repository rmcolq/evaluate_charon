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
key = {
    "NC_060925.1": "Chromosome 1",
    "NC_060926.1": "Chromosome 2",
    "NC_060927.1": "Chromosome 3",
    "NC_060928.1": "Chromosome 4",
    "NC_060929.1": "Chromosome 5",
    "NC_060930.1": "Chromosome 6",
    "NC_060931.1": "Chromosome 7",
    "NC_060932.1": "Chromosome 8",
    "NC_060933.1": "Chromosome 9",
    "NC_060934.1": "Chromosome 10",
    "NC_060935.1": "Chromosome 11",
    "NC_060936.1": "Chromosome 12",
    "NC_060937.1": "Chromosome 13",
    "NC_060938.1": "Chromosome 14",
    "NC_060939.1": "Chromosome 15",
    "NC_060940.1": "Chromosome 16",
    "NC_060941.1": "Chromosome 17",
    "NC_060942.1": "Chromosome 18",
    "NC_060943.1": "Chromosome 19",
    "NC_060944.1": "Chromosome 20",
    "NC_060945.1": "Chromosome 21",
    "NC_060946.1": "Chromosome 22",
    "NC_060947.1": "Chromosome X",
    "NC_060948.1": "Chromosome Y",
    "JAGYVL020000058.1": "Mitochondrion"
}

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
            ref_name = None
            if ref in ebv_accs:
                ref_name = "EBV"
            elif ref in other_accs:
                ref_name = "other"
            elif ref in key:
                ref_name = key[ref]
            entry={"read_id": read_id, "ref": ref, "ref_name":ref_name, "pos":pos, "ref_start":ref_coords[0], "ref_end":ref_coords[-1], "mapped_length":int(mapped_length), "mismatches": int(mismatches), "identity": 1-(float(mismatches)/float(mapped_length)), "divergence":float(divergence)}
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
        columns = ["read_id","ref","ref_name","pos","ref_start","ref_end","mapped_length","mismatches","identity","divergence","seq_length","gc_ratio","5mer_ratio","A","C","G","T","mapped_prop"]
        df = pd.DataFrame(columns=columns)
    sys.stderr.write("Found " + str(df.shape)  + " entries\n")
    return df

def load_blast_info(blast_results):
    processed_ids = set()
    sys.stderr.write("LOAD BLAST INFO FROM: " + blast_results + "\n")
    if (os.path.getsize(blast_results) > 0):
        default_ = {"read_id":None, "taxids":[], "names":[], "human_accs":[], "blast_human":False, "pident":0, "top_hit":None}
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
                        details[qseqid]["blast_human"] = True
                        details[qseqid]["human_accs"].append(f"{sacc}:{sstart}-{send}")
                        details[qseqid]["pident"] = float(pident)
        for qseqid in details:
            details[qseqid]["taxids"] = ";".join(list(set(details[qseqid]["taxids"])))
            details[qseqid]["names"] = ";".join(list(set(details[qseqid]["names"])))
            details[qseqid]["human_accs"] = ";".join(list(set(details[qseqid]["human_accs"])))
        df = pd.DataFrame(details.values())
    else:
        columns = ["read_id", "taxids", "names", "human_accs", "blast_human", "pident"]
        df = pd.DataFrame(columns=columns)
    df.set_index("read_id")
    sys.stderr.write("Found " + str(df.shape)  + " entries\n")
    return df


def load_both_sam(host_sam, microbial_sam):
    host_df = load_map_info_from_sam(host_sam)
    #host_df["sam"] = "host"
    microbial_df = load_map_info_from_sam(microbial_sam)
    host_df.drop(host_df[host_df.read_id.isin(microbial_df.read_id)].index, inplace=True)
    #microbial_df["sam"] = "microbial"
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
    df["charon"] = df["classification"].fillna("unclassified")
    df['classification'] = df['classification'].fillna("")

    for column in ["mean_quality", "length", "compression"]:
        m = df[column].mean()
        sd = df[column].std()
        df[f"{column}_num_stds"] = (df[column]-m)/sd
    return df

def load_tsv_output(path, classifier):
    sys.stderr.write("LOAD TSV OUTPUT from " + path + "\n")
    entries = []
    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            entry = {"status":"C", "read_id":row["read_id"], "classification": row["classification"]}
            entries.append(entry)
    df =  pd.DataFrame(entries)
    df[classifier] = df["classification"].fillna("unclassified")

    return df

def load_output(path):
    if "charon" in path:
        return "charon", load_charon_output(path)
    else:
        classifier = path.split(".")[-2]
        return classifier, load_tsv_output(path, classifier)
    

def add_classified_counts_to_summary(df, summary, classifier):
    #1. How many host, microbial, unclassified reads were there for charon?
    g = df.groupby(["status","classification"]).count()

    if ("C","human") in g["read_id"].index:
        summary[f"num_host_{classifier}"] = g["read_id"]["C"]["human"]
    else:
        summary[f"num_host_{classifier}"] = 0

    if ("C","microbial") in g["read_id"].index:
        summary[f"num_microbial_{classifier}"] = g["read_id"]["C"]["microbial"]
    else:
        summary[f"num_microbial_{classifier}"] = 0

    if ("U","") in g["read_id"].index:
        summary[f"num_unclassified_{classifier}"] = g["read_id"]["U"][""]
    else:
        summary[f"num_unclassified_{classifier}"] = 0

    summary["total"] = summary[f"num_host_{classifier}"] + summary[f"num_microbial_{classifier}"] + summary[f"num_unclassified_{classifier}"]
    summary[f"classified_{classifier}"] = summary[f"num_host_{classifier}"] + summary[f"num_microbial_{classifier}"]

    #2. Scale these to proportions
    summary[f"prop_host_{classifier}"] = summary[f"num_host_{classifier}"]/summary["total"]
    summary[f"prop_microbial_{classifier}"] = summary[f"num_microbial_{classifier}"]/summary["total"]
    summary[f"prop_unclassified_{classifier}"] = summary[f"num_unclassified_{classifier}"]/summary["total"]

    return summary

def add_host_counts_to_summary(df, summary, classifier, prefix):
    df_host = df[df["classification"] == "human"]
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

    data_file = Path(f"{prefix}_host_data.csv")
    host_host_df.to_csv(data_file, index=False)
    return

def add_microbial_counts_to_summary(df, summary, classifier, prefix):
    df_microbial = df[df["classification"] == "microbial"]
    microbial_total = df_microbial.shape[0]

    #6. Of the microbial reads, what proportion map back to the host reference genome, or EBV?
    microbial_unmapped_df = df_microbial[df_microbial["unmapped"] == True]
    summary[f"num_microbial_unmapped_{classifier}"] = microbial_unmapped_df.shape[0]
    df_microbial = df_microbial[df_microbial["unmapped"] == False]

    microbial_ebv_df = df_microbial[df_microbial["ref"].isin(ebv_accs)]
    summary[f"num_microbial_map_ebv_{classifier}"] = microbial_ebv_df.shape[0]

    microbial_host_df = df_microbial[~df_microbial["ref"].isin(ebv_accs + other_accs)]
    summary[f"num_microbial_map_host_{classifier}"] = microbial_host_df.shape[0]

    microbial_host_verified_df = microbial_host_df[microbial_host_df["blast_human"]==True]
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

    data_file = Path(f"{prefix}_microbial_data.csv")
    microbial_host_df.to_csv(data_file, index=False)
    return microbial_host_df

def check_related_taxa(microbial_host_df, classifier, prefix):
    #7. For reads which classify as microbial and minimap to host but do not have a blast human result, what taxa does blast return
    related_taxa = set()
    microbial_host_unverified_ids = microbial_host_df[microbial_host_df["blast_human"]==False]["taxids"]
    for i in microbial_host_unverified_ids:
        related_taxa.update(i.split(";"))

    if len(related_taxa) > 0:
        taxa_file = Path(f"{prefix}_related_taxa.csv")
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
    microbial_host_verified_accs = microbial_host_df[microbial_host_df["blast_human"]==True]["human_accs"]
    for i in microbial_host_verified_accs:
        human_accs.update(i.split(";"))

    if len(human_accs) > 0:
        accs_file = Path(f"{prefix}_human_accs.csv")
        with open(accs_file, "w") as f:
            f.write(",".join(human_accs))
        sys.stderr.write(f"Found human accessions which are classified as microbial for classifier {classifier}:\n{human_accs}\n")

def add_unclassified_to_summary(df, summary):
    #8. Collect basic stats for unclassified reads
    df_unclassified = df[df["status"] == "U"]
    unclassified_total = df_unclassified.shape[0]

    for column in ["length","mean_quality","confidence",'microbial_num_hits', 'microbial_prop','microbial_unique_prop', 'human_num_hits', 'human_prop', 'human_unique_prop']:
        summary[f"mean_{column}_unclassified"] = df_unclassified[column].mean()
        summary[f"median_{column}_unclassified"] = df_unclassified[column].median()
        summary[f"max_{column}_unclassified"] = df_unclassified[column].max()
        summary[f"min_{column}_unclassified"] = df_unclassified[column].min()

def generate_summary(df, prefix, classifier):
    sys.stderr.write("GENERATE SUMMARY\n")
    summary = {}

    add_classified_counts_to_summary(df, summary, classifier)

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
        help="TSV output from classifier",
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
        required=False,
        help="SAM file of mapping results from host-extracted file",
    )
    parser.add_argument(
        "--microbial_sam",
        dest="microbial_sam",
        required=False,
        help="SAM file of mapping results from microbial-extracted file",
    )
    parser.add_argument(
        "--blast_result",
        dest="blast_result",
        required=False,
        help="TAB separated result from blastn showing top blast hits for microbial reads which map to T2T reference",
    )

    args = parser.parse_args()

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stderr.write("PROGRAM START TIME: " + time + "\n")

    full_file = Path(args.prefix + "_full.csv")
    if not full_file.is_file():
        if args.host_sam is None or args.microbial_sam is None or args.blast_result is None:
            sys.stderr.write("If the full CSV file does not exist, then --host_sam, --microbial_sam and --blast_result must be provided\n")
            sys.exit(1)
            
        mapped_df = load_both_sam(args.host_sam, args.microbial_sam)
        mapped_df.to_csv("mapped_df.csv")

        blast_df = load_blast_info(args.blast_result)
        blast_df.to_csv("blast_df.csv")

        combined_df = mapped_df.merge(blast_df, how="left")
        combined_df.to_csv("combined_df.csv")
        assert combined_df.shape[0] == mapped_df.shape[0], "The number of rows in the combined dataframe does not match the mapped dataframe."

        classifier, classifier_df = load_output(args.input)

        sys.stderr.write("COMBINE CLASSIFIER AND MAPPING DATA\n")
        classifier_df.set_index("read_id")
        old_size = classifier_df.shape[0]
        classifier_df = classifier_df.merge(combined_df, how="left")
        assert classifier_df.shape[0] == old_size, "The number of rows in the charon dataframe changed when combined with mapping dataframe."
        
        classifier_df['ref_name'] = classifier_df['ref_name'].fillna("")
        classifier_df["unmapped"] = classifier_df["mapped_length"].isna()

        classifier_df["sample_id"] = args.input.split("/")[-1].split(".")[0]
        classifier_df.to_csv(full_file, index=False)
    else:
        classifier_df = pd.read_csv(full_file, index_col=None)
        classifier = full_file.split("_full.csv")[0].split("_")[-1]
        classifier_df['classification'] = classifier_df['classification'].fillna("")

    summary = generate_summary(classifier_df, args.prefix, classifier)

    # Save to CSV
    fieldnames = ["sample_id"] + list(summary.keys())
    summary["sample_id"] = args.input.split("/")[-1].split(".")[0]

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

