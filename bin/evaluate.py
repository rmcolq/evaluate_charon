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

ebv_accs = ["NC_007605.1","NC_009334.1"]
mycoplasma_accs = ["NC_010163.1","NC_014921.1","NZ_CP009652.1","NZ_AP014657.1","NZ_LS991950.1","NZ_LR214940.1","NZ_LR214941.1","NZ_LR214942.1","NZ_LR214943.1","NZ_LR214944.1"]
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
    if len(inputStr) == 0:
        return 0
    seq = inputStr.upper()
    return (seq.count("G") + seq.count("C"))/len(seq)

def get_kmer_ratio(inputStr, k):
    kmers = set()
    max_size = math.pow(4,k)
    if len(inputStr) < k:
        return 0
    for i in range(len(inputStr) - k + 1):
        kmers.add(inputStr[i:i+k])
        if i%max_size == 0:
            if len(kmers) == max_size:
                break
    return len(kmers)/max_size

def get_mean_quality(qual):
    """Mean Phred quality of a SAM QUAL string, averaged over error probabilities."""
    if qual is None or qual == "*" or len(qual) == 0:
        return np.nan
    phred = np.frombuffer(qual.encode(), dtype=np.uint8).astype(float) - 33
    return float(-10 * np.log10(np.mean(np.power(10, -phred / 10))))

def load_contaminant_accs(path):
    """Read extra contaminant accessions (whitespace or comma separated) which should not be counted as human."""
    accs = []
    if path is None:
        return accs
    with open(path) as f:
        for line in f:
            for acc in line.replace(",", " ").split():
                if not acc.startswith("#"):
                    accs.append(acc)
    return accs

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
            elif ref in mycoplasma_accs:
                ref_name = "mycoplasma"
            elif ref in other_accs:
                ref_name = "other"
            elif ref in key:
                ref_name = key[ref]
            entry={"read_id": read_id, "ref": ref, "ref_name":ref_name, "pos":pos, "ref_start":ref_coords[0], "ref_end":ref_coords[-1], "mapped_length":int(mapped_length), "mismatches": int(mismatches), "identity": 1-(float(mismatches)/float(mapped_length)), "divergence":float(divergence)}
            entry["seq_length"] = len(x.seq)
            entry["gc_ratio"] = get_gc_ratio(x.seq)
            entry["5mer_ratio"] = get_kmer_ratio(x.seq, 5)
            entry["sam_mean_quality"] = get_mean_quality(x.qual)
            for motif in ["A","C","G","T"]:
                entry[motif] = x.gapped('seq').count(motif)
            entry["mapped_prop"] = entry["mapped_length"] / entry["seq_length"]
            map_details.append(entry)
            if len(map_details) % 10000 == 0:
                sys.stderr.write("Processed " + str(len(map_details)) + "\n")
        df = pd.DataFrame(map_details)
    else:
        columns = ["read_id","ref","ref_name","pos","ref_start","ref_end","mapped_length","mismatches","identity","divergence","seq_length","gc_ratio","5mer_ratio","sam_mean_quality","A","C","G","T","mapped_prop"]
        df = pd.DataFrame(columns=columns)
    sys.stderr.write("Found " + str(df.shape)  + " entries\n")
    return df

BLAST_COLUMNS = ["read_id", "taxids", "names", "human_accs", "blast_human", "pident", "top_hit",
                 "human_bitscore", "human_qcov", "human_pident", "nonhuman_bitscore", "nonhuman_qcov",
                 "nonhuman_pident", "nonhuman_taxid", "nonhuman_name", "blast_nonhuman_run"]

def is_human_hit(staxids):
    return "9606" in staxids.split(";")

def parse_blast_line(line):
    """Return (hit_set, fields) where hit_set is 'all' or 'nonhuman' and fields are
    qseqid, sacc, sscinames, staxids, sstart, send, evalue, pident, length, bitscore, qcovs.
    Handles the legacy 9 column format (no bitscore or coverage) too."""
    parts = line.rstrip("\n").split("\t")
    if len(parts) == 12:
        hit_set, parts = parts[0], parts[1:]
    elif len(parts) == 11:
        hit_set = "all"
    elif len(parts) == 9:
        hit_set, parts = "all", parts + ["nan", "nan"]
    else:
        raise ValueError(f"Unexpected number of BLAST columns ({len(parts)}): {line}")
    return hit_set, parts

def load_blast_info(blast_results):
    sys.stderr.write("LOAD BLAST INFO FROM: " + blast_results + "\n")
    if (os.path.getsize(blast_results) == 0):
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        sys.stderr.write("Found " + str(df.shape)  + " entries\n")
        return df

    details = {}
    seen = set()
    nonhuman_run = False
    with open(blast_results, 'r') as f:
        for line in f:
            if len(line.strip())==0:
                continue
            hit_set, parts = parse_blast_line(line)
            qseqid,sacc,sscinames,staxids,sstart,send,evalue,pident,length,bitscore,qcovs = parts
            pident, bitscore, qcovs = float(pident), float(bitscore), float(qcovs)
            if hit_set == "nonhuman":
                nonhuman_run = True
            entry = details.setdefault(qseqid, {"read_id":qseqid, "taxids":[], "names":[], "human_accs":[], "blast_human":False,
                        "pident":0, "top_hit":None, "human_bitscore":np.nan, "human_qcov":np.nan, "human_pident":np.nan,
                        "nonhuman_bitscore":np.nan, "nonhuman_qcov":np.nan, "nonhuman_pident":np.nan,
                        "nonhuman_taxid":None, "nonhuman_name":None})
            human = is_human_hit(staxids)

            # best hit by bitscore, human and non-human separately (blast output order is not by bitscore)
            if human:
                if not bitscore <= entry["human_bitscore"]:
                    entry["human_bitscore"], entry["human_qcov"], entry["human_pident"] = bitscore, qcovs, pident
            elif not bitscore <= entry["nonhuman_bitscore"]:
                entry["nonhuman_bitscore"], entry["nonhuman_qcov"], entry["nonhuman_pident"] = bitscore, qcovs, pident
                entry["nonhuman_taxid"], entry["nonhuman_name"] = staxids, sscinames

            if hit_set == "nonhuman":
                continue

            # summary of the taxa hit (first hit, then all taxa with at least 90% identity)
            if f"{qseqid}_{staxids}" in seen:
                continue
            seen.add(f"{qseqid}_{staxids}")
            if entry["top_hit"] is None:
                entry["top_hit"] = staxids
            if pident < 90 and len(entry["taxids"]) > 0:
                continue
            if staxids not in entry["taxids"]:
                entry["taxids"].append(staxids)
                entry["names"].append(sscinames)
                if human:
                    entry["blast_human"] = True
                    entry["human_accs"].append(f"{sacc}:{sstart}-{send}")
                    entry["pident"] = pident
    for entry in details.values():
        entry["taxids"] = ";".join(sorted(set(entry["taxids"])))
        entry["names"] = ";".join(sorted(set(entry["names"])))
        entry["human_accs"] = ";".join(sorted(set(entry["human_accs"])))
        entry["blast_nonhuman_run"] = nonhuman_run
    df = pd.DataFrame(details.values(), columns=BLAST_COLUMNS)
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
    return df

def load_charon_output(path):
    sys.stderr.write("LOAD CHARON OUTPUT from " + path + "\n")
    df = pd.read_csv(path, sep="\t", index_col=False, header=None)
    entries = []
    for  i, row in df.iterrows():
        entry = {"status":row[0], "read_id":row[1], "classification": row[2], "length": row[3], "num_hashes":row[4],   "mean_quality": row[5], "confidence": row[6], "compression":row[7]}
        try:
            details = row[8].split(" ")
        except AttributeError:
            # no per-category details for this read (NaN), so do not reuse the previous read's
            sys.stderr.write(f"No category details for read {row[1]}\n")
            details = []
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

def load_output(path, classifier=None):
    if classifier is None:
        # legacy behaviour: infer from the file name
        classifier = "charon" if "charon" in path else path.split(".")[-2]
    if classifier == "charon":
        return classifier, load_charon_output(path)
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
    total = max(summary["total"], 1)
    summary[f"prop_host_{classifier}"] = summary[f"num_host_{classifier}"]/total
    summary[f"prop_microbial_{classifier}"] = summary[f"num_microbial_{classifier}"]/total
    summary[f"prop_unclassified_{classifier}"] = summary[f"num_unclassified_{classifier}"]/total

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

    host_host_df = df_host[~df_host["ref"].isin(ebv_accs + mycoplasma_accs + other_accs)]
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

    microbial_host_df = df_microbial[~df_microbial["ref"].isin(ebv_accs + mycoplasma_accs + other_accs)]
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
    for i in microbial_host_unverified_ids.dropna():
        related_taxa.update([t for t in str(i).split(";") if t.isdigit() and t != "0"])

    if len(related_taxa) > 0:
        import taxoniq
        taxa_file = Path(f"{prefix}_related_taxa.csv")
        with open(taxa_file, "w") as f:
            species_ids = []
            species_names = []
            for taxid in sorted(related_taxa):
                try:
                    t = taxoniq.Taxon(int(taxid))
                    lineage = t.ranked_lineage
                except KeyError:
                    sys.stderr.write(f"Taxon {taxid} not found in taxoniq\n")
                    continue
                for s in lineage:
                    if s.rank.name == "species":
                        if s.tax_id not in species_ids:
                            species_ids.append(s.tax_id)
                            species_names.append(s.scientific_name)
            species_ids = [str(id) for id in species_ids]
            f.write(f"{','.join(sorted(related_taxa))}\n")
            f.write(f"{','.join(species_ids)}\n")
            f.write(f"{','.join(species_names)}\n")
        sys.stderr.write(f"Found microbial taxa which are closely related to human for classifier {classifier}:\n{species_names}\n")

def check_human_accs(microbial_host_df, classifier, prefix):
    #8. For reads which classify as microbial and map to host and have a blast human result, what human accessions
    human_accs = set()
    microbial_host_verified_accs = microbial_host_df[microbial_host_df["blast_human"]==True]["human_accs"]
    for i in microbial_host_verified_accs.dropna():
        human_accs.update(str(i).split(";"))

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

def add_leak_tiers(df, min_blast_cov, min_blast_pident, bitscore_ratio):
    """Assign each read a leak tier (0 = not called microbial); each tier is a subset of the one above.
    1: called microbial
    2: 1 and maps to the human (T2T) reference, i.e. not unmapped, EBV or another contaminant
    3: 2 and a BLAST human hit with at least min_blast_cov % query coverage and min_blast_pident % identity
    4: 3 and the best human bitscore is at least bitscore_ratio times the best non-human bitscore
       (or there was no non-human hit), so the read is human-specific rather than a shared region.
       Only assigned when the non-human BLAST search was run."""
    for column in BLAST_COLUMNS[1:]:
        if column not in df.columns:
            df[column] = np.nan
    microbial = df["classification"] == "microbial"
    t2 = microbial & ~df["unmapped"] & ~df["ref"].isin(ebv_accs + other_accs)
    human_hit = df["blast_human"].fillna(False).astype(bool)
    if df["human_qcov"].notna().any():
        human_hit = human_hit & (df["human_qcov"] >= min_blast_cov) & (df["human_pident"] >= min_blast_pident)
    t3 = t2 & human_hit
    nonhuman_run = df["blast_nonhuman_run"].fillna(False).astype(bool)
    better_than_nonhuman = df["nonhuman_bitscore"].isna() | (df["human_bitscore"] >= bitscore_ratio * df["nonhuman_bitscore"])
    t4 = t3 & nonhuman_run & better_than_nonhuman
    df["leak_tier"] = microbial.astype(int) + t2.astype(int) + t3.astype(int) + t4.astype(int)
    return df

def add_leak_tiers_to_summary(df, summary, classifier):
    total = max(df.shape[0], 1)
    for tier in [1, 2, 3, 4]:
        n = int((df["leak_tier"] >= tier).sum())
        summary[f"num_leak_tier{tier}_{classifier}"] = n
        summary[f"leaks_per_million_tier{tier}_{classifier}"] = n / total * 1e6

def generate_summary(df, prefix, classifier):
    sys.stderr.write("GENERATE SUMMARY\n")
    summary = {}

    add_classified_counts_to_summary(df, summary, classifier)

    add_host_counts_to_summary(df, summary, classifier, prefix)

    microbial_host_df = add_microbial_counts_to_summary(df, summary, classifier, prefix)
    check_related_taxa(microbial_host_df, classifier, prefix)
    check_human_accs(microbial_host_df, classifier, prefix)

    add_leak_tiers_to_summary(df, summary, classifier)

    if classifier == "charon":
        add_unclassified_to_summary(df, summary)

    return summary

def write_table(df, path_without_ext):
    """Write parquet if pyarrow is installed, else gzipped CSV."""
    try:
        import pyarrow  # noqa: F401
        df.to_parquet(f"{path_without_ext}.parquet", index=False)
    except ImportError:
        sys.stderr.write("pyarrow is not installed, so writing gzipped CSV instead of parquet\n")
        df.to_csv(f"{path_without_ext}.csv.gz", index=False)

def write_compact_tables(df, prefix):
    """Compact per-sample outputs for aggregating across a large run, instead of loading the full CSV.
    _calls: one row per read with the classification, so methods can be compared read by read.
    _leaks: every read called microbial, with all features, mapping, BLAST and tier columns."""
    write_table(df[["read_id", "sample_id", "classifier", "classification", "status"]].astype(str), f"{prefix}_calls")
    write_table(df[df["classification"] == "microbial"], f"{prefix}_leaks")

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

    parser.add_argument(
        "--classifier",
        dest="classifier",
        required=False,
        help="Name of the classifier which produced the input (charon, deacon, nohuman). Inferred from the file names if not given",
    )
    parser.add_argument(
        "--sample_id",
        dest="sample_id",
        required=False,
        help="Sample identifier. Taken from the input file name if not given",
    )
    parser.add_argument(
        "--contaminant_accs",
        dest="contaminant_accs",
        required=False,
        help="File of extra reference accessions (whitespace or comma separated) which are real contaminants rather than human",
    )
    parser.add_argument("--min_blast_cov", dest="min_blast_cov", type=float, default=80, help="Minimum query coverage (%%) of a BLAST human hit to count at leak tier 3")
    parser.add_argument("--min_blast_pident", dest="min_blast_pident", type=float, default=90, help="Minimum identity (%%) of a BLAST human hit to count at leak tier 3")
    parser.add_argument("--bitscore_ratio", dest="bitscore_ratio", type=float, default=1.1, help="Human bitscore must be at least this multiple of the best non-human bitscore to count at leak tier 4")

    args = parser.parse_args()
    other_accs.extend(load_contaminant_accs(args.contaminant_accs))
    sample_id = args.sample_id if args.sample_id else args.input.split("/")[-1].split(".")[0]

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

        combined_df = mapped_df.merge(blast_df, how="left", on="read_id")
        combined_df.to_csv("combined_df.csv")
        assert combined_df.shape[0] == mapped_df.shape[0], "The number of rows in the combined dataframe does not match the mapped dataframe."

        classifier, classifier_df = load_output(args.input, args.classifier)

        sys.stderr.write("COMBINE CLASSIFIER AND MAPPING DATA\n")
        old_size = classifier_df.shape[0]
        classifier_df = classifier_df.merge(combined_df, how="left", on="read_id")
        assert classifier_df.shape[0] == old_size, "The number of rows in the charon dataframe changed when combined with mapping dataframe."
        
        classifier_df['ref_name'] = classifier_df['ref_name'].fillna("")
        classifier_df["unmapped"] = classifier_df["mapped_length"].isna()
        for column in ["taxids", "names", "human_accs"]:
            classifier_df[column] = classifier_df[column].fillna("")
        classifier_df["blast_human"] = classifier_df["blast_human"].fillna(False).astype(bool)
        classifier_df["blast_nonhuman_run"] = classifier_df["blast_nonhuman_run"].fillna(False).astype(bool)

        classifier_df["sample_id"] = sample_id
        classifier_df["classifier"] = classifier
        classifier_df = add_leak_tiers(classifier_df, args.min_blast_cov, args.min_blast_pident, args.bitscore_ratio)
        classifier_df.to_csv(full_file, index=False)
    else:
        classifier_df = pd.read_csv(full_file, index_col=None)
        classifier = args.classifier if args.classifier else str(full_file).split("_full.csv")[0].split("_")[-1]
        classifier_df['classification'] = classifier_df['classification'].fillna("")
        classifier_df = add_leak_tiers(classifier_df, args.min_blast_cov, args.min_blast_pident, args.bitscore_ratio)

    write_compact_tables(classifier_df, args.prefix)

    summary = generate_summary(classifier_df, args.prefix, classifier)

    # Save to CSV
    fieldnames = ["sample_id"] + list(summary.keys())
    summary["sample_id"] = sample_id

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

