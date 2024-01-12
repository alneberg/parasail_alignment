"""Outline:
Fetch all barcodes
insert dash instead of -NNN-
For each barcode:
    run minimap2
    filter out perfect hits
    output number of perfect hits, i5, i7, both (both in correct orientation?)
    histogram over barcode lengths (one for each i5 and i7)
    
    entropy calculation over the entire read (only works if both i5 and i7 has perfect match)
    entropy calculation over the adapter + barcode + adapter (one for each i5 and i7)
"""
import argparse
import itertools
import logging
import os
import subprocess
import uuid
from collections import deque

import numpy as np
import pandas as pd
import yaml
from scipy.stats import entropy

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("anglerfish")


class Adaptor:
    """From anglerfish"""

    def __init__(self, adaptors, delim, adaptor, i7_index=None, i5_index=None):
        self.i5 = AdaptorPart(adaptors[adaptor]["i5"], adaptor, delim, i5_index)
        self.i7 = AdaptorPart(adaptors[adaptor]["i7"], adaptor, delim, i7_index)
        self.name = f"{adaptor}"
        self.delim = delim

    def get_i5_mask(self, insert_Ns=True):
        return self.i5.get_mask(insert_Ns)

    def get_i7_mask(self, insert_Ns=True):
        return self.i7.get_mask(insert_Ns)

    def get_fastastring(self, insert_Ns=True):
        fasta_i5 = f">{self.name}_i5\n{self.get_i5_mask(insert_Ns)}\n"
        fasta_i7 = f">{self.name}_i7\n{self.get_i7_mask(insert_Ns)}\n"
        return fasta_i5 + fasta_i7


class AdaptorPart:
    def __init__(self, sequence, name, delim, index):
        self.sequence = sequence
        self.name = name
        self.delim = delim
        self.index = index

    def has_index(self):
        return self.sequence.find(self.delim) > -1

    def len_before_index(self):
        return self.sequence.find(self.delim)

    def len_after_index(self):
        return len(self.sequence) - self.sequence.find(self.delim) - len(self.delim)

    def get_mask(self, insert_Ns):
        if self.has_index():
            if not insert_Ns:
                return self.sequence.replace(self.delim, "")
            else:
                return self.sequence.replace(self.delim, "N" * len(self.index))
        else:
            return self.sequence


# Fetch all adaptors
def load_adaptors(filename):
    adaptors_raw = []
    with open(filename) as f:
        adaptors_raw = yaml.safe_load(f)

    adaptors = []
    for adaptor in adaptors_raw:
        adaptors.append(
            Adaptor(adaptors_raw, "-NNN-", adaptor, i7_index=None, i5_index=None)
        )

    return adaptors


def run_minimap2(fastq_in, indexfile, output_paf, threads, minimap_B=4):
    """
    Runs Minimap2
    """
    cmd = [
        "minimap2",
        "--cs",
        "-m8",
        "-k",
        "10",
        "-w",
        "5",
        "-A",
        "6",
        "-B",
        str(minimap_B),
        "-c",
        "-t",
        str(threads),
        indexfile,
        fastq_in,
    ]

    run_log = f"{output_paf}.log"
    with open(output_paf, "ab") as ofile, open(run_log, "ab") as log_file:
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log_file)
        subprocess.run("sort", stdin=p1.stdout, stdout=ofile, check=True)


def parse_paf_lines(paf, min_qual=10, complex_identifier=False):
    """
    Read and parse one paf alignment lines.
    Returns a dict with the import values for later use
    """
    entries = {}
    with open(paf) as paf:
        for paf_line in paf:
            aln = paf_line.split()
            try:
                # TODO: objectify this
                entry = {
                    "read": aln[0],
                    "adapter": aln[5],
                    "rlen": int(aln[1]),  # read length
                    "rstart": int(aln[2]),  # start alignment on read
                    "rend": int(aln[3]),  # end alignment on read
                    "strand": aln[4],
                    "cg": aln[-2],  # cigar string
                    "cs": aln[-1],  # cs string
                    "q": int(aln[11]),  # Q score
                    "iseq": None,
                    "sample": None,
                }
                read = entry["read"]

                if complex_identifier:
                    i5_or_i7 = entry["adapter"].split("_")[-1]
                    if entry["strand"] == "+":
                        strand_str = "positive"
                    else:
                        strand_str = "negative"
                    ix = f"{read}_{i5_or_i7}_{strand_str}"
                else:
                    ix = read
            except IndexError:
                log.debug(f"Could not find all paf columns: {read}")
                continue

            if entry["q"] < min_qual:
                log.debug(f"Low quality alignment: {read}")
                continue

            if complex_identifier:
                entries[ix] = entry
            else:
                if ix in entries.keys():
                    entries[ix].append(entry)
                else:
                    entries[ix] = [entry]

    return entries


# Rolling window implementation
def window(seq, n=3):
    win = deque(maxlen=n)
    for char in seq:
        win.append(char)
        if len(win) >= n:
            yield tuple(win)


# Each sequence gives matrix of position and k-mer occurence.
def seq_to_count_matrix(dna_seq, kmer_positions, max=50, k=2):
    """max is the length of the prefix that will be considered"""
    if len(dna_seq) < max:
        max = len(dna_seq)

    matrix = np.zeros((4**k, max - k + 1))
    for i, word in enumerate(window(dna_seq[:max], n=k)):
        matrix[kmer_positions[word], i] += 1
    return matrix


def count_for_seqs(seqs, kmer_positions, insert_length, k=2):
    seqs = [seq for seq in seqs if len(seq) == insert_length]
    matrices = [seq_to_count_matrix(seq, kmer_positions, max=50, k=k) for seq in seqs]
    stacked = np.stack(matrices)
    sum_a = np.sum(stacked, axis=0)
    even_dist = np.ones(4**k) / (4**k)
    entropies = [entropy(a, even_dist) for a in sum_a.T]

    return entropies


def extract_inserts_from_df(df):
    mim_re_cs = r"^cs:Z::[1-9][0-9]*\+([a,c,t,g]*):[1-9][0-9]*$"
    return df.cs.str.extract(mim_re_cs, expand=False).str.upper().tolist()


def calculate_relative_entropy(df_good_hits, args, insert_length):
    all_kmers = itertools.product("ACTG", repeat=args.kmer_length)
    kmer_positions = dict(
        (kmer, position)
        for kmer, position in zip(all_kmers, range(4**args.kmer_length))
    )
    seqs = extract_inserts_from_df(df_good_hits)
    entropies = count_for_seqs(seqs, kmer_positions, insert_length, k=args.kmer_length)
    return entropies


def main(args):
    run_uuid = str(uuid.uuid4())
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        log.info(f"Output directory {args.outdir} already exists")
        if not args.no_overwrite:
            log.error(
                f"Output directory {args.outdir} already exists, please use --no_overwrite to continue"
            )
            exit(1)
        else:
            pass

    log.info(f" arguments {vars(args)}")
    log.info(f" run uuid {run_uuid}")

    adaptors = load_adaptors(args.adaptors_file)
    alignments = []

    for adaptor in adaptors:
        adaptor_name = adaptor.name

        # Align
        aln_path = os.path.join(args.outdir, f"{adaptor_name}.paf")
        alignments.append((adaptor, aln_path))
        if os.path.exists(aln_path) and args.no_overwrite:
            log.info(f"Skipping {adaptor_name} as alignment already exists")
            continue
        adaptor_path = os.path.join(args.outdir, f"{adaptor_name}.fasta")
        with open(adaptor_path, "w") as f:
            f.write(adaptor.get_fastastring(insert_Ns=False))

        log.info(f"Aligning {adaptor_name}")
        run_minimap2(args.fastq, adaptor_path, aln_path, args.threads, args.minimap_B)

    # Parse alignments
    # Read all paf:s into memory with pandas, should be possible to optimize
    # by handling one line at a time for each file instead.

    entries = {}
    adaptors_included = []
    for adaptor, aln_path in alignments:
        log.info(f"Parsing {adaptor.name}")
        aln_dict = parse_paf_lines(aln_path, complex_identifier=True)
        df = pd.DataFrame.from_dict(aln_dict, orient="index")
        nr_good_hits = {}
        for adaptor_end_name, adaptor_end in zip(
            ["i5", "i7"], [adaptor.i5, adaptor.i7]
        ):
            if adaptor_end.has_index():
                # Match Insert Match = mim
                # The cs string filter is quite strict, requiring 10+ perfect match before insert and 10+ perfect match after insert
                # The cg string allows for mismatches within the matching strings but no insertions or deletions
                # All cg string matches are also cs string matches (subset) but not vice versa
                mim_re_cs = r"^cs:Z::[1-9][0-9]*\+([a,c,t,g]*):[1-9][0-9]*$"  # Match Insert Match = mim
                mim_re_cg = r"^cg:Z:([0-9]*)M([0-9]*)I([0-9]*)M$"
                df_mim = df[df.cs.str.match(mim_re_cs)]

                # Extract the match lengths
                match_col_df = df_mim.cg.str.extract(mim_re_cg).rename(
                    {0: "match_1_len", 1: "insert_len", 2: "match_2_len"}, axis=1
                )
                match_col_df = match_col_df.astype(
                    {
                        "match_1_len": "int32",
                        "insert_len": "int32",
                        "match_2_len": "int32",
                    }
                )

                df_mim.loc[match_col_df.index, match_col_df.columns] = match_col_df

                # Alignment thresholds
                before_thres = round(
                    adaptor_end.len_before_index() * args.good_hit_threshold
                )
                after_thres = round(
                    adaptor_end.len_after_index() * args.good_hit_threshold
                )
                insert_thres_low = args.insert_thres_low
                insert_thres_high = args.insert_thres_high

                requirements = (
                    (df_mim["match_1_len"] >= (before_thres))
                    & (df_mim["insert_len"] >= insert_thres_low)
                    & (df_mim["insert_len"] <= insert_thres_high)
                    & (df_mim["match_2_len"] >= (after_thres))
                )
                df_good_hits = df_mim[requirements]

                median_insert_length = df_good_hits["insert_len"].median()
                insert_lengths = df_good_hits["insert_len"].value_counts()
            else:
                m_re_cs = r"^cs:Z::([1-9][0-9]*)$"
                df_good_hits = df[df.cg.str.match(m_re_cs)]
                match_col_df = df_good_hits.cg.str.extract(m_re_cs).rename(
                    {0: "match_1_len"}, axis=1
                )
                match_col_df = match_col_df.astype({"match_1_len": "int32"})

                df_good_hits.loc[
                    match_col_df.index, match_col_df.columns
                ] = match_col_df

                thres = round(
                    (adaptor_end.len_before_index() + adaptor_end.len_after_index())
                    * args.good_hit_threshold
                )
                df_good_hits = df_good_hits[df_good_hits["match_1_len"] >= thres]

                median_insert_length = None
                insert_lengths = None

            # Filter multiple hits per read and adaptor end
            df_good_hits = df_good_hits.sort_values(
                by=["read", "match_1_len"], ascending=False
            ).drop_duplicates(subset=["read", "adapter"], keep="first")

            if adaptor.name not in entries.keys():
                entries[adaptor.name] = {}
            entries[adaptor.name][adaptor_end_name] = df_good_hits

            nr_good_hits[adaptor_end_name] = len(df_good_hits)
            log.info(
                f"{adaptor.name}:{adaptor_end_name} had {len(df_good_hits)} good hits."
            )

        if min(nr_good_hits.values()) >= args.min_hits_per_adaptor:
            log.info(f"Adaptor {adaptor.name} is included in the analysis")
            adaptors_included.append(adaptor)
        else:
            log.info(f"Adaptor {adaptor.name} is excluded from the analysis")

    # Print insert length info for adaptor types included in the analysis
    for adaptor in adaptors_included:
        for adaptor_end_name, adaptor_end in zip(
            ["i5", "i7"], [adaptor.i5, adaptor.i7]
        ):
            df_good_hits = entries[adaptor.name][adaptor_end_name]
            if adaptor_end.has_index():
                median_insert_length = df_good_hits["insert_len"].median()
                if median_insert_length > args.umi_threshold:
                    entropies = calculate_relative_entropy(
                        df_good_hits, args, median_insert_length
                    )
                    log.info(
                        f"{adaptor.name}:{adaptor_end_name} had entropy {entropies}"
                    )
                insert_lengths = df_good_hits["insert_len"].value_counts()
                log.info(
                    f"{adaptor.name}:{adaptor_end_name} had {len(df_good_hits)} good hits with median insert length {median_insert_length}"
                )
                log.info(insert_lengths[sorted(insert_lengths.index)])
            else:
                median_insert_length = None
                insert_lengths = None
                log.info(
                    f"{adaptor.name}:{adaptor_end_name} had {len(df_good_hits)} good hits (no insert length since no index)"
                )


#    entropy calculation over the entire read (only works if both i5 and i7 has perfect match)
#    entropy calculation over the adapter + barcode + adapter (one for each i5 and i7)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--fastq", help="Fastq file to align", required=True)
    parser.add_argument(
        "-a", "--adaptors_file", help="YAML file with adaptors", required=True
    )
    parser.add_argument("-o", "--outdir", help="Output directory", required=True)
    parser.add_argument(
        "-t", "--threads", help="Number of threads", default=1, type=int
    )
    parser.add_argument(
        "-n",
        "--no_overwrite",
        help="Do not overwrite existing alignments",
        action="store_true",
    )
    parser.add_argument(
        "-g",
        "--good_hit_threshold",
        help="Fraction of bases before and after index insert "
        "required to match perfectly for a hit to be considered a good hit. Default=0.9",
        default=0.9,
        type=float,
    )
    parser.add_argument(
        "-i",
        "--insert_thres_low",
        help="Lower threshold for insert length, with value included",
        default=4,
        type=int,
    )
    parser.add_argument(
        "-j",
        "--insert_thres_high",
        help="Upper threshold for insert length, with value included",
        default=30,
        type=int,
    )
    parser.add_argument(
        "-B",
        "--minimap_B",
        help="Minimap2 -B parameter, mismatch penalty",
        default=4,
        type=int,
    )
    parser.add_argument(
        "-m",
        "--min_hits_per_adaptor",
        help="Minimum number of good hits for an adaptor to be included in the analysis",
        default=50,
        type=int,
    )
    parser.add_argument(
        "-u",
        "--umi_threshold",
        help="Minimum number of bases in insert to perform entropy calculation",
        default=11,
        type=float,
    )
    parser.add_argument(
        "-k",
        "--kmer_length",
        help="Length of k-mers to use for entropy calculation",
        default=2,
        type=int,
    )
    args = parser.parse_args()

    main(args)
