"""CLI interface for parasail_alignment"""
import argparse

from parasail_alignment.align_seqs import NoAlignmentFoundError, align_single_seq


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("query_seq")
    parser.add_argument("ref_seq", help="the one with the N:s")

    args = parser.parse_args()

    try:
        print(align_single_seq(args.query_seq, args.ref_seq))
    except NoAlignmentFoundError:
        print("No alignment found")
