"""CLI interface for parasail_alignment"""
import argparse

from parasail_alignment.align_seqs import NoAlignmentFoundError, align_single_seq


def main():
    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers()

    align_single_parser = subparser.add_parser(
        "single", dest=align_single_seq, help="Align a single sequence"
    )

    align_single_parser.add_argument("query_seq")
    align_single_parser.add_argument("ref_seq", help="the one with the N:s")
    align_single_parser.set_defaults(func=align_single_seq)

    align_multiple_parser = parser.add_subparser(
        "multiple", help="Align multiple sequences"
    )

    align_multiple_parser.add_argument(
        "query_seq_file", help="File with ONT reads in fastq or fastq.gz format"
    )
    align_multiple_parser.add_argument(
        "ref_seq_file",
        help="File with illumina adapter sequences (the ones with the N:s) in fasta",
    )
    align_multiple_parser.set_defaults(func=align_single_seq)

    args = parser.parse_args()
    try:
        args.func(args)
    except NoAlignmentFoundError:
        print("No alignment found")
