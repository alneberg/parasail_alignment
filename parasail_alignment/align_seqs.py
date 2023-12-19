import sys

import editdistance
import fastq
import miniFasta
import parasail


class NoAlignmentFoundError(Exception):
    pass


def find(char, string):
    """
    stolen from https://github.com/nanoporetech/sockeye

    Return iterator of indices for positions in a string
    corresponding to a target character

    :param char: Target character whose positions you want to locate in the string
    :type char: str
    :param string: String to search for the target character positions
    :type string: str
    :return: Indices in the string corresponding to the target character
    :rtype: iterator
    """
    for i in range(len(string)):
        if string[i] == char:
            yield i


COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def reverse_complementary(sequence: str):
    return "".join([COMP[base] for base in sequence[::-1]])


def setup_matrix():
    matrix = parasail.matrix_create("ACGTN", 5, -1)

    pointers = [4, 10, 16, 22, 24, 25, 26, 27]

    for i in pointers:
        matrix.pointer[0].matrix[i] = 1

    return matrix


class AlignmentResult:
    def __init__(
        self,
        query_seq,
        alignment_line,
        ref_seq,
        adapter1,
        adapter1_ed,
        barcode,
        adapter2,
        adapter2_ed,
    ):
        self.query_seq = query_seq
        self.alignment_string = alignment_line
        self.ref_seq = ref_seq
        self.adapter1 = adapter1
        self.adapter1_ed = str(adapter1_ed)
        self.barcode = barcode
        self.adapter2 = adapter2
        self.adapter2_ed = str(adapter2_ed)
        self.query_id = None
        self.ref_id = None

    def __str__(self):
        return f"Query:\t{self.query_seq}\nAlign:\t{self.alignment_string}\nRef:\t{self.ref_seq}"

    def add_ids(self, query_id, ref_id):
        self.query_id = query_id
        self.ref_id = ref_id

    def to_table_line(self, separator="\t"):
        return separator.join(
            [
                self.query_id,
                self.ref_id,
                self.query_seq,
                self.alignment_string,
                self.ref_seq,
                self.adapter1,
                self.adapter1_ed,
                self.barcode,
                self.adapter2,
                self.adapter2_ed,
            ]
        )

    def table_header_line(separator="\t"):
        return separator.join(
            [
                "query_id",
                "ref_id",
                "query_seq",
                "alignment_string",
                "ref_seq",
                "adapter1",
                "adapter1_ed",
                "barcode",
                "adapter2",
                "adapter2_ed",
            ]
        )


def align_single_seq(query_seq: str, ref_seq: str):
    """
    Aligns a single query sequence to a reference sequence using the Smith-Waterman algorithm.

    The query sequence corresponds to the ONT read, while the reference sequence corresponds to the Illumina adapter,
    including the 'N's for the barcode. The function uses the parasail library to perform the alignment and then
    extracts the barcode and adapter sequences from the alignment.

    Parameters:
    query_seq (str): The query sequence to be aligned.
    ref_seq (str): The reference sequence to which the query sequence is aligned.

    Returns:
    AlignmentResult: An object containing the aligned sequences, the alignment line, and the extracted adapter and barcode sequences.

    Raises:
    NoAlignmentFoundError: If no alignment could be found between the query and reference sequences.
    """
    matrix = setup_matrix()

    result_plus = parasail.sw_trace(query_seq, ref_seq, 2, 4, matrix)
    result_minus = parasail.sw_trace(
        reverse_complementary(query_seq), ref_seq, 2, 4, matrix
    )

    if result_plus.score > result_minus.score:
        result = result_plus
        orientation = "+"
    else:
        result = result_minus
        orientation = "-"

    ns_in_input = list(find("N", ref_seq))
    barcode_length = len(ns_in_input)
    adapter1_probe_seq = ref_seq[0 : ns_in_input[0]]
    adapter2_probe_seq = ref_seq[(ns_in_input[-1]+1) :]

    idxs = list(find("N", result.traceback.ref))

    # TODO handle ref_seq with no Ns (i.e. single index)
    if len(idxs) > 0:
        # The Ns in the probe successfully aligned to sequence
        bc_start = min(idxs)


            
        # The read1 adapter comprises the first part of the alignment
        adapter1 = result.traceback.query[0:bc_start]
        adapter1_ed = editdistance.eval(adapter1, adapter1_probe_seq)

        # The barcode + UMI sequences in the read correspond to the
        # positions of the aligned Ns in the probe sequence
        inferred_barcode = result.traceback.query[
            bc_start : (bc_start + barcode_length)
        ]

        adapter2 = result.traceback.query[(bc_start + barcode_length) :]
        adapter2_ed = editdistance.eval(adapter2, adapter2_probe_seq)

        aligned_length = len(result.traceback.query.replace('-', '')) 
        if orientation == '-':
            alignment_start = len(query_seq) - (result.end_query + 1 - aligned_length)
            alignment_end = len(query_seq) - (result.end_query + 1)
        else:
            alignment_end = result.end_query
            alignment_start = result.end_query + 1 - aligned_length

        return AlignmentResult(
            query_seq=result.traceback.query,
            alignment_line=result.traceback.comp,
            ref_seq=result.traceback.ref,
            adapter1=adapter1,
            adapter1_ed=adapter1_ed,
            barcode=inferred_barcode,
            adapter2=adapter2,
            adapter2_ed=adapter2_ed,
        )
    else:
        raise NoAlignmentFoundError("No alignment found")


def align_multiple_seq(query_file: str, ref_file: str):
    # Populate a list to enable the reuse of the iterator, otherwise it will be empty after the first iteration
    ref_seqs = [ref_seq for ref_seq in miniFasta.read(ref_file)]
    ref_seq_ids = [ref_seq.getHead().split(" ")[0][1:] for ref_seq in ref_seqs]

    print(AlignmentResult.table_header_line())

    for query_seq in fastq.read(query_file):
        for ref_seq, ref_seq_id in zip(ref_seqs, ref_seq_ids):
            try:
                result = align_single_seq(query_seq.getSeq(), ref_seq.getSeq())
                result.add_ids(
                    query_seq.getHead().split(" ")[0], ref_seq_id
                )
                print(result.to_table_line())
            except NoAlignmentFoundError:
                sys.stderr.write(
                    f"No alignment found for {query_seq.getHead().split(' ')[0]} and {ref_seq.getHead().split(' ')[0]}"
                )


def align_single_wrapper(args):
    print(align_single_seq(args.query_seq, args.ref_seq))


def align_multiple_wrapper(args):
    align_multiple_seq(args.query_seq_file, args.ref_seq_file)
