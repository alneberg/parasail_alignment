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


def setup_matrix():
    matrix = parasail.matrix_create("ACGTN", 5, -1)

    pointers = [4, 10, 16, 22, 24, 25, 26, 27]

    for i in pointers:
        matrix.pointer[0].matrix[i] = 1

    return matrix


class AlignmentResult:
    def __init__(self, query_seq, alignment_line, ref_seq, adapter1, barcode, adapter2):
        self.query_seq = query_seq
        self.alignment_string = alignment_line
        self.ref_seq = ref_seq
        self.adapter1 = adapter1
        self.barcode = barcode
        self.adapter2 = adapter2

    def __str__(self):
        return f"Query:\t{self.query_seq}\nAlign:\t{self.alignment_string}\nRef:\t{self.ref_seq}"


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

    result = parasail.sw_trace(query_seq, ref_seq, 2, 4, matrix)

    ns_in_input = list(find("N", ref_seq))
    barcode_length = len(ns_in_input)

    idxs = list(find("N", result.traceback.ref))

    # Also stolen from sockeye
    if len(idxs) > 0:
        # The Ns in the probe successfully aligned to sequence
        bc_start = min(idxs)

        # The read1 adapter comprises the first part of the alignment
        adapter1 = result.traceback.query[0:bc_start]
        # adapter1_ed = edit_distance(adapter1, adapter1_probe_seq)

        # The barcode + UMI sequences in the read correspond to the
        # positions of the aligned Ns in the probe sequence
        inferred_barcode = result.traceback.query[
            bc_start : (bc_start + barcode_length)
        ]

        adapter2 = result.traceback.query[(bc_start + barcode_length) :]

        print(
            AlignmentResult(
                query_seq=result.traceback.query,
                alignment_line=result.traceback.comp,
                ref_seq=result.traceback.ref,
                adapter1=adapter1,
                barcode=inferred_barcode,
                adapter2=adapter2,
            )
        )
    else:
        raise NoAlignmentFoundError("No alignment found")


def align_multiple_seq(query_file: str, ref_file: str):
    # Populate a list to enable the reuse of the iterator, otherwise it will be empty after the first iteration
    ref_seqs = [ref_seq for ref_seq in miniFasta.read(ref_file)]

    for query_seq in fastq.read(query_file):
        for ref_seq in ref_seqs:
            try:
                print(align_single_seq(query_seq.getSeq(), ref_seq.getSeq()))
            except NoAlignmentFoundError:
                print(
                    f"No alignment found for {query_seq.getHead().split(' ')[0]} and {ref_seq.getHead().split(' ')[0]}"
                )


def align_single_wrapper(args):
    align_single_seq(args.query_seq, args.ref_seq)


def align_multiple_wrapper(args):
    align_multiple_seq(args.query_seq_file, args.ref_seq_file)
