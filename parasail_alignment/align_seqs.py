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


def align_single_seq(query_seq, ref_seq):
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

        AlignmentResult(
            query=result.traceback.query,
            alignment_line=result.traceback.comp,
            reference=result.traceback.ref,
            adapter1=adapter1,
            barcode=inferred_barcode,
            adapter2=adapter2,
        )
    else:
        raise NoAlignmentFoundError("No alignment found")


class AlignmentResult:
    def __init__(self, query_seq, alignment_line, ref_seq, adapter1, barcode, adapter2):
        self.query_seq = query_seq
        self.alignment_string = alignment_line
        self.ref_seq = ref_seq
        self.adapter1 = adapter1
        self.barcode = barcode
        self.adapter2 = adapter2

    def __str__(self):
        return f"Query:\t{self.query_seq}\nAlignment:\t{self.alignment_string}\nReference:\t{self.ref_seq}"
