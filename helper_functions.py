def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    n, m = len(seq1), len(seq2)

    # H = score table, T = tracecback (1 = diagonal, 2 = up, 3 = left)
    H = [[0,0] * (m + 1) for _ in range(n + 1)]
    T = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        H[i][0] = H[i - 1][0] + scoring_function(seq1[i - 1], "-")
        T[i][0] = 2
    for j in range(1, m + 1):
        H[0][j] = H[0][j - 1] + scoring_function("-", seq2[j - 1])
        T[0][j] = 3
    
    raise NotImplementedError()


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

from Bio.Align import substitution_matrices

BLOSUM62 = substitution-matrices.load("BLOSUM62")
GAP_SCORE = -4

def scoring_function_blosum62(aa_i, aa_j):
    if aa_i = "-" or aa_j = "-":
        return GAP_SCORE
    return float(BLOSUM62[aa_i,aa_j])
