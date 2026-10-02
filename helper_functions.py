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
    H = [[0.0] * (m + 1) for _ in range(n + 1)]
    T = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        H[i][0] = H[i - 1][0] + scoring_function(seq1[i - 1], "-")
        T[i][0] = 2
    for j in range(1, m + 1):
        H[0][j] = H[0][j - 1] + scoring_function("-", seq2[j - 1])
        T[0][j] = 3

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = H[i - 1][j - 1] + scoring_function(seq1[i - 1], seq2[j - 1])
            up = H[i - 1][j] + scoring_function(seq1[i - 1], "-")
            left = H[i][j - 1] + scoring_function("-", seq2[j - 1])
            if diag >= up and diag>= left:
                H[i][j], T[i][j] = diag, 1
            elif up >= left:
                H[i][j], T[i][j] = up, 2
            else:
                H[i][j], T[i][j] = left, 3
    
    out1, out2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        t = T[i][j]
        if t == 1:
            out1.append(seq1[i - 1]); out2.append(seq2[j - 1]); i -= 1; j -= 1
        elif t == 2:
            out1.append(seq1[i - 1]); out2.append("-"); i -= 1
        else:
            out1.append("-"); out2.append(seq2[j - 1]); j -= 1
            
    return "".join(reversed(out1)), "".join(reversed(out2)), float(H[n][m])


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

    n, m = len(seq1), len(seq2)

    gap1 = [scoring_function(a, "-") for a in seq1]
    gap2 = [scoring_function("-", b) for b in seq2]
    pair = {(a, b): scoring_dunction(a, b) for a in set(seq1) for b in set(seq2)}

    H = [[0.0] * (m + 1) for _ in range(n + 1)]
    T = [[0] * (m + 1) for _ in range(n + 1)]

    best, best_i, best_j = 0.0, 0, 0
    
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load("BLOSUM62")
GAP_SCORE = -4

def scoring_function_blosum62(aa_i, aa_j):
    if aa_i == "-" or aa_j == "-":
        return GAP_SCORE
    return float(BLOSUM62[aa_i,aa_j])
