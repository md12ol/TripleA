seq1 = "CATGCTGACCAGTCTGGTCAAGCTCAAGGGATGAGCGGCATGGGAACCGGAACAACAACCGGATATGATGCTGGTGGCTACGGCGGTGAGCGCCAA"
seq2 = "CATCATGATCAGTCTGGTCAAGCTCAAGCGATGGGCGGCATGGGATCCGGATATGATGCTGGTGGCTACGGTGGTGAGCACCAC"

err_penalty = 5
gap_ext_penalty = 0.5
gap_open_penalty = 10

score = [[None for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
store = [[None for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]


def align(s1, s2, m, n):
    if score[m][n] is not None:  # Check storage
        return score[m][n], store[m][n]
    if m == 0 and n == 0:  # Base case
        return 0, ""
    if m == 0:  # Base case
        return 0, "2"
    if n == 0:  # Base case
        return 0, "1"
    elif s1[m - 1] == s2[n - 1]:
        val, ans = align(s1, s2, m - 1, n - 1)
        score[m][n] = val
        store[m][n] = ans + "_"
        return val, ans + "_"
    else:
        gap_in_2, ans_2 = align(s1, s2, m - 1, n)
        gap_in_1, ans_1 = align(s1, s2, m, n - 1)
        err, ans_err = align(s1, s2, m - 1, n - 1)
        if ans_2[-1] != "_":
            gap_in_2 += gap_ext_penalty
        else:
            gap_in_2 += gap_open_penalty
            pass
        if ans_1[-1] != "_":
            gap_in_1 += gap_ext_penalty
        else:
            gap_in_1 += gap_open_penalty
            pass
        err += err_penalty
        if gap_in_2 < gap_in_1 and gap_in_2 < err:
            score[m][n] = gap_in_2
            store[m][n] = ans_2 + "1"
            return gap_in_2, ans_2 + "1"
        if gap_in_1 < err:
            score[m][n] = gap_in_1
            store[m][n] = ans_1 + "2"
            return gap_in_1, ans_1 + "2"
        score[m][n] = err
        store[m][n] = ans_err + "_"
        return err, ans_err + "_"
    pass


def construct_alignment(s1, s2, answer):
    alignments = ["", ""]
    idx1, idx2 = 0, 0
    for idx, char in enumerate(answer):
        if char == "_":
            alignments[0] += s1[idx1]
            alignments[1] += s2[idx2]
            idx1 += 1
            idx2 += 1
            pass
        elif char == "2":
            alignments[0] += "_"
            alignments[1] += s2[idx2]
            idx2 += 1
            pass
        elif char == "1":
            alignments[0] += s1[idx1]
            alignments[1] += "_"
            idx1 += 1
            pass
        pass
    return alignments


def main():
    val, ans = align(seq1, seq2, len(seq1), len(seq2))
    print(val)
    alignments = construct_alignment(seq1, seq2, ans)
    print(alignments[0])
    print(alignments[1])
    pass


main()
