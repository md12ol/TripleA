

seq1 = "CATC TGATCA  CGTCTGGTCACCAG TCAAGGCGATGAGTCGCGGCATGGGAGGCAATCACCGGATATGA TGCTGGTGGCTACGGCGGTGAGCGCCAC"
# seq2 = "CATCATGATCAGTCTGGTCAAGCTCAAGCGATGGGCGGCATGGGATCCGGATATGATGCTGGTGGCTACGGTGGTGAGCACCAC"
seq2 = "CATGCTGACCAGTCTGGTCAAGCTCAAGGGATGAGCGGCATGGGAACCGGAACAACAACCGGATATGATGCTGGTGGCTACGGCGGTGAGCGCCAA"
seq3 = "CATCATGATCAGTCTGGTCAAGCTCAAGCGATGGGCGGCATGGGATCCGGATATGATGCTGGTGGCTACGGTGGTGAGCACCAC"
seq4 = "CATCATGATCAGTCGTCTGGTCA  AGTTCAAGG GATG G  GCGGCAT  G GG  A  ACCGGATATGAC GCTGGTGGCTACGGCGGTGAGCGCCAC"


def print_LCS(src1, src2, s1, s2):
    output = ["", ""]
    start_row = len(s1)
    start_col = len(s2)
    cur_row = start_row
    cur_col = start_col
    while cur_row > 0 and cur_col > 0:
        if src1[cur_row][cur_col] == cur_row - 1 and src2[cur_row][cur_col] == cur_col - 1:  # Match
            output[0] = s1[cur_row - 1] + output[0]
            output[1] = s2[cur_col - 1] + output[1]
            cur_row = cur_row - 1
            cur_col = cur_col - 1
            pass
        elif src1[cur_row][cur_col] == cur_row - 1:  # Gap in Seq2
            output[0] = s1[cur_row - 1] + output[0]
            output[1] = "_" + output[1]
            cur_row = cur_row - 1
            pass
        else:  # Gap in Seq1
            output[0] = "_" + output[0]
            output[1] = s2[cur_col - 1] + output[1]
            cur_col = cur_col - 1
            pass
        pass

    while cur_row > 0:  # Extra characters in Seq1
        output[0] = s1[cur_row - 1] + output[0]
        output[1] = "_" + output[1]
        cur_row -= 1
        pass

    while cur_col > 0:  # Extra characters in Seq2
        output[0] = "_" + output[0]
        output[1] = s2[cur_col - 1] + output[1]
        cur_col -= 1
        pass
    return output


def LCS(s1, s2):
    longest = [[0 for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    src_i = [[0 for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    src_j = [[0 for _ in range(len(s2) + 1)] for _ in range(len(s1) + 1)]
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            if s1[i - 1] == s2[j - 1]:
                longest[i][j] = longest[i-1][j-1] + 1
                src_i[i][j] = i-1
                src_j[i][j] = j-1
                pass
            elif longest[i-1][j] >= longest[i][j-1]:
                longest[i][j] = longest[i-1][j]
                src_i[i][j] = i-1
                src_j[i][j] = j
                pass
            else:
                longest[i][j] = longest[i][j - 1]
                src_i[i][j] = i
                src_j[i][j] = j - 1
                pass
            pass
        pass
    print(longest)
    output = print_LCS(src_i, src_j, s1, s2)
    print(output)
    return longest[len(s1)][len(s2)]
    pass


def main():
    sum = 0
    sum += LCS(seq1, seq2)
    sum += LCS(seq1, seq3)
    sum += LCS(seq1, seq4)
    print(sum)
    pass


main()