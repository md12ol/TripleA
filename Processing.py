import math
import os
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

inp = "./Output/"
outp = "./Figures/"
finame = "exp.dat"
samps = 50
precision = 6
col_width = 8 + precision
sda_states = 20
sda_chars = 4
alphabet = ['G', 'C', 'A', 'T']
cull_every = 5


def writeStat(data: [], out, lower_better):
    mean = float(np.mean(data))
    mean = round(mean, precision)
    std = float(np.std(data, ddof=0))
    std = round(std, precision)  # Population standard deviation
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    diff = round(diff, precision)
    if lower_better:
        maxima = float(min(data))
        pass
    else:
        maxima = float(max(data))
        pass
    maxima = round(maxima, precision)
    out.write(str(maxima).ljust(col_width))
    out.write(str(mean).ljust(col_width))
    out.write(str(std).ljust(col_width))
    out.write(str(diff).ljust(col_width))
    return mean, maxima


def make_table(many_data: [], best_run, exp_info: [], fname: str, minimizing: bool):
    with open(fname, "w") as f:
        f.write("EXP".ljust(col_width))
        f.write("Parameters".ljust(3 * col_width))
        f.write("Best Run".ljust(col_width))
        f.write("Best Fit".ljust(col_width))
        f.write("Mean".ljust(col_width))
        f.write("SD".ljust(col_width))
        f.write("95% CI".ljust(col_width))
        f.write("\n")
        for di, data in enumerate(many_data):
            f.write(str("EXP" + str(di + 1)).ljust(col_width))
            f.write(exp_info[di].ljust(3 * col_width))
            f.write(str("Run " + str(best_run[di])).ljust(col_width))
            # f.write(str("Run " + str(data[0][0])).ljust(col_width))
            writeStat(data, f, minimizing)
            f.write("\n")
            pass
        pass
    pass


def box_plot(bp, num_splits: int, split_info: []):
    color = "#2e294e"
    for whisker in bp['whiskers']:
        whisker.set(color=color, linewidth=1.5)
        pass

    for cap in bp['caps']:
        cap.set(color=color, linewidth=1.5)
        pass

    for median in bp['medians']:
        median.set(color=color, linewidth=1.5)
        pass

    for flier in bp['fliers']:
        flier.set(marker='*', markeredgecolor=color, alpha=0.5, markersize=6)
        pass

    for info in split_info:
        for idx in info[0]:
            bp['boxes'][idx].set(facecolor=info[1][idx % len(info[1])])
            pass
        pass
    pass


def calc(data):
    mean = float(np.mean(data))
    std = float(np.std(data, ddof=0))
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    print(str(mean) + "+-" + str(diff))
    pass


def combine(dir: str, start: str, end: str):
    if os.path.exists(dir + start + end):
        os.remove(dir + start + end)
        pass
    out_file = open(dir + start + end, "a")
    for idx in range(1, 31):
        with open(dir + start + str(idx).zfill(2) + end, "r") as f:
            out_file.write(f.read())
            pass
        pass
    out_file.close()
    pass


def get_data(dir_path: str):
    match_fits = []
    novelty_fits = []
    SDAs = []
    seqs = []
    with open(dir_path + finame) as f:
        lines = f.readlines()
        next_SDA = False
        SDA = []
        for line in lines:
            if line.__contains__("best fitness"):
                line = line.rstrip()
                match_fits.append(int(line.split(" is ")[1].split(" matches ")[0]))
                novelty_fits.append(int(line.split(" novelty ")[1].rstrip()))
                pass
            elif line.__contains__(str("Best Match")):
                line = line.rstrip()
                line = line.split(" ")
                seqs.append(line[-1])
                pass
            elif line.__contains__("SDA"):
                next_SDA = True
                pass
            elif line.__contains__("Fitness Values"):
                next_SDA = False
                SDAs.append(SDA)
                pass
            elif next_SDA:
                SDA.append(line)
                pass
            pass
        pass

    # run number, fitness, profileS, dnaS, edge count
    if len(match_fits) != samps:
        print("ERROR in match fits: " + dir_path)
        pass
    if len(novelty_fits) != samps:
        print("ERROR in novelty fits: " + dir_path)
        pass
    if len(SDAs) != samps:
        print("ERROR in SDAs: " + dir_path)
        pass
    if len(seqs) != samps:
        print("ERROR in seqs: " + dir_path)
        pass

    data = [[i + 1, match_fits[i], novelty_fits[i], SDAs[i], seqs[i]] for i in range(samps)]
    data.sort(key=itemgetter(2))  # Ascending
    data.reverse()
    data.sort(key=itemgetter(1))  # Ascending
    data.reverse()
    return data


def str_to_list(orig: str):
    rtn = []
    for c in orig:
        rtn.append(int(c))
        pass
    return rtn


def cmpr(seq1, seq2):
    similar_str = ''
    for idx in range(len(seq2)):
        if seq1[idx] == seq2[idx]:
            similar_str += 'X'
        else:
            similar_str += '-'
            pass
        pass
    return similar_str


def print_best_info(path: str, info: str, dat: [], true_seq: []):
    with open(path, "w") as f:
        f.write(info + "\n")
        f.write("Run number: " + str(dat[0]) + "\n")
        f.write("With Fitness: " + str(dat[1]) + "\n")
        test_seq = str_to_list(dat[3])
        similar_str = ''
        for idx in range(len(true_seq)):
            if test_seq[idx] == true_seq[idx]:
                similar_str += 'X'
            else:
                similar_str += '-'
                pass
            pass
        f.write(int_to_DNA(test_seq) + '\n')
        f.write(similar_str + '\n')
        f.write(int_to_DNA(true_seq) + '\n')
        f.writelines(dat[2])
        pass
    pass


def get_char_freqs(sequences: [], num_chars: int):
    freqs = [['G', 'C', 'A', 'T']]
    for seq in sequences:
        dat = [0 for _ in range(num_chars)]
        for c in seq:
            dat[c] += 1
            pass
        freqs.append(dat)
        pass
    return freqs


def gen_sequences(path: str):
    seq_sets = []
    one_set = []
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            if line.__contains__('SET'):
                if len(one_set) > 0:
                    seq_sets.append(one_set)
                    pass
                one_set = []
                pass
            elif line.__contains__(">"):
                pass
            else:
                one_set.append(line)
                pass
            pass
        pass
    seq_sets.append(one_set)
    return seq_sets


def int_to_DNA(vals: []):
    rtn = ''
    for val in vals:
        if val == 0:
            rtn += 'G'
        elif val == 1:
            rtn += 'C'
        elif val == 2:
            rtn += 'A'
        elif val == 3:
            rtn += 'T'
            pass
        pass
    return rtn


def DNA_to_int(seq: str):
    rtn = []
    for c in seq:
        if c == 'g' or c == 'G':
            rtn.append(0)
        elif c == 'c' or c == 'C':
            rtn.append(1)
        elif c == 'a' or c == 'A':
            rtn.append(2)
        elif c == 't' or c == 'T':
            rtn.append(3)
            pass
        pass
    return rtn


def process_sda(batch: list, sda_states: int, sda_chars: int):
    sdas, outputs, fits = [], [], []
    one_sda = []
    for line in batch:
        if line.__contains__("SDA"):
            fit_start = True
            pass
        elif fit_start:
            fit_start = False
            sda_start = True
            fits.append(int(line.rstrip().split(": ")[1]))
            pass
        elif sda_start and not line.__contains__("-"):
            sda_start = False
            line = line.rstrip().split(" ")
            outputs.append(line)
            sdas.append(get_sda(one_sda, sda_states, sda_chars))
            one_sda = []
        elif sda_start:
            one_sda.append(str(line))
            pass
        pass
    return sdas, outputs, fits


def get_sda(lines: list[str], sda_states: int, sda_chars: int):
    init_state = int(lines[0].rstrip().split("<")[0].rstrip())
    init_char = int(lines[0].rstrip().split("-")[1].strip())
    lines = lines[1:]

    sda_trans = [[] for _ in range(sda_states)]
    sda_resps = [[] for _ in range(sda_states)]
    for state in range(sda_states):
        for char in range(sda_chars):
            one_resp = []
            sda_trans[state].append(int(lines[0].split(">")[1].lstrip().split(" ")[0]))
            line = lines[0].rstrip().split(">")[1].lstrip().split(" ")
            one_resp.append(int(line[2]))
            if len(line) == 5:
                one_resp.append(int(line[3]))
                pass
            sda_resps[state].append(one_resp)
            lines = lines[1:]
            pass
        pass
    return init_state, init_char, sda_trans, sda_resps


def make_heatmap(sda_infos: list, num_states, num_chars, out_path, num_gens, run_num, seq, exp):
    count = [[0 for _ in range(num_states)] for _ in range(num_states * num_chars)]
    for sda in sda_infos:
        transitions = sda[2]
        for st, trans in enumerate(transitions):
            for ch, tr in enumerate(trans):
                row = num_chars * st + ch
                col = tr
                count[row][col] += 1
                pass
            pass
        pass
    # out_path = out_path + "Run" + str(run_num).zfill(2) + "/"
    # if not os.path.exists(out_path):
    #     os.makedirs(out_path)
    #     pass

    plt.style.use("seaborn-v0_8")
    cmap = LinearSegmentedColormap.from_list("thing", ["#2e294e", "#1b998b", "#f46036"], 500)

    f = plt.figure()
    f.set_figheight(7)
    f.set_figwidth(7)
    plot = f.add_subplot(111)
    hm = plot.imshow(count, aspect='auto', cmap=cmap)
    # hm = plt.pcolor(count, aspect='auto', cmap=cmap, norm=norm)

    ylbls = []
    ylocs = []
    for s in range(num_states):
        for idx, c in enumerate(alphabet):
            ylbls.append(c)
            ylocs.append(s * num_chars + idx)
            pass
        pass

    plot.tick_params(axis="y", direction="in", pad=15)
    plot.set_yticks([y * num_chars + 1.5 for y in range(num_states)], (s + 1 for s in range(num_states)), minor=False,
                    fontsize=12)
    plot.set_yticks(ylocs, ylbls, minor=True, fontsize=7)
    plot.set_xticks([y for y in range(num_states)], (y + 1 for y in range(num_states)), fontsize=12)
    f.suptitle("Number of SDAs with the Same Transition", fontsize=14)
    plot.set_title("Set " + str(seq) + " - Experiment " + str(exp), fontsize=12)

    plot.set_xlabel("To State", fontsize=12)
    plot.set_ylabel("From State and the Character Driving Transition", fontsize=12)
    cbar = plot.figure.colorbar(hm, ax=plot, pad=0.01)
    plot.grid(None)
    for y in range(4 * num_states):
        plot.axhline(y=y + 0.5, color="#FFFFFF", linestyle="-", linewidth="0.25")
        pass
    for x in range(num_states):
        plot.axvline(x=x + 0.5, color="#FFFFFF", linestyle="-", linewidth="0.25")
        pass
    for y in range(3, 4 * num_states, 4):
        plot.axhline(y=y + 0.5, color="#FFFFFF", linestyle="-", linewidth="1")
        pass
    f.subplots_adjust(bottom=0.07, top=0.92, left=0.10, right=1.04)
    f.savefig(out_path + "Run" + str(run_num) + "Gen" + str(num_gens).zfill(8) + ".png", dpi=300)
    plt.close()
    pass


def make_convergence_plot(folder: str, run_num: int, out_path, num_gens, cull_every, seq, exp_num, exp_desc, seq_len):
    if isinstance(num_gens, list):
        print("ERROR!")
        pass

    means = []
    bests = []
    novs = []
    with open(folder + "run" + str(run_num).zfill(2) + ".dat") as f:
        lines = f.readlines()
        for line in lines[1:]:
            the_line = line
            if line.__contains__(" 1000"):
                new_line = line.split("1000")[1]
                the_line = line.split("1000")[0] + "    1000   " + new_line
                pass
            means.append(float(the_line.split()[2]))
            bests.append(int(the_line.split()[5]))
            novs.append(int(the_line.split()[6]))
            pass
        pass

    out_path = out_path + "ConvergenceSeq" + str(seq) + "Exp" + str(exp_num) + "Run" + str(run_num).zfill(2)

    plt.style.use("seaborn-v0_8")
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    # f = plt.figure()
    # f.set_figheight(4.5)
    # f.set_figwidth(10)
    f, plot = plt.subplots()
    f.set_figheight(5)
    f.set_figwidth(9)

    xs = [i for i in range(0, num_gens + 1, 10000)]

    xlocs = [x for x in range(0, num_gens, cull_every * 10000)]
    for x in xlocs:
        plot.axvline(x, color='gray', linewidth=0.25, zorder=1)
        pass
    plot.hlines(y=seq_len, xmin=0.5, xmax=num_gens + 0.5, color="#64C5EB",
                linestyles="dashed", linewidth=1)

    plot.plot(xs, means, label="Mean Population Match Fitness", color='#7F58AF', linewidth=1.25)
    plot.plot(xs, bests, label="Best Population Match Fitness", color='#FEB326', linewidth=1.25)
    plot.set_xlabel("Mating Event", fontsize=12)
    plot.set_ylabel("Match Fitness", fontsize=12)
    ax2 = plot.twinx()
    ax2.plot(xs, novs, label="Best Population Diversity Fitness", color='#E84D8A', linewidth=1.25, alpha=0.7)
    ax2.set_ylabel("Diversity Fitness", fontsize=12)
    ax2.grid(None)

    plot.ticklabel_format(style='plain')
    plot.margins(x=0.01)
    # ax2.margins(x=0.1)
    ax2.ticklabel_format(style='plain')

    f.suptitle("Convergence Plot", fontsize=14)
    plot.set_title("Set " + str(seq) + " - Experiment " + str(exp_desc), fontsize=12)

    lines, labels = plot.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=4, facecolor='white', frameon='true',
               fontsize=10, framealpha=0.75, ncol=1, borderaxespad=0.1)
    f.tight_layout()
    f.savefig(out_path + ".png", dpi=300)
    plt.close()
    pass


def calc_best_naive(seq_set, max_len, min_len):
    val = 0
    seq_set = seq_set.copy()
    for idx in range(min_len):
        if seq_set[0][idx] == seq_set[1][idx] and seq_set[1][idx] == seq_set[2][idx]:
            val += 3
            pass
        elif seq_set[0][idx] == seq_set[1][idx]:
            val += 2
            pass
        elif seq_set[0][idx] == seq_set[2][idx]:
            val += 2
            pass
        elif seq_set[1][idx] == seq_set[2][idx]:
            val += 2
            pass
        else:
            val += 1
            pass
        pass
    lens = [len(seq) for seq in seq_set]
    if lens[0] < lens[1] and lens[0] < lens[2]:
        seq_set = seq_set[1:]
        pass
    elif lens[1] < lens[0] and lens[1] < lens[2]:
        seq_set = [seq_set[0], seq_set[2]]
        pass
    elif lens[2] < lens[0] and lens[2] < lens[1]:
        seq_set = seq_set[:2]
        pass

    smaller = min(len(seq_set[0]), len(seq_set[1]))
    for idx in range(min_len, smaller):
        if seq_set[0][idx] == seq_set[1][idx]:
            val += 2
            pass
        else:
            val += 1
            pass
        pass
    val += max_len - smaller
    return val


def info_on_sets(seq_sets):
    tot_lens = []
    avg_lens = []
    max_lens = []
    best_naive = []
    for the_set in seq_sets:
        tot = 0
        max_len = 0
        min_len = 1000
        for seq in the_set:
            tot += len(seq)
            if len(seq) > max_len:
                max_len = len(seq)
                pass
            if len(seq) < min_len:
                min_len = len(seq)
                pass
            pass
        best = calc_best_naive(the_set, max_len, min_len)
        avg_len = tot / len(the_set)
        tot_lens.append(tot)
        avg_lens.append(avg_len)
        max_lens.append(max_len)
        best_naive.append(best)
        pass
    rtn = [[tot_lens[i], avg_lens[i], max_lens[i], best_naive[i]] for i in range(len(seq_sets))]
    return rtn


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
            output[1] = "-" + output[1]
            cur_row = cur_row - 1
            pass
        else:  # Gap in Seq1
            output[0] = "-" + output[0]
            output[1] = s2[cur_col - 1] + output[1]
            cur_col = cur_col - 1
            pass
        pass

    while cur_row > 0:  # Extra characters in Seq1
        output[0] = s1[cur_row - 1] + output[0]
        output[1] = "-" + output[1]
        cur_row -= 1
        pass

    while cur_col > 0:  # Extra characters in Seq2
        output[0] = "-" + output[0]
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
                longest[i][j] = longest[i - 1][j - 1] + 1
                src_i[i][j] = i - 1
                src_j[i][j] = j - 1
                pass
            elif longest[i - 1][j] >= longest[i][j - 1]:
                longest[i][j] = longest[i - 1][j]
                src_i[i][j] = i - 1
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
    return longest[len(s1)][len(s2)], output
    pass


def calc_matches(s1, s2):
    rtn = ''
    val = 0
    for i in range(min(len(s1), len(s2))):
        if s1[i] == s2[i]:
            rtn += "X"
            val += 1
            pass
        else:
            rtn += " "
            pass
        pass
    return rtn, val


def main():
    print("START")
    folder_names = os.listdir(inp)
    seq_idxs = [1, 2, 3, 4, 5]
    exp_descriptions = ["Naive Fitness", "LCS Fitness with MaxSeqLen", "LCS Fitness with 1.25 AvgSeqLen",
                        "LCS Fitness with 1.50 AvgSeqLen"]

    groups = []
    for seq in seq_idxs:
        groups.append(["Set" + str(seq)])
        pass

    sequences = gen_sequences("./Sequences.dat")

    for idx in range(len(sequences)):
        for seq_idx in range(len(sequences[idx])):
            sequences[idx][seq_idx] = DNA_to_int(sequences[idx][seq_idx])
            print(len(sequences[idx][seq_idx]))
            print(int_to_DNA(sequences[idx][seq_idx]))
            pass
    seq_set_info = info_on_sets(sequences)
    for seq_set in sequences:
        freq_data = get_char_freqs(seq_set, num_chars=4)
        for dat in freq_data:
            print(freq_data[0])
            print(dat)
            pass
        pass
    for idx in range(len(sequences)):
        for seq_idx in range(len(sequences[idx])):
            sequences[idx][seq_idx] = int_to_DNA(sequences[idx][seq_idx])
            pass

    # for idx in range(len(sequences)):
    #     sequences[idx] = DNA_to_int(sequences[idx])
    #     print(len(sequences[idx]))
    #     print(int_to_DNA(sequences[idx]))
    #     pass
    #
    # freq_data = get_char_freqs(sequences, num_chars=4)
    # for dat in freq_data:
    #     print(freq_data[0])
    #     print(dat)
    #     pass

    exp_dirs = []
    all_dirs = []
    for eidx, dat in enumerate(groups):
        one_exp = []
        for fld in folder_names:
            if all(fld.__contains__(itm) for itm in dat):
                one_exp.append(fld)
                all_dirs.append(fld)
            pass
        exp_dirs.append(one_exp)
        pass

    # for group in exp_dirs:
    #     for exp_idx in range(0, len(group) -1, 2):
    #         tmp = group[exp_idx]
    #         group[exp_idx] = group[exp_idx + 1]
    #         group[exp_idx + 1] = tmp
    #         pass
    #     pass

    exp_lbls = []
    exp_num = 1
    for exp_idx, exp_dat in enumerate(exp_dirs):
        for fld in exp_dat:
            fields = fld.rstrip().split(",")
            if fields[-2].__contains__("0Fit"):
                lbl = ""
                pass
            else:
                lbl = ""
                if fields[-1].__contains__("1.00"):
                    lbl += "MaxSeqLen"
                    pass
                else:
                    lbl += fields[-1].split("LCS")[0] + "AvgSeqLen"
                    pass
                pass
            if lbl not in exp_lbls:
                exp_lbls.append(lbl)
                pass
            exp_num += 1
            pass
        pass

    # mode_data[group][exp][run][0] = run num
    # mode_data[group][exp][run][1] = run's match fit (primary sort)
    # mode_data[group][exp][run][2] = run's novelty fit (secondary sort)
    # mode_data[group][exp][run][3][:] = run's SDA
    # mode_data[group][exp][run][4] = run's sequence
    mode_data = [[] for _ in range(len(exp_dirs))]
    for eidx, exp in enumerate(groups):
        for fld in exp_dirs[eidx]:
            mode_data[eidx].append(get_data(inp + fld + "/"))
            pass
        pass

    # mode_stats[seq][exp] = [run's fitness vals]
    mode_stats = [[] for _ in range(len(groups))]
    make_all = False
    make_any = False
    for gidx, seq in enumerate(groups):
        for expidx, exp in enumerate(mode_data[gidx]):
            exp_fits = []
            for runidx, run in enumerate(exp):
                exp_fits.append(run[1])
                pass
            mode_stats[gidx].append(exp_fits)
            pass
        pass

    title = "Sequence Matching using Set "
    # xsp = [[i for i in range(len(all_data[0]))], [i for i in range(len(all_data[1]))]]
    # xpos = [xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1]]
    ylb = "Fitness"
    xlb = "Experiment"

    lxpos = [1.5]
    # for i in range(2, len(mode_stats[0]), 2):
    #     lxpos.append(i + 0.5)
    #     pass
    colors = ['#7F58AF', '#64C5EB', '#E84D8A', '#FEB326']
    for gidx, ginfo in enumerate(groups):
        if len(mode_stats[gidx]) > 0:
            seq_id = int(ginfo[0][3:])
            f = open(outp + "exp_table" + str(gidx + 1) + ".dat", "w")
            for idxx, gr in enumerate(exp_dirs[gidx]):
                f.write(str(idxx + 1) + "\t")
                f.write(gr)
                f.write("\n")
                pass
            f.close()

            plt.style.use("seaborn-v0_8")
            plt.rc('xtick', labelsize=12)
            plt.rc('ytick', labelsize=12)

            f = plt.figure()
            f.set_figheight(4.5)
            f.set_figwidth(8)
            plot = f.add_subplot(111)

            bp = plot.boxplot(mode_stats[gidx], patch_artist=True)
            # box_plot(bp, 1, [[[i for i in range(len(mode_stats[gidx]))], colors]])
            box_plot(bp, 1, [[[i for i in range(len(mode_stats[gidx]))], colors]])

            plot.tick_params(axis='x', which='major', pad=20)
            plot.set_xticks([x + 1 for x in range(len(mode_stats[gidx]))], exp_lbls, minor=True, fontsize=10)
            labels = ["Naive Fitness", "LCS Fitness"]
            plot.set_xticks([1.0001, 3.0001], labels, rotation=0, minor=False, fontsize=12)

            new_title = title + str(seq_id)
            plot.set_title(new_title, fontsize=14)
            # plot.set_xlabel(xlb, fontsize=13)
            plot.set_ylabel(ylb, fontsize=12)

            plot.hlines(y=seq_set_info[gidx][3], xmin=0.5, xmax=1.5, linewidth=1.5, color='r')
            plot.hlines(y=seq_set_info[gidx][0], xmin=1.5, xmax=4.5, linewidth=1.5, color='b')
            # plot.hlines(y=len(sequences[seq_id]), xmin=0.5, xmax=len(mode_stats[gidx]) + 0.5, color="#2e294e",
            #             linestyles="dashed", linewidth=1)

            # labels = ["Naive Fitness", "LCS Fitness"]
            # patches = [Patch(facecolor=colors[0], edgecolor="#2e294e",linewidth=1, label=labels[0]),
            #            Patch(facecolor=colors[1], edgecolor="#2e294e", linewidth=1, label=labels[1])]
            # plot.legend(handles=patches, labels=labels, facecolor='white', frameon='true', fontsize=12, framealpha=1,
            #     loc='lower center', ncol=2, borderaxespad=0.1)

            for x in lxpos:
                plot.axvline(x=x, color='#2e294e', linestyle='-.', linewidth=1)
                pass
            plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
            plot.grid(visible="True", axis="x", which='minor', color="white", linewidth=0.75)
            f.tight_layout()
            f.savefig(outp + "Boxplot MatchSeq" + str(seq_id), dpi=300)
            plt.close()
            pass
        pass

    all_sda_batches = []
    for seq_num, fld_names in enumerate(exp_dirs):
        one_seq_batches = []
        for exp_num, fold in enumerate(fld_names):
            best_SDA_file = inp + fold + "/" + "SDAs" + str(mode_data[seq_num][exp_num][0][0]).zfill(2) + ".dat"
            one_exp_batches = []
            num_gens = 0
            with open(best_SDA_file, "r") as f:
                lines = f.readlines()
                one_batch = []
                for line in lines:
                    if line.__contains__("Population After"):
                        if len(one_batch) > 0:
                            one_exp_batches.append([one_batch, num_gens])
                            one_batch = []
                            pass
                        num_gens = int(line.split("Population After ")[1].split(" Mating Events")[0])
                        pass
                    else:
                        one_batch.append(line)
                        pass
                    pass
                one_exp_batches.append([one_batch, num_gens])
                pass
            one_seq_batches.append(one_exp_batches)
            pass
        all_sda_batches.append(one_seq_batches)
        pass

    # for seq_idx, seq_batches in enumerate(all_sda_batches):
    #     for exp_idx, exp_batches in enumerate(seq_batches):
    #         for batch in exp_batches:
    #             sda_infos, sda_outputs, sda_fits = process_sda(batch[0], sda_states, sda_chars)
    #             num_gens = batch[1]
    #             num_gens = int(int(num_gens) / 10) * 10
    #             make_heatmap(sda_infos, sda_states, sda_chars, outp + "HeatmapSet" + str(seq_idx + 1) + "Exp" + str(
    #                 exp_idx + 1), num_gens, mode_data[seq_idx][exp_idx][0][0], seq_idx + 1, exp_descriptions[exp_idx])
    #             pass
    #         pass
    #     pass
    #
    for seq_num, fld_names in enumerate(exp_dirs):
        for exp_num, fold in enumerate(fld_names):
            if exp_num == 0:
                best = seq_set_info[seq_num][3]
                pass
            else:
                best = seq_set_info[seq_num][0]
                pass
            make_convergence_plot(inp + fold + "/", mode_data[seq_num][exp_num][0][0], outp,
                                  int(all_sda_batches[seq_num][exp_num][-1][1]), cull_every,
                                  seq_num + 1, exp_num + 1, exp_descriptions[exp_num], best)
            pass
        pass

    for set_num, one_set in enumerate(sequences):
        for exp_num, exp_dat in enumerate(mode_data[set_num]):
            best_string = exp_dat[0][4]
            with open(outp + "BestMatch" + "Set" + str(set_num + 1) + "Exp" + str(exp_num + 1)
                      + "Run" + str(exp_dat[0][0]) + ".dat", "w") as f:
                if exp_num == 0:
                    tot = 0
                    for seq in one_set:
                        f.write(best_string + "\n")
                        thing, val = calc_matches(best_string, seq)
                        f.write(thing + "\n")
                        f.write(seq + "\n")
                        f.write("Length: " + str(len(seq)) + "\n")
                        f.write("Score: " + str(val) + "\n")
                        tot += val
                        pass
                    f.write("Total: " + str(tot) + "\n")
                    f.write("Best Possible: " + str(seq_set_info[set_num][3]))
                    pass
                else:
                    f.write(best_string + "\n\n")
                    tot = 0
                    for seq in one_set:
                        val, LCSOutput = LCS(best_string, seq)
                        f.write(LCSOutput[0] + "\n")
                        # f.write(calc_matches(LCSOutput[0], LCSOutput[1]) + "\n")
                        f.write(LCSOutput[1] + "\n")
                        f.write("Length: " + str(len(seq)) + "\n")
                        f.write("Score: " + str(val) + "\n")
                        f.write("\n")
                        tot += val
                        pass
                    f.write("Total: " + str(tot) + "\n")
                    f.write("Best Possible: " + str(seq_set_info[set_num][0]))
                    pass
                pass
            pass
        pass

    desc = ["First", "Second", "Third"]
    for set_num, one_set in enumerate(sequences):
        for exp_num, exp_dat in enumerate(mode_data[set_num]):
            best_string = exp_dat[0][4]
            with open(outp + "Seqs" + "Set" + str(set_num + 1) + "Exp" + str(exp_num + 1)
                      + "Run" + str(exp_dat[0][0]) + ".dat", "w") as f:
                f.write(">SDA" + "Seqq \n")
                f.write(best_string + "\n")
                for idx, seq in enumerate(one_set):
                    f.write(">" + desc[idx] + "Seqq \n")
                    f.write(seq + "\n")
                    pass
                pass
            pass
        pass

    # for sidx, seq in enumerate(seq_idxs):
    #     best_runs = []
    #     best_of_best_exp = -1
    #     best_of_best_val = 0
    #     for didx, dat in enumerate(mode_data[sidx]):
    #         best_runs.append(dat[0][0])
    #         if dat[0][1] > best_of_best_val:
    #             best_of_best_val = dat[0][1]
    #             best_of_best_exp = didx
    #         pass
    #     make_table(mode_stats[sidx], best_runs, exp_descriptions, outp +
    #                "Seq" + str(seq) + "table" + ".dat", False)
    #     info = "Best for Sequence " + str(seq) + " is " + \
    #            "EXP" + str(best_of_best_exp + 1) + ": " + str(exp_descriptions[best_of_best_exp])
    #     print_best_info(outp + "Seq" + str(seq) + "_best.dat", info,
    #                     mode_data[sidx][best_of_best_exp][0], sequences[seq])
    #     pass

    print("END")
    pass


main()
