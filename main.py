import matplotlib.pyplot as plt

inp = "./Sequence Files/Phi_YSK.txt"
DNA_alphabet = ["A", "C", "G", "T"]
pep_alphabet = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


def load_DNA(filename, requirement):
    infos = []
    DNA_seqs = []
    pep_seqs = []
    with open(filename, "r") as f:
        take_next = False
        lines = f.readlines()
        for line in lines:
            if "SeqID" in line and requirement in line:
                infos.append(line.rstrip())
                info = line.rstrip()
                take_next = True
                pass
            if "Peptide: " in line and take_next:
                line = line.split("Peptide: ")[1].rstrip().upper()
                if check_sequence(line, info, pep_alphabet):
                    pep_seqs.append(line)
                    pass
            if "DNA: " in line and take_next:
                line = line.split("DNA: ")[1].rstrip().upper()
                if check_sequence(line, info, DNA_alphabet):
                    DNA_seqs.append(line)
                    pass
                pass
            pass
        pass
    return infos, DNA_seqs, pep_seqs


def check_sequence(sequence, info, alphabet):
    for char in sequence:
        if char not in alphabet:
            print(info)
            print("Skipped " + char + " " + sequence)
            return False
        pass
    return True


def get_all_frequencies(sequences, alphabet):
    all_frequencies = [[] for _ in range(len(alphabet))]
    for sequence in sequences:
        one_frequency = get_frequencies(sequence, alphabet)
        for idx, freq in enumerate(one_frequency):
            all_frequencies[idx].append(freq)
            pass
        pass
    return all_frequencies


def get_frequencies(sequence, alphabet):
    count = [0 for _ in range(len(alphabet))]
    for char in sequence:
        count[alphabet.index(char)] += 1
        pass
    frequency = [val/len(sequence) for val in count]
    return frequency


def main():
    requirements = "YSK"

    infos, DNA_seqs, pep_seqs = load_DNA(inp, requirements)
    DNA_frequencies = get_all_frequencies(DNA_seqs, DNA_alphabet)
    pep_freqs = get_all_frequencies(pep_seqs, pep_alphabet)

    plt.style.use("seaborn-v0_8")
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    f = plt.figure()
    f.set_figheight(4.5)
    f.set_figwidth(10)
    plot = f.add_subplot(111)

    plot.boxplot(DNA_frequencies, labels=DNA_alphabet)
    plot.set_title("Frequencies of DNA Bases in YSK Phi-Segments")
    plot.set_xlabel("DNA Base")
    plot.set_ylabel("Frequency")
    f.tight_layout()
    f.savefig("Frequencies.png")

    f = plt.figure()
    f.set_figheight(4.5)
    f.set_figwidth(10)
    plot = f.add_subplot(111)

    plot.plot(DNA_frequencies, label=len(DNA_alphabet))
    plot.set_title("Frequencies of DNA Bases in YSK Phi-Segments")
    plot.set_xlabel("DNA Base")
    plot.set_ylabel("Frequency")
    f.tight_layout()
    f.savefig("Frequencies3.png")

    f = plt.figure()
    f.set_figheight(4.5)
    f.set_figwidth(10)
    plot = f.add_subplot(111)

    plot.boxplot(pep_freqs, labels=pep_alphabet)
    plot.set_title("Frequencies of Peptides in YSK Phi-Segments")
    plot.set_xlabel("Peptide")
    plot.set_ylabel("Frequency")
    f.tight_layout()
    f.savefig("Frequencies2.png")

    f = plt.figure()
    f.set_figheight(4.5)
    f.set_figwidth(10)
    plot = f.add_subplot(111)

    plot.plot(pep_freqs, label=len(pep_alphabet))
    plot.set_title("Frequencies of Peptides in YSK Phi-Segments")
    plot.set_xlabel("Peptide")
    plot.set_ylabel("Frequency")
    f.tight_layout()
    f.savefig("Frequencies4.png")
    pass


main()
