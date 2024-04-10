import subprocess

mast_inp = "./mast_out/mast.txt"
dna_inp = "./Sequence Files/S1_YSK.fasta"
pep_inp = "./Sequence Files/S2_YSK.fasta"
motif_len = 15
outp = "./Sequence Files/Phi_YSK.txt"


def get_peptide_seq_locs(filename):
    seq_ids = []
    seq_locs = []
    seq_lengths = []
    with open(filename, "r") as f:
        lines = f.readlines()
        close = False
        seq_next = False
        for line in lines:
            if "SECTION III:" in line:
                close = True
                pass
            elif close and "|" in line:
                seq_ids.append(line.rstrip())
                seq_next = True
                pass
            elif seq_next and "DIAGRAM:" in line:
                num_K_segs = line.count("-[1]-")
                if num_K_segs < 2:
                    seq_ids = seq_ids[:-1]
                    pass
                else:
                    if "-[1]-[1]-" in line:
                        print("ERROR: Two K-Segments, one right after the other")
                        pass
                    if "[1]-" in line[:4]:
                        print("ERROR: Starting with a K-Segment")
                        pass
                    if "-[1]" in line[-4:]:
                        print("ERROR: Ending with a K-Segment")
                        pass
                    chars_in_K_segs = (num_K_segs - 1) * motif_len
                    start_idx = 0
                    for val in line.split("DIAGRAM: ")[1].split("-[1]-")[0:num_K_segs - 1]:
                        start_idx += int(val)
                        pass
                    seq_locs.append(start_idx + chars_in_K_segs)
                    seq_lengths.append(int(line.split("DIAGRAM: ")[1].split("-[1]-")[num_K_segs - 1]))
                    pass
                seq_next = False
                pass
            elif seq_next and "+" in line:
                seq_ids = seq_ids[:-1]
                seq_next = False
                pass
            pass
        pass
    return [[seq_ids[i], seq_locs[i], seq_lengths[i]] for i in range(len(seq_ids))]


def get_pep_seqs(seq_id, seq_start, seq_len):
    with open(pep_inp, "r") as f:
        lines = f.readlines()
        next_seq = False
        for line in lines:
            if seq_id in line:
                next_seq = True
                pass
            elif next_seq:
                return line[seq_start:seq_start + seq_len]
            pass
        pass
    pass


def get_DNA_seqs(peptide_seqs):
    with open(outp, "w") as f:
        for pep_seq in peptide_seqs:
            next_seq = False
            with open(dna_inp, "r") as f2:
                lines = f2.readlines()
                for line in lines:
                    if pep_seq[0] in line:
                        next_seq = True
                        pass
                    elif next_seq:
                        f.write("SeqID: " + pep_seq[0] + "\n")
                        f.write("Peptide: " + get_pep_seqs(pep_seq[0], pep_seq[1], pep_seq[2]) + "\n")
                        start = pep_seq[1] * 3
                        end = start + pep_seq[2] * 3
                        f.write("DNA: " + line[start:end] + "\n")
                        next_seq = False
                        pass
                    pass
                pass
            pass
        pass
    pass


def get_both_seqs(species_names, filename):
    with open("./Sequence Files/Phi_A_YSK.txt", "w") as f:
        for type in ["Peptide: ", "DNA: "]:
            with open(filename, "r") as f2:
                lines = f2.readlines()
                get_next = False
                for line in lines:
                    if "SeqID" in line:
                        get_next = False
                        if any(name in line for name in species_names):
                            f.write(">" + line.split("SeqID: ")[1])
                            get_next = True
                            pass
                        pass
                    elif get_next and type in line:
                        f.write(line.split(type)[1])
                        pass
                    pass
                pass
            pass
            f.write("\n\n")
        pass
    pass


def main():
    mast_command = "mast './Sequence Files/YSK_KSeg.txt' './Sequence Files/S2_YSK.fasta'"
    # subprocess.run(mast_command, shell=True, executable="/bin/bash")
    try:
        result = subprocess.check_output(mast_command, shell=True, executable="/bin/bash", stderr=subprocess.STDOUT)
        pass

    except subprocess.CalledProcessError as cpe:
        result = cpe.output
        pass

    finally:
        for line in result.splitlines():
            print(line.decode())
            pass
        pass

    seq_locs = get_peptide_seq_locs(mast_inp)
    get_DNA_seqs(seq_locs)

    species_names = ["Ahalleri", "Alyrata", "Athaliana"]
    get_both_seqs(species_names, outp)
    pass


main()
