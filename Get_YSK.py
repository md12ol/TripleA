
inp = "./Sequence Files/S1_ErrorsRemoved.fasta"
outp = "./S1_Temp.txt"


def get_seqs(filename, type):
    thing = []
    with open(filename, "r") as f:
        with open(outp, "w") as g:
            lines = f.readlines()
            get_next = False
            for line in lines:
                if get_next:
                    if line in thing:
                        print("ERROR\n" + line)
                    g.write(line)
                    thing.append(line)
                    get_next = False
                if line.__contains__(type):
                    if line in thing:
                        print("ERROR\n" + line)
                    g.write(line)
                    thing.append(line)
                    get_next = True
                    pass
                pass
            pass
        pass
    pass


def main():
    get_seqs(inp, "YSK")
    pass


main()