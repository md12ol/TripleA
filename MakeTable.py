import random

prog_name = "./cmake-build-release---beluga/TripleA"


def main():
    popsize = ["50"]
    first = "4 20"
    # second = "50 10000000 2 2 0 10"
    second = "5 10000000 2 2 0 10"
    set = ["1", "2", "3", "4", "5"]
    tourn_size = ["7"]
    # third = "0 0.5 1.0 0.25 1 5 1"
    third = "0 0.5 1.0 0.25 1 5"
    start_run = [str(i * 2 + 1) for i in range(25)]
    fit_fun = ["0", "1"]
    LCS_size = ["1", "1.25", "1.5"]

    with open("./table.dat", "w") as f:
        for ps in popsize:
            for seq in set:
                for ts in tourn_size:
                    for sr in start_run:
                        for ff in fit_fun:
                            if ff == "0":
                                f.writelines(prog_name + " " + ps + " " + first + " "
                                             + str(random.randint(1000, 9999)) + " " + second + " " + seq + " "
                                             + ts + " " + third + " " + sr + " " + ff + " 1" + "\n")
                            else:
                                for ls in LCS_size:
                                    f.writelines(prog_name + " " + ps + " " + first + " "
                                                 + str(random.randint(1000, 9999)) + " " + second + " " + seq + " "
                                                 + ts + " " + third + " " + sr + " " + ff + " " + ls + "\n")
                                    pass
                                pass
                            pass
                        pass
                    pass
                pass
            pass
        pass
    print("DONE")
    pass


main()
