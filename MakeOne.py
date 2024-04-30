import os

tld = "./Output/"


def get_one_best(filename, low, high):
    with open(filename, "r") as f:
        one_best = f.readlines()
        match_val = 0
        nov_val = 0
        for line in one_best:
            if line.__contains__("Best Fitness: "):
                match_val = float(line.split("Best Fitness: ")[1].split(" of ")[0])
                nov_val = float(line.rstrip().split("novelty score of ")[1])
                pass
            if line.__contains__("Best Run: "):
                test_val = float(line.rstrip().split("Best Run: ")[1]) - 1 + low
                one_best[1] = "Best Run: " + str(int(test_val)) + "\n"
                if test_val < low or test_val > high:
                    print("ERROR get_one_best(): " + filename)
                    pass
                pass
            pass
        pass
    return [match_val, nov_val, one_best]


def get_the_best(fold, prefix, ends):
    data = []
    for end in ends:
        val = int(end)
        low = val
        high = val + 1
        dat = get_one_best(fold + "/" + prefix + end + ".dat", low, high)
        data.append(dat)
        pass

    best_match = 0
    best_novelty = 0
    best = []
    for dat in data:
        if dat[0] > best_match:
            best_match = dat[0]
            best_novelty = dat[1]
            best = dat[2]
            pass
        elif dat[0] == best_match and dat[1] < best_novelty:
            best_match = dat[0]
            best_novelty = dat[1]
            best = dat[2]
            pass
        pass
    return best


def get_exp(filename, start, low, high):
    idx = start
    one_run = []
    all_runs = []
    with open(filename, "r") as f:
        lines = f.readlines()
        first = True
        for line in lines:
            if not first and line.__contains__("The best fitness is "):
                all_runs.append(one_run)
                one_run = []
                idx += 1
                pass
            one_run.append(line)
            first = False
            pass
        all_runs.append(one_run)
        pass
    if len(all_runs) != high - low + 1:
        print("ERROR get_exp(): " + filename)
        pass
    return all_runs


def get_all_exp(fold, prefix, ends):
    all_exp = []

    for end in ends:
        val = int(end)
        low = val
        high = val + 1
        one_exp = get_exp(fold + "/" + prefix + end + ".dat", int(end), low, high)
        all_exp.extend(one_exp)
        pass
    return all_exp


def main():
    best_idxs = [str(i * 2 + 1).zfill(2) for i in range(25)]

    folds = os.listdir(tld)
    for fold in folds:
        if fold.__contains__("TripleA"):
            with open(tld + fold + "/best.dat", "w") as f:
                f.writelines(get_the_best(tld + fold, "best", best_idxs))
            pass
            all_exp = get_all_exp(tld + fold, "exp", best_idxs)
            with open(tld + fold + "/exp.dat", "w") as f:
                for exp in all_exp:
                    f.writelines(exp)
                    pass
                pass
            pass
        pass
    print("DONE!")
    pass


main()
