import os,sys
import subprocess
import argparse
import glob
from collections import defaultdict

import matplotlib.pyplot as plt

def toLatexTable(elems, sample):

    print("\\begin{table}")
    print("\\caption{{Runtimes for dataset {}}}".format(sample))
    print("\\label{{tab:testset:{}}}".format(sample))
    print("\\begin{tabular}{" + "l"*len(elems[0])+"}" )

    inSerialBlock = False
    for elem in elems:
        if elem[0].upper().startswith("METHOD"):
            print("&".join(["\\textbf{Method/Threads}"] + ["{:.0f}".format(x) for x in elem[1:]]) + "\\\\")
        else:

            if elem[0].endswith("[s]") and elem[0].upper().startswith("SERIAL"):
                inSerialBlock = True
            elif elem[0].endswith("[s]"):
                inSerialBlock = False

            if not inSerialBlock:
                print("&".join([elem[0]] + ["{:.2f}".format(x) for x in elem[1:]]) + "\\\\")
            else:
                remElems = ["-" for x in elem[2:]]

                elem = [elem[0]] + ["{:.2f}".format(elem[1])] + remElems
                print("&".join(elem) + "\\\\")


    print("\\end{tabular}")
    print("\\end{table}")

if __name__ == '__main__':
    def readable_dir(prospective_dir):
        try:
            os.makedirs(prospective_dir)
        except:
            pass

        if not os.path.isdir(prospective_dir):
            raise Exception("readable_dir:\"{0}\" is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            return prospective_dir
        else:
            raise Exception("readable_dir:{0} is not a readable dir".format(prospective_dir))


    parser = argparse.ArgumentParser(description='Process some gtf.')
    parser.add_argument("-o", "--output", type=readable_dir, required=True)
    parser.add_argument('-s', '--sample', type=str, default='small_t7.1000')
    args = parser.parse_args()

    method2color = {
        'TSX': (213,94,0),
        'OMP': (0,114,178),
        'PTHREAD': (0, 158,115),
        'SERIAL': (230,159,0),
        'CAS': (86,180,233)
    }

    for x in method2color:
        method2color[x] = tuple([(y/255) for y in method2color[x]])

    allRuns = defaultdict(lambda: defaultdict(list))

    seenThreads = set()

    for fname in glob.glob(args.output + "/" + args.sample + "*.time"):


        aname = fname.split(".")
        retry = aname[-2]
        threads = int(aname[-3])
        method = aname[-4]
        dataset = aname[-5]
        
        with open(fname) as fin:

            hours = 0
            mins = 0
            secs = 0

            wasParsed=False

            for line in fin:
                line = line.strip()

                if "Elapsed (wall clock) time (h:mm:ss or m:ss)" in line:
                    
                    aline = line.split("):")

                    wcTime = aline[1]

                    #2:54.28
                    awcTime = wcTime.split(":")

                    

                    if len(awcTime) == 2:
                        mins = int(awcTime[0])
                        secs = float(awcTime[1])
                    elif len(awcTime) == 3:
                        hours = int(awcTime[0])
                        mins = int(awcTime[1])
                        secs = float(awcTime[2])

                    wasParsed=True

                elif "Maximum resident set size" in line:

                    aline =line.split("):")

                    memSize = int(aline[1]) / 4096

        seenThreads.add(threads)

        print(wasParsed, fname, retry, threads, method, dataset, memSize, hours, mins, secs)

        if wasParsed:
            allRuns[method][threads].append((memSize, (hours, mins, secs), retry))


    plt.figure()

    origThreads = seenThreads
    seenThreads = ["Method"] + sorted(seenThreads)
    
    print(seenThreads)

    latexTable = []
    latexTable.append(seenThreads)

    xoffsets=0.25

    xoffset = 2*xoffsets / len(allRuns)

    for midx, method in enumerate(allRuns):

        threads = sorted([x for x in allRuns[method]])

        times = []

        for thread in threads:
            data = allRuns[method][thread]

            threadtimes = []
            for run in data:

                runtime = ((run[1][0] * 60 + run[1][1]) * 60) + run[1][2]

                threadtimes.append(runtime)

            times.append(sum(threadtimes)/len(threadtimes))

        

        speedUps = ["Speed-up", 1]
        efficiency = ["Efficiency", 1]


        for i in range(1, len(times)):
            if times[i] == 0:
                speedUps.append(0)
            else:
                speedUps.append(times[0]/times[i])

            efficiency.append(speedUps[-1] / seenThreads[i+1])


        latexTable.append([method + " [s]"] + times)
        latexTable.append(speedUps)
        latexTable.append(efficiency)

        threads = [x-xoffsets+midx*xoffset for x in threads]


        plt.scatter(threads, times, label=method, color=method2color[method])

    toLatexTable(latexTable, args.sample)

    plt.title("Runtime for sample {} with method and number of threads".format(args.sample))
    plt.xlabel("threads")
    plt.ylabel("time[s]")
    plt.yscale('log')

    origThreads = sorted([int(x) for x in origThreads])
    plt.xticks(origThreads)

    plt.legend()
    plt.savefig(args.output + "/" + args.sample.replace(".", "_") + "_times.png")

                




    