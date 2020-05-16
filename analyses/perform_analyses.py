import os,sys
import subprocess
import argparse



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
    parser.add_argument('-b', '--sample-base', type=str, default="small_t7")
    parser.add_argument('-s', '--samples', nargs='+', type=str, default=["1000", "3000", "5000"])
    parser.add_argument('-m', '--methods', nargs='+', type=str, default=["TSX", "CAS", "OMP"])
    parser.add_argument('-t', '--threads', nargs='+', type=int, default=[1,2,4,8])
    parser.add_argument('-r', '--repeats', type=int, default=3)
    parser.add_argument('-mt', '--max-time', type=int, default=30)
    parser.add_argument('-e', '--executable', type=argparse.FileType('r'), required=True)

    parser.add_argument('-l', '--L', type=int, default=26)

    args = parser.parse_args()

    repeats = args.repeats

    threads = args.threads
    samples = args.samples
    methods = args.methods

    #threads = [1,2,4]
    #methods = ["OMPPERF", "OMP", "TSX", "CAS", "SERIAL", "SERIALPERF"]

    for sample in samples:

        for method in methods:

            for threadCount in threads:

                for rep in range(0, repeats):
                    outfile = os.path.join(args.output, ".".join([str(x) for x in [args.sample_base, sample, method, threadCount, rep, "time"]]))
                    outerr = os.path.join(args.output, ".".join([str(x) for x in [args.sample_base, sample, method, threadCount, "err"]]))



                    
                    #print(args.executable.name)

                    samplefile = ".".join([args.sample_base, sample, 'fastq'])

                    cmd=["/usr/bin/time", "--verbose", args.executable.name, '--input='+samplefile, '--mode='+method, '--threads='+str(threadCount), "--l=" + str(args.L)]
                    cmdStr = " ".join(cmd)


                    if os.path.exists(outfile):
                        print(cmdStr)
                        print(sample, method, threadCount, rep)
                        print("Outfile already exists: ", outfile)
                        continue

                    runSucceeded = False

                    while not runSucceeded:
                        print(sample, method, threadCount, rep)
                        print(cmdStr)
                        try:
                            output = subprocess.check_output(cmdStr, shell=True, stderr=subprocess.STDOUT) #, timeout=args.max_time * 60

                            with open(outfile, 'bw') as fout:
                                fout.write(output)

                            runSucceeded = True
                            break
                        
                        except subprocess.CalledProcessError as err:
                            print("CalledProcessError", err.returncode)

                            with open(outerr, 'bw') as fout:
                                fout.write(err.output)
                            pass
                    





