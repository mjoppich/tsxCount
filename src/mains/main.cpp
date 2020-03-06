
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>
#include "testExecution.h"
#include <stddef.h>
#include <immintrin.h>


#include <argp.h>
#include <stdbool.h>
#include <tsxcount/TSXHashMapTSXPerf.h>
#include <tsxcount/TSXHashMapOMP.h>
#include <tsxcount/TSXHashMapTSXSmall.h>
#include <tsxcount/TSXHashMapCAS.h>
#include <tsxcount/TSXHashMapOMPPerf.h>
#include <tsxcount/TSXHashMapPerf.h>

const char *argp_program_version = "tsxCount 1.0";
const char *argp_program_bug_address = "joppich@bio.ifi.lmu.de";
static char doc[] = "Count k-mers.";
static char args_doc[] = "[FILENAME]...";
static struct argp_option options[] = {
        { "k", 'k', "K", 0, "parameter k"},
        { "s", 's', "STORAGE", 0, "parameter s"},
        { "l", 'l', "L", 0, "parameter l"},
        { "input", 'i', "INPUT_FASTA", 0, "input string"},
        { "check", 'c', 0, OPTION_ARG_OPTIONAL, "check counts"},
        { "threads", 't', "THREADS", OPTION_ARG_OPTIONAL, "Number of threads."},
        { "mode", 'm', "MODE", OPTION_ARG_OPTIONAL, "counting mode"},
        { 0 }
};
enum tsx_mode{ SERIAL, PTHREAD, OMP, CAS, TRANSACTIONS, OMPPERF, SERIALPERF };
const char * TSXModeStrings[] = { "SERIAL", "PTHREAD", "OMP", "CAS", "TRANSACTIONS/TSX", "OMPPERF", "SERIALPERF" };

struct arguments {
    uint16_t k,l,storagebits;
    uint8_t threads;
    std::string input_path;
    bool check;
    tsx_mode mode;
};

tsx_mode strToMode(char* pArg)
{
    std::string argStr(pArg);
    std::transform(argStr.begin(), argStr.end(),argStr.begin(), ::toupper);

    if (argStr == "SERIAL")
    {
        return tsx_mode::SERIAL;
    } else if (argStr == "SERIALPERF")
    {
        return tsx_mode::SERIALPERF;
    } else if (argStr == "PTHREAD")
    {
        return tsx_mode::PTHREAD;
    } else if (argStr == "OMP")
    {
        return tsx_mode::OMP;
    } else if (argStr == "OMPPERF")
    {
        return tsx_mode::OMPPERF;

    } else if (argStr == "CAS")
    {
        return tsx_mode::CAS;
    } else if (argStr == "TRANSACTIONS")
    {
        return tsx_mode::TRANSACTIONS;
    }  else if (argStr == "TSX")
    {
        return tsx_mode::TRANSACTIONS;
    }

    return tsx_mode::TRANSACTIONS;
}

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = (struct arguments*) state->input;

    switch (key) {
        case 'k': arguments->k = atoi(arg); break;
        case 'l': arguments->l = atoi(arg); break;
        case 's': arguments->storagebits = atoi(arg); break;
        case 't': arguments->threads= atoi(arg); break;
        case 'i': arguments->input_path = (arg != NULL) ? std::string(arg) : ""; break;
        case 'c': arguments->check = true; break;
        case 'm': arguments->mode = strToMode(arg);
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


void countKMers(TSXHashMap* pMap, struct arguments* pARGP)
{

    const size_t iMaxCount = 2048*4*16;

    uint8_t threads = pMap->getThreads();
    omp_set_num_threads(pMap->getThreads());

    std::cout << "Running on " << (int) omp_get_max_threads() << " threads" << std::endl;
    omp_set_dynamic(0);

    /*
    std::string sFileName = "/mnt/d/owncloud/data/tsx/usmall_t7.fastq";
    sFileName = "/mnt/d/owncloud/data/tsx/small2_t7.fastq";
    sFileName = "/mnt/d/owncloud/data/tsx/small_t7.3000.fastq";
    sFileName = "/mnt/d/owncloud/data/tsx/small_t7.5000.fastq";
    sFileName = "/mnt/d/owncloud/data/tsx/small_t7.8000.fastq";
    sFileName = "/mnt/d/owncloud/data/tsx/small_t7.fastq";
    */
    std::string sFileName = pARGP->input_path;

    FASTXreader<FASTQEntry>* pReader = new FASTXreader<FASTQEntry>(&sFileName);

    uint32_t iK = pMap->getK();


#pragma omp parallel num_threads(threads)
    {
#pragma omp master
        {
            bool exitNow = false;
            bool verbose = false;

            while (pReader->hasNext())
            {

                std::vector<FASTQEntry>* pEntries = pReader->getEntries(40);
                uint64_t iPosOfInterest = 0;

#pragma omp task firstprivate(pEntries) shared(iPosOfInterest)
                {

#pragma omp critical
                    {
                        std::cout << "In thread " << omp_get_thread_num() << " reading entries " << pEntries->size() << std::endl;
                    }

                    MemoryPool<FIELDTYPE>* pPool = pMap->getMemoryPool();

                    for (size_t i = 0; i < pEntries->size(); ++i)
                    {
                        FASTQEntry* pEntry = &(pEntries->at(i));

                        std::string sSeq = pEntry->getSequence();
                        std::vector<std::string> allKmers = createKMers(sSeq, iK, pPool);

                        uint32_t iAddedKmers = 0;
                        uint32_t iTotalKmers = allKmers.size();

                        for (auto kmerStr : allKmers)
                        {
                            TSX::tsx_kmer_t oKmer = TSXSeqUtils::fromSequence(kmerStr, pPool);
                            bool vadd = false;
                            std::vector<uint64_t> vAllPos;

                            pMap->addKmer(oKmer, verbose); //vadd
                            //iAddedKmers += 1;
                        }
                    }


                    delete pEntries;
                }



            }


        }

    }



    std::cout << "Added a total of " << pMap->getKmerCount() << " different kmers" << std::endl;

    if (pARGP->check)
    {
        std::map<std::string, uint32_t> kmer2c = loadReferences(sFileName + "." + std::to_string(iK) + ".count");


        // checking counts
        std::cout << "Checking kmer counts against manual hashmap ..." << std::endl;
        std::cout << "Manual counts: " << kmer2c.size() << std::endl;


#pragma omp parallel
        {
#pragma omp single
            {
                for (auto& elem: kmer2c)
                {
#pragma omp task
                    {
                        // Do something with x, e.g.
                        UBigInt tkmer = TSXSeqUtils::fromSequenceD(elem.first, pMap->getMemoryPool());
                        evaluate(pMap, tkmer, elem.second, &(elem.first));
                    }
                }
            }
        }
    }




    std::cout << "Kmer count check completed." << std::endl;


}

int main(int argc, char *argv[])
{

    struct arguments arguments;

    arguments.check = false;
    arguments.threads = omp_get_max_threads();
    arguments.k = 14;
    arguments.l = 26;
    arguments.storagebits = 4;

    argp_parse(&argp, argc, argv, 0, NULL, &arguments);

    tsx_mode runMode = arguments.mode;
    TSXHashMap* pMap = NULL;

    std::cout << "Running with parameters " << std::endl;
    std::cerr << "K=" << (int)arguments.k << std::endl;
    std::cerr << "L=" << (int)arguments.l << std::endl;
    std::cerr << "StorageBits=" << (int)arguments.storagebits << std::endl;
    std::cerr << "Check=" << (arguments.check ? "Yes" : "No") << std::endl;
    std::cerr << "Input=" << arguments.input_path << std::endl;
    std::cerr << "Threads=" << (int)arguments.threads << std::endl;
    std::cerr << "Mode=" << TSXModeStrings[arguments.mode] << std::endl;

    switch (runMode)
    {
        case SERIAL:
            std::cerr << "Creating TSXHashMap SERIAL" << std::endl;
            pMap = new TSXHashMap(arguments.l, arguments.storagebits, arguments.k);

            if (arguments.threads != 1)
            {
                std::cerr << "Requesting to run SERIAL with Threads != 1 => EXIT(0)" << std::endl;
                return 0;
            }
            break;

        case SERIALPERF:
            std::cerr << "Creating TSXHashMap SERIALPERF" << std::endl;
            pMap = new TSXHashMapPerf(arguments.l, arguments.storagebits, arguments.k);

            if (arguments.threads != 1)
            {
                std::cerr << "Requesting to run SERIALPERF with Threads != 1 => EXIT(0)" << std::endl;
                return 0;
            }
            break;


        case PTHREAD:
            std::cerr << "Creating TSXHashMap PTHREAD" << std::endl;
            pMap = new TSXHashMapPThread(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;

        case OMP:
            std::cerr << "Creating TSXHashMap OMP" << std::endl;
            pMap = new TSXHashMapOMP(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
        case OMPPERF:
            std::cerr << "Creating TSXHashMap OMPPERF" << std::endl;
            pMap = new TSXHashMapOMPPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
        case CAS:
            std::cerr << "Creating TSXHashMap CAS" << std::endl;
            pMap = new TSXHashMapCAS(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;

        case TRANSACTIONS:
            std::cerr << "Creating TSXHashMap TRANSACTIONS/TSX" << std::endl;
            pMap = new TSXHashMapTSXPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;

        default:
            std::cerr << "Creating TSXHashMap DEFAULT(TRANSACTIONS/TSX)" << std::endl;

            pMap = new TSXHashMapTSXPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
    }

    countKMers(pMap, &arguments);

    pMap->print_stats();

    TSXHashMapTSXPerf* pPerfTSX = NULL;
    TSXHashMapCAS* pCAS = NULL;
    if (pPerfTSX = dynamic_cast<TSXHashMapTSXPerf*>(pMap))
    {

        std::cerr << "Used fields: " << pPerfTSX->getUsedPositions() << std::endl;
        std::cerr << "adds: " << pPerfTSX->iAddCount << std::endl;
        std::cerr << "add calls: " << pPerfTSX->iAddKmerCount << std::endl;
        std::cerr << "aborts calls: " << pPerfTSX->iAborts << std::endl;
    } else if (pCAS = dynamic_cast<TSXHashMapCAS*>(pMap))
    {

        std::cerr << "Used fields: " << pCAS->getUsedPositions() << std::endl;
        std::cerr << "adds: " << pCAS->iAddCount << std::endl;
        std::cerr << "add calls: " << pCAS->iAddKmerCount << std::endl;
        std::cerr << "aborts calls: " << pCAS->iAborts << std::endl;
    }


    return 0;


}