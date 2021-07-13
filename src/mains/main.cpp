
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>
#include "testExecution.h"
#include <stddef.h>
#include <immintrin.h>
#include <omp.h>


#include <argp.h>
#include <stdbool.h>
#include <tsxcount/TSXHashMapTSXPerf.h>
#include <tsxcount/TSXHashMapCAS.h>
#include <tsxcount/TSXHashMapOMPPerf.h>
#include <tsxcount/TSXHashMapPerf.h>
#include <tsxcount/TSXHashMapPThreadPerf.h>
#include <tsxcount/TSXHashMapTSXPerf_EXP.h>

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
        { "checkabort", 'a', 0, OPTION_ARG_OPTIONAL, "abort if check count raises error"},
        { "threads", 't', "THREADS", OPTION_ARG_OPTIONAL, "Number of threads."},
        { "mode", 'm', "MODE", OPTION_ARG_OPTIONAL, "counting mode"},
        { 0 }
};
enum tsx_mode{ SERIAL, PTHREAD,OMP, CAS, TRANSACTIONAL, EXPERIMENTAL };
const char * TSXModeStrings[] = { "SERIAL", "PTHREAD", "OMP", "CAS", "TSX", "EXPERIMENTAL" };

struct arguments {
    uint16_t k,l,storagebits;
    uint8_t threads;
    std::string input_path;
    bool check;
    bool checkabort;
    tsx_mode mode;
};

tsx_mode strToMode(char* pArg)
{
    std::string argStr(pArg);
    std::transform(argStr.begin(), argStr.end(),argStr.begin(), ::toupper);

    if (argStr == "SERIAL") {
        return tsx_mode::SERIAL;
    } else if (argStr == "PTHREAD")
    {
        return tsx_mode::PTHREAD;
    } else if (argStr == "OMP")
    {
        return tsx_mode::OMP;
    } else if (argStr == "CAS")
    {
        return tsx_mode::CAS;
    } else if (argStr == "TSX")
    {
        return tsx_mode::TRANSACTIONAL;
    } else if (argStr == "EXPERIMENTAL")
    {
        return tsx_mode::EXPERIMENTAL;
    }

    return tsx_mode::TRANSACTIONAL;
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
        case 'a': arguments->checkabort = true; break;
        case 'm': arguments->mode = strToMode(arg);
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


void countKMers(TSXHashMap* pMap, struct arguments* pARGP, uint8_t threads)
{

//    uint8_t threads = pMap->getThreads();
    std::cout << "Threads should be set to " << (int) threads << std::endl;
    //std::cout << "Running on " << (int) omp_get_max_threads() << " threads (set to " <<  << ")" << std::endl;
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
    std::set<std::string> oTestSet;

    FASTXreader<FASTQEntry>* pReader = new FASTXreader<FASTQEntry>(&sFileName);

    uint32_t iK = pMap->getK();

    std::map<std::string, int> m;

//omp_set_num_threads(pMap->getThreads());
omp_set_num_threads(threads);

#pragma omp parallel
    {
#pragma omp single
	{
	    std::cout << "Running on " << (int) omp_get_max_threads() << " threads (set to " << omp_get_num_threads() << ")" << std::endl;

            bool exitNow = false;
            bool verbose = false;

            while (pReader->hasNext())
            {

                std::vector<FASTQEntry>* pEntries = pReader->getEntries(40);
                uint64_t iPosOfInterest = 0;

		//std::cout << "master " << omp_get_thread_num() << std::endl;

#pragma omp task firstprivate(pEntries) shared(iPosOfInterest)
                {

                   uint8_t taskThreads = omp_get_num_threads();

#pragma omp critical
                    {
                        std::cout << "In thread " << omp_get_thread_num() << " of " << (int) taskThreads << " reading entries " << pEntries->size() << std::endl;
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
                            //bool vadd = false;
                            //std::vector<uint64_t> vAllPos;

                            //#pragma omp critical
                            //{
                            //    m[kmerStr] = m[kmerStr] + 1;
                            //}

                            //bool strInSet = true;
                            //#pragma omp critical
                            //{
                            //std::set<std::string>::iterator oKit = oTestSet.find(kmerStr);
                            //strInSet = oKit != oTestSet.cend();
                            //oTestSet.insert(kmerStr);
                            //}

                            //uint32_t iKmerCount1 = pMap->getKmerCount(oKmer, false).toUInt();

                            pMap->addKmer(oKmer, verbose); //vadd

                            //uint32_t iKmerCount2 = pMap->getKmerCount(oKmer, false).toUInt();

                            //if (iKmerCount2 == 0)
                            //{
                            //    std::cout << "no increment happened! " << kmerStr << std::endl;
                            //}
                            //iAddedKmers += 1;
                        }
                    }


                    delete pEntries;
                }



            }


        }
        #pragma omp taskwait
        #pragma omp barrier


    }



    std::cout << "Added a total of " << pMap->getKmerCount() << " different kmers" << std::endl;
    uint64_t iRefCount = 0;
    if (pARGP->check)
    {
        std::string sRefFilename = sFileName + "." + std::to_string(iK) + ".count";
        uint64_t totalerrors = 0;
        
        MemoryPool<FIELDTYPE>* pDirectPool = new DirectMemoryPool<FIELDTYPE>();

        UBigInt queriedPositions = UBigInt(0, pDirectPool);
        queriedPositions.resize(pMap->getMaxElements());

        #pragma omp parallel
        {
        #pragma omp single
            {

                std::cout << "Checking kmer counts against manual hashmap ..." << std::endl;

                std::cerr << "Loading reference file: " << sRefFilename << std::endl;
                std::ifstream file(sRefFilename);


                if (file.is_open()) {
                    std::string line;



                    std::map<std::string, int>* pKmer2count = new std::map<std::string, int>();


                    while (getline(file, line)) {

                        std::vector<std::string> vElems = split(line);

                        uint32_t iCount = std::atoi(vElems[1].c_str());
                        std::string kmer = vElems[0];

                        pKmer2count->insert(std::pair<std::string, int>(kmer, iCount));


                        if (pKmer2count->size() >= 100000)
                        {
                            std::cout << "Going to check " << pKmer2count->size() << " kmers" << std::endl;
                            iRefCount += pKmer2count->size();

                            #pragma omp task firstprivate(pKmer2count)
                            {

                                std::map<std::string, int>::iterator oIt;
                                for (oIt = pKmer2count->begin(); oIt != pKmer2count->end(); ++oIt)
                                {

                                    std::string kmerStr = oIt->first;
                                    int kmerCount = oIt->second;

                                    //std::cout << kmerStr << " " << kmerCount << std::endl;

                                    UBigInt tkmer = TSXSeqUtils::fromSequenceD(kmerStr, pMap->getMemoryPool());
                                    KmerCountDebug oRes = pMap->getKmerCountDebug(tkmer);

                                    uint8_t errorCount = evaluate(pMap, tkmer, kmerCount, &(kmerStr), true, NULL);

                                    if (errorCount)
                                    {
                                        if (pARGP->checkabort)
                                        {
                                            exit(200);
                                        }
                                    }

                                    if (errorCount)
                                    {
                                        #pragma omp critical
                                        {
                                            totalerrors += errorCount;
                                            queriedPositions.setBit(oRes.iFirstPos, 1);

                                        }
                                        
                                    } else {
                                        #pragma omp critical
                                        {
                                            queriedPositions.setBit(oRes.iFirstPos, 1);
                                        }
                                    }
                                }
                                std::cout << "Checked " << pKmer2count->size() << " kmers" << std::endl;


                                delete pKmer2count;
                            }


                            pKmer2count = new std::map<std::string, int>();
                        }


                    }
                    file.close();

                    // any remaining kmers are handled here!

                    iRefCount += pKmer2count->size();

                    #pragma omp task firstprivate(pKmer2count)
                    {

                        std::map<std::string, int>::iterator oIt;
                        for (oIt = pKmer2count->begin(); oIt != pKmer2count->end(); ++oIt)
                        {

                            std::string kmerStr = oIt->first;
                            int kmerCount = oIt->second;

                            //std::cout << kmerStr << " " << kmerCount << std::endl;

                            UBigInt tkmer = TSXSeqUtils::fromSequenceD(kmerStr, pMap->getMemoryPool());
                            KmerCountDebug oRes = pMap->getKmerCountDebug(tkmer);
                            uint8_t errorCount = evaluate(pMap, tkmer, kmerCount, &(kmerStr), true, NULL);

                            if (errorCount)
                            {
                                #pragma omp critical
                                {
                                    totalerrors += errorCount;
                                    queriedPositions.setBit(oRes.iFirstPos, 1);

                                }
                                
                            } else {
                                #pragma omp critical
                                {
                                    queriedPositions.setBit(oRes.iFirstPos, 1);
                                }
                            }
                        }
                        std::cout << "Checked " << pKmer2count->size() << " kmers" << std::endl;


                        delete pKmer2count;
                    }

                }
            }
        }

    	std::cout << "total errors" << (int) totalerrors << std::endl;
        std::cout << "Kmer count check completed." << std::endl;

        // checking counts
        std::cout << "Reference kmer count: " << (int) iRefCount << std::endl;
        std::cout << "queried kmer count: " << (int) queriedPositions.sumBits() << std::endl;
        std::cout << "tsxCount kmer count: " << (int) pMap->getKmerCount() << std::endl;

        std::cout << "calculating quered positions" << std::endl;
        UBigInt& oKmerStarts = pMap->getKmerStartsRef();
        std::cout << "got kmer starts" << std::endl;
        queriedPositions.bitXor(oKmerStarts);
        std::cout << "queried (kmerstarts) kmer count: " << (int) oKmerStarts.sumBits() << std::endl;
        std::cout << "queried (Xor) kmer count: " << (int) queriedPositions.sumBits() << std::endl;

        if (queriedPositions.sumBits() > 0)
        {  
            std::vector<uint64_t> vOnePos = queriedPositions.onePositions();
            std::cout << "queried one pos count: " << (int) vOnePos.size() << std::endl;

            for (size_t j = 0; j < vOnePos.size(); ++j)
            {
                uint64_t iCurPos = vOnePos[j];
                //pMap->getCountAtPosition(iCurPos);
            }
        }

    }





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
            pMap = new TSXHashMapPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);

            if (arguments.threads != 1)
            {
                std::cerr << "Requesting to run SERIAL with Threads != 1 => EXIT(0)" << std::endl;
		std::cerr << "CONTINUEING AT OWN RISK!" << std::endl;
                //return 0;
            }
            break;



        case PTHREAD:
            std::cerr << "Creating TSXHashMap PTHREAD" << std::endl;
            pMap = new TSXHashMapPThreadPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
        case OMP:
            std::cerr << "Creating TSXHashMap OMP" << std::endl;
            pMap = new TSXHashMapOMPPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;

        case CAS:
            std::cerr << "Creating TSXHashMap CAS" << std::endl;
            pMap = new TSXHashMapCAS(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
        case TRANSACTIONAL:
            std::cerr << "Creating TSXHashMap TRANSACTIONS/TSX" << std::endl;
            pMap = new TSXHashMapTSXPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
        case EXPERIMENTAL:
            std::cerr << "Creating TSXHashMap TRANSACTIONS/TSX" << std::endl;
            pMap = new TSXHashMapTSXPerfExperimental(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
        default:
            std::cerr << "Creating TSXHashMap DEFAULT(TRANSACTIONS/TSX)" << std::endl;

            pMap = new TSXHashMapTSXPerf(arguments.l, arguments.storagebits, arguments.k, arguments.threads);
            break;
    }

    countKMers(pMap, &arguments, arguments.threads);

    pMap->print_stats();

    TSXHashMapTSXPerf* pPerfTSX = NULL;
    TSXHashMapCAS* pCAS = NULL;
    if (pPerfTSX = dynamic_cast<TSXHashMapTSXPerf*>(pMap))
    {

        std::cerr << "addkmer-inserted count" << pPerfTSX->iInsertedCount << std::endl;
        std::cerr << "Used fields: " << pPerfTSX->getUsedPositions() << std::endl;
        std::cerr << "adds: " << pPerfTSX->iAddCount << std::endl;
        std::cerr << "add calls: " << pPerfTSX->iAddKmerCount << std::endl;
        std::cerr << "aborts calls: " << pPerfTSX->iAborts << std::endl;

        pPerfTSX->print_counts();

    } else if (pCAS = dynamic_cast<TSXHashMapCAS*>(pMap))
    {

        std::cerr << "Used fields: " << pCAS->getUsedPositions() << std::endl;
        std::cerr << "adds: " << pCAS->iAddCount << std::endl;
        std::cerr << "add calls: " << pCAS->iAddKmerCount << std::endl;
        std::cerr << "aborts calls: " << pCAS->iAborts << std::endl;
    }


    return 0;


}
