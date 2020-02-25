
#include <iostream>


#include <tsxcount/TSXHashMap.h>
#include <tsxcount/TSXTypes.h>
#include <utils/CLParser.h>
#include <fastxutils/FastXReader.h>
#include <utils/SequenceUtils.h>
#include <tsxcount/TSXHashMapTSXSmall.h>
#include "testExecution.h"
#include <stddef.h>
#include <immintrin.h>


#include <argp.h>
#include <stdbool.h>

const char *argp_program_version = "tsxCount 1.0";
const char *argp_program_bug_address = "joppich@bio.ifi.lmu.de";
static char doc[] = "Count k-mers.";
static char args_doc[] = "[FILENAME]...";
static struct argp_option options[] = {
        { "k", 'k', 0, 0, "parameter k"},
        { "l", 'l', 0, 0, "parameter l"},
        { "input", 'i', 0, 0, "input string"},
        { "check", 'c', 0, OPTION_ARG_OPTIONAL, "check counts"},
        { "threads", 't', 0, OPTION_ARG_OPTIONAL, "Number of threads."},

        { 0 }
};

struct arguments {
    uint32_t k,l;
    uint8_t threads;
    char* input_path;
    bool check;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = (struct arguments*) state->input;
    switch (key) {
        case 'k': arguments->k = atoi(arg); break;
        case 'l': arguments->l = atoi(arg); break;
        case 't': arguments->threads= atoi(arg); break;
        case 'i': arguments->input_path = arg; break;
        case 'c': arguments->check = true; break;
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };


int main(int argc, char *argv[])
{

    struct arguments arguments;

    arguments.check = false;
    arguments.threads = omp_get_max_threads();
    arguments.k = 14;
    arguments.l = 26;

    argp_parse(&argp, argc, argv, 0, 0, &arguments);


    uint32_t itK = 14;
    uint32_t iThreads = 4;

    TSXHashMapTSX* ptMap = new TSXHashMapTSX(26, 4, itK, iThreads);

    //testHashMapOld(ptMap, true);
    testHashMap(ptMap, iThreads != 1);


    std::cerr << "Used fields: " << ptMap->getUsedPositions() << std::endl;
    std::cerr << "adds: " << ptMap->iAddCount << std::endl;
    std::cerr << "add calls: " << ptMap->iAddKmerCount << std::endl;
    std::cerr << "aborts calls: " << ptMap->iAborts << std::endl;

    return 1;


}