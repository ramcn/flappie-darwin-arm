/*
MIT License

Copyright (c) 2018 Yatish Turakhia, Gill Bejerano and William Dally

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#include <iostream>
#include <string>
#include <cstring>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <atomic>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <map>
#include <thread>
#include <mutex>
#include <algorithm> 
#include "fasta.h"
#include "ntcoding.h"
#include "seed_pos_table.h"
#include "gact.h"
#include "ConfigFile.h"


#include <dirent.h>
#include <glob.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>

#include "decode.h"
#include "fast5_interface.h"
#include "layers.h"
#include "networks.h"
#include "flappie_common.h"
#include "flappie_licence.h"
#include "flappie_output.h"
#include "flappie_stdlib.h"
#include "flappie_structures.h"
#include "util.h"
#include "version.h"

#if !defined(FLAPPIE_VERSION)
#    define FLAPPIE_VERSION "unknown"
#endif

// Doesn't play nice with other headers, include last
#include <argp.h>
#define MAX_THREADS 4


static char doc[] = "Flappie basecaller -- basecall from raw signal";
static char args_doc[] = "fast5 [fast5 ...]";
static struct argp_option options[] = {
    {"format", 'f', "format", 0, "Format to output reads (FASTA or SAM)"},
    {"limit", 'l', "nreads", 0, "Maximum number of reads to call (0 is unlimited)"},
    {"model", 'm', "name", 0, "Model to use (\"help\" to list)"},
    {"output", 'o', "filename", 0, "Write to file rather than stdout"},
    {"prefix", 'p', "string", 0, "Prefix to append to name of each read"},
    {"temperature", 7, "factor", 0, "Temperature for weights"},
    {"trim", 't', "start:end", 0, "Number of samples to trim, as start:end"},
    {"trace", 'T', "filename", 0, "Dump trace to HDF5 file"},
    {"licence", 10, 0, 0, "Print licensing information"},
    {"license", 11, 0, OPTION_ALIAS, "Print licensing information"},
    {"segmentation", 3, "chunk:percentile", 0, "Chunk size and percentile for variance based segmentation"},
    {"hdf5-compression", 12, "level", 0,
     "Gzip compression level for HDF5 output (0:off, 1: quickest, 9: best)"},
    {"hdf5-chunk", 13, "size", 0, "Chunk size for HDF5 output"},

    {"uuid", 14, 0, 0, "Output UUID"},
    {"no-uuid", 15, 0, OPTION_ALIAS, "Output read file"},
    {0}
};


#define DEFAULT_MODEL FLAPPIE_MODEL_R941_NATIVE

struct arguments {
    int compression_level;
    int compression_chunk_size;
    char * trace;
    enum flappie_outformat_type outformat;
    int limit;
    enum model_type model;
    FILE * output;
    char * prefix;
    float temperature;
    int trim_start;
    int trim_end;
    int varseg_chunk;
    float varseg_thresh;
    char ** files;
    bool uuid;
};

static struct arguments args = {
     1,
     200,
     NULL,
    FLAPPIE_OUTFORMAT_FASTQ,
     0,
     DEFAULT_MODEL,
     NULL,
    "",
     1.0f,
     200,
    10,
     100,
    0.0f,
     NULL,
     true
};



void fprint_flappie_models(FILE * fh, enum model_type default_model){
    if(NULL == fh){
        return;
    }

    for(size_t mdl=0 ; mdl < flappie_nmodel ; mdl++){
        fprintf(fh, "%10s : %s  %s\n", flappie_model_string(mdl), flappie_model_description(mdl),
                                      (default_model == mdl) ? "(default)" : "");
    }
}


static error_t parse_arg(int key, char * arg, struct  argp_state * state){
    int ret = 0;
    char * next_tok = NULL;

    switch(key){
    case 'f':
        args.outformat = get_outformat(arg);
        if(FLAPPIE_OUTFORMAT_INVALID == args.outformat){
            errx(EXIT_FAILURE, "Unrecognised output format \"%s\".", arg);
        }
        break;
    case 'l':
        args.limit = atoi(arg);
        assert(args.limit > 0);
        break;
    case 'm':
        if(0 == strcasecmp(arg, "help")){
            fprint_flappie_models(stdout, DEFAULT_MODEL);
            exit(EXIT_SUCCESS);
        }
        args.model = get_flappie_model_type(arg);
        if(FLAPPIE_MODEL_INVALID == args.model){
            fprintf(stdout, "Invalid Flappie model \"%s\".\n", arg);
            fprint_flappie_models(stdout, DEFAULT_MODEL);
            exit(EXIT_FAILURE);
        }
        break;
    case 'o':
        args.output = fopen(arg, "w");
        if(NULL == args.output){
            errx(EXIT_FAILURE, "Failed to open \"%s\" for output.", arg);
        }
        break;
    case 'p':
        args.prefix = arg;
        break;
    case 't':
        args.trim_start = atoi(strtok(arg, ":"));
        next_tok = strtok(NULL, ":");
        if(NULL != next_tok){
            args.trim_end = atoi(next_tok);
        } else {
            args.trim_end = args.trim_start;
        }
        assert(args.trim_start >= 0);
        assert(args.trim_end >= 0);
        break;
    case 'T':
        args.trace = arg;
        break;
    case 3:
        args.varseg_chunk = atoi(strtok(arg, ":"));
        next_tok = strtok(NULL, ":");
        if(NULL == next_tok){
            errx(EXIT_FAILURE, "--segmentation should be of form chunk:percentile");
        }
        args.varseg_thresh = atof(next_tok) / 100.0;
        assert(args.varseg_chunk >= 0);
        assert(args.varseg_thresh > 0.0 && args.varseg_thresh < 1.0);
        break;
    case 7:
	args.temperature = atof(arg);
	assert(isfinite(args.temperature) && args.temperature > 0.0f);
        break;
    case 10:
    case 11:
        ret = fputs(flappie_licence_text, stdout);
        exit((EOF != ret) ? EXIT_SUCCESS : EXIT_FAILURE);
        break;
    case 12:
        args.compression_level = atoi(arg);
        assert(args.compression_level >= 0 && args.compression_level <= 9);
        break;
    case 13:
        args.compression_chunk_size = atoi(arg);
        assert(args.compression_chunk_size > 0);
        break;
    case 14:
        args.uuid = true;
        break;
    case 15:
        args.uuid = false;
        break;
    case ARGP_KEY_NO_ARGS:
        argp_usage (state);
        break;

    case ARGP_KEY_ARG:
        args.files = &state->argv[state->next - 1];
        state->next = state->argc;
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


static struct argp argp = {options, parse_arg, args_doc, doc};


static struct _raw_basecall_info calculate_post(char * filename, enum model_type model){
    RETURN_NULL_IF(NULL == filename, (struct _raw_basecall_info){0});

    raw_table rt = read_raw(filename, true);
    RETURN_NULL_IF(NULL == rt.raw, (struct _raw_basecall_info){0});

    rt = trim_and_segment_raw(rt, args.trim_start, args.trim_end, args.varseg_chunk, args.varseg_thresh);
    RETURN_NULL_IF(NULL == rt.raw, (struct _raw_basecall_info){0});

    medmad_normalise_array(rt.raw + rt.start, rt.end - rt.start);

    flappie_matrix trans_weights = flipflop_transitions(rt, args.temperature, model);
    if (NULL == trans_weights) {
        free(rt.raw);
        free(rt.uuid);
        return (struct _raw_basecall_info){0};
    }

    const size_t nbase = nbase_from_flipflop_nparam(trans_weights->nr);
    const size_t nblock = trans_weights->nc;
    int * path = calloc(nblock + 2, sizeof(int));
    int * path_idx = calloc(nblock + 2, sizeof(int));
    float * qpath = calloc(nblock + 2, sizeof(float));
    int * pos = calloc(nblock + 1, sizeof(int));

    float score = NAN;

    flappie_matrix posterior = transpost_crf_flipflop(trans_weights, true);
    score = decode_crf_flipflop(posterior, false, path, qpath);
    size_t path_nidx = change_positions(path, nblock, path_idx);

    char * basecall = calloc(path_nidx + 1, sizeof(char));
    char * quality = calloc(path_nidx + 1, sizeof(char));
    for(size_t i=0 ; i < path_nidx ; i++){
        const size_t idx = path_idx[i];
        basecall[i] = base_lookup[path[idx] % nbase];
        quality[i] = phredf(expf(qpath[idx]));
    }

    exp_activation_inplace(posterior);
    flappie_imatrix trace = trace_from_posterior(posterior);
    posterior = free_flappie_matrix(posterior);
    free(qpath);
    free(path_idx);
    free(path);
    trans_weights = free_flappie_matrix(trans_weights);
    const size_t basecall_length = strlen(basecall);

    return (struct _raw_basecall_info) {
    	.score = score, 
        .rt = rt,
        .basecall = basecall,
        .quality = quality,
        .basecall_length = basecall_length,
        .trace = trace,
        .pos = pos,
        .nblock = nblock};
}



// GACT scoring
int gact_sub_mat[25];
int gap_open;
int gap_extend;

// D-SOFT parameters
std::string seed_shape;
std::string seed_shape_str;
uint32_t bin_size;
int dsoft_threshold;
int num_seeds;
int seed_occurence_multiple;
int max_candidates;
int num_nz_bins;

#ifdef FPGA
// Function prototypes
bool init();
void cleanup();
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name);
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name);
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name);
static void device_info_string( cl_device_id device, cl_device_info param, const char* name);
static void display_device_info( cl_device_id device );
void checkErr(cl_int err, const char * name);
cl_platform_id platform = NULL; 
scoped_array<cl_device_id> devices; // num_devices elements
cl_context context = NULL;
scoped_array<cl_command_queue> queues; // num_devices elements
cl_program program = NULL;
scoped_array<cl_kernel> kernels; // num_devices elements
scoped_array<cl_event> kernel_events;
unsigned int num_devices=NUM_DEVICES, max_num_devices=MAX_NUM_DEVICES;
#endif
int forward_kernel_counter=0;
int reverse_kernel_counter=0;


// GACT first tile
int first_tile_size;
int first_tile_score_threshold;

//GACT extend 
int tile_size;
int tile_overlap;

//Multi-threading
int num_references;


static std::string reference_string[MAX_THREADS];
static std::string query_string;

uint32_t reference_length[MAX_THREADS];
std::string reference_filenames[MAX_THREADS];
uint32_t query_length;

std::vector<long long int> reference_lengths[MAX_THREADS];
std::vector<std::string> reference_seqs[MAX_THREADS];

std::vector<long long int> reads_lengths;
uint32_t read_length;
std::vector<std::string> reads_seqs;
std::vector<std::string> rev_reads_seqs;

std::vector<std::vector<std::string> > reference_descrips[MAX_THREADS];
std::vector<long long int> reference_fileposs[MAX_THREADS];

std::vector<std::vector<std::string> > read_descrip;
std::vector<long long int> reads_fileposs;

char* reference_char[MAX_THREADS];
char** reads_char;
char** rev_reads_char;

std::map<int, uint32_t> chr_id_to_start_bin[MAX_THREADS];
std::map<uint32_t, int> bin_to_chr_id[MAX_THREADS];

SeedPosTable *sa[MAX_THREADS];

std::mutex io_lock;

std::atomic<int> num_aligned(0);

std::string RevComp(std::string seq) {
    std::string rc = "";
    for (int i = seq.size()-1; i >= 0; i--) {
        if (seq[i] != 'a' && seq[i] != 'A' &&
                seq[i] != 'c' && seq[i] != 'C' &&
                seq[i] != 'g' && seq[i] != 'G' &&
                seq[i] != 't' && seq[i] != 'T' &&
                seq[i] != 'n' && seq[i] != 'N') {
            std::cerr<<"Bad Nt char: "<< seq[i] <<std::endl;
            exit(1);
        }
        else {
            switch (seq[i]) {
                case 'a': rc += 't';
                          break;
                case 'A': rc += 'T';
                          break;
                case 'c': rc += 'g';
                          break;
                case 'C': rc += 'G';
                          break;
                case 'g': rc += 'c';
                          break;
                case 'G': rc += 'C';
                          break;
                case 't': rc += 'a';
                          break;
                case 'T': rc += 'A';
                          break;
                case 'n': rc += 'n';
                          break;
                case 'N': rc += 'N';
                          break;
            }
        }
    }
    return rc;
}


bool CompareAlignments(Alignment a1, Alignment a2) {
    return (a1.score > a2.score);
}

int mode;
int l1l2enable;

// called as align_threads.push_back(std::thread(AlignRead, thread_id, num_reads, num_references));
//void AlignRead (int start_read_num, int last_read_num, int num_references) {
void AlignRead (int k, int num_reads, int num_references) {

    uint32_t log_bin_size = (uint32_t) (log2(bin_size));
    int num_bins = 1 + (reference_length[k] >> log_bin_size);
    uint64_t* candidate_hit_offset;

    uint32_t* nz_bins_array = new uint32_t[num_nz_bins];
    uint64_t* bin_count_offset_array = new uint64_t[num_bins];
    candidate_hit_offset = new uint64_t[max_candidates];

    std::cout << "Aligned read called from thread " << k << " " << std::endl;
    for(int i=0; i < num_bins; i++) {
        bin_count_offset_array[i] = 0;
    }

    //for (int k = start_read_num; k < last_read_num; k+=num_threads) {

        int len = read_length;
        vector<Alignment> alignments;
        alignments.clear();
        //int thread_id = k % num_threads;

        // Forward reads
        
        int num_candidates;
           num_candidates = sa[k]->DSOFT(reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);


        // TODO: Need to apply more filtering here to avoid duplicate and
        // unnecessary alignments
        for (int i = 0; i < num_candidates; i++) {
            uint32_t candidate_hit = (candidate_hit_offset[i] >> 32);
            uint32_t last_hit_offset = ((candidate_hit_offset[i] << 32) >> 32);

            int chr_id = bin_to_chr_id[k][candidate_hit/bin_size];
            std::string chrom = reference_descrips[k][chr_id][0];
            uint32_t start_bin = chr_id_to_start_bin[k][chr_id];

            uint32_t ref_pos = candidate_hit - (start_bin*bin_size);
            uint32_t query_pos = last_hit_offset;
            uint32_t ref_len = reference_lengths[k][chr_id];
            uint32_t query_len = read_length;
            char* ref_start = reference_char[k] + (start_bin*bin_size);
            char strand = '+';

            Alignment align = GACT(ref_start, reads_char[k], chrom,  gact_sub_mat, gap_open, gap_extend, tile_size, tile_overlap, ref_pos, query_pos, ref_len, query_len, strand, first_tile_score_threshold, mode, k);
            alignments.push_back(align);
        }

        // Reverse complement reads
             num_candidates = sa[k]->DSOFT(rev_reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);

        // TODO: Need to apply more filtering here to avoid duplicate and
        // unnecessary alignments
        for (int i = 0; i < num_candidates; i++) {
            uint32_t candidate_hit = (candidate_hit_offset[i] >> 32);
            uint32_t last_hit_offset = ((candidate_hit_offset[i] << 32) >> 32);

            int chr_id = bin_to_chr_id[k][candidate_hit/bin_size];
            std::string chrom = reference_descrips[k][chr_id][0];
            uint32_t start_bin = chr_id_to_start_bin[k][chr_id];

            uint32_t ref_pos = candidate_hit - (start_bin*bin_size);
            uint32_t query_pos = last_hit_offset;
            uint32_t ref_len = reference_lengths[k][chr_id];
            uint32_t query_len = read_length;
            char* ref_start = reference_char[k] + (start_bin*bin_size);
            char strand = '-';
            

            Alignment align = GACT(ref_start, rev_reads_char[k], chrom, gact_sub_mat, gap_open, gap_extend, tile_size, tile_overlap, ref_pos, query_pos, ref_len, query_len, strand, first_tile_score_threshold, mode, k);
            alignments.push_back(align);
        }

        std::stable_sort(alignments.begin(), alignments.end(), CompareAlignments);
        int num_alignments = alignments.size();
        int* flags = (int*) calloc(num_alignments, sizeof(int));

        for (int m=0; m < num_alignments; m++) {
            if (flags[m] < 0) {
                continue;
            }
            if (alignments[m].aligned_query_len == 0) {
                flags[m] = 0;
                continue;
            }
            uint32_t s1 = alignments[m].query_start; 
            uint32_t e1 = s1+alignments[m].aligned_query_len; 
            for (int n=m+1; n < num_alignments; n++) {
                if (flags[n] < 0) {
                    continue;
                }
                uint32_t s2 = alignments[n].query_start; 
                uint32_t e2 = s2+alignments[n].aligned_query_len; 
                uint32_t s = std::max(s1 , s2);
                uint32_t e = std::min(e1, e2);
                uint32_t overlap = 0;
                if (s < e) {
                    overlap = e-s;
                }
                if (2*overlap >= alignments[n].aligned_query_len) {
                    flags[n] = -1;
                }
            }
        }

        io_lock.lock();
        for (int m=0; m < num_alignments; m++) {
            if (flags[m] >= 0) {
                Alignment align = alignments[m];
                std::string ext = ".maf";
    		std::ofstream out(reference_filenames[k].append(ext));
    		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(out.rdbuf());

                std::cout << "a score=" << align.score << std::endl;
                std::cout << "s\t" << align.ref_name << "\t" << 1+align.ref_start << "\t" << align.aligned_ref_len << "\t+\t" << align.ref_len << "\t" << align.aligned_ref_str << std::endl;
                std::cout << "s\t" << align.query_name << "\t" << 1+align.query_start << "\t" << align.aligned_query_len << "\t" << align.strand << "\t" << align.query_len << "\t" << align.aligned_query_str << std::endl;
                std::cout << std::endl;
    		std::cout.rdbuf(coutbuf); //reset to standard output again
            }
        }
        io_lock.unlock();
        int n = ++num_aligned;
        if (n % 100 == 0) {
            io_lock.lock();
            std::cerr << n << " reads aligned calling forward kernel " << forward_kernel_counter << " times and reverse kernel " << 
									  reverse_kernel_counter << " times\n";
            io_lock.unlock();
        }
    // } end of num reads

    delete[] bin_count_offset_array;
    delete[] nz_bins_array;
    delete[] candidate_hit_offset;
}


int main(int argc, char *argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: ./darwin_ref_guided <READS.fasta> <REFERENCE1>.fasta <REFERECE2>.fasta ..."<< endl;
        exit(1);
    }
    struct timeval start, end_time;


    gettimeofday(&start, NULL);
    std::cerr << "\nReading configuration file ..." << std::endl;
    std::ifstream config_file("params.cfg");
    if (!config_file) {
        std::cerr << "Configuration file <params.cfg> not found! Exiting." << std::endl;
        exit(1);
    }
    ConfigFile cfg("params.cfg");

    // GACT scoring
    int sub_N = cfg.Value("GACT_scoring", "sub_N");
    for (int i = 0; i < 25; i++) {
        gact_sub_mat[i] = sub_N;
    }
    gact_sub_mat[0] = cfg.Value("GACT_scoring", "sub_AA");
    gact_sub_mat[1] = gact_sub_mat[5]  = cfg.Value("GACT_scoring", "sub_AC");
    gact_sub_mat[2] = gact_sub_mat[10] = cfg.Value("GACT_scoring", "sub_AG");
    gact_sub_mat[3] = gact_sub_mat[15] = cfg.Value("GACT_scoring", "sub_AT");
    gact_sub_mat[6] = cfg.Value("GACT_scoring", "sub_CC");
    gact_sub_mat[7] = gact_sub_mat[11] = cfg.Value("GACT_scoring", "sub_CG");
    gact_sub_mat[8] = gact_sub_mat[16] = cfg.Value("GACT_scoring", "sub_CT");
    gact_sub_mat[12] = cfg.Value("GACT_scoring", "sub_GG");
    gact_sub_mat[13] = gact_sub_mat[17] = cfg.Value("GACT_scoring", "sub_GT");
    gact_sub_mat[18] = cfg.Value("GACT_scoring", "sub_TT");
    gap_open        = cfg.Value("GACT_scoring", "gap_open");
    gap_extend      = cfg.Value("GACT_scoring", "gap_extend");

    // D-SOFT parameters
    seed_shape_str          = (std::string) cfg.Value("DSOFT_params", "seed_shape");
    bin_size                = cfg.Value("DSOFT_params", "bin_size");
    dsoft_threshold         = cfg.Value("DSOFT_params", "threshold");
    num_seeds               = cfg.Value("DSOFT_params", "num_seeds");
    seed_occurence_multiple = cfg.Value("DSOFT_params", "seed_occurence_multiple");
    max_candidates          = cfg.Value("DSOFT_params", "max_candidates");
    num_nz_bins             = cfg.Value("DSOFT_params", "num_nz_bins");

    // GACT first tile
    first_tile_size            = cfg.Value("GACT_first_tile", "first_tile_size");
    first_tile_score_threshold = cfg.Value("GACT_first_tile", "first_tile_score_threshold");
    std::cerr << "Running with configuration" << " bin size=" << bin_size << " threshold=" << dsoft_threshold << " ft_size=" << first_tile_size << " ft_threshold=" << first_tile_score_threshold << endl;

    // GACT extend
    tile_size    = cfg.Value("GACT_extend", "tile_size");
    tile_overlap = cfg.Value("GACT_extend", "tile_overlap");

    // Multi-threading
    num_references = cfg.Value("Multithreading", "num_threads");

    seed_shape = seed_shape_str.c_str();

    for(int j=0;j<num_references; j++) {
    	reference_filenames[j].assign(argv[2+j]);
    }

    std::string reads_dirname(argv[1]);


#ifdef FPGA
    if(!init()) {
        std::cerr << "Device initialization failed \n";
        exit(1);
    }
#endif

    gettimeofday(&end_time, NULL);
    long useconds = end_time.tv_usec - start.tv_usec;
    long seconds = end_time.tv_sec - start.tv_sec;
    long mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (reading configuration file): " << mseconds <<" msec" << std::endl;

  for(int k=0; k < num_references; k++) {
    // LOAD REFERENCE
    std::cerr << "\nLoading reference genome" << k << "...\n"<< std::endl;
    gettimeofday(&start, NULL);

    ParseFastaFile(reference_filenames[k], reference_descrips[k], reference_seqs[k], reference_lengths[k], reference_fileposs[k]);

    reference_string[k] = "";

    int curr_bin = 0;

    for (size_t i=0; i < reference_seqs[k].size(); i++) {
        chr_id_to_start_bin[k][i] =  curr_bin;
        reference_string[k] += reference_seqs[k][i];
        for (size_t j = 0; j < (reference_seqs[k][i].length() / bin_size); j++) {
            bin_to_chr_id[k][curr_bin++] = i;
        }
        if (reference_seqs[k][i].length() % bin_size > 0) {
            reference_string[k] += std::string((bin_size - (reference_seqs[k][i].length() % bin_size)), 'N');
            bin_to_chr_id[k][curr_bin++] = i;
        }
    }

    reference_length[k] = reference_string[k].length();
    reference_char[k] = (char*) reference_string[k].c_str();
    std::cerr << "Reference length (after padding): " << (unsigned int) reference_length[k] << std::endl;

    gettimeofday(&end_time, NULL);
    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (loading reference genome): " << mseconds <<" msec" << std::endl;

    // CONSTRUCT SEED POSITION TABLE
    gettimeofday(&start, NULL);
       std::cerr << "\nConstructing seed position table ...\n";
       sa[k] = new SeedPosTable(reference_char[k], reference_length[k], seed_shape, seed_occurence_multiple, bin_size);
    gettimeofday(&end_time, NULL);
    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cerr << "Time elapsed (seed position table construction): " << mseconds <<" msec" << std::endl;

  }
  
    int reads_started = 0;

        //  Iterate through all files and directories on command line.
    argp_parse(&argp, argc, argv, 0, 0, NULL);
    if(NULL == args.output){
        args.output = stdout;
    }

        glob_t globbuf;
        {
            // Find all files matching commandline argument using system glob
            const size_t rootlen = strlen(reads_dirname.c_str());
            char * globpath = calloc(rootlen + 9, sizeof(char));
            memcpy(globpath, reads_dirname.c_str(), rootlen * sizeof(char));
            {
                DIR * dirp = opendir(reads_dirname.c_str());
                if(NULL != dirp){
                    // If filename is a directory, add wildcard to find all fast5 files within it
                    memcpy(globpath + rootlen, "/*.fast5", 8 * sizeof(char));
                    closedir(dirp);
                }
            }
            int globret = glob(globpath, GLOB_NOSORT, NULL, &globbuf);
            free(globpath);
            if(0 != globret){
                if(GLOB_NOMATCH == globret){
                    warnx("File or directory \"%s\" does not exist or no fast5 files found.", reads_dirname.c_str());
                }
                globfree(&globbuf);
            }
        }

        for(size_t fn2=0 ; fn2 < globbuf.gl_pathc ; fn2++){
            reads_started += 1;

            char * filename = globbuf.gl_pathv[fn2];
	    //fprintf(stderr,"basecalling %s\n",filename);
            struct _raw_basecall_info res = calculate_post(filename, args.model);
            if(NULL == res.basecall){
                warnx("No basecall returned for %s", filename);
                continue;
            }

            //ParseFastaFile(reads_filename, reads_descrips, reads_seqs, reads_lengths, reads_fileposs);
            read_length = strlen(res.basecall);		
            std::string reads_seqs(res.basecall);
            //int num_reads = reads_seqs.size();
            int num_reads = 1;
            std::string rev_read = RevComp(reads_seqs);

            reads_char = new char*[num_reads];
            rev_reads_char = new char*[num_reads];
            for (int i =0; i < num_reads; i++) {
                reads_char[i] = (char*) reads_seqs.c_str();
                rev_reads_char[i] = (char*) rev_read.c_str();
            }
            // RUN D-SOFT TO MAP READS
            gettimeofday(&start, NULL);
            std::cerr << "\nFinding candidate bin locations for each read: " << std::endl;
            std::vector<std::thread> align_threads;
            for (int k = 0; k < num_references; k++) {
                align_threads.push_back(std::thread(AlignRead, k, num_reads, num_references));
            }
            std::cerr << "Using " << align_threads.size() << " threads ...\n";
            for (auto& th : align_threads) th.join();
            std::cerr << "Synchronizing threads. " << num_aligned << " reads aligned\n";
            gettimeofday(&end_time, NULL);
            useconds = end_time.tv_usec - start.tv_usec;
            seconds = end_time.tv_sec - start.tv_sec;
            mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
            std::cerr << "Time elapsed (seed table querying): " << mseconds <<" msec" << std::endl;

            //fprintf_format(args.outformat, args.output, res.rt.uuid, basename(filename), args.uuid, args.prefix, res);

            //write_summary(hdf5out, args.uuid ? res.rt.uuid : basename(filename), res,
             //             args.compression_chunk_size, args.compression_level);


            free_raw_basecall_info(&res);
        }
        globfree(&globbuf);




#ifdef FPGA
    cleanup();
#endif
    return 0;
}


#ifdef FPGA
/////// OPENCL HELPER FUNCTIONS ///////

bool init() {
  cl_int status;

  if(!setCwdToExeDir()) {
    return false;
  }

  // Get the OpenCL platform.
  platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
  if(platform == NULL) {
    printf("ERROR: Unable to find Intel(R) FPGA OpenCL platform.\n");
    return false;
  }

  // Query the available OpenCL devices.
  devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &max_num_devices));
  printf("Number of devices is %d\n", max_num_devices);

  // Create the context.
  //context = clCreateContext(NULL, num_devices, devices, NULL, NULL, &status);
  context = clCreateContext(NULL, num_devices, &devices[0], NULL, NULL, &status);
  checkError(status, "Failed to create context");

  queues.reset(num_devices*num_reference);
  kernels.reset(num_devices*num_threads);
  kernel_events.reset(num_devices*num_threads*20);

  // Create command queue.
  for(int i = 0; i<num_devices ; i++) {
    for (int j = 0; j<num_threads ; j++){
          int id = (i*num_devices)+j;
          //queues[id] = clCreateCommandQueue(context, devices[i], CL_QUEUE_PROFILING_ENABLE, &status);
          queues[id] = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &status);
          checkError(status, "Failed to create command queue");
    }
  }


  // Create the program.
  std::string binary_file = getBoardBinaryFile(PRECOMPILED_BINARY, devices[0]);
  printf("Using AOCX: %s\n", binary_file.c_str());
  program = createProgramFromBinary(context, binary_file.c_str(), devices, num_devices);

  // Build the program that was just created.
  status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
  checkError(status, "Failed to build program");

  for(int i= 0; i<num_devices ; i++) {
  for (int j = 0; j<num_threads ; j++){
        // Create the kernel - name passed in here must match kernel name in the
        // original CL file, that was compiled into an AOCX file using the AOC tool
        const char * kernel_name = "xl";  // Kernel name, as defined in the CL file
        kernels[(i*num_devices)+j] = clCreateKernel(program, kernel_name, &status);
        checkError(status, "Failed to create kernel");
  }
  }

  return true;
}


void cleanup() {
  for(int i = 0; i < num_devices; i++) {
    for(unsigned j = 0; j < num_threads; j++) {
      int id = (i*num_devices)+j;
      if(kernels[id]) {
        clReleaseKernel(kernels[id]);
      }
      if(queues[id]) {
        clReleaseCommandQueue(queues[id]);
      }
    }
  }
  if(program) {
    clReleaseProgram(program);
  }
  if(context) {
    clReleaseContext(context);
  }
}

// Helper functions to display parameters returned by OpenCL queries
static void device_info_ulong( cl_device_id device, cl_device_info param, const char* name) {
   cl_ulong a;
   clGetDeviceInfo(device, param, sizeof(cl_ulong), &a, NULL);
   printf("%-40s = %lu\n", name, a);
}
static void device_info_uint( cl_device_id device, cl_device_info param, const char* name) {
   cl_uint a;
   clGetDeviceInfo(device, param, sizeof(cl_uint), &a, NULL);
   printf("%-40s = %u\n", name, a);
}
static void device_info_bool( cl_device_id device, cl_device_info param, const char* name) {
   cl_bool a;
   clGetDeviceInfo(device, param, sizeof(cl_bool), &a, NULL);
   printf("%-40s = %s\n", name, (a?"true":"false"));
}
static void device_info_string( cl_device_id device, cl_device_info param, const char* name) {
   char a[STRING_BUFFER_LEN];
   clGetDeviceInfo(device, param, STRING_BUFFER_LEN, &a, NULL);
   printf("%-40s = %s\n", name, a);
}


// Query and display OpenCL information on device and runtime environment
static void display_device_info( cl_device_id device ) {

   printf("Querying device for info:\n");
   printf("========================\n");
   device_info_string(device, CL_DEVICE_NAME, "CL_DEVICE_NAME");
   device_info_string(device, CL_DEVICE_VENDOR, "CL_DEVICE_VENDOR");
   device_info_uint(device, CL_DEVICE_VENDOR_ID, "CL_DEVICE_VENDOR_ID");
   device_info_string(device, CL_DEVICE_VERSION, "CL_DEVICE_VERSION");
   device_info_string(device, CL_DRIVER_VERSION, "CL_DRIVER_VERSION");
   device_info_uint(device, CL_DEVICE_ADDRESS_BITS, "CL_DEVICE_ADDRESS_BITS");
   device_info_bool(device, CL_DEVICE_AVAILABLE, "CL_DEVICE_AVAILABLE");
   device_info_bool(device, CL_DEVICE_ENDIAN_LITTLE, "CL_DEVICE_ENDIAN_LITTLE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, "CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE");
   device_info_ulong(device, CL_DEVICE_GLOBAL_MEM_SIZE, "CL_DEVICE_GLOBAL_MEM_SIZE");
   device_info_bool(device, CL_DEVICE_IMAGE_SUPPORT, "CL_DEVICE_IMAGE_SUPPORT");
   device_info_ulong(device, CL_DEVICE_LOCAL_MEM_SIZE, "CL_DEVICE_LOCAL_MEM_SIZE");
   device_info_ulong(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, "CL_DEVICE_MAX_CLOCK_FREQUENCY");
   device_info_ulong(device, CL_DEVICE_MAX_COMPUTE_UNITS, "CL_DEVICE_MAX_COMPUTE_UNITS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_ARGS, "CL_DEVICE_MAX_CONSTANT_ARGS");
   device_info_ulong(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE");
   device_info_uint(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MEM_BASE_ADDR_ALIGN, "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS");
   device_info_uint(device, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE, "CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT");
   device_info_uint(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, "CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE");

   {
      cl_command_queue_properties ccp;
      clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(cl_command_queue_properties), &ccp, NULL);
      printf("%-40s = %s\n", "Command queue out of order? ", ((ccp & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)?"true":"false"));
      printf("%-40s = %s\n", "Command queue profiling enabled? ", ((ccp & CL_QUEUE_PROFILING_ENABLE)?"true":"false"));
   }
}

/* Error checking */
void checkErr(cl_int err, const char * name)
{
        if (err != CL_SUCCESS) {
                printf("\n ERROR (%d): %s\n",err,name);
                exit(EXIT_FAILURE);
        }
}

#endif // end OpenCL helper

