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

#include "gact.h"
#include <math.h>


enum states {Z, D, I, M};
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};

extern int INF;

int dequeue(int *BT_states, int *front, int *queuesize) {
   int data = BT_states[(*front)++];
	
   if(*front == MAX_TILE_SIZE) {
      *front = 0;
   }
	
   (*queuesize)--;
   return data;  
}


std::queue<int> CpuXL(char* ref_str, long long int ref_tile_length, char* query_str, long long int query_tile_length, int* sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first_tile, int early_terminate, int* BT_states, int *queuesize, char strand, int *rear, int *front);

#ifdef FPGA
std::queue<int> FpgaXL(char* ref_str, long long int ref_tile_length, char* query_str, long long int query_tile_length, int* sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first_tile, int early_terminate, int* BT_states, int *queuesize, char strand, int *rear, int *front, int thread_id); 
#endif

Alignment GACT (char* ref_str, char* query_str, std::string ref_name, int* sub_mat, int gap_open, int gap_extend, int tile_size, int tile_overlap, int ref_pos, int query_pos, uint32_t ref_length, uint32_t query_length, char strand, int first_tile_score_threshold, int mode, int thread_id) {
    std::queue<int> BT_states_std;

    std::string aligned_ref_str = "";
    std::string aligned_query_str = "";

    Alignment alignment;
    alignment.ref_name = ref_name;
    //alignment.query_name = query_name;
    alignment.aligned_ref_str = "";
    alignment.aligned_query_str = "";
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;
    alignment.aligned_ref_len = 0;
    alignment.aligned_query_len = 0;
    alignment.ref_len = ref_length;
    alignment.query_len = query_length;
    alignment.strand = strand;
    alignment.score = 0;
    alignment.flag = 0;
    
    int ref_tile_length = tile_size;
    int query_tile_length = tile_size;
    int rev_ref_pos = ref_pos;
    int rev_query_pos = query_pos;
   
    int max_ref_pos = 0;
    int max_query_pos = 0;
    int i = 0;
    int j = 0;
    
    int first_tile_score = 0;
    bool first_tile = true;

    
    while ((ref_pos > 0) && (query_pos > 0) && (((i > 0) && (j > 0)) || first_tile)) {
        //change the tile length if elements less than that of the tile size
        ref_tile_length = (ref_pos > tile_size) ? tile_size : ref_pos;
        query_tile_length = (query_pos > tile_size) ? tile_size : query_pos;

    	int *BT_states = (int *) malloc(sizeof(int)*MAX_TILE_SIZE);
    	int *front = (int *) malloc(sizeof(int));
    	int *rear = (int *) malloc(sizeof(int));
    	int *queuesize = (int *) malloc(sizeof(int));
    	*front = 0; *rear=-1; *queuesize=0; 

        if(mode==1) { 
       	       BT_states_std =  CpuXL (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, ref_tile_length, query_tile_length, false, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front);
        }
#ifdef FPGA
        else if(mode==2) { 
       	       BT_states_std = FpgaXL (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, ref_tile_length, query_tile_length, false, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, thread_id);
        }
#endif
        else {
        	BT_states_std = AlignWithBT (ref_str+(ref_pos-ref_tile_length), ref_tile_length, 
				query_str+(query_pos-query_tile_length), query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, false, 
				first_tile, (tile_size - tile_overlap));
	}

        i = 0;
        j = 0;
        int tile_score = BT_states_std.front();
        BT_states_std.pop();

        
        if (first_tile) {
            ref_pos = ref_pos - ref_tile_length + BT_states_std.front();
            max_ref_pos = BT_states_std.front(); BT_states_std.pop();
            query_pos = query_pos - query_tile_length + BT_states_std.front();
            max_query_pos = BT_states_std.front(); BT_states_std.pop();
            rev_ref_pos = ref_pos;
            rev_query_pos = query_pos;
            first_tile_score = tile_score;
        }

        int num_tb = BT_states_std.size();
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = num_tb-1;
        int query_buf_curr = num_tb-1;

        while (!BT_states_std.empty()) {
            first_tile = false;
            int state = BT_states_std.front(); BT_states_std.pop();
            if (state == M) {
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
                j += 1;
            }
            if (state == I) {
                ref_buf[ref_buf_curr--] = ref_str[ref_pos - j - 1];
                query_buf[query_buf_curr--] = '-';
                j += 1;
            }
            if (state == D) {
                ref_buf[ref_buf_curr--] = '-';
                query_buf[query_buf_curr--] = query_str[query_pos - i - 1];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = std::string(ref_buf, num_tb) + aligned_ref_str;
            aligned_query_str = std::string(query_buf, num_tb) + aligned_query_str;
        }

    	free(BT_states);
    	free(front);
    	free(rear);
    	free(queuesize);
        free(ref_buf);
        free(query_buf);
        ref_pos -= (j);
        query_pos -= (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }
    
    alignment.ref_start = ref_pos;
    alignment.query_start = query_pos;

    ref_pos = rev_ref_pos;
    query_pos = rev_query_pos;
    
    i =  tile_size;
    j = tile_size;
    
    //starts with the first tile
    while ((ref_pos < ref_length) && (query_pos < query_length) && (((i > 0) && (j > 0)) || first_tile)) {
        ref_tile_length = (ref_pos + tile_size < ref_length) ? tile_size : ref_length - ref_pos;
        query_tile_length = (query_pos + tile_size < query_length) ? tile_size : query_length - query_pos;
    	int *BT_states = (int *) malloc(sizeof(int)*MAX_TILE_SIZE);
    	int *front = (int *) malloc(sizeof(int));
    	int *rear = (int *) malloc(sizeof(int));
    	int *queuesize = (int *) malloc(sizeof(int));
    	*front = 0; *rear=-1; *queuesize=0; 

        if(mode == 1) { 
                BT_states_std = CpuXL(ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, ref_tile_length, query_tile_length, true, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front);
        }
#ifdef FPGA
        else if(mode == 2) { 
                BT_states_std = FpgaXL (ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, ref_tile_length, query_tile_length, true, 
				first_tile, (tile_size - tile_overlap), BT_states, queuesize, strand, rear, front, thread_id);
	}
#endif
        else {
        	BT_states_std = AlignWithBT (ref_str+ref_pos, ref_tile_length, query_str+query_pos, query_tile_length, 
				sub_mat, gap_open, gap_extend, query_tile_length, ref_tile_length, true, 
				first_tile, (tile_size - tile_overlap));
	}
        i = 0;
        j = 0;
        int tile_score = BT_states_std.front(); BT_states_std.pop(); 

        if (first_tile) {
            ref_pos = ref_pos + ref_tile_length - BT_states_std.front();
            max_ref_pos = BT_states_std.front(); BT_states_std.pop();
            query_pos = query_pos + query_tile_length - BT_states_std.front();
            max_query_pos = BT_states_std.front(); BT_states_std.pop();
        }

        int num_tb = BT_states_std.size();;
        char* ref_buf = (char*) malloc(num_tb);
        char* query_buf = (char*) malloc(num_tb);
        int ref_buf_curr = 0;
        int query_buf_curr = 0;
        
        while (!BT_states_std.empty()) {
            first_tile = false;
            int state = BT_states_std.front(); BT_states_std.pop();
            if (state == M) {
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
                j += 1;
            }
            if (state == I) {
                ref_buf[ref_buf_curr++] = ref_str[ref_pos + j];
                query_buf[query_buf_curr++] = '-';
                j += 1;
            }
            if (state == D) {
                ref_buf[ref_buf_curr++] = '-';
                query_buf[query_buf_curr++] = query_str[query_pos + i];
                i += 1;
            }
        }
        if (num_tb > 0) {
            aligned_ref_str = aligned_ref_str + std::string(ref_buf, num_tb);
            aligned_query_str = aligned_query_str + std::string(query_buf, num_tb);
        }

    	free(BT_states);
    	free(front);
    	free(rear);
    	free(queuesize);
        free(ref_buf);
        free(query_buf);
        ref_pos += (j);
        query_pos += (i);
        alignment.aligned_ref_len += j;
        alignment.aligned_query_len += i;
    }

    int total_score = 0;
    bool open = true;
    for (uint32_t j = 0; j < aligned_ref_str.length(); j++) {
        char ref_nt = aligned_ref_str[j];
        char query_nt = aligned_query_str[j];
        if (ref_nt == '-' || query_nt == '-') {
            total_score += (open) ? gap_open : gap_extend;
            open = false;
        }
        else {
            total_score += sub_mat[5*NtChar2Int(query_nt) + NtChar2Int(ref_nt)];
            open = true;
        }
    }
    alignment.aligned_ref_str = aligned_ref_str;
    alignment.aligned_query_str = aligned_query_str;
    alignment.score = total_score;
    return alignment;
}


void xl_cpu_score( const char *  a, const short i_start, const char *  b, const short j_start,
                                const short open_extend_gap, const short extend_gap, const short match, const short mismatch,
                                short *maxScore, short *  globalH, short tile_num);

void xl_cpu_direction (short m, const short n, const short gap_open, const short gap_extend,
                                        short *  dir_matrix, short * h_matrix, 
                                        short *h_matrix_rd, short *m_matrix_rd, short *i_matrix_rd, short *d_matrix_rd,
                                        short *h_matrix_wr, short *m_matrix_wr, short *i_matrix_wr, short *d_matrix_wr);

std::queue<int> CpuXL(char* ref_str, long long int ref_tile_length, char* query_str, long long int query_tile_length, int* sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first_tile, int early_terminate, int* BT_states_arg, int *queuesize, char strand, int *rear, int *front) {

    char * a, * a_header, *b, * b_header;
    short m, n, ext_m, ext_n, nbb, open_extend_gap, num_active_devices, * scores[MAX_NUM_DEVICES], score=0;
    int i, ii, j, jj, flag=0, d;
    short match = 1, mismatch = -1;
    short *maxScore=(short *)malloc(3*sizeof(short));
    int pos_score;
    short h_matrix[MAX_TILE_SIZE+1][MAX_TILE_SIZE+1];
    short h_matrix_wr[MAX_TILE_SIZE + 1];
    short m_matrix_wr[MAX_TILE_SIZE + 1];
    short i_matrix_wr[MAX_TILE_SIZE + 1];
    short d_matrix_wr[MAX_TILE_SIZE + 1];
    short h_matrix_rd[MAX_TILE_SIZE + 1];
    short m_matrix_rd[MAX_TILE_SIZE + 1];
    short i_matrix_rd[MAX_TILE_SIZE + 1];
    short d_matrix_rd[MAX_TILE_SIZE + 1];
    short dir_matrix[MAX_TILE_SIZE+1][MAX_TILE_SIZE+1];
    for (int i = 0; i < m + 1; i++) {
      h_matrix_rd[i] = 0; m_matrix_rd[i] = 0; i_matrix_rd[i] = -INF; d_matrix_rd[i] = -INF;
      h_matrix_wr[i] = 0; m_matrix_wr[i] = 0; i_matrix_wr[i] = -INF; d_matrix_wr[i] = -INF;
    }
    for (int i = 0; i < n + 1; i++) {
      dir_matrix[i][0] = ZERO_OP;
    }
    for (int j = 0; j < m + 1; j++) {
      dir_matrix[0][j] = ZERO_OP;
    }

    int reverse_int = reverse, first_tile_int = first_tile;
    
    a = ref_str;
    b = query_str;
    m = ref_tile_length;
    n = query_tile_length;
    ext_m = m;
    ext_n = ceil((double) n / TILE_SIZE) * TILE_SIZE;
    nbb = ext_n / TILE_SIZE;

    for(jj=0; jj<nbb ; jj++) {
  	xl_cpu_score( a, m, b, n, gap_open,  gap_extend,  match, mismatch, maxScore, (short *)h_matrix, jj);
  	//xl_cpu_score( a, m, b, n, OPEN_GAP,  EXTEND_GAP,  MATCH, MISMATCH, maxScore, (short *)h_matrix, jj);
    }

    xl_cpu_direction(m, n, gap_open, gap_extend, (short *)dir_matrix, (short *)h_matrix,
						 h_matrix_rd, m_matrix_rd, i_matrix_rd, d_matrix_rd,
                                                 h_matrix_wr, m_matrix_wr, i_matrix_wr, d_matrix_wr);
    std::queue<int> BT_states;

    int i_curr=ref_pos, j_curr=query_pos;
    int i_steps = 0, j_steps = 0;

    int open = 0;
    if (first_tile) {
      i_curr = maxScore[1];
      j_curr = maxScore[2];
      BT_states.push(maxScore[0]);
      BT_states.push(i_curr);
      BT_states.push(j_curr);
    }
    else {
      BT_states.push(pos_score);
    }

    int state = dir_matrix[i_curr][j_curr] % 4;

    while (state != Z) {
      if ((i_steps >= early_terminate) || (j_steps >= early_terminate)) {
          break;
      }
      BT_states.push(state);
      if (state == M) {
          state = (dir_matrix[i_curr-1][j_curr-1] % 4);
          i_curr--;
          j_curr--;
          i_steps++;
          j_steps++;
      }
      else if (state == I) {
          state = (dir_matrix[i_curr][j_curr] & (2 << INSERT_OP)) ? M : I;
          i_curr--;
          i_steps++;
      }
      else if (state == D) {
          state = (dir_matrix[i_curr][j_curr] & (2 << DELETE_OP)) ? M : D;
          j_curr--;
          j_steps++;
      }
    };

    return BT_states;

}



#ifdef FPGA


std::queue<int> FpgaXL(char* a, long long int m, char* b, long long int n, int* sub_mat, int gap_open, int gap_extend, int ref_pos, int query_pos, bool reverse, bool first, int early_terminate, int* BT_states_arg, int *queuesize, char strand, int *rear, int *front, int thread_id) {

    cl_int status;
    cl_mem cl_a[MAX_NUM_DEVICES], cl_b[MAX_NUM_DEVICES], cl_scores[MAX_NUM_DEVICES], cl_prev_maxRow[MAX_NUM_DEVICES], cl_next_maxRow[MAX_NUM_DEVICES];
    cl_mem cl_bt_states[MAX_NUM_DEVICES], cl_queuesize[MAX_NUM_DEVICES], cl_rear[MAX_NUM_DEVICES];
    cl_mem cl_prev_lastCol[MAX_NUM_DEVICES], cl_next_lastCol[MAX_NUM_DEVICES];
    cl_mem cl_sub_mat[MAX_NUM_DEVICES];
    cl_mem cl_dir_matrix[MAX_NUM_DEVICES], cl_max_score[MAX_NUM_DEVICES], cl_max_i[MAX_NUM_DEVICES], cl_max_j[MAX_NUM_DEVICES], cl_pos_score[MAX_NUM_DEVICES];

    size_t wgSize[3] = {1, 1, 1};
    size_t gSize[3] = {1, 1, 1};

    int d = thread_id;
    int score;
    int nbb=1, jj=0;

    int dir_matrix[MAX_TILE_SIZE+1][MAX_TILE_SIZE+1];
    int max_score = 0, max_i = 0, max_j = 0, pos_score = 0;

    for (int i = 0; i < m+1; i++)
           dir_matrix[i][0] = ZERO_OP;
    for (int j = 0; j < n + 1; j++)
           dir_matrix[0][j] = ZERO_OP;


                cl_a[d] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, m* sizeof(char), a, &status);
                checkErr(status,"clCreateBuffer cl_a");
                cl_b[d] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, n* sizeof(char), b, &status);
                checkErr(status,"clCreateBuffer cl_b");
                cl_prev_lastCol[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
                checkErr(status,"clCreateBuffer cl_prev_lastCol");
                cl_next_lastCol[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
                checkErr(status,"clCreateBuffer cl_next_lastCol");
                cl_prev_maxRow[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
                checkErr(status,"clCreateBuffer cl_prev_maxRow");
                cl_next_maxRow[d] = clCreateBuffer(context, CL_MEM_READ_WRITE, m*sizeof(int), NULL, &status);
                checkErr(status,"clCreateBuffer cl_next_maxRow");
                cl_scores[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), &score, &status);
                checkErr(status,"clCreateBuffer cl_scores");
                cl_bt_states[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*MAX_TILE_SIZE, BT_states_arg, &status);
                checkErr(status,"clCreateBuffer cl_bt_states");
                cl_queuesize[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), queuesize, &status);
                checkErr(status,"clCreateBuffer cl_queuesize");
                cl_rear[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), rear, &status);
                checkErr(status,"clCreateBuffer cl_rear");
                cl_sub_mat[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*25, sub_mat, &status);
                checkErr(status,"clCreateBuffer cl_sub_mat");
                cl_dir_matrix[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int)*(MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1), (int *)dir_matrix, &status);
                checkErr(status,"clCreateBuffer cl_dir_matric");
                cl_max_score[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), &max_score, &status);
                checkErr(status,"clCreateBuffer cl_max_score");
                cl_max_i[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), &max_i, &status);
                checkErr(status,"clCreateBuffer max_i");
                cl_max_j[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), &max_j, &status);
                checkErr(status,"clCreateBuffer max_j");
                cl_pos_score[d] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(int), &pos_score, &status);
                checkErr(status,"clCreateBuffer cl_pos_score");


                status = clSetKernelArg(kernels[d], 0, sizeof(cl_mem), &cl_a[d]);
                checkError(status, "Failed to set kernels[d] arg 0");

                status = clSetKernelArg(kernels[d], 1, sizeof(int), &m);
                checkError(status, "Failed to set kernels[d] arg 1");

                status = clSetKernelArg(kernels[d], 2, sizeof(cl_mem), &cl_b[d]);
                checkError(status, "Failed to set kernels[d] arg 2");

                status = clSetKernelArg(kernels[d], 3, sizeof(int), &n);
                checkError(status, "Failed to set kernels[d] arg 3");

                status = clSetKernelArg(kernels[d], 4, sizeof(int), &gap_open);
                checkError(status, "Failed to set kernels[d] arg 4");

                status = clSetKernelArg(kernels[d], 5, sizeof(int), &gap_extend);
                checkError(status, "Failed to set kernels[d] arg 5");

                status = clSetKernelArg(kernels[d], 6, sizeof(int), &jj);
                checkError(status, "Failed to set kernels[d] arg 6");

                status = clSetKernelArg(kernels[d], 7, sizeof(int), &query_pos);
                checkError(status, "Failed to set kernels[d] arg 7");
                status = clSetKernelArg(kernels[d], 8, sizeof(int), &ref_pos);
                checkError(status, "Failed to set kernels[d] arg 8");
                int reverse_int = reverse, first_tile_int = first;
                status = clSetKernelArg(kernels[d], 9, sizeof(int), &reverse_int);
                checkError(status, "Failed to set kernels[d] arg 9");
                status = clSetKernelArg(kernels[d], 10, sizeof(int), &first_tile_int);
                checkError(status, "Failed to set kernels[d] arg 10");
                status = clSetKernelArg(kernels[d], 11, sizeof(int), &early_terminate);
                checkError(status, "Failed to set kernels[d] arg 11");
                status = clSetKernelArg(kernels[d], 12, sizeof(cl_mem), &cl_sub_mat[d]);
                checkError(status, "Failed to set kernels[d] arg 12");

                status = clSetKernelArg(kernels[d], 13, sizeof(cl_mem), &cl_dir_matrix[d]);
                checkError(status, "Failed to set kernels[d] arg 13");
                status = clSetKernelArg(kernels[d], 14, sizeof(cl_mem), &cl_max_score[d]);
                checkError(status, "Failed to set kernels[d] arg 14");
                status = clSetKernelArg(kernels[d], 15, sizeof(cl_mem), &cl_max_i[d]);
                checkError(status, "Failed to set kernels[d] arg 15");
                status = clSetKernelArg(kernels[d], 16, sizeof(cl_mem), &cl_max_j[d]);
                checkError(status, "Failed to set kernels[d] arg 16");
                status = clSetKernelArg(kernels[d], 17, sizeof(cl_mem), &cl_pos_score[d]);
                checkError(status, "Failed to set kernels[d] arg 17");


                status = clEnqueueNDRangeKernel(queues[d], kernels[d], 1, NULL, gSize, wgSize, 0, NULL, NULL);
                checkError(status, "Failed to launch kernel");

                status = clFinish(queues[d]);
                checkError(status, "Failed to finish");

                status = clEnqueueReadBuffer(queues[d], cl_dir_matrix[d], CL_TRUE, 0, sizeof(int)*(MAX_TILE_SIZE+1)*(MAX_TILE_SIZE+1), (int *)dir_matrix, 0, NULL, NULL);
                checkErr(status,"clEnqueueReadBuffer: Couldn't read cl_dir_matrix buffer");
                status = clEnqueueReadBuffer(queues[d], cl_max_score[d], CL_TRUE, 0, sizeof(int), &max_score, 0, NULL, NULL);
                checkErr(status,"clEnqueueReadBuffer: Couldn't read cl_max_score buffer");
                status = clEnqueueReadBuffer(queues[d], cl_max_i[d], CL_TRUE, 0, sizeof(int), &max_i, 0, NULL, NULL);
                checkErr(status,"clEnqueueReadBuffer: Couldn't read cl_max_i buffer");
                status = clEnqueueReadBuffer(queues[d], cl_max_j[d], CL_TRUE, 0, sizeof(int), &max_j, 0, NULL, NULL);
                checkErr(status,"clEnqueueReadBuffer: Couldn't read cl_max_j buffer");
                status = clEnqueueReadBuffer(queues[d], cl_pos_score[d], CL_TRUE, 0, sizeof(int), &pos_score, 0, NULL, NULL);
                checkErr(status,"clEnqueueReadBuffer: Couldn't read cl_pos_score buffer");

                clReleaseMemObject(cl_a[d]);
                clReleaseMemObject(cl_b[d]);
                clReleaseMemObject(cl_prev_maxRow[d]);
                clReleaseMemObject(cl_next_maxRow[d]);
                clReleaseMemObject(cl_prev_lastCol[d]);
                clReleaseMemObject(cl_next_lastCol[d]);
                clReleaseMemObject(cl_scores[d]);
                clReleaseMemObject(cl_bt_states[d]);
                clReleaseMemObject(cl_queuesize[d]);
                clReleaseMemObject(cl_rear[d]);
                clReleaseMemObject(cl_dir_matrix[d]);
                clReleaseMemObject(cl_max_score[d]);
                clReleaseMemObject(cl_max_j[d]);
                clReleaseMemObject(cl_max_i[d]);
                clReleaseMemObject(cl_pos_score[d]);

  std::queue<int> BT_states;

  int i_curr=ref_pos, j_curr=query_pos;
  int i_steps = 0, j_steps = 0;

  int open = 0;
  if (first) {
      i_curr = max_i;
      j_curr = max_j;
      BT_states.push(max_score);
      BT_states.push(i_curr);
      BT_states.push(j_curr);
  }
  else {
      BT_states.push(pos_score);
  }

  int state = dir_matrix[i_curr][j_curr] % 4;

  while (state != Z) {
      if ((i_steps >= early_terminate) || (j_steps >= early_terminate)) {
          break;
      }
      BT_states.push(state);
      if (state == M) {
          state = (dir_matrix[i_curr-1][j_curr-1] % 4);
          i_curr--;
          j_curr--;
          i_steps++;
          j_steps++;
      }
      else if (state == I) {
          state = (dir_matrix[i_curr][j_curr] & (2 << INSERT_OP)) ? M : I;
          i_curr--;
          i_steps++;
      }
      else if (state == D) {
          state = (dir_matrix[i_curr][j_curr] & (2 << DELETE_OP)) ? M : D;
          j_curr--;
          j_steps++;
      }
  };

  return BT_states;

}


#endif

