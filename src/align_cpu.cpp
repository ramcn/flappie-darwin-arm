#include <queue>
#define BLOCK_WIDTH 128 // update this value in host/src/arguments.h in case of modification
#define MAX_TILE_SIZE 512

extern int INF;
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
enum states {Z, D, I, M};

#define A_NT 0
#define C_NT 1
#define G_NT 2
#define T_NT 3
#define N_NT 4

int max(int a, int b) {
	if(a>b) return a;
	else return b;
}

extern unsigned int NtChar2Int (char nt);


void xl_cpu ( const char *  a, const int m, 
					 const char *  b, const int nbb,
					const int gap_open, const int gap_extend, const int match, const int mismatch,
					 const int *  prev_lastCol,   int *  next_lastCol,
					 const int *  prev_maxRow,  int *  next_maxRow,
					 int *  maxScore, const int jj,
					 int *  BT_states_local,  int *  qsize,  int *  rear,
					const int query_pos, const int ref_pos, const int reverse, const int first, const int early_terminate,
					 int *  sub_mat,
                			 int * dir_matrix_arg, int * max_score_arg, int * max_i_arg, int * max_j_arg, int * pos_score_arg)
{
	// auxiliar buffers
	int row1[BLOCK_WIDTH]={0}, maxCol1[BLOCK_WIDTH]={0};
	int score=0, auxLastCol=0;

        // MICRO
  	int h_matrix_wr[MAX_TILE_SIZE + 1], m_matrix_wr[MAX_TILE_SIZE + 1], i_matrix_wr[MAX_TILE_SIZE + 1], d_matrix_wr[MAX_TILE_SIZE + 1];
  	int h_matrix_rd[MAX_TILE_SIZE + 1], m_matrix_rd[MAX_TILE_SIZE + 1], i_matrix_rd[MAX_TILE_SIZE + 1], d_matrix_rd[MAX_TILE_SIZE + 1];
  	int dir_matrix[MAX_TILE_SIZE+1][BLOCK_WIDTH+1];

  	int local_sub_mat[25];

  	for (int i = 0; i < 25; i++) {
       		local_sub_mat[i] = 0;
  	}
  	local_sub_mat[0] = 1; local_sub_mat[1] = -1; local_sub_mat[2] = -1; local_sub_mat[3] = -1; local_sub_mat[5] = -1; local_sub_mat[6] = 1; 
	local_sub_mat[7] = -1; local_sub_mat[8] = -1; local_sub_mat[10] = -1; local_sub_mat[11] = -1; local_sub_mat[12] = 1; 
	local_sub_mat[13] = -1; local_sub_mat[15] = -1; local_sub_mat[16] = -1; local_sub_mat[17] = -1; local_sub_mat[18] = 1;


	// private buffer for query sequence 
	char private_b1[BLOCK_WIDTH];

	// pointer to query sequence 
	const char * ptr_b = b + jj*BLOCK_WIDTH;

	// copy sequence b ty private memory
	for(int i = 0; i < BLOCK_WIDTH; i++)
		private_b1[i] = ptr_b[i];

	// MICRO
  	for (int i = 0; i < BLOCK_WIDTH + 1; i++) {
    		h_matrix_rd[i] = 0; m_matrix_rd[i] = 0; i_matrix_rd[i] = -INF; d_matrix_rd[i] = -INF;
		h_matrix_wr[i] = 0; m_matrix_wr[i] = 0; i_matrix_wr[i] = -INF; d_matrix_wr[i] = -INF;
  	}
  	for (int i = 0; i < m+1; i++) 
      		dir_matrix[i][0] = ZERO_OP;
  	for (int j = 0; j < BLOCK_WIDTH + 1; j++) 
      		dir_matrix[0][j] = ZERO_OP;

	// MICRO
  	int max_score = 0;
  	int pos_score = 0;
  	int max_i = 0;
  	int max_j = 0;

	for(int i = 1; i < m+1; i++){
		// MICRO
      		for (int k = 1; k < MAX_TILE_SIZE + 1; k++) {
          		m_matrix_rd[k] = m_matrix_wr[k]; h_matrix_rd[k] = h_matrix_wr[k];
          		i_matrix_rd[k] = i_matrix_wr[k]; d_matrix_rd[k] = d_matrix_wr[k];
      		}

		#pragma unroll
		for (int j=1; j < BLOCK_WIDTH+1 ; j++){

          		int ref_nt = (reverse) ? NtChar2Int(a[m-i]) : NtChar2Int(a[i-1]);
          		int query_nt = (reverse) ? NtChar2Int(b[BLOCK_WIDTH-j]) : NtChar2Int(b[j-1]);
          		int match;
          		//case of unknown nucleotide in either reference or query
          		if (ref_nt == N_NT || query_nt == N_NT) {
              			match = -INF; // Force N's to align with gaps
          		} else {
              			//value from the W matrix for the match/mismatch penalty/point
              			match = sub_mat[query_nt*5 + ref_nt];
          		}
          		//columnwise calculations
          		if (m_matrix_rd[j-1] > i_matrix_rd[j-1] && m_matrix_rd[j-1] > d_matrix_rd[j-1]) {
              			m_matrix_wr[j] = m_matrix_rd[j-1] + match;
         		 } else if (i_matrix_rd[j-1] > d_matrix_rd[j-1]) {
              			m_matrix_wr[j] = i_matrix_rd[j-1] + match;
          		} else {
              			m_matrix_wr[j] = d_matrix_rd[j-1] + match;
          		}
          		if (m_matrix_wr[j] < 0) {
              			m_matrix_wr[j] = 0;
          		}
          		int ins_open   = m_matrix_rd[j] + gap_open;
          		int ins_extend = i_matrix_rd[j] + gap_extend;
          		int del_open   = m_matrix_wr[j-1] + gap_open;
          		int del_extend = d_matrix_wr[j-1] + gap_extend;
          		i_matrix_wr[j] = (ins_open > ins_extend) ? ins_open : ins_extend;
          		d_matrix_wr[j] = (del_open > del_extend) ? del_open : del_extend;
         	 	int max1 = m_matrix_wr[j] > i_matrix_wr[j] ? m_matrix_wr[j] : i_matrix_wr[j];
          		int max2 = d_matrix_wr[j] > 0 ? d_matrix_wr[j] : 0;
          		h_matrix_wr[j] = max1 > max2 ? max1 : max2;
          		(dir_matrix)[i][j] = ((m_matrix_wr[j] >= i_matrix_wr[j]) ? ((m_matrix_wr[j] >= d_matrix_wr[j]) ? MATCH_OP : DELETE_OP) :
										 ((i_matrix_wr[j] >= d_matrix_wr[j]) ? INSERT_OP : DELETE_OP));
          		if ((m_matrix_wr[j] <= 0) && (i_matrix_wr[j] <= 0) && (d_matrix_wr[j] <= 0)) {
              			(dir_matrix)[i][j] = ZERO_OP;
          		}
          		(dir_matrix)[i][j] += (ins_open >= ins_extend) ? (2 << INSERT_OP) : 0;
          		(dir_matrix)[i][j] += (del_open >= del_extend) ? (2 << DELETE_OP) : 0;
          		if (h_matrix_wr[j] >= max_score) {
              			max_score = h_matrix_wr[j];
              			max_i = i;
              			max_j = j;
          		}
          		if ((j == query_pos) && (i == ref_pos)) {
              			pos_score = h_matrix_wr[j];
          		}

		} // end pragma unroll

	}
 
	for (int i = 0; i < m+1; i++) 
  	for (int j = 0; j < BLOCK_WIDTH + 1; j++) 
      		dir_matrix_arg[i*(BLOCK_WIDTH+1) + j] = dir_matrix[i][j]; // bug 2 block_width+1. otherwise hang. 

  *max_score_arg = max_score; 
  *max_i_arg = max_i;
  *max_j_arg = max_j;
  *pos_score_arg = pos_score;

}

std::queue<int> AlignWithBTCpuXL(char* a, int m, char* b, int nbb, int jj, int* sub_mat, int gap_open, int gap_extend, int query_pos, int ref_pos, bool reverse, bool first, int early_terminate, int* BT_states_local, int *qsize, int *rear, int *prev_lastCol, int *next_lastCol, int *prev_maxRow, int *next_maxRow) {
    int ext_n = nbb * BLOCK_WIDTH;
    int score;
    int match = 1, mismatch = -1;
    int max_score, max_i, max_j, pos_score; 
    int dir_matrix[MAX_TILE_SIZE+1][BLOCK_WIDTH+1];

    int reverse_int = reverse, first_tile_int = first;

  	for (int i = 0; i < m+1; i++) 
      		dir_matrix[i][0] = ZERO_OP;
  	for (int j = 0; j < BLOCK_WIDTH + 1; j++) 
      		dir_matrix[0][j] = ZERO_OP;

   xl_cpu( a, m, b, nbb, gap_open,  gap_extend,  match, mismatch,  prev_lastCol,  next_lastCol,  prev_maxRow,  next_maxRow,
                &score, jj, BT_states_local, qsize, rear, query_pos, ref_pos, reverse_int, first_tile_int, early_terminate, sub_mat,
		(int *)dir_matrix, &max_score, &max_i, &max_j, &pos_score);

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
