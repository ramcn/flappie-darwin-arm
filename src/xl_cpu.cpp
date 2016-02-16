#include <queue>
#include <stdio.h>
#include <align.h>

enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
enum states {Z, D, I, M};

#define A_NT 0
#define C_NT 1
#define G_NT 2
#define T_NT 3
#define N_NT 4


void xl_cpu_score( const char *  a, const short i_start, const char *  b, const short j_start,
				const short open_extend_gap, const short extend_gap, const short match, const short mismatch,
				short *maxScore, short *globalH, short tile_num){
	// matrix H
	short h[TILE_SIZE][TILE_SIZE];	
	// auxiliar buffers
	short row[TILE_SIZE]={0};
	short maxCol[TILE_SIZE]={0};
	//max score of current alignment
	short score= 0; 
	// max score position
	short score_pos_i, score_pos_j; 

	// private buffer for reference sequence 
	char private_b[TILE_SIZE];

	// poshorter to query sequence
	 const char * ptr_a = a + i_start;
	// poshorter to reference sequence
	 const char * ptr_b = b + j_start;

	// copy sequence b to private memory
	for(short i = 0; i < TILE_SIZE; i++)
		private_b[i] = ptr_b[i];

	for(short i = 0; i < TILE_SIZE; i++){ // this loop could be from 0 to m-1, m being a kernel argument

		// copy resiude to local memory
		char a_i = ptr_a[i];
		// tmp var's for current row processing
		short maxRow_i = 0;
		short score_i = 0;
		short previous = 0;
		short score_i_pos_i = 0, score_i_pos_j = 0;

		for (short j=0; j < TILE_SIZE; j++){ 
			//calcuate the diagonal value
			short current = row[j] + (a_i==private_b[j] ? match : -mismatch);
			// calculate current max value
			current = max(current, maxRow_i);
			current = max(current, maxCol[j]);
			current = max(current, (short) 0);
			// update max score and its position
			score_i_pos_i = (current >= score_i ? i : score_i_pos_i);
			score_i_pos_j = (current >= score_i ? j : score_i_pos_j);
			score_i = max(score_i, current);
			// update maxRow and maxCol
			short aux1 = maxRow_i - extend_gap;
			short aux2 = maxCol[j] - extend_gap;
			short aux3 = current - open_extend_gap;
			maxRow_i = max(aux1, aux3);
			maxCol[j] =  max(aux2, aux3);	
			// update row
			row[j] = previous;
			previous = current;
			// update h
			h[i][j] = current ;
		}
		// update max score and its position
                score_pos_i = (score_i >= score ? score_i_pos_i : score_pos_i);
                score_pos_j = (score_i >= score ? score_i_pos_j : score_pos_j);
                score = max(score, score_i);
	}

        // save H in global memory
        for(short i = 0; i < TILE_SIZE; i++)
                for (short j=0; j < TILE_SIZE; j++)
                        globalH[i*TILE_SIZE+(j+(tile_num*TILE_SIZE))] = h[i][j];

        maxScore[0] = score;
        maxScore[1] = score_pos_i;
        maxScore[2] = score_pos_j;
}


void xl_cpu_direction (short m, const short n, const short gap_open, const short gap_extend, 
                                        short *  dir_matrix, short * h_matrix,
                                        short *h_matrix_rd, short *m_matrix_rd, short *i_matrix_rd, short *d_matrix_rd,
                                        short *h_matrix_wr, short *m_matrix_wr, short *i_matrix_wr, short *d_matrix_wr)
{
  for (short i = 1; i < n + 1; i++) {
      for (short k = 1; k < MAX_TILE_SIZE + 1; k++) {
          m_matrix_rd[k] = m_matrix_wr[k];
          h_matrix_rd[k] = h_matrix_wr[k];
          i_matrix_rd[k] = i_matrix_wr[k];
          d_matrix_rd[k] = d_matrix_wr[k];
      }
      for (short j = 1; j < m + 1; j++) {
          if (m_matrix_rd[j-1] > i_matrix_rd[j-1] && m_matrix_rd[j-1] > d_matrix_rd[j-1]) {
              m_matrix_wr[j] = m_matrix_rd[j-1] + h_matrix[i*TILE_SIZE+j];
          } else if (i_matrix_rd[j-1] > d_matrix_rd[j-1]) {
              m_matrix_wr[j] = i_matrix_rd[j-1] + h_matrix[i*TILE_SIZE+j];
          } else {
              m_matrix_wr[j] = d_matrix_rd[j-1] + h_matrix[i*TILE_SIZE+j];
          }
          if (m_matrix_wr[j] < 0) {
              m_matrix_wr[j] = 0;
          }
          short ins_open   = m_matrix_rd[j] + gap_open;
          short ins_extend = i_matrix_rd[j] + gap_extend;
          short del_open   = m_matrix_wr[j-1] + gap_open;
          short del_extend = d_matrix_wr[j-1] + gap_extend;
          i_matrix_wr[j] = (ins_open > ins_extend) ? ins_open : ins_extend;
          d_matrix_wr[j] = (del_open > del_extend) ? del_open : del_extend;
          short max1 = m_matrix_wr[j] > i_matrix_wr[j] ? m_matrix_wr[j] : i_matrix_wr[j];
          short max2 = d_matrix_wr[j] > 0 ? d_matrix_wr[j] : 0;
          h_matrix_wr[j] = max1 > max2 ? max1 : max2;
          (dir_matrix)[i*(MAX_TILE_SIZE+1)+j] = ((m_matrix_wr[j] >= i_matrix_wr[j]) ? ((m_matrix_wr[j] >= d_matrix_wr[j]) ? MATCH_OP : DELETE_OP) : ((i_matrix_wr[j] >= d_matrix_wr[j]) ? INSERT_OP : DELETE_OP));

          if ((m_matrix_wr[j] <= 0) && (i_matrix_wr[j] <= 0) && (d_matrix_wr[j] <= 0)) {
              (dir_matrix)[i*(MAX_TILE_SIZE+1)+j] = ZERO_OP;
          }

          (dir_matrix)[i*(MAX_TILE_SIZE+1)+j] += (ins_open >= ins_extend) ? (2 << INSERT_OP) : 0;
          (dir_matrix)[i*(MAX_TILE_SIZE+1)+j] += (del_open >= del_extend) ? (2 << DELETE_OP) : 0;
        }
  }
}

