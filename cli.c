#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include "ksw2.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}

int main(int argc, char *argv[])
{
	int8_t a = 6, b = 18, q = 30, e = 5;
	int c, w = -1;

	while ((c = getopt(argc, argv, "w:A:B:O:E:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'A') a = atoi(optarg);
		else if (c == 'B') b = atoi(optarg);
		else if (c == 'O') q = atoi(optarg);
		else if (c == 'E') e = atoi(optarg);
	}
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: ksw2-test [options] <DNA-target> <DNA-query>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -A INT        match score [%d]\n", a);
		fprintf(stderr, "  -B INT        mismatch penalty [%d]\n", b);
		fprintf(stderr, "  -O INT        gap open penalty [%d]\n", q);
		fprintf(stderr, "  -E INT        gap extension penalty [%d]\n", e);
		return 1;
	} else {
		int tlen, qlen, i, j, k, l, score, m_cigar = 0, n_cigar = 0, blen;
		uint32_t *cigar = 0;
		uint8_t *tseq, *qseq;
		int8_t mat[25], *pso = 0;
		char *out[3];

		// prepare for the alignment
		ksw_gen_simple_mat(5, mat, a, -b);
		tlen = strlen(argv[optind]);
		qlen = strlen(argv[optind+1]);
		if (optind + 2 < argc) {
			j = strlen(argv[optind+2]);
			assert(j == tlen);
		}
		tseq = (uint8_t*)calloc(tlen * 2 + qlen, 1);
		qseq = tseq + tlen;
		pso = (int8_t*)qseq + qlen;
		for (j = 0; j < tlen; ++j) tseq[j] = seq_nt4_table[(uint8_t)argv[optind][j]];
		for (j = 0; j < qlen; ++j) qseq[j] = seq_nt4_table[(uint8_t)argv[optind + 1][j]];
		if (optind + 2 < argc) {
			for (j = 0; j < tlen; ++j) {
				int c = (int)argv[optind + 2][j] - '0';
				assert(c >= 0 && c <= 127);
				pso[j] = c;
			}
		}
		if (w < 0) w = tlen > qlen? tlen : qlen;

		// perform alignment
		score = ksw_ggd(0, qlen, qseq, tlen, tseq, 5, mat, q, e, w, pso, &m_cigar, &n_cigar, &cigar);

		// print alignment
		for (k = blen = 0; k < n_cigar; ++k) blen += cigar[k] >> 4;
		out[0] = (char*)calloc((blen + 1) * 3, 1);
		out[1] = out[0] + blen + 1;
		out[2] = out[1] + blen + 1;
		for (k = i = j = l = 0; k < n_cigar; ++k) {
			int t, op = cigar[k] & 0xf, len = cigar[k] >> 4;
			if (op == 0) {
				for (t = 0; t < len; ++t) {
					out[0][l] = argv[optind][i + t];
					out[2][l] = argv[optind+1][j + t];
					out[1][l++] = tseq[i + t] == qseq[j + t] && tseq[i + t] < 4? '|' : ' ';
				}
				i += len, j += len;
			} else if (op == 1) {
				for (t = 0; t < len; ++t)
					out[0][l] = '-', out[2][l] = argv[optind+1][j + t], out[1][l++] = ' ';
				j += len;
			} else {
				for (t = 0; t < len; ++t)
					out[0][l] = argv[optind][i + t], out[2][l] = '-', out[1][l++] = ' ';
				i += len;
			}
		}
		printf("Score: %d\n%s\n%s\n%s\n", score, out[0], out[1], out[2]);
		free(tseq);
		free(out[0]);
	}
	return 0;
}
