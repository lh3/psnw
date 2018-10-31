#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
#include "ksw2.h"

#define BUF_SIZE 0x10000

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
	int8_t a = 6, b = 18, gapo = 30, gape = 5;
	int c, w = -1;

	while ((c = getopt(argc, argv, "w:A:B:O:E:")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'A') a = atoi(optarg);
		else if (c == 'B') b = atoi(optarg);
		else if (c == 'O') gapo = atoi(optarg);
		else if (c == 'E') gape = atoi(optarg);
	}
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: psnw [options] <input>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -A INT        match score [%d]\n", a);
		fprintf(stderr, "  -B INT        mismatch penalty [%d]\n", b);
		fprintf(stderr, "  -O INT        gap open penalty [%d]\n", gapo);
		fprintf(stderr, "  -E INT        gap extension penalty [%d]\n", gape);
		fprintf(stderr, "Note: input format: name, target sequence, query, gap open for ins\n");
		fprintf(stderr, "  and gap extension for del. TAB delimited.\n");
		return 1;
	} else {
		int m_cigar = 0, n_cigar = 0;
		uint32_t *cigar = 0;
		int8_t mat[25];
		char *buf, *cigar_str, *out[3];
		FILE *fp;

		buf = (char*)malloc(BUF_SIZE * 5);
		cigar_str = buf + BUF_SIZE;
		out[0] = cigar_str + BUF_SIZE, out[1] = out[0] + BUF_SIZE, out[2] = out[1] + BUF_SIZE;
		ksw_gen_simple_mat(5, mat, a, -b);

		fp = strcmp(argv[optind], "-") == 0? stdin : fopen(argv[optind], "r");
		assert(fp);

		while (fgets(buf, BUF_SIZE - 1, fp) != NULL) {
			char *name, *p, *q;
			int8_t *pso = 0, *pse = 0;
			uint8_t *tseq = 0, *qseq = 0;
			int i, j, k, l, bw, score, tlen, qlen, nm;

			// parse
			for (k = 0, p = q = buf;; ++q) {
				if (*q == 0 || *q == '\n' || isspace(*q)) {
					int oq = *q;
					*q = 0;
					if (k == 0) {
						name = p;
					} else if (k == 1) {
						tseq = (uint8_t*)p;
						tlen = q - p;
						for (i = 0; i < tlen; ++i)
							tseq[i] = seq_nt4_table[tseq[i]];
					} else if (k == 2) {
						qseq = (uint8_t*)p;
						qlen = q - p;
						for (i = 0; i < qlen; ++i)
							qseq[i] = seq_nt4_table[qseq[i]];
					} else if (k == 3) {
						assert(tlen == q - p);
						pso = (int8_t*)p;
						for (i = 0; i < tlen; ++i) {
							c = p[i] - '0';
							assert(c >= 0 && c <= 127);
							pso[i] = c;
						}
					} else if (k == 4) {
						assert(tlen == q - p);
						pse = (int8_t*)p;
						for (i = 0; i < tlen; ++i) {
							c = p[i] - '0';
							assert(c >= 0 && c <= 127);
							pse[i] = c;
						}
					}
					++k, p = q + 1;
					if (oq == 0 || oq == '\n') break;
				}
			}

			// perform alignment
			bw = w > 0? w : tlen > qlen? tlen : qlen;
			score = ksw_ggd(0, qlen, qseq, tlen, tseq, 5, mat, gapo, gape, bw, pso, pse, &m_cigar, &n_cigar, &cigar);

			// output
			printf(">%s\t", name);
			for (k = 0; k < n_cigar; ++k)
				printf("%d%c", cigar[k]>>4, "MID"[cigar[k]&0xf]);
			printf("\tMD:Z:");
			for (k = i = j = l = nm = 0; k < n_cigar; ++k) {
				int t, op = cigar[k] & 0xf, len = cigar[k] >> 4;
				if (op == 0) {
					for (t = 0; t < len; ++t) {
						if (tseq[i + t] == qseq[j + t] && tseq[i + t] < 4) {
							++l;
						} else {
							printf("%d%c", l, "ACGTN"[tseq[i + t]]);
							l = 0, ++nm;
						}
					}
					i += len, j += len;
				} else if (op == 1) { // insertion
					j += len, nm += len;
				} else { // op == 2
					printf("%d^", l);
					for (t = 0; t < len; ++t)
						putchar("ACGTN"[tseq[i + t]]);
					l = 0;
					i += len, nm += len;
				}
			}
			printf("%d\tNM:i:%d\tAS:i:%d\n", l, nm, score);

			// print base alignment
			for (k = i = j = l = 0; k < n_cigar; ++k) {
				int t, op = cigar[k] & 0xf, len = cigar[k] >> 4;
				if (op == 0) {
					for (t = 0; t < len; ++t) {
						out[0][l] = "ACGTN"[tseq[i + t]];
						out[2][l] = "ACGTN"[qseq[j + t]];
						out[1][l++] = tseq[i + t] == qseq[j + t] && tseq[i + t] < 4? '|' : ' ';
					}
					i += len, j += len;
				} else if (op == 1) { // insertion
					for (t = 0; t < len; ++t)
						out[0][l] = '-', out[2][l] = "ACGTN"[qseq[j + t]], out[1][l++] = ' ';
					j += len;
				} else { // op == 2
					for (t = 0; t < len; ++t)
						out[0][l] = "ACGTN"[tseq[i + t]], out[2][l] = '-', out[1][l++] = ' ';
					i += len;
				}
			}
			printf("%s\n%s\n%s\n", out[0], out[1], out[2]);
		}

		free(cigar); free(buf);
		fclose(fp);
	}
	return 0;
}
