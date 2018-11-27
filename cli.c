#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <ctype.h>
#include "ksw2.h"

#define BUF_SIZE 0x10000

#define MALLOC(type, len) ((type*)malloc((len) * sizeof(type)))
#define CALLOC(type, len) ((type*)calloc((len), sizeof(type)))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#define EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		REALLOC((a), (m)); \
	} while (0)

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

typedef struct {
	int offset;
	uint32_t *cigar;
} shifted_cigar_t;

static shifted_cigar_t *shift_cigar(int n_cigar, const uint32_t *cigar, int tlen, const uint8_t *tseq, int qlen, const uint8_t *qseq, int *n_)
{
	int n = 1, m = 8;
	int gk = -1, glen = -1, gx = -1, gy = -1;
	int k, x, y;
	shifted_cigar_t *cs;

	cs = CALLOC(shifted_cigar_t, m);
	cs[0].cigar = MALLOC(uint32_t, n_cigar);
	memcpy(cs[0].cigar, cigar, 4 * n_cigar);

	for (k = x = y = 0; k < n_cigar; ++k) {
		int op = cigar[k] & 0xf, len = cigar[k] >> 4;
		if (op == 0) x += len, y += len;
		else {
			if (len > glen) gk = k, glen = len, gx = x, gy = y; // find the longest gap
			if (op == 1) y += len;
			else if (op == 2) x += len;
		}
	}

	if (gk > 0 && gk < n_cigar - 1) {
		int o;
		for (o = 1, x = gx, y = gy; x < tlen - 1 && y < qlen - 1 && (cigar[gk+1]>>4) - o > 1; ++o, ++x, ++y) { // move the gap to the right
			if (tseq[x] == qseq[y]) {
				if (n == m) EXPAND(cs, m);
				cs[n].cigar = MALLOC(uint32_t, n_cigar);
				memcpy(cs[n].cigar, cigar, 4 * n_cigar);
				cs[n].cigar[gk-1] += o<<4;
				cs[n].cigar[gk+1] -= o<<4;
				cs[n++].offset = o;
			} else break;
		}
		if ((cigar[gk]&0xf) == 1) gy += cigar[gk]>>4;
		else gx += cigar[gk]>>4;
		for (o = -1, x = gx, y = gy; x > 0 && y > 0 && (cigar[gk-1]>>4) - o > 1; --o, --x, --y) { // move the gap to the left
			if (tseq[x - 1] == qseq[y - 1]) {
				if (n == m) EXPAND(cs, m);
				cs[n].cigar = MALLOC(uint32_t, n_cigar);
				memcpy(cs[n].cigar, cigar, 4 * n_cigar);
				cs[n].cigar[gk-1] -= (-o)<<4;
				cs[n].cigar[gk+1] += (-o)<<4;
				cs[n++].offset = o;
			} else break;
		}
	}
	*n_ = n;
	return cs;
}

int main(int argc, char *argv[])
{
	int8_t a = 6, b = 18, gapo = 30, gape = 5;
	int c, w = -1, print_aln = 0;

	while ((c = getopt(argc, argv, "w:A:B:O:E:a")) >= 0) {
		if (c == 'w') w = atoi(optarg);
		else if (c == 'A') a = atoi(optarg);
		else if (c == 'B') b = atoi(optarg);
		else if (c == 'O') gapo = atoi(optarg);
		else if (c == 'E') gape = atoi(optarg);
		else if (c == 'a') print_aln = 1;
	}
	if (argc - optind < 1) {
		fprintf(stderr, "Usage: psnw [options] <input>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -A INT        match score [%d]\n", a);
		fprintf(stderr, "  -B INT        mismatch penalty [%d]\n", b);
		fprintf(stderr, "  -O INT        gap open penalty [%d]\n", gapo);
		fprintf(stderr, "  -E INT        gap extension penalty [%d]\n", gape);
		fprintf(stderr, "  -w INT        band width [inf]\n");
		fprintf(stderr, "  -a            print detailed base alignment\n");
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

		printf("CC\tCC     comment\n");
		printf("CC\tSQ     name    #records  AS\n");
		printf("CC\tCG     cigar   MD        NM\n");
		printf("CC\tO[123] alignment-strings\n");
		printf("CC\t//     end-of-report\n");
		printf("CC\n");

		while (fgets(buf, BUF_SIZE - 1, fp) != NULL) {
			char *name, *p, *q;
			int8_t *pso = 0, *pse = 0;
			uint8_t *tseq = 0, *qseq = 0;
			int i, j, k, l, bw, score, tlen, qlen, nm, n_cs;
			shifted_cigar_t *cs;

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
			cs = shift_cigar(n_cigar, cigar, tlen, tseq, qlen, qseq, &n_cs);
			printf("SQ\t%s\t%d\tAS:i:%d\n", name, n_cs, score);
			for (c = 0; c < n_cs; ++c) {
				uint32_t *cigar = cs[c].cigar; // WARNING: this shadows "cigar" in the upper level
				printf("CG\t");
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
				printf("%d\tNM:i:%d\n", l, nm);

				// print base alignment
				if (!print_aln) goto end_print;
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
				printf("O1\t%s\nO2\t%s\nO3\t%s\n", out[0], out[1], out[2]);
end_print:
				free(cigar);
			}
			printf("//\n");
			free(cs);
		}

		free(cigar); free(buf);
		fclose(fp);
	}
	return 0;
}
