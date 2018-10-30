#include <stdio.h> // for debugging only
#include <assert.h>
#include "ksw2.h"

typedef struct { int32_t h, e; } eh_t;

int ksw_ggd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w,
		  	int8_t *pso, int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{
	eh_t *eh;
	int8_t *qp; // query profile
	int32_t i, j, k, gapoe = gapo + gape, score, n_col, *off = 0;
	uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

	assert(m_cigar_ && n_cigar_ && cigar_);
	// allocate memory
	if (w < 0) w = tlen > qlen? tlen : qlen;
	n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
	qp = (int8_t*)malloc(qlen * m);
	eh = (eh_t*)calloc(qlen + 1, 8);
	if (m_cigar_ && n_cigar_ && cigar_) {
		*n_cigar_ = 0;
		z = (uint8_t*)malloc((size_t)n_col * tlen);
		off = (int32_t*)calloc(tlen, 4);
	}

	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}

	// fill the first row
	eh[0].h = 0, eh[0].e = -gapoe - gapoe;
	for (j = 1; j <= qlen && j <= w; ++j)
		eh[j].h = -(gapoe + gape * (j - 1)), eh[j].e = -(gapoe + gapoe + gape * j);
	for (; j <= qlen; ++j) eh[j].h = eh[j].e = KSW_NEG_INF; // everything is -inf outside the band

	// DP loop
	for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
		int32_t f = KSW_NEG_INF, h1, st, en, os = pso? -pso[i] : 0, gapoe2 = gapoe + os;
		int8_t *q = &qp[target[i] * qlen];
		uint8_t *zi = &z[(long)i * n_col];
		st = i > w? i - w : 0;
		en = i + w + 1 < qlen? i + w + 1 : qlen;
		h1 = st > 0? KSW_NEG_INF : -(gapoe + gape * i);
		f  = st > 0? KSW_NEG_INF : -(gapoe + gapoe + os + gape * i);
		off[i] = st;
		for (j = st; j < en; ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1) + S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int32_t h = p->h, e = p->e;
			uint8_t d; // direction
			p->h = h1;
			h += q[j];
			d = h >= e? 0 : 1;
			h = h >= e? h : e;
			d = h >= f? d : 2;
			h = h >= f? h : f;
			h1 = h;
			h -= gapoe2;
			e -= gape;
			d |= e > h? 0x08 : 0;
			e  = e > h? e    : h;
			p->e = e;
			f -= gape;
			d |= f > h? 0x10 : 0; // if we want to halve the memory, use one bit only, instead of two
			f  = f > h? f    : h;
			zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
		}
		eh[en].h = h1, eh[en].e = KSW_NEG_INF;
	}

	// backtrack
	score = eh[qlen].h;
	free(qp); free(eh);
	if (m_cigar_ && n_cigar_ && cigar_) {
		ksw_backtrack(0, 0, 0, z, off, 0, n_col, tlen-1, qlen-1, m_cigar_, n_cigar_, cigar_);
		free(z); free(off);
	}
	return score;
}
