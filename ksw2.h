#ifndef KSW2_H_
#define KSW2_H_

#include <stdint.h>
#include <stdlib.h>

#define KSW_NEG_INF -0x40000000

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Global alignment
 *
 * (first 10 parameters identical to ksw_extz_sse())
 * @param m_cigar   (modified) max CIGAR length; feed 0 if cigar==0
 * @param n_cigar   (out) number of CIGAR elements
 * @param cigar     (out) BAM-encoded CIGAR; caller need to deallocate with kfree(km, )
 *
 * @return          score of the alignment
 */
int ksw_ggd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t gapo, int8_t gape, int w,
		  	int8_t *pso, int *m_cigar_, int *n_cigar_, uint32_t **cigar_);

#ifdef __cplusplus
}
#endif

/************************************
 *** Private macros and functions ***
 ************************************/

static inline uint32_t *ksw_push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = (uint32_t*)realloc(cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

// In the backtrack matrix, value p[] has the following structure:
//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
static inline void ksw_backtrack(int is_rot, int is_rev, int min_intron_len, const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0,
								 int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{ // p[] - lower 3 bits: which type gets the max; bit
	int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
	uint32_t *cigar = *cigar_, tmp;
	while (i >= 0 && j >= 0) { // at the beginning of the loop, _state_ tells us which state to check
		int force_state = -1;
		if (is_rot) {
			r = i + j;
			if (i < off[r]) force_state = 2;
			if (off_end && i > off_end[r]) force_state = 1;
			tmp = force_state < 0? p[(size_t)r * n_col + i - off[r]] : 0;
		} else {
			if (j < off[i]) force_state = 2;
			if (off_end && j > off_end[i]) force_state = 1;
			tmp = force_state < 0? p[(size_t)i * n_col + j - off[i]] : 0;
		}
		if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
		else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
		if (state == 0) state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
		if (force_state >= 0) state = force_state;
		if (state == 0) cigar = ksw_push_cigar(&n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
		else if (state == 1 || (state == 3 && min_intron_len <= 0)) cigar = ksw_push_cigar(&n_cigar, &m_cigar, cigar, 2, 1), --i; // deletion
		else if (state == 3 && min_intron_len > 0) cigar = ksw_push_cigar(&n_cigar, &m_cigar, cigar, 3, 1), --i; // intron
		else cigar = ksw_push_cigar(&n_cigar, &m_cigar, cigar, 1, 1), --j; // insertion
	}
	if (i >= 0) cigar = ksw_push_cigar(&n_cigar, &m_cigar, cigar, min_intron_len > 0 && i >= min_intron_len? 3 : 2, i + 1); // first deletion
	if (j >= 0) cigar = ksw_push_cigar(&n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
	if (!is_rev)
		for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
			tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
	*m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

#endif
