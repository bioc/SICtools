#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "bam.h"
#include "bam2bcf.h"
#include "kaln.h"
#include "kprobaln.h"
#include "khash.h"
KHASH_SET_INIT_STR(rg)

#include "ksort.h"
KSORT_INIT_GENERIC(uint32_t)

#define MINUS_CONST 0x10000000
#define INDEL_WINDOW_SIZE 50

void *bcf_call_add_rg(void *_hash, const char *hdtext, const char *list)
{
	const char *s, *p, *q, *r, *t;
	khash_t(rg) *hash;
	if (list == 0 || hdtext == 0) return _hash;
	if (_hash == 0) _hash = kh_init(rg);
	hash = (khash_t(rg)*)_hash;
	if ((s = strstr(hdtext, "@RG\t")) == 0) return hash;
	do {
		t = strstr(s + 4, "@RG\t"); // the next @RG
		if ((p = strstr(s, "\tID:")) != 0) p += 4;
		if ((q = strstr(s, "\tPL:")) != 0) q += 4;
		if (p && q && (t == 0 || (p < t && q < t))) { // ID and PL are both present
			int lp, lq;
			char *x;
			for (r = p; *r && *r != '\t' && *r != '\n'; ++r); lp = r - p;
			for (r = q; *r && *r != '\t' && *r != '\n'; ++r); lq = r - q;
			x = calloc((lp > lq? lp : lq) + 1, 1);
			for (r = q; *r && *r != '\t' && *r != '\n'; ++r) x[r-q] = *r;
			if (strstr(list, x)) { // insert ID to the hash table
				khint_t k;
				int ret;
				for (r = p; *r && *r != '\t' && *r != '\n'; ++r) x[r-p] = *r;
				x[r-p] = 0;
				k = kh_get(rg, hash, x);
				if (k == kh_end(hash)) k = kh_put(rg, hash, x, &ret);
				else free(x);
			} else free(x);
		}
		s = t;
	} while (s);
	return hash;
}

void bcf_call_del_rghash(void *_hash)
{
	khint_t k;
	khash_t(rg) *hash = (khash_t(rg)*)_hash;
	if (hash == 0) return;
	for (k = kh_begin(hash); k < kh_end(hash); ++k)
		if (kh_exist(hash, k))
			free((char*)kh_key(hash, k));
	kh_destroy(rg, hash);
}

/*static int tpos2qpos(const bam1_core_t *c, const uint32_t *cigar, int32_t tpos, int is_left, int32_t *_tpos)
{
	int k, x = c->pos, y = 0, last_y = 0;
 *_tpos = c->pos;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		int l = cigar[k] >> BAM_CIGAR_SHIFT;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			if (c->pos > tpos) return y;
			if (x + l > tpos) {
 *_tpos = tpos;
				return y + (tpos - x);
			}
			x += l; y += l;
			last_y = y;
		} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
			if (x + l > tpos) {
 *_tpos = is_left? x : x + l;
				return y;
			}
			x += l;
		}
	}
 *_tpos = x;
	return last_y;
}
 */
// FIXME: check if the inserted sequence is consistent with the homopolymer run
// l is the relative gap length and l_run is the length of the homopolymer on the reference
static inline int est_seqQ(const bcf_callaux_t *bca, int l, int l_run)
{
	int q, qh;
	q = bca->openQ + bca->extQ * (abs(l) - 1);
	qh = l_run >= 3? (int)(bca->tandemQ * (double)abs(l) / l_run + .499) : 1000;
	return q < qh? q : qh;
}

static inline int est_indelreg(int pos, const char *ref, int l, char *ins4)
{
	int i, j, max = 0, max_i = pos, score = 0;
	l = abs(l);
	for (i = pos + 1, j = 0; ref[i]; ++i, ++j) {
		if (ins4) score += (toupper(ref[i]) != "ACGTN"[(int)ins4[j%l]])? -10 : 1;
		else score += (toupper(ref[i]) != toupper(ref[pos+1+j%l]))? -10 : 1;
		if (score < 0) break;
		if (max < score) max = score, max_i = i;
	}
	return max_i - pos;
}

/*
 *  @n:  number of samples
 */
int bcf_call_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref,
		const void *rghash)
{
	int i, s, j, k, t, n_types, *types, max_rd_len, left, right, max_ins, *score1, *score2, max_ref2;
	int N, l_run, ref_type;
	char *inscns = 0, *ref2, *query, **ref_sample;
	khash_t(rg) *hash = (khash_t(rg)*)rghash;
	if (ref == 0 || bca == 0) return -1;
	// mark filtered reads
	if (rghash) {
		N = 0;
		for (s = N = 0; s < n; ++s) {
			for (i = 0; i < n_plp[s]; ++i) {
				bam_pileup1_t *p = plp[s] + i;
				const uint8_t *rg = bam_aux_get(p->b, "RG");
				p->aux = 1; // filtered by default
				if (rg) {
					khint_t k = kh_get(rg, hash, (const char*)(rg + 1));
					if (k != kh_end(hash)) p->aux = 0, ++N; // not filtered
				}
			}
		}
		if (N == 0) return -1; // no reads left
	}
	// determine if there is a gap
	for (s = N = 0; s < n; ++s) {
		for (i = 0; i < n_plp[s]; ++i)
			if (plp[s][i].indel != 0) break;
		if (i < n_plp[s]) break;
	}
	if (s == n) return -1; // there is no indel at this position.
	for (s = N = 0; s < n; ++s) N += n_plp[s]; // N is the total number of reads
	{ // find out how many types of indels are present
		bca->max_support = bca->max_frac = 0;
		int m, n_alt = 0, n_tot = 0, indel_support_ok = 0;
		uint32_t *aux;
		aux = calloc(N + 1, 4);
		m = max_rd_len = 0;
		aux[m++] = MINUS_CONST; // zero indel is always a type
		for (s = 0; s < n; ++s) {
			int na = 0, nt = 0;
			for (i = 0; i < n_plp[s]; ++i) {
				const bam_pileup1_t *p = plp[s] + i;
				if (rghash == 0 || p->aux == 0) {
					++nt;
					if (p->indel != 0) {
						++na;
						aux[m++] = MINUS_CONST + p->indel;
					}
				}
				j = bam_cigar2qlen(&p->b->core, bam1_cigar(p->b));
				if (j > max_rd_len) max_rd_len = j;
			}
			float frac = (float)na/nt;
			if ( !indel_support_ok && na >= bca->min_support && frac >= bca->min_frac )
				indel_support_ok = 1;
			if ( na > bca->max_support && frac > 0 ) bca->max_support = na, bca->max_frac = frac;
			n_alt += na;
			n_tot += nt;
		}
		// To prevent long stretches of N's to be mistaken for indels (sometimes thousands of bases),
		//  check the number of N's in the sequence and skip places where half or more reference bases are Ns.
		int nN=0; for (i=pos; i-pos<max_rd_len && ref[i]; i++) if ( ref[i]=='N' ) nN++;
		if ( nN*2>i ) { free(aux); return -1; }

		ks_introsort(uint32_t, m, aux);
		// squeeze out identical types
		for (i = 1, n_types = 1; i < m; ++i)
			if (aux[i] != aux[i-1]) ++n_types;
		// Taking totals makes it hard to call rare indels
		if ( !bca->per_sample_flt )
			indel_support_ok = ( (float)n_alt / n_tot < bca->min_frac || n_alt < bca->min_support ) ? 0 : 1;
		if ( n_types == 1 || !indel_support_ok ) { // then skip
			free(aux); return -1;
		}
		if (n_types >= 64) {
			free(aux);
			if (bam_verbose >= 2) 
				fprintf(stderr, "[%s] excessive INDEL alleles at position %d. Skip the position.\n", __func__, pos + 1);
			return -1;
		}
		types = (int*)calloc(n_types, sizeof(int));
		t = 0;
		types[t++] = aux[0] - MINUS_CONST; 
		for (i = 1; i < m; ++i)
			if (aux[i] != aux[i-1])
				types[t++] = aux[i] - MINUS_CONST;
		free(aux);
		for (t = 0; t < n_types; ++t)
			if (types[t] == 0) break;
		ref_type = t; // the index of the reference type (0)
	}
	{ // calculate left and right boundary
		left = pos > INDEL_WINDOW_SIZE? pos - INDEL_WINDOW_SIZE : 0;
		right = pos + INDEL_WINDOW_SIZE;
		if (types[0] < 0) right -= types[0];
		// in case the alignments stand out the reference
		for (i = pos; i < right; ++i)
			if (ref[i] == 0) break;
		right = i;
	}
	/* The following block fixes a long-existing flaw in the INDEL
	 * calling model: the interference of nearby SNPs. However, it also
	 * reduces the power because sometimes, substitutions caused by
	 * indels are not distinguishable from true mutations. Multiple
	 * sequence realignment helps to increase the power.
	 *
	 * Masks mismatches present in at least 70% of the reads with 'N'.
	 */
	{ // construct per-sample consensus
		int L = right - left + 1, max_i, max2_i;
		uint32_t *cns, max, max2;
		char *ref0, *r;
		ref_sample = calloc(n, sizeof(void*));
		cns = calloc(L, 4);
		ref0 = calloc(L, 1);
		for (i = 0; i < right - left; ++i)
			ref0[i] = bam_nt16_table[(int)ref[i+left]];
		for (s = 0; s < n; ++s) {
			r = ref_sample[s] = calloc(L, 1);
			memset(cns, 0, sizeof(int) * L);
			// collect ref and non-ref counts
			for (i = 0; i < n_plp[s]; ++i) {
				bam_pileup1_t *p = plp[s] + i;
				bam1_t *b = p->b;
				uint32_t *cigar = bam1_cigar(b);
				uint8_t *seq = bam1_seq(b);
				int x = b->core.pos, y = 0;
				for (k = 0; k < b->core.n_cigar; ++k) {
					int op = cigar[k]&0xf;
					int j, l = cigar[k]>>4;
					if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
						for (j = 0; j < l; ++j)
							if (x + j >= left && x + j < right)
								cns[x+j-left] += (bam1_seqi(seq, y+j) == ref0[x+j-left])? 1 : 0x10000;
						x += l; y += l;
					} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) x += l;
					else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
				}
			}
			// determine the consensus
			for (i = 0; i < right - left; ++i) r[i] = ref0[i];
			max = max2 = 0; max_i = max2_i = -1;
			for (i = 0; i < right - left; ++i) {
				if (cns[i]>>16 >= max>>16) max2 = max, max2_i = max_i, max = cns[i], max_i = i;
				else if (cns[i]>>16 >= max2>>16) max2 = cns[i], max2_i = i;
			}
			if ((double)(max&0xffff) / ((max&0xffff) + (max>>16)) >= 0.7) max_i = -1;
			if ((double)(max2&0xffff) / ((max2&0xffff) + (max2>>16)) >= 0.7) max2_i = -1;
			if (max_i >= 0) r[max_i] = 15;
			if (max2_i >= 0) r[max2_i] = 15;
			//for (i = 0; i < right - left; ++i) fputc("=ACMGRSVTWYHKDBN"[(int)r[i]], stderr); fputc('\n', stderr);
		}
		free(ref0); free(cns);
	}
	{ // the length of the homopolymer run around the current position
		int c = bam_nt16_table[(int)ref[pos + 1]];
		if (c == 15) l_run = 1;
		else {
			for (i = pos + 2; ref[i]; ++i)
				if (bam_nt16_table[(int)ref[i]] != c) break;
			l_run = i;
			for (i = pos; i >= 0; --i)
				if (bam_nt16_table[(int)ref[i]] != c) break;
			l_run -= i + 1;
		}
	}
	// construct the consensus sequence
	max_ins = types[n_types - 1];   // max_ins is at least 0
	if (max_ins > 0) {
		int *inscns_aux = calloc(5 * n_types * max_ins, sizeof(int));
		// count the number of occurrences of each base at each position for each type of insertion
		for (t = 0; t < n_types; ++t) {
			if (types[t] > 0) {
				for (s = 0; s < n; ++s) {
					for (i = 0; i < n_plp[s]; ++i) {
						bam_pileup1_t *p = plp[s] + i;
						if (p->indel == types[t]) {
							uint8_t *seq = bam1_seq(p->b);
							for (k = 1; k <= p->indel; ++k) {
								int c = bam_nt16_nt4_table[bam1_seqi(seq, p->qpos + k)];
								assert(c<5);
								++inscns_aux[(t*max_ins+(k-1))*5 + c];
							}
						}
					}
				}
			}
		}
		// use the majority rule to construct the consensus
		inscns = calloc(n_types * max_ins, 1);
		for (t = 0; t < n_types; ++t) {
			for (j = 0; j < types[t]; ++j) {
				int max = 0, max_k = -1, *ia = &inscns_aux[(t*max_ins+j)*5];
				for (k = 0; k < 5; ++k)
					if (ia[k] > max)
						max = ia[k], max_k = k;
				inscns[t*max_ins + j] = max? max_k : 4;
				if ( max_k==4 ) { types[t] = 0; break; } // discard insertions which contain N's
			}
		}
		free(inscns_aux);
	}
	// compute the likelihood given each type of indel for each read
	max_ref2 = right - left + 2 + 2 * (max_ins > -types[0]? max_ins : -types[0]);
	ref2  = calloc(max_ref2, 1);
	query = calloc(right - left + max_rd_len + max_ins + 2, 1);
	score1 = calloc(N * n_types, sizeof(int));
	score2 = calloc(N * n_types, sizeof(int));
	bca->indelreg = 0;

	int ir;
	for (t = 0; t < n_types; ++t){
		if (types[t] == 0) ir = 0;
		else if (types[t] > 0) ir = est_indelreg(pos, ref, types[t], &inscns[t*max_ins]);
		else ir = est_indelreg(pos, ref, -types[t], 0);
		if (ir > bca->indelreg) bca->indelreg = ir;
		//		fprintf(stderr, "%d, %d, %d, %d\n", pos+1, types[t], ir,l_run);
		printf("%d,",bca->indelreg);
	}
	putchar('\t');

	return 1;

}
