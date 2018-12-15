#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int st, en, x;
} intv1_t;

typedef struct {
	int n, m;
	intv1_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)
KHASH_SET_INIT_INT64(64)

typedef kh_reg_t reghash_t;

reghash_t *reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t str = {0,0,0};

	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		int i, st = -1, en = -1, x = -1;
		char *p, *q, *ctg = 0;
		reglist_t *r;
		khint_t k;
		for (p = q = str.s, i = 0;; ++q) {
			if (*q == 0 || *q == '\n') {
				int c = *q;
				*q = 0;
				if (i == 0) ctg = p;
				else if (i == 1) st = atoi(p);
				else if (i == 2) en = atoi(p);
				else if (i == 3) x = atoi(p);
				++i, p = q + 1;
				if (c == 0) break;
			}
		}
		if (i < 3 || st < 0 || en <= 0 || x <= 0 || x > 127 - 33) continue;
		k = kh_get(reg, h, ctg);
		if (k == kh_end(h)) {
			int ret;
			char *s;
			s = strdup(ctg);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		r = &kh_val(h, k);
		if (r->n == r->m) {
			r->m = r->m? r->m<<1 : 4;
			r->a = realloc(r->a, r->m * sizeof(intv1_t));
		}
		r->a[r->n].st = st;
		r->a[r->n].en = en;
		r->a[r->n++].x = x;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str.s);
	return h;
}

void reg_destroy(reghash_t *h)
{
	khint_t k;
	if (h == 0) return;
	for (k = 0; k < kh_end(h); ++k) {
		if (!kh_exist(h, k)) continue;
		free(kh_val(h, k).a);
		free((char*)kh_key(h, k));
	}
	kh_destroy(reg, h);
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	reghash_t *h;
	gzFile fp;
	kseq_t *ks;
	int c;

	while ((c = ketopt(&o, argc, argv, 1, "", 0)) >= 0) {
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: gen-fq <in.fa> <in.bed>\n");
		return 1;
	}

	fp = strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(0, "r");
	ks = kseq_init(fp);
	h = reg_read(argv[o.ind+1]);
	while (kseq_read(ks) >= 0) {
		khint_t k;
		if (ks->qual.m < ks->seq.m) {
			ks->qual.m = ks->seq.m;
			ks->qual.s = (char*)realloc(ks->qual.s, ks->qual.m);
		}
		ks->qual.l = ks->seq.l;
		memset(ks->qual.s, 33, ks->qual.l);
		ks->qual.s[ks->qual.l] = 0;
		k = kh_get(reg, h, ks->name.s);
		if (k != kh_end(h)) {
			reglist_t *p = &kh_val(h, k);
			int i, j;
			for (i = 0; i < p->n; ++i)
				for (j = p->a[i].st; j < p->a[i].en; ++j)
					ks->qual.s[j] += p->a[i].x;
		}
		printf("@%s\n", ks->name.s);
		puts(ks->seq.s);
		puts("+");
		puts(ks->qual.s);
	}
	reg_destroy(h);
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}
