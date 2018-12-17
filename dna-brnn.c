#include <zlib.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ketopt.h"
#include "dna-io.h"
#include "kann.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

kann_t *dbr_model_gen(int n_lbl, int n_layer, int n_neuron, float h_dropout, int is_tied)
{
	kad_node_t *s[2], *t, *w, *b, *y, *par[256]; // for very unreasonably deep models, this may overflow
	int i, k, offset;
	memset(par, 0, sizeof(kad_node_p) * 256);
	for (k = 0; k < 2; ++k) {
		s[k] = kad_feed(2, 1, 4), s[k]->ext_flag = KANN_F_IN, s[k]->ext_label = k + 1;
		offset = 0;
		for (i = 0; i < n_layer; ++i) {
			if (is_tied) {
				kad_node_t *h0;
				h0 = kann_new_leaf2(&offset, par, KAD_CONST, 0.0f, 2, 1, n_neuron);
				s[k] = kann_layer_gru2(&offset, par, s[k], h0, KANN_RNN_NORM);
			} else s[k] = kann_layer_gru(s[k], n_neuron, KANN_RNN_NORM);
			if (h_dropout > 0.0f) s[k] = kann_layer_dropout(s[k], h_dropout);
		}
		s[k] = kad_stack(1, &s[k]); // first and second pivot
	}
	s[1] = kad_reverse(s[1], 0);
	t = kad_concat(2, 2, s[0], s[1]), w = kann_new_weight(n_lbl, n_neuron * 2);
//	t = kad_avg(2, s), w= kann_new_weight(2, n_neuron);
	b = kann_new_bias(n_lbl);
	t = kad_softmax(kad_add(kad_cmul(t, w), b)), t->ext_flag = KANN_F_OUT;
	y = kad_feed(2, 1, n_lbl), y->ext_flag = KANN_F_TRUTH;
	y = kad_stack(1, &y); // third pivot
	t = kad_ce_multi(t, y), t->ext_flag = KANN_F_COST;
	return kann_new(t, 0);
}

int dbr_get_n_lbl(const kann_t *ann)
{
	int out_id;
	kad_node_t *p;
	out_id = kann_find(ann, KANN_F_OUT, 0);
	assert(out_id >= 0);
	p = ann->v[out_id];
	assert(p->n_d == 2);
	return p->d[1];
}

void dbr_train(kann_t *ann, dn_seqs_t *dr, int ulen, float lr, int m_epoch, int mbs, int n_threads, int batch_len, const char *fn)
{
	float **x[2], **y, *r, grad_clip = 10.0f, min_cost = 1e30f;
	kann_t *ua;
	int epoch, u, k, n_var, n_lbl;

	n_lbl = dbr_get_n_lbl(ann);
	x[0] = (float**)calloc(ulen, sizeof(float*));
	x[1] = (float**)calloc(ulen, sizeof(float*));
	y    = (float**)calloc(ulen, sizeof(float*));
	for (u = 0; u < ulen; ++u) {
		x[0][u] = (float*)calloc(4 * mbs, sizeof(float));
		x[1][u] = (float*)calloc(4 * mbs, sizeof(float));
		y[u]    = (float*)calloc(n_lbl * mbs, sizeof(float));
	}
	n_var = kann_size_var(ann);
	r = (float*)calloc(n_var, sizeof(float));

	ua = kann_unroll(ann, ulen, ulen, ulen);
	kann_mt(ua, n_threads, mbs);
	kann_set_batch_size(ua, mbs);
	kann_switch(ua, 1);
	kann_feed_bind(ua, KANN_F_IN,    1, x[0]);
	kann_feed_bind(ua, KANN_F_IN,    2, x[1]);
	kann_feed_bind(ua, KANN_F_TRUTH, 0, y);
	for (epoch = 0; epoch < m_epoch; ++epoch) {
		double cost = 0.0;
		int i, j, b, tot = 0, ctot = 0, n_cerr = 0;
		for (i = 0; i < batch_len; i += mbs * ulen) {
			for (u = 0; u < ulen; ++u) {
				memset(x[0][u], 0, 4 * mbs * sizeof(float));
				memset(x[1][u], 0, 4 * mbs * sizeof(float));
				memset(y[u],    0, n_lbl * mbs * sizeof(float));
			}
			for (b = 0; b < mbs; ++b) {
				for (;;) {
					k = dn_select_seq(dr, kann_drand());
					if (dr->len[k] >= ulen) break;
				}
				j = (int)((dr->len[k] - ulen) * kad_drand(0));
				for (u = 0; u < ulen; ++u) {
					int c = (uint8_t)dr->seq[k][j + u];
					int a = dr->lbl[k][j + u];
					if (c >= 4) continue;
					x[0][u][b * 4 + c] = 1.0f;
					x[1][ulen - 1 - u][b * 4 + (3 - c)] = 1.0f;
					y[u][b * n_lbl + a] = 1.0f;
				}
			}
			cost += kann_cost(ua, 0, 1) * ulen * mbs;
			n_cerr += kann_class_error(ua, &b);
			tot += ulen * mbs, ctot += b;
			if (grad_clip > 0.0f) kann_grad_clip(grad_clip, n_var, ua->g);
			kann_RMSprop(n_var, lr, 0, 0.9f, ua->g, ua->x, r);
		}
		fprintf(stderr, "epoch: %d; running cost: %g (class error: %.2f%%)\n", epoch+1, cost / tot, 100.0 * n_cerr / ctot);
		if (fn && cost / tot < min_cost) kann_save(fn, ann);
		if (cost / tot < min_cost) min_cost = cost / tot;
	}
	kann_delete_unrolled(ua);

	for (u = 0; u < ulen; ++u) { free(x[0][u]); free(x[1][u]); free(y[u]); }
	free(r); free(y); free(x[0]); free(x[1]);
}

uint8_t *dbr_predict(kann_t *ua, const char *str)
{
	float **x[2], *z;
	int mbs = -1, ulen;
	int step, n_lbl, l_seq, i, j, b, u, *st;
	uint8_t *lbl;
	kad_node_t *out;

	out = ua->v[kann_find(ua, KANN_F_OUT, 0)];
	assert(out->n_d == 2);
	for (i = 0; i < ua->n; ++i)
		if (ua->v[i]->ext_flag & KANN_F_IN) {
			mbs = ua->v[i]->d[0];
			break;
		}
	assert(mbs > 0);
	assert(out->d[0] % mbs == 0);
	ulen = out->d[0] / mbs;
	n_lbl = out->d[1];
	step = ulen / 2;

	x[0] = (float**)calloc(ulen, sizeof(float*));
	x[1] = (float**)calloc(ulen, sizeof(float*));
	for (u = 0; u < ulen; ++u) {
		x[0][u] = (float*)calloc(4 * mbs, sizeof(float));
		x[1][u] = (float*)calloc(4 * mbs, sizeof(float));
	}
	kann_feed_bind(ua, KANN_F_IN, 1, x[0]);
	kann_feed_bind(ua, KANN_F_IN, 2, x[1]);

	l_seq = strlen(str);
	lbl = (uint8_t*)calloc(l_seq, 1);
	z = (float*)calloc(l_seq, sizeof(float));
	st = (int*)calloc(mbs, sizeof(int));

	for (i = b = 0; i < l_seq; i += step) {
		for (j = i; j < l_seq && j < i + ulen; ++j) {
			int c = seq_nt4_table[(uint8_t)str[j]];
			int u = j - i;
			if (c >= 4) continue;
			x[0][u][b * 4 + c] = 1.0f;
			x[1][ulen - 1 - u][b * 4 + (3 - c)] = 1.0f;
		}
		st[b++] = i;
		if (b == mbs || i + step >= l_seq) {
			int k;
			kann_eval_out(ua);
			for (k = 0; k < b; ++k) {
				for (j = st[k]; j < l_seq && j < st[k] + ulen; ++j) {
					int u = j - st[k], a, max_a;
					float *y = &out->x[(u * mbs + k) * n_lbl], max;
					max_a = 0, max = y[0];
					for (a = 1; a < n_lbl; ++a)
						if (y[a] > max) max = y[a], max_a = a;
					if (max > z[j]) z[j] = max, lbl[j] = max_a;
				}
			}
			for (u = 0; u < ulen; ++u) {
				memset(x[0][u], 0, 4 * mbs * sizeof(float));
				memset(x[1][u], 0, 4 * mbs * sizeof(float));
			}
			b = 0;
		}
	}

	for (u = 0; u < ulen; ++u) { free(x[0][u]); free(x[1][u]); }
	free(x[0]); free(x[1]); free(z); free(st);
	return lbl;
}

int main(int argc, char *argv[])
{
	kann_t *ann = 0;
	int c, n_layer = 1, n_neuron = 128, ulen = 100, to_apply = 0, to_eval = 0;
	int batch_len = 1000000, mbs = 64, m_epoch = 50, n_threads = 1, is_tied = 1, seed = 11;
	float h_dropout = 0.0f, lr = 0.001f;
	char *fn_out = 0, *fn_in = 0;
	ketopt_t o = KETOPT_INIT;

	while ((c = ketopt(&o, argc, argv, 1, "Au:l:n:m:B:o:i:t:Tb:Ed:s:", 0)) >= 0) {
		if (c == 'u') ulen = atoi(o.arg);
		else if (c == 'l') n_layer = atoi(o.arg);
		else if (c == 'n') n_neuron = atoi(o.arg);
		else if (c == 'r') lr = atof(o.arg);
		else if (c == 's') seed = atoi(o.arg);
		else if (c == 'm') m_epoch = atoi(o.arg);
		else if (c == 'd') h_dropout = atof(o.arg);
		else if (c == 'B') mbs = atoi(o.arg);
		else if (c == 'o') fn_out = o.arg;
		else if (c == 'i') fn_in = o.arg;
		else if (c == 'A') to_apply = 1;
		else if (c == 'E') to_eval = to_apply = 1;
		else if (c == 't') n_threads = atoi(o.arg);
		else if (c == 'T') is_tied = 0; // for debugging only; weights should be tiled for DNA sequences
		else if (c == 'b') {
			double x;
			char *s;
			x = strtod(o.arg, &s);
			if (*s == 'g' || *s == 'G') x *= 1e9;
			else if (*s == 'm' || *s == 'M') x *= 1e6;
			else if (*s == 'k' || *s == 'K') x *= 1e3;
			batch_len = (int)(x + .499);
		}
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: dna-brnn [options] <seq.fq>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  General:\n");
		fprintf(stderr, "    -i FILE    read a trained model from FILE []\n");
		fprintf(stderr, "    -o FILE    write model to FILE []\n");
		fprintf(stderr, "    -u INT     window length [%d]\n", ulen);
		fprintf(stderr, "    -B INT     mini-batch size [%d]\n", mbs);
		fprintf(stderr, "    -t INT     number of threads [%d]\n", n_threads);
		fprintf(stderr, "    -A         apply a trained model (req. -i)\n");
		fprintf(stderr, "    -E         apply and evaluate a trained model (req. -i)\n");
		fprintf(stderr, "  Model construction:\n");
		fprintf(stderr, "    -l INT     number of GRU layers [%d]\n", n_layer);
		fprintf(stderr, "    -n INT     number of hidden neurons [%d]\n", n_neuron);
		fprintf(stderr, "    -d FLOAT   dropout rate [%g]\n", h_dropout);
		fprintf(stderr, "  Training:\n");
		fprintf(stderr, "    -r FLOAT   learning rate [%g]\n", lr);
		fprintf(stderr, "    -m INT     number of epochs [%d]\n", m_epoch);
		fprintf(stderr, "    -b INT     batch size [%d]\n", batch_len);
		fprintf(stderr, "    -s INT     PRNG seed [%d]\n", seed);
		return 1;
	}

	fprintf(stderr, "[M::%s] command line: ", __func__);
	for (c = 0; c < argc; ++c) {
		if (c) fprintf(stderr, " ");
		fprintf(stderr, "%s", argv[c]);
	}
	fputc('\n', stderr);

	if (fn_in) ann = kann_load(fn_in);
	if (!to_apply) {
		dn_seqs_t *dr;
		dr = dn_read(argv[o.ind]);
		if (ann == 0) ann = dbr_model_gen(dr->n_lbl, n_layer, n_neuron, h_dropout, is_tied);
		dbr_train(ann, dr, ulen, lr, m_epoch, mbs, n_threads, batch_len, fn_out);
	} else if (ann) {
		gzFile fp;
		kseq_t *ks;
		int n_lbl, *cnt;
		kann_t *ua;

		n_lbl = dbr_get_n_lbl(ann);
		cnt = (int*)calloc(n_lbl * n_lbl, sizeof(int));
		fp = strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(0, "r");
		ks = kseq_init(fp);

		ua = kann_unroll(ann, ulen, ulen, ulen);
		kann_mt(ua, n_threads, mbs);
		kann_set_batch_size(ua, mbs);
		kann_switch(ua, 0);
		while (kseq_read(ks) >= 0) {
			uint8_t *lbl;
			int i;
			lbl = dbr_predict(ua, ks->seq.s);
			if (to_eval && ks->qual.l > 0) {
				for (i = 0; i < ks->seq.l; ++i) {
					int c = ks->qual.s[i] - 33;
					if (c < 0 || c >= n_lbl) continue;
					++cnt[c * n_lbl + lbl[i]];
				}
			}
			printf("@%s\n", ks->name.s);
			puts(ks->seq.s);
			printf("+\n");
			for (i = 0; i < ks->seq.l; ++i) lbl[i] += 33;
			fwrite(lbl, 1, ks->seq.l, stdout);
			putchar('\n');
			free(lbl);
		}
		kann_delete_unrolled(ua);

		kseq_destroy(ks);
		gzclose(fp);

		if (to_eval) {
			int a, b;
			for (a = 0; a < n_lbl; ++a) {
				fprintf(stderr, "[M::%s] true label %d:", __func__, a);
				for (b = 0; b < n_lbl; ++b)
					fprintf(stderr, " %11d", cnt[a * n_lbl + b]);
				fputc('\n', stderr);
			}
		}
		free(cnt);
	}
	return 0;
}
