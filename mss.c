#include <stdlib.h>
#include "mss.h"

#define kvec_t(type) struct { size_t n, m; type *a; }

#define kv_push(type, v, x) do { \
        if ((v).n == (v).m) { \
            (v).m = (v).m? (v).m<<1 : 2; \
            (v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
        } \
        (v).a[(v).n++] = (x); \
    } while (0)

#define kv_pushp(type, v, p) do { \
        if ((v).n == (v).m) { \
            (v).m = (v).m? (v).m<<1 : 2; \
            (v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
        } \
        *(p) = &(v).a[(v).n++]; \
    } while (0)

typedef struct {
    int st, en;
    MSS_FLOAT L, R;
    int pre;
} msseg_aux_t;

typedef kvec_t(msseg_t) msseg_v;
typedef kvec_t(msseg_aux_t) msseg_aux_v;

static void add_segs(msseg_v *ret, msseg_aux_v *seg, int min_sc)
{
    int i;
    for (i = 0; i < seg->n; ++i) {
        msseg_aux_t *p = &seg->a[i];
        if (p->R - p->L >= min_sc) {
            msseg_t *q;
            kv_pushp(msseg_t, *ret, &q);
            q->st = p->st, q->en = p->en, q->sc = p->R - p->L;
        }
    }
    seg->n = 0;
}

msseg_t *mss_find_all(int n, const MSS_FLOAT *S, MSS_FLOAT min_sc, int *n_seg)
{
    int i, j;
    MSS_FLOAT L;
    msseg_v ret = {0,0,0};
    msseg_aux_v seg = {0,0,0};
    msseg_aux_t t;

    for (i = L = 0; i < n;) {
        if (S[i] > 0) {
            int k;
            MSS_FLOAT R = L + S[i];
            for (k = i + 1; k < n && S[k] > 0.; ++k)
                R += S[k];
            t.st = i, t.en = k, t.L = L, t.R = R;
            while (1) {
                msseg_aux_t *p;
                for (j = seg.n - 1; j >= 0;) {
                    p = &seg.a[j];
                    if (p->L < t.L) break;
                    j = p->pre >= 0? p->pre : j - 1;
                }
                if (j >= 0 && seg.a[j].R < t.R) {
                    p = &seg.a[j];
                    t.st = p->st, t.L = p->L, t.pre = p->pre;
                    seg.n = j;
                } else {
                    if (j < 0) add_segs(&ret, &seg, min_sc);
                    t.pre = j;
                    kv_push(msseg_aux_t, seg, t);
                    break;
                }
            }
            L = R, i = k;
        } else L += S[i++];
    }
    add_segs(&ret, &seg, min_sc);
    free(seg.a);
    ret.a = (msseg_t*)realloc(ret.a, ret.n * sizeof(msseg_t));
    *n_seg = ret.n;
    return ret.a;
}
