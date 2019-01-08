#ifndef MSS_H
#define MSS_H

#define MSS_FLOAT double

typedef struct {
    int st, en;
    MSS_FLOAT sc;
} msseg_t;

msseg_t *mss_find_all(int n, const MSS_FLOAT *S, MSS_FLOAT min_sc, int *n_seg);

#endif
