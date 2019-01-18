/* Classify protein sequences
 *
 * Copyright 2014 Peter Meinicke, Robin Martinjak
 *
 * This file is part of libuproc.
 *
 * libuproc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * libuproc is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libuproc.  If not, see <http://www.gnu.org/licenses/>.
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>
#include <string.h>

#include "uproc/common.h"
#include "uproc/error.h"
#include "uproc/bst.h"
#include "uproc/list.h"
#include "uproc/protclass.h"

struct uproc_protclass_s
{
    enum uproc_protclass_mode mode;
    const uproc_substmat *substmat;
    const uproc_ecurve *fwd;
    const uproc_ecurve *rev;
    uproc_protfilter *filter;
    void *filter_arg;
    struct uproc_protclass_trace
    {
        uproc_protclass_trace_cb *cb;
        void *cb_arg;
    } trace;
};

/*********************
 * score computation *
 *********************/

struct sc
{
    size_t index;
    double total, dist[UPROC_WORD_LEN];
};

static void
reverse_array(void *p, size_t n, size_t sz)
{
    unsigned char *s = p, tmp;
    size_t i, k, i1, i2;

    for (i = 0; i < n / 2; i++) {
        for (k = 0; k < sz; k++) {
            i1 = sz * i + k;
            i2 = sz * (n - i - 1) + k;
            tmp = s[i1];
            s[i1] = s[i2];
            s[i2] = tmp;
        }
    }
}

static void
sc_init(struct sc *s)
{
    size_t i;
    s->index = -1;
    s->total = 0.0;
    for (i = 0; i < UPROC_WORD_LEN; i++) {
        s->dist[i] = -INFINITY;
    }
}

static void
sc_add(struct sc *score, size_t index, double dist[static UPROC_SUFFIX_LEN],
       bool reverse)
{
    size_t i, diff;
    double tmp[UPROC_WORD_LEN];

    for (i = 0; i < UPROC_PREFIX_LEN; i++) {
        tmp[i] = -INFINITY;
    }
    memcpy(tmp + UPROC_PREFIX_LEN, dist, sizeof *dist * UPROC_SUFFIX_LEN);
    if (reverse) {
        reverse_array(tmp, UPROC_WORD_LEN, sizeof *tmp);
    }

    if (score->index != (size_t) -1) {
        diff = index - score->index;
        if (diff > UPROC_WORD_LEN) {
            diff = UPROC_WORD_LEN;
        }
        for (i = 0; i < diff; i++) {
            if (isfinite(score->dist[i])) {
                score->total += score->dist[i];
                score->dist[i] = -INFINITY;
            }
        }
    }
    else {
        diff = 0;
    }

    for (i = 0; i + diff < UPROC_WORD_LEN; i++) {
#define MAX(a, b) (a > b ? a : b)
        score->dist[i] = MAX(score->dist[i + diff], tmp[i]);
    }
    for (; i < UPROC_WORD_LEN; i++) {
        score->dist[i] = tmp[i];
    }
    score->index = index;
}

static double
sc_finalize(struct sc *score)
{
    size_t i;
    for (i = 0; i < UPROC_WORD_LEN; i++) {
        if (isfinite(score->dist[i])) {
            score->total += score->dist[i];
        }
    }
    return score->total;
}

static int
scores_add(uproc_bst *scores, uproc_family family, size_t index,
           double dist[static UPROC_SUFFIX_LEN], bool reverse)
{
    struct sc sc;
    union uproc_bst_key key = { .uint = family };
    sc_init(&sc);
    (void) uproc_bst_get(scores, key, &sc);
    sc_add(&sc, index, dist, reverse);
    return uproc_bst_update(scores, key, &sc);
}


static int
scores_add_word(const uproc_protclass *pc, uproc_bst *scores,
                const struct uproc_word *word,
                size_t index, bool reverse, const uproc_ecurve *ecurve,
                const uproc_substmat *substmat)
{
    int res;
    struct uproc_word
        lower_nb = UPROC_WORD_INITIALIZER,
        upper_nb = UPROC_WORD_INITIALIZER;
    uproc_family lower_family, upper_family;
    double dist[UPROC_SUFFIX_LEN];

    if (!ecurve) {
        return 0;
    }
    uproc_ecurve_lookup(ecurve, word, &lower_nb, &lower_family, &upper_nb,
                        &upper_family);
    uproc_substmat_align_suffixes(substmat, word->suffix, lower_nb.suffix,
                                  dist);
    if (pc->trace.cb) {
        pc->trace.cb(&lower_nb, lower_family, index, reverse, dist,
                     pc->trace.cb_arg);
    }
    res = scores_add(scores, lower_family, index, dist, reverse);
    if (res || !uproc_word_cmp(&lower_nb, &upper_nb)) {
        return res;
    }
    uproc_substmat_align_suffixes(substmat, word->suffix, upper_nb.suffix,
                                  dist);
    if (pc->trace.cb) {
        pc->trace.cb(&upper_nb, upper_family, index, reverse, dist,
                     pc->trace.cb_arg);
    }
    res = scores_add(scores, upper_family, index, dist, reverse);
    return res;
}

static int
scores_compute(const struct uproc_protclass_s *pc, const char *seq,
               uproc_bst *scores)
{
    int res;
    uproc_worditer *iter;
    size_t index;
    struct uproc_word
        fwd_word = UPROC_WORD_INITIALIZER,
        rev_word = UPROC_WORD_INITIALIZER;

    iter = uproc_worditer_create(seq, uproc_ecurve_alphabet(pc->fwd));
    if (!iter) {
        return -1;
    }

    while (res = uproc_worditer_next(iter, &index, &fwd_word, &rev_word),
           !res)
    {
        res = scores_add_word(pc, scores, &fwd_word, index, false, pc->fwd,
                              pc->substmat);
        if (res) {
            break;
        }
        res = scores_add_word(pc, scores, &rev_word, index, true, pc->rev,
                              pc->substmat);
        if (res) {
            break;
        }
    }
    uproc_worditer_destroy(iter);
    return res == -1 ? -1 : 0;
}


/****************
 * finalization *
 ****************/

static int
scores_finalize(const struct uproc_protclass_s *pc, const char *seq,
                uproc_bst *score_tree, uproc_list *results)
{
    int res = 0;
    uproc_bstiter *iter;
    union uproc_bst_key key;
    struct sc value;
    size_t seq_len = strlen(seq);
    struct uproc_protresult pred, pred_max = { .score = -INFINITY };

    iter = uproc_bstiter_create(score_tree);
    if (!iter) {
        return -1;
    }
    while (!uproc_bstiter_next(iter, &key, &value)) {
        uproc_family family = key.uint;
        double score = sc_finalize(&value);
        if (pc->filter &&
            !pc->filter(seq, seq_len, family, score, pc->filter_arg)) {
            continue;
        }
        pred.score = score;
        pred.family = family;
        if (pc->mode == UPROC_PROTCLASS_MAX) {
            if (!uproc_list_size(results)) {
                pred_max = pred;
                res = uproc_list_append(results, &pred);
                if (res) {
                    break;
                }
            }
            else if (pred.score > pred_max.score) {
                pred_max = pred;
                uproc_list_set(results, 0, &pred_max);
            }
        }
        else {
            uproc_list_append(results, &pred);
        }
    }
    uproc_bstiter_destroy(iter);
    return res;
}


/**********************
 * exported functions *
 **********************/

uproc_protclass *
uproc_protclass_create(enum uproc_protclass_mode mode, const uproc_ecurve *fwd,
                       const uproc_ecurve *rev, const uproc_substmat *substmat,
                       uproc_protfilter *filter, void *filter_arg)
{
    struct uproc_protclass_s *pc;
    if (!(fwd || rev)) {
        uproc_error_msg(UPROC_EINVAL,
                        "protein classifier requires at least one ecurve");
        return NULL;
    }
    pc = malloc(sizeof *pc);
    if (!pc) {
        uproc_error(UPROC_ENOMEM);
        return NULL;
    }
    *pc = (struct uproc_protclass_s) {
        .mode = mode,
        .substmat = substmat,
        .fwd = fwd,
        .rev = rev,
        .filter = filter,
        .filter_arg = filter_arg,
        .trace = {
            .cb = NULL,
            .cb_arg = NULL,
        },
    };
    return pc;
}

void
uproc_protclass_destroy(uproc_protclass *pc)
{
    free(pc);
}


static void
map_list_protresult_free(void *value, void *opaque)
{
    (void) opaque;
    uproc_protresult_free(value);
}


int
uproc_protclass_classify(const uproc_protclass *pc, const char *seq,
                         uproc_list **results)
{
    int res;
    uproc_bst *scores;

    if (!*results) {
        *results = uproc_list_create(sizeof (struct uproc_protresult));
        if (!*results) {
            return -1;
        }
    }
    else {
        uproc_list_map(*results, map_list_protresult_free, NULL);
        uproc_list_clear(*results);
    }

    scores = uproc_bst_create(UPROC_BST_UINT, sizeof (struct sc));
    if (!scores) {
        return -1;
    }
    res = scores_compute(pc, seq, scores);
    if (res || uproc_bst_isempty(scores)) {
        goto error;
    }
    res = scores_finalize(pc, seq, scores, *results);
error:
    uproc_bst_destroy(scores);
    return res;
}

void
uproc_protclass_set_trace(uproc_protclass *pc, uproc_protclass_trace_cb *cb,
                          void *cb_arg)
{
    pc->trace.cb = cb;
    pc->trace.cb_arg = cb_arg;
}

void
uproc_protresult_init(struct uproc_protresult *results)
{
    *results = (struct uproc_protresult) UPROC_PROTRESULT_INITIALIZER;
}

void
uproc_protresult_free(struct uproc_protresult *results)
{
    (void) results;
}

int
uproc_protresult_copy(struct uproc_protresult *dest,
                      const struct uproc_protresult *src)
{
    *dest = *src;
    return 0;
}
