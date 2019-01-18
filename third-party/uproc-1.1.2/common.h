/* Copyright 2014 Peter Meinicke, Robin Martinjak
 *
 * This file is part of uproc.
 *
 * uproc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * uproc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with uproc.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COMMON_H
#define COMMON_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_GETOPT_LONG
#define _GNU_SOURCE
#include <getopt.h>
#else
#include <unistd.h>
#endif

#include <uproc.h>


/* Stringify a macro value */
#define STR1(x) #x
#define STR(x) STR1(x)


/* Open file for reading or stdin if `path` is "-" */
uproc_io_stream *open_read(const char *path);

/* Open file for writing or stdout if `path` is "-" */
uproc_io_stream *open_write(const char *path, enum uproc_io_type type);


/* Print the program version and copyright note */
void print_version(const char *progname);


/* Parse int from string */
int parse_int(const char *arg, int *x);

/* Parse int and check whether it is 0, 2 or 3 */
int parse_prot_thresh_level(const char *arg, int *x);

/* Parse int and check whether it is 0, 1 or 2 */
int parse_orf_thresh_level(const char *arg, int *x);


/* Error handler that prints the message and exits the program */
void errhandler_bail(enum uproc_error_code num, const char *msg,
                     const char *loc);


/* Trim to the first word containing neither a comma nor any whitespace */
void trim_header(char *s);


/* Create dir (or fail silently) */
void make_dir(const char *path);


/* Struct representing the "database" */
struct database
{
    uproc_ecurve *fwd, *rev;
    uproc_idmap *idmap;
    uproc_matrix *prot_thresh;
};

#define DATABASE_INITIALIZER { 0, 0, 0, 0 }

int database_load(struct database *db, const char *path, int prot_thresh_level,
                  enum uproc_ecurve_format format);
void database_free(struct database *db);


/* Struct representing the "model" */
struct model
{
    uproc_substmat *substmat;
    uproc_matrix *codon_scores, *orf_thresh;
};

#define MODEL_INITIALIZER { 0, 0, 0 }

/* Load model */
int model_load(struct model *m, const char *path, int orf_thresh_level);

void model_free(struct model *m);

/* Create classifiers
 *
 * `dc` may be NULL if no DNA classifier is needed
 * */
int create_classifiers(uproc_protclass **pc, uproc_dnaclass **dc,
                       const struct database *db, const struct model *model,
                       bool short_read_mode);


#if defined(TIMEIT) && HAVE_CLOCK_GETTIME
#include <time.h>
typedef struct {
    int running;
    double total;
    struct timespec start;
} timeit;

#define TIMEIT_INITIALIZER { 0, 0, { 0, 0 } }

void timeit_start(timeit *t);
void timeit_stop(timeit *t);
void timeit_print(timeit *t, const char *s);
#else
typedef int timeit;
#define TIMEIT_INITIALIZER 0
#define timeit_start(...) ((void)0)
#define timeit_stop(...) ((void)0)
#define timeit_print(...) ((void)0)
#endif


void progress(uproc_io_stream *stream, const char *new_label, double percent);
#endif
