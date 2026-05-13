#define _POSIX_C_SOURCE 200809L

/*
 * OptJV3 2026
 *
 * Genetic Algorithm optimizer for JARVIS3 parameter search.
 *
 * Design:
 *   - black-box optimization over the JARVIS3 CLI parameter space
 *   - thread-parallel candidate evaluation (one JARVIS3 compression per worker)
 *   - configurable search space for:
 *       * number of context models (-cm)
 *       * number of repeat models (-rm)
 *       * global neural mixer parameters (-hs, -lr, -sd)
 *   - robust temp-file isolation for concurrent runs
 *   - optional decompression verification
 *   - optional restart from a known-good parameter string
 *
 * Repeat-model order used everywhere in this program:
 *   nr:ctx:beta:limit:gamma:ir:weight:cache
 *
 * Compile:
 *   gcc -O3 -march=native -flto=auto -std=c11 -pthread -Wall -Wextra -o GA ga.c -lm
 */

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

#define HARD_MAX_CMODELS 16
#define HARD_MAX_RMODELS 8
#define MAX_ARGV_ITEMS   256
#define MAX_CLI_CHARS    8192
#define CACHE_LOAD_NUM   7
#define CACHE_LOAD_DEN   10

#define INVALID_FITNESS 1.0e300

/* Parser-accurate hard limits from JARVIS3 */
#define J3_CM_CTX_MIN      1
#define J3_CM_CTX_MAX      14
#define J3_CM_DEN_MIN      1
#define J3_CM_DEN_MAX      5000
#define J3_CM_IR_MIN       0
#define J3_CM_IR_MAX       2
#define J3_CM_EDITS_MIN    0
#define J3_CM_EDITS_MAX    20
#define J3_CM_EDEN_MIN     1
#define J3_CM_EDEN_MAX     50000
#define J3_CM_EIR_MIN      0
#define J3_CM_EIR_MAX      1

#define J3_RM_NR_MIN       1
#define J3_RM_NR_MAX       100000
#define J3_RM_CTX_MIN      1
#define J3_RM_CTX_MAX      14
#define J3_RM_LIMIT_MIN    1
#define J3_RM_LIMIT_MAX    20
#define J3_RM_IR_MIN       0
#define J3_RM_IR_MAX       2
#define J3_RM_CACHE_MIN    1
#define J3_RM_CACHE_MAX    15

typedef enum {
  OBJ_BYTES = 0,
  OBJ_BPS   = 1,
  OBJ_BYTES_PLUS_TIME = 2
} objective_t;

typedef struct { int lo, hi; } irange_t;
typedef struct { double lo, hi; } drange_t;
typedef struct { long long lo, hi; } lrange_t;

typedef struct {
  int enabled;
  int ctx;
  int den;
  int ir;
  double gamma;
  int edits;
  int eden;
  int eir;
  double egamma;
} cmodel_gene_t;

typedef struct {
  int enabled;
  int nr;
  int ctx;
  double beta;
  int limit;
  double gamma;
  int ir;
  double weight;
  long long cache;
} rmodel_gene_t;

typedef struct {
  int hs;
  double lr;
  int seed;

  int n_active_c;
  int n_active_r;

  cmodel_gene_t c[HARD_MAX_CMODELS];
  rmodel_gene_t r[HARD_MAX_RMODELS];

  double fitness;
  size_t comp_bytes;
  double elapsed_s;
  int valid;

  char key[MAX_CLI_CHARS];
  char param_string[MAX_CLI_CHARS];
} individual_t;

typedef struct {
  char jarvis_path[PATH_MAX];
  char input_path[PATH_MAX];
  char workdir[PATH_MAX];
  char best_out[PATH_MAX];
  char history_out[PATH_MAX];

  size_t input_bytes;

  int population;
  int generations;
  int threads;
  int elite_count;
  int tournament_size;

  double crossover_rate;
  double mutation_rate;
  double toggle_rate;
  double blend_alpha;

  int max_cmodels;
  int max_rmodels;
  int min_cmodels;
  int min_rmodels;

  int verify;
  int keep_temps;
  int quiet;

  int optimize_hs;
  int optimize_lr;
  int optimize_seed;

  int fixed_hs;
  double fixed_lr;
  int fixed_seed;

  irange_t hs_range;
  drange_t lr_range;
  irange_t seed_range;

  irange_t cm_ctx_range;
  irange_t cm_den_range;
  irange_t cm_ir_range;
  drange_t cm_gamma_range;
  irange_t cm_edits_range;
  irange_t cm_eden_range;
  irange_t cm_eir_range;
  drange_t cm_egamma_range;

  irange_t rm_nr_range;
  irange_t rm_ctx_range;
  drange_t rm_beta_range;
  irange_t rm_limit_range;
  drange_t rm_gamma_range;
  irange_t rm_ir_range;
  drange_t rm_weight_range;
  lrange_t rm_cache_range;

  objective_t objective;
  double time_weight;

  uint64_t master_seed;

  int have_restart;
  char restart_from[MAX_CLI_CHARS];
} ga_config_t;

typedef struct {
  double fitness;
  size_t comp_bytes;
  double elapsed_s;
  int valid;
} eval_result_t;

typedef struct {
  char *key;
  eval_result_t value;
  int used;
} cache_entry_t;

typedef struct {
  cache_entry_t *table;
  size_t cap;
  size_t count;
  pthread_mutex_t mutex;
} fitness_cache_t;

typedef struct {
  const ga_config_t *cfg;
  individual_t *pop;
  int *job_indices;
  int job_count;
  int next_job;
  pthread_mutex_t next_mutex;
  fitness_cache_t *cache;
} worker_ctx_t;

typedef struct {
  uint64_t s;
} rng_t;

/* ------------------------------------------------------------------------- */
/* Utilities                                                                 */
/* ------------------------------------------------------------------------- */

static void die(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  exit(EXIT_FAILURE);
}

static void *xcalloc(size_t n, size_t s)
{
  void *p = calloc(n, s);
  if (!p) die("Out of memory (%zu bytes)", n * s);
  return p;
}

static char *xstrdup(const char *s)
{
  char *p = strdup(s ? s : "");
  if (!p) die("Out of memory duplicating string");
  return p;
}

static int file_exists_exec(const char *path)
{
  return access(path, X_OK) == 0;
}

static int file_exists_read(const char *path)
{
  return access(path, R_OK) == 0;
}

static size_t get_file_size_or_die(const char *path)
{
  struct stat st;
  if (stat(path, &st) != 0)
    die("Cannot stat '%s': %s", path, strerror(errno));
  if (st.st_size < 0)
    die("Negative file size for '%s'", path);
  return (size_t) st.st_size;
}

static int ensure_dir(const char *path)
{
  struct stat st;
  if (stat(path, &st) == 0) {
    if (!S_ISDIR(st.st_mode))
      return -1;
    return 0;
  }
  if (mkdir(path, 0755) != 0)
    return -1;
  return 0;
}

static int files_equal(const char *a, const char *b)
{
  FILE *fa = fopen(a, "rb");
  FILE *fb = fopen(b, "rb");
  unsigned char ba[1 << 15];
  unsigned char bb[1 << 15];
  size_t na, nb;
  int eq = 1;

  if (!fa || !fb) {
    if (fa) fclose(fa);
    if (fb) fclose(fb);
    return 0;
  }

  while (1) {
    na = fread(ba, 1, sizeof(ba), fa);
    nb = fread(bb, 1, sizeof(bb), fb);
    if (na != nb || memcmp(ba, bb, na) != 0) {
      eq = 0;
      break;
    }
    if (na == 0)
      break;
  }

  fclose(fa);
  fclose(fb);
  return eq;
}

static double now_sec(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
}

static void path_join(char *dst, size_t cap, const char *a, const char *b)
{
  if (snprintf(dst, cap, "%s/%s", a, b) >= (int) cap)
    die("Path too long while joining '%s' and '%s'", a, b);
}

static void remove_if_exists(const char *path)
{
  if (unlink(path) != 0 && errno != ENOENT) {
    // ignore cleanup failures 
  }
}

static void rmdir_if_exists(const char *path)
{
  if (rmdir(path) != 0 && errno != ENOENT) {
    // ignore cleanup failures 
  }
}

static void clamp_irange(irange_t *r, int lo, int hi)
{
  if (r->lo < lo) r->lo = lo;
  if (r->hi > hi) r->hi = hi;
  if (r->lo > r->hi) r->lo = r->hi = lo;
}

static void clamp_drange(drange_t *r, double lo, double hi)
{
  if (r->lo < lo) r->lo = lo;
  if (r->hi > hi) r->hi = hi;
  if (r->lo > r->hi) r->lo = r->hi = lo;
}

static double clampd(double x, double lo, double hi)
{
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static int clampi(int x, int lo, int hi)
{
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static long long clampll(long long x, long long lo, long long hi)
{
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static double quantize_lr(double lr)
{
  if (lr <= 0.0) return 0.0;
  return floor(lr * 65534.0) / 65534.0;
}

static int appendf(char *dst, size_t cap, const char *fmt, ...)
{
  size_t used = strlen(dst);
  int written;
  va_list ap;
  if (used >= cap) return -1;
  va_start(ap, fmt);
  written = vsnprintf(dst + used, cap - used, fmt, ap);
  va_end(ap);
  if (written < 0 || used + (size_t) written >= cap)
    return -1;
  return 0;
}

static double compute_bps(size_t comp_bytes, size_t input_bytes)
{
  if (input_bytes == 0) return 0.0;
  return (double) comp_bytes * 8.0 / (double) input_bytes;
}

static int format_cmodel_spec(char *dst, size_t cap, const cmodel_gene_t *g)
{
  int n = snprintf(dst, cap, "%d:%d:%d:%.3f/%d:%d:%d:%.3f",
                   g->ctx, g->den, g->ir, g->gamma,
                   g->edits, g->eden, g->eir, g->egamma);
  return (n < 0 || (size_t) n >= cap) ? -1 : 0;
}

static int format_rmodel_spec(char *dst, size_t cap, const rmodel_gene_t *g)
{
  int n = snprintf(dst, cap, "%d:%d:%.3f:%d:%.3f:%d:%.3f:%lld",
                   g->nr, g->ctx, g->beta, g->limit,
                   g->gamma, g->ir, g->weight, g->cache);
  return (n < 0 || (size_t) n >= cap) ? -1 : 0;
}

static void build_param_string(const ga_config_t *cfg, individual_t *ind)
{
  int i;
  ind->param_string[0] = '\0';

  if (appendf(ind->param_string, sizeof(ind->param_string),
              "-hs %d -lr %.3f -sd %d",
              ind->hs, quantize_lr(ind->lr), ind->seed) != 0)
    die("Parameter string too long");

  for (i = 0; i < cfg->max_cmodels; ++i) {
    char cm[256];
    const cmodel_gene_t *g = &ind->c[i];
    if (!g->enabled) continue;
    if (format_cmodel_spec(cm, sizeof(cm), g) != 0)
      die("C-model string too long");
    if (appendf(ind->param_string, sizeof(ind->param_string), " -cm %s", cm) != 0)
      die("Parameter string too long");
  }

  for (i = 0; i < cfg->max_rmodels; ++i) {
    char rm[256];
    const rmodel_gene_t *g = &ind->r[i];
    if (!g->enabled) continue;
    if (format_rmodel_spec(rm, sizeof(rm), g) != 0)
      die("R-model string too long");
    if (appendf(ind->param_string, sizeof(ind->param_string), " -rm %s", rm) != 0)
      die("Parameter string too long");
  }

  snprintf(ind->key, sizeof(ind->key), "%s", ind->param_string);
}

static void print_full_command_line(FILE *fp, const ga_config_t *cfg, const individual_t *ind)
{
  fprintf(fp, "%s %s %s\n", cfg->jarvis_path, ind->param_string, cfg->input_path);
}

static int validate_cmodel_gene(const cmodel_gene_t *g)
{
  if (!g->enabled) return 1;

  if (g->ctx < J3_CM_CTX_MIN || g->ctx > J3_CM_CTX_MAX) return 0;
  if (g->den < J3_CM_DEN_MIN || g->den > J3_CM_DEN_MAX) return 0;
  if (g->ir  < J3_CM_IR_MIN  || g->ir  > J3_CM_IR_MAX)  return 0;
  if (g->gamma <= 0.0 || g->gamma >= 1.0) return 0;
  if (g->edits < J3_CM_EDITS_MIN || g->edits > J3_CM_EDITS_MAX) return 0;
  if (g->eden < J3_CM_EDEN_MIN || g->eden > J3_CM_EDEN_MAX) return 0;
  if (g->eir  < J3_CM_EIR_MIN  || g->eir  > J3_CM_EIR_MAX)  return 0;
  if (g->egamma <= 0.0 || g->egamma >= 1.0) return 0;

  return 1;
}

static int validate_rmodel_gene(const rmodel_gene_t *g)
{
  if (!g->enabled) return 1;

  if (g->nr    < J3_RM_NR_MIN    || g->nr    > J3_RM_NR_MAX)    return 0;
  if (g->ctx   < J3_RM_CTX_MIN   || g->ctx   > J3_RM_CTX_MAX)   return 0;
  if (!(g->beta   > 0.0 && g->beta   < 1.0)) return 0;
  if (g->limit < J3_RM_LIMIT_MIN || g->limit > J3_RM_LIMIT_MAX) return 0;
  if (!(g->gamma  > 0.0 && g->gamma  < 1.0)) return 0;
  if (g->ir    < J3_RM_IR_MIN    || g->ir    > J3_RM_IR_MAX)    return 0;
  if (!(g->weight > 0.0 && g->weight < 1.0)) return 0;
  if (g->cache < J3_RM_CACHE_MIN || g->cache > J3_RM_CACHE_MAX) return 0;

  return 1;
}

static int validate_individual_for_jarvis(const ga_config_t *cfg, const individual_t *ind)
{
  int i;

  (void) cfg;

  for (i = 0; i < HARD_MAX_CMODELS; ++i)
    if (!validate_cmodel_gene(&ind->c[i])) return 0;

  for (i = 0; i < HARD_MAX_RMODELS; ++i)
    if (!validate_rmodel_gene(&ind->r[i])) return 0;

  return 1;
}

/* ------------------------------------------------------------------------- */
/* RNG                                                                       */
/* ------------------------------------------------------------------------- */

static uint64_t splitmix64_next(uint64_t *x)
{
  uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static void rng_seed(rng_t *r, uint64_t seed)
{
  if (seed == 0) seed = 1;
  r->s = seed;
}

static uint64_t rng_u64(rng_t *r)
{
  return splitmix64_next(&r->s);
}

static double rng_unit(rng_t *r)
{
  return (rng_u64(r) >> 11) * (1.0 / 9007199254740992.0);
}

static int rng_bool(rng_t *r, double p)
{
  return rng_unit(r) < p;
}

static int rng_int_closed(rng_t *r, int lo, int hi)
{
  uint64_t span;
  if (hi <= lo) return lo;
  span = (uint64_t) (hi - lo + 1);
  return lo + (int) (rng_u64(r) % span);
}

static long long rng_ll_closed(rng_t *r, long long lo, long long hi)
{
  uint64_t span;
  if (hi <= lo) return lo;
  span = (uint64_t) (hi - lo + 1);
  return lo + (long long) (rng_u64(r) % span);
}

static double rng_double_closed(rng_t *r, double lo, double hi)
{
  if (hi <= lo) return lo;
  return lo + (hi - lo) * rng_unit(r);
}

static int mutate_int_local(rng_t *r, int value, irange_t range)
{
  int width = range.hi - range.lo + 1;
  int step = width / 10;
  if (step < 1) step = 1;
  if (rng_bool(r, 0.5))
    return rng_int_closed(r, range.lo, range.hi);
  return clampi(value + rng_int_closed(r, -step, step), range.lo, range.hi);
}

static long long mutate_ll_local(rng_t *r, long long value, lrange_t range)
{
  long long width = range.hi - range.lo + 1;
  long long step = width / 10;
  if (step < 1) step = 1;
  if (rng_bool(r, 0.5))
    return rng_ll_closed(r, range.lo, range.hi);
  return clampll(value + rng_ll_closed(r, -step, step), range.lo, range.hi);
}

static double mutate_double_local(rng_t *r, double value, drange_t range)
{
  double width = range.hi - range.lo;
  double step = width * 0.12;
  if (step < 1e-12) step = width > 0.0 ? width : 0.001;
  if (rng_bool(r, 0.45))
    return rng_double_closed(r, range.lo, range.hi);
  return clampd(value + (rng_double_closed(r, -1.0, 1.0) * step), range.lo, range.hi);
}

static int cross_int(rng_t *r, int a, int b, irange_t range)
{
  if (rng_bool(r, 0.5))
    return a;
  if (rng_bool(r, 0.5))
    return b;
  return clampi((a + b) / 2, range.lo, range.hi);
}

static long long cross_ll(rng_t *r, long long a, long long b, lrange_t range)
{
  if (rng_bool(r, 0.5))
    return a;
  if (rng_bool(r, 0.5))
    return b;
  return clampll((a + b) / 2, range.lo, range.hi);
}

static double cross_double(rng_t *r, double a, double b, drange_t range, double alpha)
{
  double lo = a < b ? a : b;
  double hi = a > b ? a : b;
  double span = hi - lo;
  if (rng_bool(r, 0.4))
    return rng_bool(r, 0.5) ? a : b;
  return clampd(rng_double_closed(r, lo - alpha * span, hi + alpha * span), range.lo, range.hi);
}

/* ------------------------------------------------------------------------- */
/* Cache                                                                     */
/* ------------------------------------------------------------------------- */

static uint64_t fnv1a64(const char *s)
{
  uint64_t h = 1469598103934665603ULL;
  while (*s) {
    h ^= (unsigned char) *s++;
    h *= 1099511628211ULL;
  }
  return h;
}

static size_t next_pow2(size_t x)
{
  size_t p = 1;
  while (p < x) p <<= 1;
  return p;
}

static void cache_init(fitness_cache_t *c, size_t expected_entries)
{
  c->cap = next_pow2(expected_entries * 2 + 32);
  c->count = 0;
  c->table = (cache_entry_t *) xcalloc(c->cap, sizeof(cache_entry_t));
  pthread_mutex_init(&c->mutex, NULL);
}

static void cache_destroy(fitness_cache_t *c)
{
  size_t i;
  if (!c) return;
  for (i = 0; i < c->cap; ++i)
    free(c->table[i].key);
  free(c->table);
  pthread_mutex_destroy(&c->mutex);
}

static int cache_lookup_unsafe(const fitness_cache_t *c, const char *key, eval_result_t *out)
{
  size_t mask = c->cap - 1;
  size_t idx = (size_t) fnv1a64(key) & mask;
  size_t start = idx;

  while (c->table[idx].used) {
    if (strcmp(c->table[idx].key, key) == 0) {
      if (out) *out = c->table[idx].value;
      return 1;
    }
    idx = (idx + 1) & mask;
    if (idx == start) break;
  }
  return 0;
}

static void cache_grow_unsafe(fitness_cache_t *c)
{
  size_t old_cap = c->cap;
  cache_entry_t *old_table = c->table;
  size_t i;

  c->cap <<= 1;
  c->table = (cache_entry_t *) xcalloc(c->cap, sizeof(cache_entry_t));
  c->count = 0;

  for (i = 0; i < old_cap; ++i) {
    if (old_table[i].used) {
      size_t mask = c->cap - 1;
      size_t idx = (size_t) fnv1a64(old_table[i].key) & mask;
      while (c->table[idx].used)
        idx = (idx + 1) & mask;
      c->table[idx] = old_table[i];
      c->count++;
    }
  }
  free(old_table);
}

static void cache_insert_unsafe(fitness_cache_t *c, const char *key, eval_result_t val)
{
  size_t mask;
  size_t idx;

  if ((c->count + 1) * CACHE_LOAD_DEN > c->cap * CACHE_LOAD_NUM)
    cache_grow_unsafe(c);

  mask = c->cap - 1;
  idx = (size_t) fnv1a64(key) & mask;
  while (c->table[idx].used) {
    if (strcmp(c->table[idx].key, key) == 0) {
      c->table[idx].value = val;
      return;
    }
    idx = (idx + 1) & mask;
  }

  c->table[idx].used = 1;
  c->table[idx].key = xstrdup(key);
  c->table[idx].value = val;
  c->count++;
}

static int cache_lookup(fitness_cache_t *c, const char *key, eval_result_t *out)
{
  int found;
  pthread_mutex_lock(&c->mutex);
  found = cache_lookup_unsafe(c, key, out);
  pthread_mutex_unlock(&c->mutex);
  return found;
}

static void cache_insert(fitness_cache_t *c, const char *key, eval_result_t val)
{
  pthread_mutex_lock(&c->mutex);
  cache_insert_unsafe(c, key, val);
  pthread_mutex_unlock(&c->mutex);
}

/* ------------------------------------------------------------------------- */
/* Bounds parsing                                                            */
/* ------------------------------------------------------------------------- */

static void parse_key_value_pair(const char *token, char *key, size_t kcap, char *value, size_t vcap)
{
  const char *eq = strchr(token, '=');
  size_t klen;
  if (!eq)
    die("Invalid bounds token '%s' (expected key=min:max)", token);
  klen = (size_t) (eq - token);
  if (klen == 0 || klen >= kcap)
    die("Invalid bounds key in '%s'", token);
  memcpy(key, token, klen);
  key[klen] = '\0';
  snprintf(value, vcap, "%s", eq + 1);
}

static void parse_int_bounds_value(const char *s, int *lo, int *hi)
{
  char *tmp = xstrdup(s);
  char *colon = strchr(tmp, ':');
  if (!colon) die("Invalid integer bounds '%s' (expected lo:hi)", s);
  *colon = '\0';
  *lo = atoi(tmp);
  *hi = atoi(colon + 1);
  if (*lo > *hi)
    die("Invalid integer bounds '%s' (lo > hi)", s);
  free(tmp);
}

static void parse_ll_bounds_value(const char *s, long long *lo, long long *hi)
{
  char *tmp = xstrdup(s);
  char *colon = strchr(tmp, ':');
  if (!colon) die("Invalid integer bounds '%s' (expected lo:hi)", s);
  *colon = '\0';
  *lo = atoll(tmp);
  *hi = atoll(colon + 1);
  if (*lo > *hi)
    die("Invalid integer bounds '%s' (lo > hi)", s);
  free(tmp);
}

static void parse_double_bounds_value(const char *s, double *lo, double *hi)
{
  char *tmp = xstrdup(s);
  char *colon = strchr(tmp, ':');
  if (!colon) die("Invalid floating bounds '%s' (expected lo:hi)", s);
  *colon = '\0';
  *lo = atof(tmp);
  *hi = atof(colon + 1);
  if (*lo > *hi)
    die("Invalid floating bounds '%s' (lo > hi)", s);
  free(tmp);
}

static void parse_global_bounds(ga_config_t *cfg, const char *s)
{
  char *tmp = xstrdup(s);
  char *save = NULL;
  char *tok = strtok_r(tmp, ",", &save);

  while (tok) {
    char key[64], val[128];
    parse_key_value_pair(tok, key, sizeof(key), val, sizeof(val));
    if (strcmp(key, "hs") == 0) {
      parse_int_bounds_value(val, &cfg->hs_range.lo, &cfg->hs_range.hi);
    } else if (strcmp(key, "lr") == 0) {
      parse_double_bounds_value(val, &cfg->lr_range.lo, &cfg->lr_range.hi);
    } else if (strcmp(key, "seed") == 0) {
      parse_int_bounds_value(val, &cfg->seed_range.lo, &cfg->seed_range.hi);
    } else {
      die("Unknown global bound key '%s'", key);
    }
    tok = strtok_r(NULL, ",", &save);
  }
  free(tmp);
}

static void parse_cm_bounds(ga_config_t *cfg, const char *s)
{
  char *tmp = xstrdup(s);
  char *save = NULL;
  char *tok = strtok_r(tmp, ",", &save);

  while (tok) {
    char key[64], val[128];
    parse_key_value_pair(tok, key, sizeof(key), val, sizeof(val));
    if (strcmp(key, "ctx") == 0) parse_int_bounds_value(val, &cfg->cm_ctx_range.lo, &cfg->cm_ctx_range.hi);
    else if (strcmp(key, "den") == 0) parse_int_bounds_value(val, &cfg->cm_den_range.lo, &cfg->cm_den_range.hi);
    else if (strcmp(key, "ir") == 0) parse_int_bounds_value(val, &cfg->cm_ir_range.lo, &cfg->cm_ir_range.hi);
    else if (strcmp(key, "gamma") == 0) parse_double_bounds_value(val, &cfg->cm_gamma_range.lo, &cfg->cm_gamma_range.hi);
    else if (strcmp(key, "edits") == 0) parse_int_bounds_value(val, &cfg->cm_edits_range.lo, &cfg->cm_edits_range.hi);
    else if (strcmp(key, "eden") == 0) parse_int_bounds_value(val, &cfg->cm_eden_range.lo, &cfg->cm_eden_range.hi);
    else if (strcmp(key, "eir") == 0) parse_int_bounds_value(val, &cfg->cm_eir_range.lo, &cfg->cm_eir_range.hi);
    else if (strcmp(key, "egamma") == 0) parse_double_bounds_value(val, &cfg->cm_egamma_range.lo, &cfg->cm_egamma_range.hi);
    else die("Unknown C-model bound key '%s'", key);
    tok = strtok_r(NULL, ",", &save);
  }
  free(tmp);
}

static void parse_rm_bounds(ga_config_t *cfg, const char *s)
{
  char *tmp = xstrdup(s);
  char *save = NULL;
  char *tok = strtok_r(tmp, ",", &save);

  while (tok) {
    char key[64], val[128];
    parse_key_value_pair(tok, key, sizeof(key), val, sizeof(val));
    if (strcmp(key, "nr") == 0) parse_int_bounds_value(val, &cfg->rm_nr_range.lo, &cfg->rm_nr_range.hi);
    else if (strcmp(key, "ctx") == 0) parse_int_bounds_value(val, &cfg->rm_ctx_range.lo, &cfg->rm_ctx_range.hi);
    else if (strcmp(key, "beta") == 0) parse_double_bounds_value(val, &cfg->rm_beta_range.lo, &cfg->rm_beta_range.hi);
    else if (strcmp(key, "limit") == 0) parse_int_bounds_value(val, &cfg->rm_limit_range.lo, &cfg->rm_limit_range.hi);
    else if (strcmp(key, "gamma") == 0) parse_double_bounds_value(val, &cfg->rm_gamma_range.lo, &cfg->rm_gamma_range.hi);
    else if (strcmp(key, "ir") == 0) parse_int_bounds_value(val, &cfg->rm_ir_range.lo, &cfg->rm_ir_range.hi);
    else if (strcmp(key, "weight") == 0) parse_double_bounds_value(val, &cfg->rm_weight_range.lo, &cfg->rm_weight_range.hi);
    else if (strcmp(key, "cache") == 0) parse_ll_bounds_value(val, &cfg->rm_cache_range.lo, &cfg->rm_cache_range.hi);
    else die("Unknown R-model bound key '%s'", key);
    tok = strtok_r(NULL, ",", &save);
  }
  free(tmp);
}

/* ------------------------------------------------------------------------- */
/* Candidate generation and repair                                           */
/* ------------------------------------------------------------------------- */

static void set_default_config(ga_config_t *cfg)
{
  memset(cfg, 0, sizeof(*cfg));

  cfg->population = 48;
  cfg->generations = 30;
  cfg->threads = 4;
  cfg->elite_count = 4;
  cfg->tournament_size = 3;
  cfg->crossover_rate = 0.85;
  cfg->mutation_rate = 0.12;
  cfg->toggle_rate = 0.08;
  cfg->blend_alpha = 0.35;

  cfg->max_cmodels = 4;
  cfg->max_rmodels = 2;
  cfg->min_cmodels = 0;
  cfg->min_rmodels = 0;

  cfg->verify = 0;
  cfg->keep_temps = 0;
  cfg->quiet = 0;

  cfg->optimize_hs = 1;
  cfg->optimize_lr = 1;
  cfg->optimize_seed = 0;

  cfg->fixed_hs = 42;
  cfg->fixed_lr = 0.03;
  cfg->fixed_seed = 17;

  cfg->hs_range.lo = 8;
  cfg->hs_range.hi = 512;
  cfg->lr_range.lo = 0.0;
  cfg->lr_range.hi = 0.30;
  cfg->seed_range.lo = 1;
  cfg->seed_range.hi = 599999;

  cfg->cm_ctx_range.lo = J3_CM_CTX_MIN;       cfg->cm_ctx_range.hi = J3_CM_CTX_MAX;
  cfg->cm_den_range.lo = J3_CM_DEN_MIN;       cfg->cm_den_range.hi = J3_CM_DEN_MAX;
  cfg->cm_ir_range.lo = J3_CM_IR_MIN;         cfg->cm_ir_range.hi = J3_CM_IR_MAX;
  cfg->cm_gamma_range.lo = 0.0001;            cfg->cm_gamma_range.hi = 0.999999;
  cfg->cm_edits_range.lo = J3_CM_EDITS_MIN;   cfg->cm_edits_range.hi = J3_CM_EDITS_MAX;
  cfg->cm_eden_range.lo = J3_CM_EDEN_MIN;     cfg->cm_eden_range.hi = J3_CM_EDEN_MAX;
  cfg->cm_eir_range.lo = J3_CM_EIR_MIN;       cfg->cm_eir_range.hi = J3_CM_EIR_MAX;
  cfg->cm_egamma_range.lo = 0.0001;           cfg->cm_egamma_range.hi = 0.999999;

  cfg->rm_nr_range.lo = J3_RM_NR_MIN;         cfg->rm_nr_range.hi = J3_RM_NR_MAX;
  cfg->rm_ctx_range.lo = J3_RM_CTX_MIN;       cfg->rm_ctx_range.hi = J3_RM_CTX_MAX;
  cfg->rm_beta_range.lo = 0.0001;             cfg->rm_beta_range.hi = 0.999999;
  cfg->rm_limit_range.lo = J3_RM_LIMIT_MIN;   cfg->rm_limit_range.hi = J3_RM_LIMIT_MAX;
  cfg->rm_gamma_range.lo = 0.0001;            cfg->rm_gamma_range.hi = 0.999999;
  cfg->rm_ir_range.lo = J3_RM_IR_MIN;         cfg->rm_ir_range.hi = J3_RM_IR_MAX;
  cfg->rm_weight_range.lo = 0.0001;           cfg->rm_weight_range.hi = 0.999999;
  cfg->rm_cache_range.lo = J3_RM_CACHE_MIN;   cfg->rm_cache_range.hi = J3_RM_CACHE_MAX;

  cfg->objective = OBJ_BYTES;
  cfg->time_weight = 1.0;

  cfg->master_seed = (uint64_t) time(NULL) ^ 0x5bf03635f123abULL;

  cfg->have_restart = 0;
  cfg->restart_from[0] = '\0';

  snprintf(cfg->workdir, sizeof(cfg->workdir), "/tmp");
}

static void sanitize_config(ga_config_t *cfg)
{
  if (cfg->population < 2) cfg->population = 2;
  if (cfg->generations < 1) cfg->generations = 1;
  if (cfg->threads < 1) cfg->threads = 1;
  if (cfg->elite_count < 1) cfg->elite_count = 1;
  if (cfg->elite_count > cfg->population) cfg->elite_count = cfg->population;
  if (cfg->tournament_size < 2) cfg->tournament_size = 2;
  if (cfg->tournament_size > cfg->population) cfg->tournament_size = cfg->population;

  if (cfg->max_cmodels < 0) cfg->max_cmodels = 0;
  if (cfg->max_rmodels < 0) cfg->max_rmodels = 0;
  if (cfg->max_cmodels > HARD_MAX_CMODELS) die("--max-cmodels exceeds hard limit %d", HARD_MAX_CMODELS);
  if (cfg->max_rmodels > HARD_MAX_RMODELS) die("--max-rmodels exceeds hard limit %d", HARD_MAX_RMODELS);

  if (cfg->min_cmodels < 0) cfg->min_cmodels = 0;
  if (cfg->min_rmodels < 0) cfg->min_rmodels = 0;
  if (cfg->min_cmodels > cfg->max_cmodels) cfg->min_cmodels = cfg->max_cmodels;
  if (cfg->min_rmodels > cfg->max_rmodels) cfg->min_rmodels = cfg->max_rmodels;

  if (cfg->max_cmodels == 0 && cfg->max_rmodels == 0)
    die("At least one model slot must be available (C or R)");

  cfg->crossover_rate = clampd(cfg->crossover_rate, 0.0, 1.0);
  cfg->mutation_rate = clampd(cfg->mutation_rate, 0.0, 1.0);
  cfg->toggle_rate = clampd(cfg->toggle_rate, 0.0, 1.0);
  cfg->blend_alpha = clampd(cfg->blend_alpha, 0.0, 2.0);

  clamp_irange(&cfg->hs_range, 1, 999999);
  clamp_drange(&cfg->lr_range, 0.0, 1.0);
  clamp_irange(&cfg->seed_range, 1, 599999);

  clamp_irange(&cfg->cm_ctx_range, J3_CM_CTX_MIN, J3_CM_CTX_MAX);
  clamp_irange(&cfg->cm_den_range, J3_CM_DEN_MIN, J3_CM_DEN_MAX);
  clamp_irange(&cfg->cm_ir_range, J3_CM_IR_MIN, J3_CM_IR_MAX);
  clamp_drange(&cfg->cm_gamma_range, 0.001, 0.999999);
  clamp_irange(&cfg->cm_edits_range, J3_CM_EDITS_MIN, J3_CM_EDITS_MAX);
  clamp_irange(&cfg->cm_eden_range, J3_CM_EDEN_MIN, J3_CM_EDEN_MAX);
  clamp_irange(&cfg->cm_eir_range, J3_CM_EIR_MIN, J3_CM_EIR_MAX);
  clamp_drange(&cfg->cm_egamma_range, 0.001, 0.999999);

  clamp_irange(&cfg->rm_nr_range, J3_RM_NR_MIN, J3_RM_NR_MAX);
  clamp_irange(&cfg->rm_ctx_range, J3_RM_CTX_MIN, J3_RM_CTX_MAX);
  clamp_drange(&cfg->rm_beta_range, 0.001, 0.999999);
  clamp_irange(&cfg->rm_limit_range, J3_RM_LIMIT_MIN, J3_RM_LIMIT_MAX);
  clamp_drange(&cfg->rm_gamma_range, 0.001, 0.999999);
  clamp_irange(&cfg->rm_ir_range, J3_RM_IR_MIN, J3_RM_IR_MAX);
  clamp_drange(&cfg->rm_weight_range, 0.001, 0.999999);
  if (cfg->rm_cache_range.lo < J3_RM_CACHE_MIN) cfg->rm_cache_range.lo = J3_RM_CACHE_MIN;
  if (cfg->rm_cache_range.hi > J3_RM_CACHE_MAX) cfg->rm_cache_range.hi = J3_RM_CACHE_MAX;
  if (cfg->rm_cache_range.hi < cfg->rm_cache_range.lo) cfg->rm_cache_range.hi = cfg->rm_cache_range.lo;

  cfg->fixed_hs = clampi(cfg->fixed_hs, cfg->hs_range.lo, cfg->hs_range.hi);
  cfg->fixed_lr = clampd(cfg->fixed_lr, cfg->lr_range.lo, cfg->lr_range.hi);
  cfg->fixed_lr = quantize_lr(cfg->fixed_lr);
  cfg->fixed_seed = clampi(cfg->fixed_seed, cfg->seed_range.lo, cfg->seed_range.hi);

  if (ensure_dir(cfg->workdir) != 0)
    die("Cannot access/create workdir '%s': %s", cfg->workdir, strerror(errno));
}

static void randomize_cgene(const ga_config_t *cfg, rng_t *rng, cmodel_gene_t *g, int can_enable)
{
  g->enabled = can_enable ? rng_bool(rng, 0.55) : 0;
  g->ctx = rng_int_closed(rng, cfg->cm_ctx_range.lo, cfg->cm_ctx_range.hi);
  g->den = rng_int_closed(rng, cfg->cm_den_range.lo, cfg->cm_den_range.hi);
  g->ir = rng_int_closed(rng, cfg->cm_ir_range.lo, cfg->cm_ir_range.hi);
  g->gamma = rng_double_closed(rng, cfg->cm_gamma_range.lo, cfg->cm_gamma_range.hi);
  g->edits = rng_int_closed(rng, cfg->cm_edits_range.lo, cfg->cm_edits_range.hi);
  g->eden = rng_int_closed(rng, cfg->cm_eden_range.lo, cfg->cm_eden_range.hi);
  g->eir = rng_int_closed(rng, cfg->cm_eir_range.lo, cfg->cm_eir_range.hi);
  g->egamma = rng_double_closed(rng, cfg->cm_egamma_range.lo, cfg->cm_egamma_range.hi);
}

static void randomize_rgene(const ga_config_t *cfg, rng_t *rng, rmodel_gene_t *g, int can_enable)
{
  g->enabled = can_enable ? rng_bool(rng, 0.70) : 0;
  g->nr = rng_int_closed(rng, cfg->rm_nr_range.lo, cfg->rm_nr_range.hi);
  g->ctx = rng_int_closed(rng, cfg->rm_ctx_range.lo, cfg->rm_ctx_range.hi);
  g->beta = rng_double_closed(rng, cfg->rm_beta_range.lo, cfg->rm_beta_range.hi);
  g->limit = rng_int_closed(rng, cfg->rm_limit_range.lo, cfg->rm_limit_range.hi);
  g->gamma = rng_double_closed(rng, cfg->rm_gamma_range.lo, cfg->rm_gamma_range.hi);
  g->ir = rng_int_closed(rng, cfg->rm_ir_range.lo, cfg->rm_ir_range.hi);
  g->weight = rng_double_closed(rng, cfg->rm_weight_range.lo, cfg->rm_weight_range.hi);
  g->cache = rng_ll_closed(rng, cfg->rm_cache_range.lo, cfg->rm_cache_range.hi);
}

static void count_active_models(const ga_config_t *cfg, individual_t *ind)
{
  int i;
  ind->n_active_c = 0;
  ind->n_active_r = 0;
  for (i = 0; i < cfg->max_cmodels; ++i)
    ind->n_active_c += ind->c[i].enabled ? 1 : 0;
  for (i = 0; i < cfg->max_rmodels; ++i)
    ind->n_active_r += ind->r[i].enabled ? 1 : 0;
}

static void repair_individual(const ga_config_t *cfg, rng_t *rng, individual_t *ind)
{
  int i;

  if (cfg->optimize_hs) ind->hs = clampi(ind->hs, cfg->hs_range.lo, cfg->hs_range.hi);
  else ind->hs = cfg->fixed_hs;

  if (cfg->optimize_lr) ind->lr = clampd(ind->lr, cfg->lr_range.lo, cfg->lr_range.hi);
  else ind->lr = cfg->fixed_lr;
  ind->lr = quantize_lr(ind->lr);

  if (cfg->optimize_seed) ind->seed = clampi(ind->seed, cfg->seed_range.lo, cfg->seed_range.hi);
  else ind->seed = cfg->fixed_seed;

  for (i = 0; i < HARD_MAX_CMODELS; ++i) {
    cmodel_gene_t *g = &ind->c[i];
    if (i >= cfg->max_cmodels) g->enabled = 0;
    g->ctx = clampi(g->ctx, cfg->cm_ctx_range.lo, cfg->cm_ctx_range.hi);
    g->den = clampi(g->den, cfg->cm_den_range.lo, cfg->cm_den_range.hi);
    g->ir = clampi(g->ir, cfg->cm_ir_range.lo, cfg->cm_ir_range.hi);
    g->gamma = clampd(g->gamma, cfg->cm_gamma_range.lo, cfg->cm_gamma_range.hi);
    g->edits = clampi(g->edits, cfg->cm_edits_range.lo, cfg->cm_edits_range.hi);
    g->eden = clampi(g->eden, cfg->cm_eden_range.lo, cfg->cm_eden_range.hi);
    g->eir = clampi(g->eir, cfg->cm_eir_range.lo, cfg->cm_eir_range.hi);
    g->egamma = clampd(g->egamma, cfg->cm_egamma_range.lo, cfg->cm_egamma_range.hi);
    g->enabled = g->enabled ? 1 : 0;
  }

  for (i = 0; i < HARD_MAX_RMODELS; ++i) {
    rmodel_gene_t *g = &ind->r[i];
    if (i >= cfg->max_rmodels) g->enabled = 0;
    g->nr = clampi(g->nr, cfg->rm_nr_range.lo, cfg->rm_nr_range.hi);
    g->ctx = clampi(g->ctx, cfg->rm_ctx_range.lo, cfg->rm_ctx_range.hi);
    g->beta = clampd(g->beta, cfg->rm_beta_range.lo, cfg->rm_beta_range.hi);
    g->limit = clampi(g->limit, cfg->rm_limit_range.lo, cfg->rm_limit_range.hi);
    g->gamma = clampd(g->gamma, cfg->rm_gamma_range.lo, cfg->rm_gamma_range.hi);
    g->ir = clampi(g->ir, cfg->rm_ir_range.lo, cfg->rm_ir_range.hi);
    g->weight = clampd(g->weight, cfg->rm_weight_range.lo, cfg->rm_weight_range.hi);
    g->cache = clampll(g->cache, cfg->rm_cache_range.lo, cfg->rm_cache_range.hi);
    g->enabled = g->enabled ? 1 : 0;
  }

  count_active_models(cfg, ind);

  while (ind->n_active_c < cfg->min_cmodels && cfg->max_cmodels > 0) {
    int idx = rng_int_closed(rng, 0, cfg->max_cmodels - 1);
    if (!ind->c[idx].enabled) {
      ind->c[idx].enabled = 1;
      ind->n_active_c++;
    }
  }

  while (ind->n_active_r < cfg->min_rmodels && cfg->max_rmodels > 0) {
    int idx = rng_int_closed(rng, 0, cfg->max_rmodels - 1);
    if (!ind->r[idx].enabled) {
      ind->r[idx].enabled = 1;
      ind->n_active_r++;
    }
  }

  if (ind->n_active_c + ind->n_active_r == 0) {
    if (cfg->max_rmodels > 0) {
      ind->r[rng_int_closed(rng, 0, cfg->max_rmodels - 1)].enabled = 1;
    } else if (cfg->max_cmodels > 0) {
      ind->c[rng_int_closed(rng, 0, cfg->max_cmodels - 1)].enabled = 1;
    }
    count_active_models(cfg, ind);
  }
}

static void random_individual(const ga_config_t *cfg, rng_t *rng, individual_t *ind)
{
  int i;
  memset(ind, 0, sizeof(*ind));

  ind->hs = cfg->optimize_hs ? rng_int_closed(rng, cfg->hs_range.lo, cfg->hs_range.hi) : cfg->fixed_hs;
  ind->lr = cfg->optimize_lr ? rng_double_closed(rng, cfg->lr_range.lo, cfg->lr_range.hi) : cfg->fixed_lr;
  ind->seed = cfg->optimize_seed ? rng_int_closed(rng, cfg->seed_range.lo, cfg->seed_range.hi) : cfg->fixed_seed;

  for (i = 0; i < HARD_MAX_CMODELS; ++i)
    randomize_cgene(cfg, rng, &ind->c[i], i < cfg->max_cmodels);

  for (i = 0; i < HARD_MAX_RMODELS; ++i)
    randomize_rgene(cfg, rng, &ind->r[i], i < cfg->max_rmodels);

  repair_individual(cfg, rng, ind);

  ind->fitness = INVALID_FITNESS;
  ind->comp_bytes = 0;
  ind->elapsed_s = 0.0;
  ind->valid = 0;
}

static void mutate_individual(const ga_config_t *cfg, rng_t *rng, individual_t *ind)
{
  int i;
  if (cfg->optimize_hs && rng_bool(rng, cfg->mutation_rate))
    ind->hs = mutate_int_local(rng, ind->hs, cfg->hs_range);

  if (cfg->optimize_lr && rng_bool(rng, cfg->mutation_rate))
    ind->lr = mutate_double_local(rng, ind->lr, cfg->lr_range);

  if (cfg->optimize_seed && rng_bool(rng, cfg->mutation_rate))
    ind->seed = mutate_int_local(rng, ind->seed, cfg->seed_range);

  for (i = 0; i < cfg->max_cmodels; ++i) {
    cmodel_gene_t *g = &ind->c[i];
    if (rng_bool(rng, cfg->toggle_rate))
      g->enabled = !g->enabled;
    if (rng_bool(rng, cfg->mutation_rate)) g->ctx = mutate_int_local(rng, g->ctx, cfg->cm_ctx_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->den = mutate_int_local(rng, g->den, cfg->cm_den_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->ir = mutate_int_local(rng, g->ir, cfg->cm_ir_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->gamma = mutate_double_local(rng, g->gamma, cfg->cm_gamma_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->edits = mutate_int_local(rng, g->edits, cfg->cm_edits_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->eden = mutate_int_local(rng, g->eden, cfg->cm_eden_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->eir = mutate_int_local(rng, g->eir, cfg->cm_eir_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->egamma = mutate_double_local(rng, g->egamma, cfg->cm_egamma_range);
  }

  for (i = 0; i < cfg->max_rmodels; ++i) {
    rmodel_gene_t *g = &ind->r[i];
    if (rng_bool(rng, cfg->toggle_rate))
      g->enabled = !g->enabled;
    if (rng_bool(rng, cfg->mutation_rate)) g->nr = mutate_int_local(rng, g->nr, cfg->rm_nr_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->ctx = mutate_int_local(rng, g->ctx, cfg->rm_ctx_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->beta = mutate_double_local(rng, g->beta, cfg->rm_beta_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->limit = mutate_int_local(rng, g->limit, cfg->rm_limit_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->gamma = mutate_double_local(rng, g->gamma, cfg->rm_gamma_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->ir = mutate_int_local(rng, g->ir, cfg->rm_ir_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->weight = mutate_double_local(rng, g->weight, cfg->rm_weight_range);
    if (rng_bool(rng, cfg->mutation_rate)) g->cache = mutate_ll_local(rng, g->cache, cfg->rm_cache_range);
  }

  repair_individual(cfg, rng, ind);
  ind->fitness = INVALID_FITNESS;
  ind->valid = 0;
}

static void crossover_individual(const ga_config_t *cfg, rng_t *rng,
                                 const individual_t *a, const individual_t *b,
                                 individual_t *child)
{
  int i;
  memset(child, 0, sizeof(*child));

  child->hs = cfg->optimize_hs ? cross_int(rng, a->hs, b->hs, cfg->hs_range) : cfg->fixed_hs;
  child->lr = cfg->optimize_lr ? cross_double(rng, a->lr, b->lr, cfg->lr_range, cfg->blend_alpha) : cfg->fixed_lr;
  child->seed = cfg->optimize_seed ? cross_int(rng, a->seed, b->seed, cfg->seed_range) : cfg->fixed_seed;

  for (i = 0; i < cfg->max_cmodels; ++i) {
    child->c[i].enabled = rng_bool(rng, 0.5) ? a->c[i].enabled : b->c[i].enabled;
    child->c[i].ctx = cross_int(rng, a->c[i].ctx, b->c[i].ctx, cfg->cm_ctx_range);
    child->c[i].den = cross_int(rng, a->c[i].den, b->c[i].den, cfg->cm_den_range);
    child->c[i].ir = cross_int(rng, a->c[i].ir, b->c[i].ir, cfg->cm_ir_range);
    child->c[i].gamma = cross_double(rng, a->c[i].gamma, b->c[i].gamma, cfg->cm_gamma_range, cfg->blend_alpha);
    child->c[i].edits = cross_int(rng, a->c[i].edits, b->c[i].edits, cfg->cm_edits_range);
    child->c[i].eden = cross_int(rng, a->c[i].eden, b->c[i].eden, cfg->cm_eden_range);
    child->c[i].eir = cross_int(rng, a->c[i].eir, b->c[i].eir, cfg->cm_eir_range);
    child->c[i].egamma = cross_double(rng, a->c[i].egamma, b->c[i].egamma, cfg->cm_egamma_range, cfg->blend_alpha);
  }

  for (i = 0; i < cfg->max_rmodels; ++i) {
    child->r[i].enabled = rng_bool(rng, 0.5) ? a->r[i].enabled : b->r[i].enabled;
    child->r[i].nr = cross_int(rng, a->r[i].nr, b->r[i].nr, cfg->rm_nr_range);
    child->r[i].ctx = cross_int(rng, a->r[i].ctx, b->r[i].ctx, cfg->rm_ctx_range);
    child->r[i].beta = cross_double(rng, a->r[i].beta, b->r[i].beta, cfg->rm_beta_range, cfg->blend_alpha);
    child->r[i].limit = cross_int(rng, a->r[i].limit, b->r[i].limit, cfg->rm_limit_range);
    child->r[i].gamma = cross_double(rng, a->r[i].gamma, b->r[i].gamma, cfg->rm_gamma_range, cfg->blend_alpha);
    child->r[i].ir = cross_int(rng, a->r[i].ir, b->r[i].ir, cfg->rm_ir_range);
    child->r[i].weight = cross_double(rng, a->r[i].weight, b->r[i].weight, cfg->rm_weight_range, cfg->blend_alpha);
    child->r[i].cache = cross_ll(rng, a->r[i].cache, b->r[i].cache, cfg->rm_cache_range);
  }

  repair_individual(cfg, rng, child);
  mutate_individual(cfg, rng, child);
}

static const individual_t *tournament_select(const individual_t *pop, int n, int tsize, rng_t *rng)
{
  int i;
  int best = rng_int_closed(rng, 0, n - 1);
  for (i = 1; i < tsize; ++i) {
    int c = rng_int_closed(rng, 0, n - 1);
    if (pop[c].fitness < pop[best].fitness)
      best = c;
  }
  return &pop[best];
}

/* ------------------------------------------------------------------------- */
/* Restart-from parsing                                                      */
/* ------------------------------------------------------------------------- */

static int parse_cmodel_spec(const char *s, cmodel_gene_t *g)
{
  int n = sscanf(s, "%d:%d:%d:%lf/%d:%d:%d:%lf",
                 &g->ctx, &g->den, &g->ir, &g->gamma,
                 &g->edits, &g->eden, &g->eir, &g->egamma);
  if (n != 8) return 0;
  g->enabled = 1;
  return 1;
}

static int parse_rmodel_spec(const char *s, rmodel_gene_t *g)
{
  int n;
  memset(g, 0, sizeof(*g));

  n = sscanf(s, "%d:%d:%lf:%d:%lf:%d:%lf:%lld",
             &g->nr, &g->ctx, &g->beta, &g->limit,
             &g->gamma, &g->ir, &g->weight, &g->cache);

  if (n != 8) return 0;
  g->enabled = 1;
  return 1;
}

static void init_restart_base(const ga_config_t *cfg, individual_t *ind)
{
  int i;
  memset(ind, 0, sizeof(*ind));

  ind->hs = cfg->fixed_hs;
  ind->lr = cfg->fixed_lr;
  ind->seed = cfg->fixed_seed;

  for (i = 0; i < HARD_MAX_CMODELS; ++i) {
    ind->c[i].enabled = 0;
    ind->c[i].ctx = cfg->cm_ctx_range.lo;
    ind->c[i].den = cfg->cm_den_range.lo;
    ind->c[i].ir = cfg->cm_ir_range.lo;
    ind->c[i].gamma = cfg->cm_gamma_range.lo;
    ind->c[i].edits = cfg->cm_edits_range.lo;
    ind->c[i].eden = cfg->cm_eden_range.lo;
    ind->c[i].eir = cfg->cm_eir_range.lo;
    ind->c[i].egamma = cfg->cm_egamma_range.lo;
  }

  for (i = 0; i < HARD_MAX_RMODELS; ++i) {
    ind->r[i].enabled = 0;
    ind->r[i].nr = cfg->rm_nr_range.lo;
    ind->r[i].ctx = cfg->rm_ctx_range.lo;
    ind->r[i].beta = cfg->rm_beta_range.lo;
    ind->r[i].limit = cfg->rm_limit_range.lo;
    ind->r[i].gamma = cfg->rm_gamma_range.lo;
    ind->r[i].ir = cfg->rm_ir_range.lo;
    ind->r[i].weight = cfg->rm_weight_range.lo;
    ind->r[i].cache = cfg->rm_cache_range.lo;
  }

  ind->fitness = INVALID_FITNESS;
  ind->comp_bytes = 0;
  ind->elapsed_s = 0.0;
  ind->valid = 0;
}

static void parse_restart_individual(const ga_config_t *cfg, const char *s, individual_t *ind)
{
  char *tmp = xstrdup(s);
  char *save = NULL;
  char *tok = NULL;
  int cslot = 0;
  int rslot = 0;
  rng_t rr;

  init_restart_base(cfg, ind);

  tok = strtok_r(tmp, " \t\r\n", &save);
  while (tok) {
    if (strcmp(tok, "-hs") == 0) {
      tok = strtok_r(NULL, " \t\r\n", &save);
      if (!tok) die("Missing value after -hs in --restart-from");
      ind->hs = atoi(tok);
    } else if (strcmp(tok, "-lr") == 0) {
      tok = strtok_r(NULL, " \t\r\n", &save);
      if (!tok) die("Missing value after -lr in --restart-from");
      ind->lr = atof(tok);
    } else if (strcmp(tok, "-sd") == 0) {
      tok = strtok_r(NULL, " \t\r\n", &save);
      if (!tok) die("Missing value after -sd in --restart-from");
      ind->seed = atoi(tok);
    } else if (strcmp(tok, "-cm") == 0) {
      tok = strtok_r(NULL, " \t\r\n", &save);
      if (!tok) die("Missing value after -cm in --restart-from");
      if (cslot >= cfg->max_cmodels)
        die("Too many -cm models in --restart-from (increase --max-cmodels)");
      if (!parse_cmodel_spec(tok, &ind->c[cslot]))
        die("Invalid -cm spec in --restart-from: '%s'", tok);
      cslot++;
    } else if (strcmp(tok, "-rm") == 0) {
      tok = strtok_r(NULL, " \t\r\n", &save);
      if (!tok) die("Missing value after -rm in --restart-from");
      if (rslot >= cfg->max_rmodels)
        die("Too many -rm models in --restart-from (increase --max-rmodels)");
      if (!parse_rmodel_spec(tok, &ind->r[rslot]))
        die("Invalid -rm spec in --restart-from: '%s'", tok);
      rslot++;
    } else if (tok[0] == '-') {
      die("Unknown option '%s' in --restart-from", tok);
    }

    tok = strtok_r(NULL, " \t\r\n", &save);
  }

  rng_seed(&rr, cfg->master_seed ^ 0x9e3779b97f4a7c15ULL);
  repair_individual(cfg, &rr, ind);
  build_param_string(cfg, ind);

  free(tmp);
}

/* ------------------------------------------------------------------------- */
/* Process execution                                                         */
/* ------------------------------------------------------------------------- */

static int run_command(char *const argv[],
                       const char *stdout_path,
                       const char *stderr_path,
                       double *elapsed_s)
{
  pid_t pid;
  int status;
  double t0 = now_sec();

  pid = fork();
  if (pid < 0)
    return -1;

  if (pid == 0) {
    int fdout = open(stdout_path, O_CREAT | O_WRONLY | O_TRUNC, 0644);
    int fderr = open(stderr_path, O_CREAT | O_WRONLY | O_TRUNC, 0644);
    if (fdout >= 0) {
      dup2(fdout, STDOUT_FILENO);
      close(fdout);
    }
    if (fderr >= 0) {
      dup2(fderr, STDERR_FILENO);
      close(fderr);
    }
    execv(argv[0], argv);
    _exit(127);
  }

  if (waitpid(pid, &status, 0) < 0)
    return -1;

  if (elapsed_s)
    *elapsed_s = now_sec() - t0;

  if (WIFEXITED(status) && WEXITSTATUS(status) == 0)
    return 0;
  return -1;
}

static int build_compress_argv(const ga_config_t *cfg, const individual_t *ind,
                               const char *out_path, char ***argv_out)
{
  char **argv = (char **) xcalloc(MAX_ARGV_ITEMS, sizeof(char *));
  char buf_hs[64], buf_lr[64], buf_seed[64];
  int argc = 0;
  int i;

  argv[argc++] = xstrdup(cfg->jarvis_path);
  argv[argc++] = xstrdup("-o");
  argv[argc++] = xstrdup(out_path);

  snprintf(buf_hs, sizeof(buf_hs), "%d", ind->hs);
  snprintf(buf_lr, sizeof(buf_lr), "%.3f", quantize_lr(ind->lr));
  snprintf(buf_seed, sizeof(buf_seed), "%d", ind->seed);

  argv[argc++] = xstrdup("-hs");
  argv[argc++] = xstrdup(buf_hs);
  argv[argc++] = xstrdup("-lr");
  argv[argc++] = xstrdup(buf_lr);
  argv[argc++] = xstrdup("-sd");
  argv[argc++] = xstrdup(buf_seed);

  for (i = 0; i < cfg->max_cmodels; ++i) {
    char cm[256];
    const cmodel_gene_t *g = &ind->c[i];
    if (!g->enabled) continue;
    if (format_cmodel_spec(cm, sizeof(cm), g) != 0)
      die("C-model string too long");
    argv[argc++] = xstrdup("-cm");
    argv[argc++] = xstrdup(cm);
    if (argc >= MAX_ARGV_ITEMS - 4)
      die("Too many command-line items for compression argv");
  }

  for (i = 0; i < cfg->max_rmodels; ++i) {
    char rm[256];
    const rmodel_gene_t *g = &ind->r[i];
    if (!g->enabled) continue;
    if (format_rmodel_spec(rm, sizeof(rm), g) != 0)
      die("R-model string too long");
    argv[argc++] = xstrdup("-rm");
    argv[argc++] = xstrdup(rm);
    if (argc >= MAX_ARGV_ITEMS - 4)
      die("Too many command-line items for compression argv");
  }

  argv[argc++] = xstrdup(cfg->input_path);
  argv[argc] = NULL;

  *argv_out = argv;
  return argc;
}

static int build_decompress_argv(const ga_config_t *cfg, const char *compressed_path,
                                 const char *out_path, char ***argv_out)
{
  char **argv = (char **) xcalloc(8, sizeof(char *));
  int argc = 0;
  argv[argc++] = xstrdup(cfg->jarvis_path);
  argv[argc++] = xstrdup("-d");
  argv[argc++] = xstrdup("-o");
  argv[argc++] = xstrdup(out_path);
  argv[argc++] = xstrdup(compressed_path);
  argv[argc] = NULL;
  *argv_out = argv;
  return argc;
}

static void free_argv(char **argv)
{
  int i;
  if (!argv) return;
  for (i = 0; argv[i]; ++i)
    free(argv[i]);
  free(argv);
}

static eval_result_t compute_fitness(const ga_config_t *cfg, size_t comp_bytes, double elapsed_s, int valid)
{
  eval_result_t r;
  r.comp_bytes = comp_bytes;
  r.elapsed_s = elapsed_s;
  r.valid = valid;
  if (!valid) {
    r.fitness = INVALID_FITNESS;
    return r;
  }

  switch (cfg->objective) {
    case OBJ_BYTES:
      r.fitness = (double) comp_bytes;
      break;
    case OBJ_BPS:
      if (cfg->input_bytes == 0) r.fitness = INVALID_FITNESS;
      else r.fitness = (double) comp_bytes * 8.0 / (double) cfg->input_bytes;
      break;
    case OBJ_BYTES_PLUS_TIME:
      r.fitness = (double) comp_bytes + cfg->time_weight * elapsed_s;
      break;
    default:
      r.fitness = (double) comp_bytes;
      break;
  }
  return r;
}

static eval_result_t evaluate_candidate(const ga_config_t *cfg, individual_t *ind)
{
  char template_dir[PATH_MAX];
  const char suffix[] = "/j3ga_XXXXXX";
  char run_dir[PATH_MAX];
  char out_path[PATH_MAX], dec_path[PATH_MAX];
  char comp_out_log[PATH_MAX], comp_err_log[PATH_MAX];
  char dec_out_log[PATH_MAX], dec_err_log[PATH_MAX];
  char **argv = NULL;
  char **dargv = NULL;
  double elapsed = 0.0;
  eval_result_t res;
  size_t wlen;

  build_param_string(cfg, ind);

  if (!validate_individual_for_jarvis(cfg, ind)) {
    if (!cfg->quiet)
      fprintf(stderr, "[warn] candidate rejected before JARVIS: %s\n", ind->param_string);
    return compute_fitness(cfg, 0, 0.0, 0);
  }

  wlen = strlen(cfg->workdir);
  if (wlen + sizeof(suffix) > sizeof(template_dir)) {
    res = compute_fitness(cfg, 0, 0.0, 0);
    return res;
  }

  memcpy(template_dir, cfg->workdir, wlen);
  memcpy(template_dir + wlen, suffix, sizeof(suffix));

  if (!mkdtemp(template_dir)) {
    res = compute_fitness(cfg, 0, 0.0, 0);
    return res;
  }
  snprintf(run_dir, sizeof(run_dir), "%s", template_dir);

  path_join(out_path, sizeof(out_path), run_dir, "cand.jc");
  path_join(dec_path, sizeof(dec_path), run_dir, "cand.dec");
  path_join(comp_out_log, sizeof(comp_out_log), run_dir, "compress.stdout");
  path_join(comp_err_log, sizeof(comp_err_log), run_dir, "compress.stderr");
  path_join(dec_out_log, sizeof(dec_out_log), run_dir, "decompress.stdout");
  path_join(dec_err_log, sizeof(dec_err_log), run_dir, "decompress.stderr");

  build_compress_argv(cfg, ind, out_path, &argv);

  if (run_command(argv, comp_out_log, comp_err_log, &elapsed) != 0 || access(out_path, R_OK) != 0) {
    if (!cfg->quiet) {
      fprintf(stderr, "[warn] compression failed for: %s\n", ind->param_string);
      fprintf(stderr, "       logs: %s  %s\n", comp_out_log, comp_err_log);
    }
    free_argv(argv);
    if (!cfg->keep_temps) {
      remove_if_exists(comp_out_log);
      remove_if_exists(comp_err_log);
      rmdir_if_exists(run_dir);
    }
    return compute_fitness(cfg, 0, elapsed, 0);
  }

  free_argv(argv);
  argv = NULL;

  if (cfg->verify) {
    build_decompress_argv(cfg, out_path, dec_path, &dargv);
    if (run_command(dargv, dec_out_log, dec_err_log, NULL) != 0 ||
        access(dec_path, R_OK) != 0 ||
        !files_equal(cfg->input_path, dec_path)) {
      if (!cfg->quiet) {
        fprintf(stderr, "[warn] verification failed for: %s\n", ind->param_string);
        fprintf(stderr, "       logs: %s  %s\n", dec_out_log, dec_err_log);
      }
      free_argv(dargv);
      if (!cfg->keep_temps) {
        remove_if_exists(out_path);
        remove_if_exists(dec_path);
        remove_if_exists(comp_out_log);
        remove_if_exists(comp_err_log);
        remove_if_exists(dec_out_log);
        remove_if_exists(dec_err_log);
        rmdir_if_exists(run_dir);
      }
      return compute_fitness(cfg, 0, elapsed, 0);
    }
    free_argv(dargv);
  }

  res = compute_fitness(cfg, get_file_size_or_die(out_path), elapsed, 1);

  if (!cfg->keep_temps) {
    remove_if_exists(out_path);
    remove_if_exists(dec_path);
    remove_if_exists(comp_out_log);
    remove_if_exists(comp_err_log);
    remove_if_exists(dec_out_log);
    remove_if_exists(dec_err_log);
    rmdir_if_exists(run_dir);
  }

  return res;
}

/* ------------------------------------------------------------------------- */
/* Parallel evaluation                                                       */
/* ------------------------------------------------------------------------- */

static void *worker_main(void *arg)
{
  worker_ctx_t *ctx = (worker_ctx_t *) arg;

  while (1) {
    int job_slot;
    int pop_idx;
    eval_result_t res;

    pthread_mutex_lock(&ctx->next_mutex);
    job_slot = ctx->next_job++;
    pthread_mutex_unlock(&ctx->next_mutex);

    if (job_slot >= ctx->job_count)
      break;

    pop_idx = ctx->job_indices[job_slot];
    res = evaluate_candidate(ctx->cfg, &ctx->pop[pop_idx]);

    ctx->pop[pop_idx].fitness = res.fitness;
    ctx->pop[pop_idx].comp_bytes = res.comp_bytes;
    ctx->pop[pop_idx].elapsed_s = res.elapsed_s;
    ctx->pop[pop_idx].valid = res.valid;

    cache_insert(ctx->cache, ctx->pop[pop_idx].key, res);
  }

  return NULL;
}

static void evaluate_population(const ga_config_t *cfg, individual_t *pop, int n, fitness_cache_t *cache)
{
  int i, j;
  int *jobs = (int *) xcalloc((size_t) n, sizeof(int));
  int *dup_of = (int *) xcalloc((size_t) n, sizeof(int));
  int job_count = 0;
  worker_ctx_t ctx;
  pthread_t *threads;
  int thread_count;

  for (i = 0; i < n; ++i)
    dup_of[i] = -1;

  for (i = 0; i < n; ++i) {
    eval_result_t cached;
    build_param_string(cfg, &pop[i]);
    if (cache_lookup(cache, pop[i].key, &cached)) {
      pop[i].fitness = cached.fitness;
      pop[i].comp_bytes = cached.comp_bytes;
      pop[i].elapsed_s = cached.elapsed_s;
      pop[i].valid = cached.valid;
      continue;
    }

    for (j = 0; j < i; ++j) {
      if (strcmp(pop[i].key, pop[j].key) == 0) {
        dup_of[i] = j;
        break;
      }
    }

    if (dup_of[i] < 0)
      jobs[job_count++] = i;
  }

  memset(&ctx, 0, sizeof(ctx));
  ctx.cfg = cfg;
  ctx.pop = pop;
  ctx.job_indices = jobs;
  ctx.job_count = job_count;
  ctx.next_job = 0;
  ctx.cache = cache;
  pthread_mutex_init(&ctx.next_mutex, NULL);

  thread_count = cfg->threads;
  if (thread_count > job_count) thread_count = job_count;
  if (thread_count < 1) thread_count = 1;

  threads = (pthread_t *) xcalloc((size_t) thread_count, sizeof(pthread_t));
  for (i = 0; i < thread_count; ++i) {
    if (pthread_create(&threads[i], NULL, worker_main, &ctx) != 0)
      die("pthread_create failed");
  }
  for (i = 0; i < thread_count; ++i)
    pthread_join(threads[i], NULL);

  for (i = 0; i < n; ++i) {
    if (dup_of[i] >= 0) {
      pop[i].fitness = pop[dup_of[i]].fitness;
      pop[i].comp_bytes = pop[dup_of[i]].comp_bytes;
      pop[i].elapsed_s = pop[dup_of[i]].elapsed_s;
      pop[i].valid = pop[dup_of[i]].valid;
    }
  }

  pthread_mutex_destroy(&ctx.next_mutex);
  free(threads);
  free(jobs);
  free(dup_of);
}

/* ------------------------------------------------------------------------- */
/* Sorting and reporting                                                     */
/* ------------------------------------------------------------------------- */

static int cmp_individual(const void *a, const void *b)
{
  const individual_t *ia = (const individual_t *) a;
  const individual_t *ib = (const individual_t *) b;

  if (ia->fitness < ib->fitness) return -1;
  if (ia->fitness > ib->fitness) return 1;
  if (ia->comp_bytes < ib->comp_bytes) return -1;
  if (ia->comp_bytes > ib->comp_bytes) return 1;
  if (ia->elapsed_s < ib->elapsed_s) return -1;
  if (ia->elapsed_s > ib->elapsed_s) return 1;
  return strcmp(ia->key, ib->key);
}

static const char *objective_name(objective_t obj)
{
  switch (obj) {
    case OBJ_BYTES: return "bytes";
    case OBJ_BPS: return "bps";
    case OBJ_BYTES_PLUS_TIME: return "bytes+time";
    default: return "bytes";
  }
}

static void print_individual(FILE *fp, const ga_config_t *cfg, const individual_t *ind, const char *prefix)
{
  fprintf(fp, "%sfitness=%.2f bytes=%zu bps=%.4f elapsed=%.2f valid=%d objective=%s\n",
          prefix, ind->fitness, ind->comp_bytes,
          compute_bps(ind->comp_bytes, cfg->input_bytes),
          ind->elapsed_s, ind->valid,
          objective_name(cfg->objective));
  fprintf(fp, "%sparams: %s\n", prefix, ind->param_string);
  fprintf(fp, "%scommand: %s %s %s\n",
          prefix, cfg->jarvis_path, ind->param_string, cfg->input_path);
}

static void write_best_file(const ga_config_t *cfg, const individual_t *best)
{
  FILE *f;
  if (cfg->best_out[0] == '\0')
    return;
  f = fopen(cfg->best_out, "w");
  if (!f) {
    fprintf(stderr, "Warning: cannot open '%s' for writing: %s\n",
            cfg->best_out, strerror(errno));
    return;
  }
  print_individual(f, cfg, best, "");
  fclose(f);
}

static void append_history(const ga_config_t *cfg, int gen, const individual_t *best)
{
  FILE *f;
  int need_header = 0;
  struct stat st;

  if (cfg->history_out[0] == '\0')
    return;

  if (stat(cfg->history_out, &st) != 0 || st.st_size == 0)
    need_header = 1;

  f = fopen(cfg->history_out, "a");
  if (!f) {
    fprintf(stderr, "Warning: cannot open '%s' for appending: %s\n",
            cfg->history_out, strerror(errno));
    return;
  }

  if (need_header)
    fprintf(f, "generation,fitness,bytes,bps,elapsed_s,valid,objective,params\n");

  fprintf(f, "%d,%.2f,%zu,%.4f,%.2f,%d,%s,\"%s\"\n",
          gen, best->fitness, best->comp_bytes,
          compute_bps(best->comp_bytes, cfg->input_bytes),
          best->elapsed_s, best->valid,
          objective_name(cfg->objective), best->param_string);
  fclose(f);
}

/* ------------------------------------------------------------------------- */
/* CLI parsing                                                               */
/* ------------------------------------------------------------------------- */

static void print_usage(FILE *fp)
{
  fprintf(fp,
    "Usage:\n"
    "  GA --jarvis PATH --input FILE [options]\n\n"
    "Required:\n"
    "  --jarvis PATH           Path to the JARVIS3 executable\n"
    "  --input FILE            DNA input file to optimize on\n\n"
    "GA options:\n"
    "  --population N          Population size (default: 48)\n"
    "  --generations N         Number of generations (default: 30)\n"
    "  --threads N             Parallel worker threads (default: 4)\n"
    "  --elite N               Elite count (default: 4)\n"
    "  --tournament N          Tournament size (default: 3)\n"
    "  --crossover X           Crossover rate in [0,1] (default: 0.85)\n"
    "  --mutation X            Mutation rate in [0,1] (default: 0.12)\n"
    "  --toggle X              Model on/off toggle rate (default: 0.08)\n"
    "  --blend X               Blend alpha for real crossover (default: 0.35)\n"
    "  --objective NAME        bytes | bps | bytes+time (default: bytes)\n"
    "  --time-weight X         Weight for bytes+time objective (default: 1.0)\n"
    "  --seed U64              Master GA RNG seed\n"
    "  --restart-from STR      Seed population around a known parameter string\n\n"
    "Model counts:\n"
    "  --max-cmodels N         Maximum active/searchable -cm slots (default: 4)\n"
    "  --max-rmodels N         Maximum active/searchable -rm slots (default: 2)\n"
    "  --min-cmodels N         Minimum enabled -cm models (default: 0)\n"
    "  --min-rmodels N         Minimum enabled -rm models (default: 0)\n"
    "                          Set min=max to force an exact count.\n\n"
    "Global parameter control:\n"
    "  --optimize-hs 0|1       Search hidden size (default: 1)\n"
    "  --optimize-lr 0|1       Search learning rate (default: 1)\n"
    "  --optimize-seed 0|1     Search JARVIS seed (default: 0)\n"
    "  --fixed-hs N            Fixed -hs when optimize-hs=0\n"
    "  --fixed-lr X            Fixed -lr when optimize-lr=0\n"
    "  --fixed-seed N          Fixed -sd when optimize-seed=0\n\n"
    "Bounds (comma-separated key=min:max lists):\n"
    "  --global-bounds STR     hs=8:512,lr=0.0:0.2,seed=1:599999\n"
    "  --cm-bounds STR         ctx=1:14,den=1:5000,ir=0:2,gamma=0.0:0.999,\n"
    "                          edits=0:256,eden=1:50000,eir=0:1,egamma=0.0:0.999\n"
    "  --rm-bounds STR         nr=1:100000,ctx=1:14,beta=0.001:0.999,limit=1:21,\n"
    "                          gamma=0.001:0.999,ir=0:2,weight=0.001:0.999,\n"
    "                          cache=1:50\n\n"
    "Execution:\n"
    "  --verify                Run JARVIS3 decompression and compare with input\n"
    "  --workdir DIR           Temporary workspace (default: /tmp)\n"
    "  --keep-temps            Keep per-candidate temp directories/logs\n"
    "  --best-out FILE         Save best result summary\n"
    "  --history-out FILE      Append per-generation CSV history\n"
    "  --quiet                 Suppress warning chatter from failed candidates\n"
    "  --help                  Show this help\n");
}

static void parse_args(ga_config_t *cfg, int argc, char **argv)
{
  enum {
    OPT_JARVIS = 1000,
    OPT_INPUT,
    OPT_POPULATION,
    OPT_GENERATIONS,
    OPT_THREADS,
    OPT_ELITE,
    OPT_TOURNAMENT,
    OPT_CROSSOVER,
    OPT_MUTATION,
    OPT_TOGGLE,
    OPT_BLEND,
    OPT_OBJECTIVE,
    OPT_TIME_WEIGHT,
    OPT_MAX_CMODELS,
    OPT_MAX_RMODELS,
    OPT_MIN_CMODELS,
    OPT_MIN_RMODELS,
    OPT_OPT_HS,
    OPT_OPT_LR,
    OPT_OPT_SEED,
    OPT_FIXED_HS,
    OPT_FIXED_LR,
    OPT_FIXED_SEED,
    OPT_GLOBAL_BOUNDS,
    OPT_CM_BOUNDS,
    OPT_RM_BOUNDS,
    OPT_VERIFY,
    OPT_KEEP_TEMPS,
    OPT_WORKDIR,
    OPT_BEST_OUT,
    OPT_HISTORY_OUT,
    OPT_QUIET,
    OPT_RESTART_FROM,
    OPT_SEED
  };

  static struct option long_opts[] = {
    {"jarvis", required_argument, 0, OPT_JARVIS},
    {"input", required_argument, 0, OPT_INPUT},
    {"population", required_argument, 0, OPT_POPULATION},
    {"generations", required_argument, 0, OPT_GENERATIONS},
    {"threads", required_argument, 0, OPT_THREADS},
    {"elite", required_argument, 0, OPT_ELITE},
    {"tournament", required_argument, 0, OPT_TOURNAMENT},
    {"crossover", required_argument, 0, OPT_CROSSOVER},
    {"mutation", required_argument, 0, OPT_MUTATION},
    {"toggle", required_argument, 0, OPT_TOGGLE},
    {"blend", required_argument, 0, OPT_BLEND},
    {"objective", required_argument, 0, OPT_OBJECTIVE},
    {"time-weight", required_argument, 0, OPT_TIME_WEIGHT},
    {"max-cmodels", required_argument, 0, OPT_MAX_CMODELS},
    {"max-rmodels", required_argument, 0, OPT_MAX_RMODELS},
    {"min-cmodels", required_argument, 0, OPT_MIN_CMODELS},
    {"min-rmodels", required_argument, 0, OPT_MIN_RMODELS},
    {"optimize-hs", required_argument, 0, OPT_OPT_HS},
    {"optimize-lr", required_argument, 0, OPT_OPT_LR},
    {"optimize-seed", required_argument, 0, OPT_OPT_SEED},
    {"fixed-hs", required_argument, 0, OPT_FIXED_HS},
    {"fixed-lr", required_argument, 0, OPT_FIXED_LR},
    {"fixed-seed", required_argument, 0, OPT_FIXED_SEED},
    {"global-bounds", required_argument, 0, OPT_GLOBAL_BOUNDS},
    {"cm-bounds", required_argument, 0, OPT_CM_BOUNDS},
    {"rm-bounds", required_argument, 0, OPT_RM_BOUNDS},
    {"verify", no_argument, 0, OPT_VERIFY},
    {"keep-temps", no_argument, 0, OPT_KEEP_TEMPS},
    {"workdir", required_argument, 0, OPT_WORKDIR},
    {"best-out", required_argument, 0, OPT_BEST_OUT},
    {"history-out", required_argument, 0, OPT_HISTORY_OUT},
    {"quiet", no_argument, 0, OPT_QUIET},
    {"restart-from", required_argument, 0, OPT_RESTART_FROM},
    {"seed", required_argument, 0, OPT_SEED},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, "h", long_opts, NULL)) != -1) {
    switch (c) {
      case OPT_JARVIS: snprintf(cfg->jarvis_path, sizeof(cfg->jarvis_path), "%s", optarg); break;
      case OPT_INPUT: snprintf(cfg->input_path, sizeof(cfg->input_path), "%s", optarg); break;
      case OPT_POPULATION: cfg->population = atoi(optarg); break;
      case OPT_GENERATIONS: cfg->generations = atoi(optarg); break;
      case OPT_THREADS: cfg->threads = atoi(optarg); break;
      case OPT_ELITE: cfg->elite_count = atoi(optarg); break;
      case OPT_TOURNAMENT: cfg->tournament_size = atoi(optarg); break;
      case OPT_CROSSOVER: cfg->crossover_rate = atof(optarg); break;
      case OPT_MUTATION: cfg->mutation_rate = atof(optarg); break;
      case OPT_TOGGLE: cfg->toggle_rate = atof(optarg); break;
      case OPT_BLEND: cfg->blend_alpha = atof(optarg); break;
      case OPT_OBJECTIVE:
        if (strcmp(optarg, "bytes") == 0) cfg->objective = OBJ_BYTES;
        else if (strcmp(optarg, "bps") == 0) cfg->objective = OBJ_BPS;
        else if (strcmp(optarg, "bytes+time") == 0) cfg->objective = OBJ_BYTES_PLUS_TIME;
        else die("Unknown objective '%s'", optarg);
        break;
      case OPT_TIME_WEIGHT: cfg->time_weight = atof(optarg); break;
      case OPT_MAX_CMODELS: cfg->max_cmodels = atoi(optarg); break;
      case OPT_MAX_RMODELS: cfg->max_rmodels = atoi(optarg); break;
      case OPT_MIN_CMODELS: cfg->min_cmodels = atoi(optarg); break;
      case OPT_MIN_RMODELS: cfg->min_rmodels = atoi(optarg); break;
      case OPT_OPT_HS: cfg->optimize_hs = atoi(optarg) ? 1 : 0; break;
      case OPT_OPT_LR: cfg->optimize_lr = atoi(optarg) ? 1 : 0; break;
      case OPT_OPT_SEED: cfg->optimize_seed = atoi(optarg) ? 1 : 0; break;
      case OPT_FIXED_HS: cfg->fixed_hs = atoi(optarg); break;
      case OPT_FIXED_LR: cfg->fixed_lr = atof(optarg); break;
      case OPT_FIXED_SEED: cfg->fixed_seed = atoi(optarg); break;
      case OPT_GLOBAL_BOUNDS: parse_global_bounds(cfg, optarg); break;
      case OPT_CM_BOUNDS: parse_cm_bounds(cfg, optarg); break;
      case OPT_RM_BOUNDS: parse_rm_bounds(cfg, optarg); break;
      case OPT_VERIFY: cfg->verify = 1; break;
      case OPT_KEEP_TEMPS: cfg->keep_temps = 1; break;
      case OPT_WORKDIR: snprintf(cfg->workdir, sizeof(cfg->workdir), "%s", optarg); break;
      case OPT_BEST_OUT: snprintf(cfg->best_out, sizeof(cfg->best_out), "%s", optarg); break;
      case OPT_HISTORY_OUT: snprintf(cfg->history_out, sizeof(cfg->history_out), "%s", optarg); break;
      case OPT_QUIET: cfg->quiet = 1; break;
      case OPT_RESTART_FROM:
        cfg->have_restart = 1;
        snprintf(cfg->restart_from, sizeof(cfg->restart_from), "%s", optarg);
        break;
      case OPT_SEED: cfg->master_seed = (uint64_t) strtoull(optarg, NULL, 10); break;
      case 'h':
        print_usage(stdout);
        exit(EXIT_SUCCESS);
      default:
        print_usage(stderr);
        exit(EXIT_FAILURE);
    }
  }

  if (cfg->jarvis_path[0] == '\0' || cfg->input_path[0] == '\0') {
    print_usage(stderr);
    fprintf(stderr, "\n");
    die("Both --jarvis and --input are required");
  }
}

/* ------------------------------------------------------------------------- */
/* Main GA loop                                                              */
/* ------------------------------------------------------------------------- */

static void initialize_population(const ga_config_t *cfg, individual_t *pop, int n, rng_t *rng)
{
  int i;

  if (cfg->have_restart) {
    individual_t seed;
    parse_restart_individual(cfg, cfg->restart_from, &seed);

    pop[0] = seed;

    for (i = 1; i < n; ++i) {
      if (i < (n * 3) / 4) {
        pop[i] = seed;
        mutate_individual(cfg, rng, &pop[i]);
        if (rng_bool(rng, 0.35))
          mutate_individual(cfg, rng, &pop[i]);
      } else {
        random_individual(cfg, rng, &pop[i]);
      }
    }
    return;
  }

  for (i = 0; i < n; ++i)
    random_individual(cfg, rng, &pop[i]);
}

static void breed_next_population(const ga_config_t *cfg,
                                  individual_t *cur,
                                  individual_t *next,
                                  int n,
                                  rng_t *rng)
{
  int i;
  qsort(cur, (size_t) n, sizeof(individual_t), cmp_individual);

  for (i = 0; i < cfg->elite_count; ++i)
    next[i] = cur[i];

  for (; i < n; ++i) {
    const individual_t *p1 = tournament_select(cur, n, cfg->tournament_size, rng);
    const individual_t *p2 = tournament_select(cur, n, cfg->tournament_size, rng);

    if (rng_bool(rng, cfg->crossover_rate))
      crossover_individual(cfg, rng, p1, p2, &next[i]);
    else
      next[i] = rng_bool(rng, 0.5) ? *p1 : *p2;

    mutate_individual(cfg, rng, &next[i]);
  }
}

int main(int argc, char **argv)
{
  ga_config_t cfg;
  rng_t rng;
  fitness_cache_t cache;
  individual_t *pop = NULL;
  individual_t *next = NULL;
  int gen;

  set_default_config(&cfg);
  parse_args(&cfg, argc, argv);
  sanitize_config(&cfg);

  if (!file_exists_exec(cfg.jarvis_path))
    die("JARVIS executable not found or not executable: '%s'", cfg.jarvis_path);
  if (!file_exists_read(cfg.input_path))
    die("Input file not readable: '%s'", cfg.input_path);

  cfg.input_bytes = get_file_size_or_die(cfg.input_path);

  rng_seed(&rng, cfg.master_seed);

  pop = (individual_t *) xcalloc((size_t) cfg.population, sizeof(individual_t));
  next = (individual_t *) xcalloc((size_t) cfg.population, sizeof(individual_t));

  cache_init(&cache, (size_t) cfg.population * (size_t) (cfg.generations + 2));

  initialize_population(&cfg, pop, cfg.population, &rng);
  evaluate_population(&cfg, pop, cfg.population, &cache);
  qsort(pop, (size_t) cfg.population, sizeof(individual_t), cmp_individual);

  printf("\033[33m=>\033[0m Generation \033[32m%d\033[0m of %d best: fitness=%.2f bytes=\033[32m%zu\033[0m bps=\033[32m%.4f\033[0m time=\033[32m%.2f\033[0m valid=%d\n",
         0, cfg.generations, pop[0].fitness, pop[0].comp_bytes,
         compute_bps(pop[0].comp_bytes, cfg.input_bytes),
         pop[0].elapsed_s, pop[0].valid);
  print_full_command_line(stdout, &cfg, &pop[0]);
  append_history(&cfg, 0, &pop[0]);

  for (gen = 1; gen <= cfg.generations; ++gen) {
    breed_next_population(&cfg, pop, next, cfg.population, &rng);
    evaluate_population(&cfg, next, cfg.population, &cache);
    qsort(next, (size_t) cfg.population, sizeof(individual_t), cmp_individual);

    printf("\033[33m=>\033[0m Generation \033[32m%d\033[0m of %d best: fitness=%.2f bytes=\033[32m%zu\033[0m bps=\033[32m%.4f\033[0m time=\033[32m%.2f\033[0m valid=%d\n",
           gen, cfg.generations, next[0].fitness, next[0].comp_bytes,
           compute_bps(next[0].comp_bytes, cfg.input_bytes),
           next[0].elapsed_s, next[0].valid);
    print_full_command_line(stdout, &cfg, &next[0]);
    append_history(&cfg, gen, &next[0]);

    {
      individual_t *tmp = pop;
      pop = next;
      next = tmp;
    }
  }

  qsort(pop, (size_t) cfg.population, sizeof(individual_t), cmp_individual);

  printf("\nBest solution found\n");
  printf("===================\n");
  print_individual(stdout, &cfg, &pop[0], "");
  write_best_file(&cfg, &pop[0]);

  cache_destroy(&cache);
  free(pop);
  free(next);
  return 0;
}
