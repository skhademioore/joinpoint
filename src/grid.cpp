#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// Exhaustive grid search for segmented (joinpoint) regression.
//
// Columns of the design are ordered as p0 fixed base columns followed by the
// m candidate hinge columns (x - tau_j)_+ for sorted candidates tau_1 < ... < tau_m.
// The search works on cross products only:
//   G  = X' W X       ((p0 + m) x (p0 + m)) for all base and candidate columns,
//   Bt = Y' W X       (nd x (p0 + m)) for nd response vectors,
//   yy = diag(Y' W Y) (nd).
// A depth first enumeration adds one joinpoint at a time and extends the
// Cholesky factor of the selected Gram block by one row, so each visited
// combination costs O(p^2 + nd * p). The residual sum of squares of a model is
// yy - ||z||^2 where z = L^{-1} X' W y.
//
// Constraints are given per candidate:
//   first_ok[j]  candidate j may be the first joinpoint,
//   last_ok[j]   candidate j may be the last joinpoint,
//   next_min[j]  smallest candidate index (0 based) allowed after j.

struct GridState {
  const double *G;
  const double *Bt;           // nd x q cross products X' W Y (transposed)
  const double *yy;
  int q, p0, m, nd, kmax;
  const int *first_ok;
  const int *last_ok;
  const int *next_min;
  std::vector<double> L;      // pmax x pmax lower triangle, row major
  std::vector<double> Z;      // pmax x nd
  std::vector<double> ss;     // (kmax + 1) x nd partial sums of squares
  std::vector<int> sel;       // selected Gram column indices
  std::vector<int> path;      // selected candidate indices (0 based)
  std::vector<double> best;   // nd x (kmax + 1) minimum SSE
  std::vector<int> argbest;   // per k: nd x k candidate indices (1 based)
  std::vector<int> offs;      // offsets of each k block in argbest
  std::vector<double> lrow;
  std::vector<double> tmp;    // nd scratch vector
  double tol;
  long long visited;
  bool profile;               // record profiles (first response only)
  std::vector<double> prof;   // per k: m x k minimum SSE with tau_j fixed
  std::vector<long long> poff;
};

static inline double gram(const GridState &s, int a, int b) {
  return s.G[a + (long long)b * s.q];
}

static void record(GridState &s, int k) {
  const double *ssk = &s.ss[(long long)k * s.nd];
  if (s.profile && k > 0) {
    double sse = s.yy[0] - ssk[0];
    if (sse < 0) sse = 0;
    for (int r = 0; r < k; r++) {
      double &cur = s.prof[s.poff[k] + s.path[r] + (long long)r * s.m];
      if (sse < cur) cur = sse;
    }
  }
  const int nd = s.nd;
  const double *__restrict__ yy = s.yy;
  double *__restrict__ best = &s.best[(long long)k * nd];
  for (int d = 0; d < nd; d++) {
    const double sse = yy[d] - ssk[d];
    if (sse < best[d]) {
      best[d] = sse < 0 ? 0 : sse;
      if (k > 0) {
        int *dst = &s.argbest[s.offs[k] + (long long)d * k];
        for (int r = 0; r < k; r++) dst[r] = s.path[r] + 1;
      }
    }
  }
}

// Adds Gram column c as row p of the factor; returns false when the column is
// numerically dependent on the selected ones.
static bool extend(GridState &s, int p, int c, int k) {
  int pmax = s.p0 + s.kmax;
  for (int r = 0; r < p; r++) {
    double v = gram(s, s.sel[r], c);
    const double *Lr = &s.L[(long long)r * pmax];
    for (int t = 0; t < r; t++) v -= Lr[t] * s.lrow[t];
    s.lrow[r] = v / Lr[r];
  }
  double gcc = gram(s, c, c);
  double d2 = gcc;
  for (int t = 0; t < p; t++) d2 -= s.lrow[t] * s.lrow[t];
  if (!(d2 > s.tol * gcc) || !std::isfinite(d2)) return false;
  double dg = std::sqrt(d2);
  double *Lp = &s.L[(long long)p * pmax];
  for (int t = 0; t < p; t++) Lp[t] = s.lrow[t];
  Lp[p] = dg;
  s.sel[p] = c;
  const int nd = s.nd;
  double *__restrict__ tmp = s.tmp.data();
  const double *__restrict__ bc = &s.Bt[(long long)c * nd];
  for (int d = 0; d < nd; d++) tmp[d] = bc[d];
  for (int t = 0; t < p; t++) {
    const double lt = Lp[t];
    const double *__restrict__ Zt = &s.Z[(long long)t * nd];
    for (int d = 0; d < nd; d++) tmp[d] -= lt * Zt[d];
  }
  const double inv = 1.0 / dg;
  double *__restrict__ Zp = &s.Z[(long long)p * nd];
  const double *__restrict__ ssprev = &s.ss[(long long)(k - 1) * nd];
  double *__restrict__ ssk = &s.ss[(long long)k * nd];
  for (int d = 0; d < nd; d++) {
    const double v = tmp[d] * inv;
    Zp[d] = v;
    ssk[d] = ssprev[d] + v * v;
  }
  return true;
}

static void descend(GridState &s, int k, int from) {
  // k joinpoints already placed; try candidate j >= from as joinpoint k + 1.
  int p = s.p0 + k;
  for (int j = from; j < s.m; j++) {
    if (k == 0 && !s.first_ok[j]) continue;
    if ((++s.visited & 4095) == 0) Rcpp::checkUserInterrupt();
    if (!extend(s, p, s.p0 + j, k + 1)) continue;
    s.path[k] = j;
    if (s.last_ok[j]) record(s, k + 1);
    if (k + 1 < s.kmax && s.next_min[j] < s.m) descend(s, k + 1, s.next_min[j]);
  }
}

// [[Rcpp::export(rng = false)]]
List jp_grid_cpp(NumericMatrix G, NumericMatrix Bt, NumericVector yy, int p0,
                 int kmax, IntegerVector first_ok, IntegerVector last_ok,
                 IntegerVector next_min, bool profile = false,
                 double tol = 1e-10) {
  GridState s;
  s.profile = profile;
  s.q = G.nrow();
  s.p0 = p0;
  s.m = s.q - p0;
  s.nd = Bt.nrow();
  if (kmax > s.m) kmax = s.m;
  if (kmax < 0) kmax = 0;
  s.kmax = kmax;
  s.G = G.begin();
  s.Bt = Bt.begin();
  s.yy = yy.begin();
  s.first_ok = first_ok.begin();
  s.last_ok = last_ok.begin();
  s.next_min = next_min.begin();
  s.tol = tol;
  s.visited = 0;
  int pmax = p0 + kmax;
  s.L.assign((size_t)pmax * pmax, 0.0);
  s.Z.assign((size_t)pmax * s.nd, 0.0);
  s.ss.assign((size_t)(kmax + 1) * s.nd, 0.0);
  s.sel.assign(pmax, 0);
  s.path.assign(kmax > 0 ? kmax : 1, 0);
  s.lrow.assign(pmax, 0.0);
  s.tmp.assign(s.nd, 0.0);
  s.best.assign((size_t)s.nd * (kmax + 1), std::numeric_limits<double>::infinity());
  s.offs.assign(kmax + 2, 0);
  for (int k = 1; k <= kmax + 1; k++) s.offs[k] = s.offs[k - 1] + (k - 1) * s.nd;
  s.argbest.assign((size_t)s.offs[kmax + 1] + 1, NA_INTEGER);
  s.poff.assign(kmax + 2, 0);
  for (int k = 1; k <= kmax + 1; k++) s.poff[k] = s.poff[k - 1] + (long long)s.m * (k - 1);
  if (profile) {
    s.prof.assign((size_t)s.poff[kmax + 1], std::numeric_limits<double>::infinity());
  }

  // Base model: factor the fixed columns one by one.
  for (int c = 0; c < p0; c++) {
    double *Lc = &s.L[(long long)c * pmax];
    for (int r = 0; r < c; r++) {
      double v = gram(s, r, c);
      const double *Lr = &s.L[(long long)r * pmax];
      for (int t = 0; t < r; t++) v -= Lr[t] * Lc[t];
      Lc[r] = v / Lr[r];
    }
    double d2 = gram(s, c, c);
    for (int t = 0; t < c; t++) d2 -= Lc[t] * Lc[t];
    if (!(d2 > tol * gram(s, c, c))) stop("base design columns are linearly dependent");
    Lc[c] = std::sqrt(d2);
    s.sel[c] = c;
    for (int d = 0; d < s.nd; d++) {
      double v = s.Bt[d + (long long)c * s.nd];
      for (int t = 0; t < c; t++) v -= Lc[t] * s.Z[(long long)t * s.nd + d];
      s.Z[(long long)c * s.nd + d] = v / Lc[c];
    }
  }
  for (int d = 0; d < s.nd; d++) {
    double acc = 0;
    for (int t = 0; t < p0; t++) {
      double z = s.Z[(long long)t * s.nd + d];
      acc += z * z;
    }
    s.ss[d] = acc;
  }
  record(s, 0);
  if (kmax > 0) descend(s, 0, 0);

  NumericMatrix best(s.nd, kmax + 1);
  std::copy(s.best.begin(), s.best.end(), best.begin());
  List idx(kmax + 1);
  idx[0] = IntegerMatrix(s.nd, 0);
  for (int k = 1; k <= kmax; k++) {
    IntegerMatrix a(s.nd, k);
    for (int d = 0; d < s.nd; d++)
      for (int r = 0; r < k; r++) a(d, r) = s.argbest[s.offs[k] + (long long)d * k + r];
    idx[k] = a;
  }
  List prof(kmax + 1);
  for (int k = 1; k <= kmax && profile; k++) {
    NumericMatrix pk(s.m, k);
    std::copy(s.prof.begin() + s.poff[k], s.prof.begin() + s.poff[k + 1], pk.begin());
    prof[k] = pk;
  }
  return List::create(_["sse"] = best, _["idx"] = idx, _["profile"] = prof);
}
