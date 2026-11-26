#include "chebpol.h"
#ifdef HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#endif

static const double STALKER_EPS = 1e3 * 0x1p-26; // 1e3 * sqrt(DBL_EPSILON)
typedef struct {
  double dmin, dplus, vmin, vplus;
  int pos;
} r_block;
typedef struct {
  double kq, vmin, vplus;
} s_block;
typedef struct {
  const double *det, *pmin, *pplus;
  bool uniform;
} precomp;

#ifdef HAVE_GSL
static double r_fun(double r, void *P) {
  r_block *p = (r_block *)P;
  const double dmin = p->dmin, dplus = p->dplus, vmin = p->vmin, vplus = p->vplus;
  int pos = p->pos;
  if (pos)
    return vmin * pow(dplus, r) - vplus * pow(dmin, r) - r * pow(dplus, r - 1) * (vmin * dplus + vplus * dmin);
  else
    return vmin * pow(dplus, r) - vplus * pow(dmin, r) + r * pow(dmin, r - 1) * (vmin * dplus + vplus * dmin);
}

static double find_mono_r(double dmin, double dplus, double vmin, double vplus) {
  double D2 = (vmin * dplus * dplus - vplus * dmin * dmin) / (vplus * dmin + dplus * vmin);
  double r;
  if (D2 > 2 * dplus || D2 < -2 * dmin) {
    r = 2.0;
  } else if (fabs(dplus - dmin) < STALKER_EPS * dplus) {
    // uniform grid
    r = fabs((vplus - vmin) / (vplus + vmin));
  } else {
    // Find the largest r which keeps it monotonic, i.e. derivative 0 in one
    // of the end points
    r_block rb = {dmin, dplus, vmin, vplus, 1};
    if (r_fun(1, &rb) * r_fun(2, &rb) > 0)
      rb.pos = 0;
    gsl_function F = {.function = &r_fun, .params = &rb};
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(s, &F, 1, 2);
    int status, iter = 0;
    double lo, hi;
    do {
      iter++;
      gsl_root_fsolver_iterate(s);
      lo = gsl_root_fsolver_x_lower(s);
      hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(lo, hi, 0, 1e-8);
    } while (status == GSL_CONTINUE && iter < 50);
    // Err on the low side, this will ensure that we do not
    // overshoot due to numerical inaccuracy in the solution
    r = lo;
    if (r < 1)
      r = 1;
    // If the solver did not converge we keep the current bracket: we may be inside an
    // OpenMP region, so we cannot signal an R error from here.
    gsl_root_fsolver_free(s);
  }
  return r;
}

static double r_fun_new(double r, void *P) {
  s_block *p = (s_block *)P;
  const double kq = p->kq, vmin = p->vmin, vplus = p->vplus;
  return vplus * pow(-kq, r) - vmin - (vplus * kq - vmin) * r;
}

static double solver(double kmin, double kplus, double vmin, double vplus) {
  double bn = (vplus * kmin * kmin - vmin * kplus * kplus);
  double ncn = vplus * kmin - vmin * kplus;
  double r;
  if (bn / ncn <= 2 * kmin || bn / ncn >= 2 * kplus) {
    r = 2.0;
  } else if (fabs(kplus + kmin) < STALKER_EPS * kplus) {
    // uniform grid
    r = fabs((vplus - vmin) / (vplus + vmin));
    if (r < 1)
      r = 1;
    else if (ISNAN(r) || r > 2)
      r = 2;
  } else {
    // Find the largest r which keeps it monotonic, i.e. derivative 0 in one
    // of the end points
    s_block rb;
    int eq1 = (vplus > 0 && ncn > 0) || (vplus < 0 && ncn < 0);
    if (eq1) {
      rb.kq = kmin / kplus;
      rb.vmin = vmin;
      rb.vplus = vplus;
    } else {
      rb.kq = kplus / kmin;
      rb.vmin = vplus;
      rb.vplus = vmin;
    }
    gsl_function F = {.function = &r_fun_new, .params = &rb};
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(s, &F, 1, 2);
    int status, iter = 0;
    do {
      iter++;
      gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      double lo = gsl_root_fsolver_x_lower(s);
      double hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(lo, hi, 0, 1e-8);
    } while (status == GSL_CONTINUE && iter < 50);
    // If the solver did not converge we keep its current estimate.
    gsl_root_fsolver_free(s);
  }
  return r;
}
#endif

// hyperbolic stalker
static double stalk_hyp(double x, double vmin, double vplus, double dmin, double dplus) {
  if (fabs(x - dplus) < 1e3 * DBL_EPSILON * dplus)
    return vplus;
  if (fabs(x + dmin) < 1e3 * DBL_EPSILON * dmin)
    return vmin;
  // hyperbolic function
  // a+bx-a/(cx+1)
  // The special cases are detected relative to the size of the terms that cancel, so that they do not depend on
  // the scale of the values or of the grid
  if (fabs(vplus * dmin + vmin * dplus) <= 1e3 * DBL_EPSILON * (fabs(vplus) * dmin + fabs(vmin) * dplus))
    return x * (vplus - vmin) / (dplus + dmin); // linear
  if (vplus * vmin < 0) {
    // monotonic, set b=0
    const double c = (vmin * dplus + vplus * dmin) / (dmin * dplus * (vmin - vplus));
    const double a = vplus + vplus / (c * dplus);
    return a - a / (1 + c * x);
  } else {
    // non-monotonic
    const double D = dplus * dplus * vmin - dmin * dmin * vplus;

    if (fabs(D) <= 1e3 * DBL_EPSILON * (dplus * dplus * fabs(vmin) + dmin * dmin * fabs(vplus))) {
      // if c==0, it's a parabola y = a*x^2
      // we have (0,0), (-dmin,vmin), and (dplus,vplus)
      // vmin = a*dmin^2
      // vplus = a*dplus^2
      const double a = vmin / (dmin * dmin);
      return a * x * x;
    }
    const double c = D / (dmin * dplus * (dmin * vplus + dplus * vmin));
    const double b = vplus * (1 + c * dplus) / (c * dplus * dplus);
    const double a = -b / c;
    return a + b * x - a / (1 + c * x);
  }
}

// Variable degree stalker spline
static double stalk_1(double x, double vmin, double v0, double vplus, double dmin, double dplus, double r, double iD,
                      double powmin, double powplus, bool uniform) {
  // Do linear end points

  if (ISNAN(vmin))
    return v0 + x / dplus * (vplus - v0);
  if (ISNAN(vplus))
    return v0 + x / dmin * (v0 - vmin);
  vmin -= v0;
  vplus -= v0;
  if (r == 0)
    return v0 + stalk_hyp(x, vmin, vplus, dmin, dplus);
  if (!ISNAN(r)) {
    // use precomputed degree
    const double b = (vplus * powmin - vmin * powplus) * iD;
    const double c = (vplus * dmin + vmin * dplus) * iD;
    return v0 + b * x + c * pow(fabs(x), r);
  } else {
    // Compute r for this basis function
#ifdef HAVE_GSL
    if (uniform) {
#endif
      // Uniform grid
      // This works when dplus==dmin, otherwise we must solve an equation for r.
      x /= dplus;
      const double b = 0.5 * (vplus - vmin);
      const double c = 0.5 * (vplus + vmin);

      if (c == 0.0)
        return v0 + b * x;

      // Find the degree. For uniform grids in normalized coordinates, this is quite easy
      double ac = fabs(c), ab = fabs(b);
      if (vplus * vmin <= 0) {
        // monotonic
        r = (ab < 2.0 * ac) ? ab / ac : 2.0;
      } else {
        // non-monotonic
        // pretend sign change of vplus for computation of r, that is, interchange |b| and |c|
        r = (ac < 2.0 * ab) ? ac / ab : 2.0;
      }

      if (r == 2.0) {
        return v0 + b * x + c * x * x;
      } else if (r != 1.0) {
        return v0 + b * x + c * pow(fabs(x), r);
      } else {
        return v0 + b * x + c * fabs(x);
      }
#ifdef HAVE_GSL
    } else {

      // We must solve for r:
      double b, c;
      if (vmin * vplus < 0) {
        // monotonic
        r = find_mono_r(dmin, dplus, vmin, vplus);
      } else {
        // non-monotonic, should we keep f'(0) = 0? No, won't work.
        // We can use r=2, but if we just precisely is non-monotonic
        // we should use a smaller r to avoid overshoot. But which?
        // Keep it in some way "symmetric" with the monotonic case?  How?
        // Change sign on one of the values, compute the r which
        // fits these monotonic points.
        r = find_mono_r(dmin, dplus, -vmin, vplus);
      }
      if (ISNAN(r)) {
        c = 0;
        b = vplus / dplus;
      } else {
        c = (dplus * vmin + dmin * vplus) / (pow(dmin, r) * dplus + pow(dplus, r) * dmin);
        b = (vplus - c * pow(dplus, r)) / dplus;
      }

      return v0 + b * x + c * pow(fabs(x), r);
    }
#endif
  }
}

static void make_stalker(const int ndim, const int *dims, const double **grid, const double *val, double *b, double *c,
                         double *r, int depth, int *point) {
  // Loop over every grid point
  // In each grid point, find the function living there, i.e. b, c and r for each dimension.
  // Recurse over the dimensions to visit every grid point
  if (depth < ndim) {
    for (int i = 0; i < dims[depth]; i++) {
      point[depth] = i;
      make_stalker(ndim, dims, grid, val, b, c, r, depth + 1, point);
    }
    return;
  }

  // What is the index of the point?  (We could do this in the recursion as well).
  R_xlen_t index = 0, stride = 1;

  for (int i = 0; i < ndim; i++) {
    index += point[i] * stride;
    stride *= dims[i];
  }

  // Inner loop, we have a grid point in 'point'
  // loop over the dimensions, collecting b, c and r
  double *bb = &b[ndim * index], *cc = &c[ndim * index], *rr = &r[ndim * index];
  stride = 1;

  for (int i = 0; i < ndim; i++) {
    const double *gr = grid[i];
    int idx = point[i];

    // Check border, we do it linearly there
    if (idx == 0) {
      // Nothing below us, we are linear upwards
      bb[i] = (val[index + stride] - val[index]) / (gr[idx + 1] - gr[idx]);
      cc[i] = rr[i] = 0.0;

    } else if (idx == dims[i] - 1) {
      // Nothing above us, linear downwards
      bb[i] = (val[index - stride] - val[index]) / (gr[idx - 1] - gr[idx]);
      cc[i] = rr[i] = 0.0;
    } else {
      // interior grid point
      double v0 = val[index];
      double vplus = val[index + stride] - v0;
      double vmin = val[index - stride] - v0;
      double kplus = gr[idx + 1] - gr[idx];
      double kmin = gr[idx - 1] - gr[idx];
      bool uniform = fabs(kplus + kmin) < STALKER_EPS * kplus;
      double r;

      // If non-monotonic, simulate monotonic by changing sign on vmin when computing r
      if (vplus * vmin <= 0) {
        // monotonic
        // relative tests, independent of the scale of the values
        if (vplus == 0 && vmin == 0) {
          // constant zero
          r = 2.0;
        } else if (fabs(kmin * vplus - kplus * vmin) <= STALKER_EPS * (fabs(kmin * vplus) + fabs(kplus * vmin))) {
          // Linear
          r = 2.0;
        } else {
          // general
          if (uniform) {
            r = fabs((vplus - vmin) / (vplus + vmin));
            if (r < 1)
              r = 1;
            if (ISNAN(r) || r > 2.0)
              r = 2.0;
          } else {
#ifdef HAVE_GSL
            r = solver(kmin, kplus, vmin, vplus);
#else
            r = 2.0;
#endif
          }
        }
      } else {
        // non-monotonic
        // At the outset, r = 2.0
        // But if we change sign on vmin, and the resulting monotonic triple
        // can't be handled by r=2, we use the r from there
#ifdef HAVE_GSL
        r = solver(kmin, kplus, -vmin, vplus);
#else
        r = 2.0;
#endif
      }
      rr[i] = r;
      bb[i] = (vplus - vmin * pow(-kplus / kmin, r)) / (kplus - kmin * pow(-kplus / kmin, r));
      cc[i] = (kplus * vmin - kmin * vplus) / (kplus * pow(-kmin, r) - kmin * pow(kplus, r));
    }
    stride *= dims[i];
  }
}
static double eval_stalker(const double *x, const int ndim, const int *dims, const double **grid, const double *val,
                           blend_type blend, const double *b, const double *c, const double *r) {
  // loop over the corners of the hypercube surrounding x.
  // For each corner, evaluate the function there.
  // Make weighted sum.

  R_xlen_t stride = 1;
  // the blended weights of the lower and upper end of the interval in each dimension
  double wlo[ndim], wup[ndim];
  int gbase[ndim];
  R_xlen_t llcorner = 0;
  // Find the lower left corner for each dimension, store it in gbase
  for (int g = 0; g < ndim; g++) {
    const double *gr = grid[g];
    int gp = find_interval(gr, dims[g], x[g]);
    // gp is below x
    gbase[g] = gp;
    double weight = 1 - (x[g] - gr[gbase[g]]) / (gr[gbase[g] + 1] - gr[gbase[g]]);
    // This ensures extrapolation is linear
    if (weight < 0)
      weight = 0;
    if (weight > 1)
      weight = 1;
    wlo[g] = blend_fun(weight, blend);
    wup[g] = blend_fun(1 - weight, blend);
    llcorner += stride * gp;
    stride *= dims[g];
  }

  double V = 0.0;
  double zx[ndim];
  double sumw = 0.0;
  for (int i = 0; i < (1 << ndim); i++) {
    // i is a corner represented as ndim bits
    // bit j is 1 if above in dimension j, 0 if below
    // Find the position and weight of the corner
    R_xlen_t stride = 1;
    R_xlen_t corner = llcorner;
    double cw = 1.0;
    for (int g = 0; g < ndim; g++) {
      const double *gr = grid[g];
      if (i & (1 << g)) {
        cw *= wup[g];
        corner += stride;
        // zero based in this corner
        zx[g] = x[g] - gr[gbase[g] + 1];
      } else {
        zx[g] = x[g] - gr[gbase[g]];
        cw *= wlo[g];
      }
      stride *= dims[g];
    }
    if (cw == 0.0)
      continue;
    // Find b, c, and r for the corner. These have stride 'ndim'
    const double *bb = &b[ndim * corner], *cc = &c[ndim * corner], *rr = &r[ndim * corner];
    // These are for a zero-based function, i.e. based on the corner being 0,
    // and f(0)=0
    // Compute the stalker in x
    // The zero-based x is in zx
    double fval = val[corner];
    for (int j = 0; j < ndim; j++) {
      if (cc[j] == 0) {
        fval += bb[j] * zx[j];
      } else {
        fval += bb[j] * zx[j] + cc[j] * pow(fabs(zx[j]), rr[j]);
      }
    }
    sumw += cw;
    V += cw * fval;
  }
  return V / sumw;
}

// Find the parameters for the hyperbolic stalker in each grid point
// use normalized coordinates as in the formulas, f(0) = 0
static void make_hstalker(const int ndim, const int *dims, const double **grid, const double *val, double *a, double *b,
                          double *c, double *d, int depth, int *point) {
  // Loop over every grid point
  // In each grid point, find the function living there, i.e. an a,
  // and b,c,d for each dimension.
  // Recurse over the dimensions to visit every grid point
  if (depth < ndim) {
    for (int i = 0; i < dims[depth]; i++) {
      point[depth] = i;
      make_hstalker(ndim, dims, grid, val, a, b, c, d, depth + 1, point);
    }
    return;
  }

  // What is the index of the point?  (We could do this in the recursion as well).
  R_xlen_t index = 0, stride = 1;

  for (int i = 0; i < ndim; i++) {
    index += point[i] * stride;
    stride *= dims[i];
  }

  // Inner loop, we have a grid point in 'point'
  // loop over the dimensions, collecting b, c, and d

  double *bb = &b[ndim * index], *cc = &c[ndim * index], *dd = &d[ndim * index];
  a[index] = 0.0;
  stride = 1;
  for (int i = 0; i < ndim; i++) {
    const double *gr = grid[i];
    int idx = point[i];

    // Check border, we do it linearly there
    if (idx == 0) {
      // Nothing below us, we are linear upwards
      bb[i] = (val[index + stride] - val[index]) / (gr[idx + 1] - gr[idx]);
      cc[i] = dd[i] = 0.0;

    } else if (idx == dims[i] - 1) {
      // Nothing above us, linear downwards
      bb[i] = (val[index - stride] - val[index]) / (gr[idx - 1] - gr[idx]);
      cc[i] = dd[i] = 0.0;

    } else {
      // interior grid point
      double v0 = val[index];
      double vplus = val[index + stride] - v0;
      double vmin = val[index - stride] - v0;
      double kplus = gr[idx + 1] - gr[idx];
      double kmin = gr[idx - 1] - gr[idx];
      // Work in normalized coordinates, x / kplus and f / vscale, so that the special cases
      // and the intermediate products do not depend on the scale of the values or of the grid.
      const double vscale = fmax(fabs(vplus), fabs(vmin));
      if (vscale == 0) {
        // constant zero
        bb[i] = cc[i] = dd[i] = 0.0;
      } else {
        const double p = vplus / vscale;
        const double m = vmin / vscale;
        const double k = kmin / kplus;
        const double D = k * p - m;
        const double D2 = k * k * p - m;
        // The special cases are detected relative to the size of the terms that cancel
        if (fabs(p) <= 4 * DBL_EPSILON || fabs(m) <= 4 * DBL_EPSILON) {
          // Flat towards one neighbor. The hyperbola degenerates to zero with its pole on the
          // other neighbor, where evaluating it would give 0 / 0.
          bb[i] = cc[i] = dd[i] = 0.0;
        } else if (p * m <= 0) {
          // Monotonic
          if (fabs(D) <= STALKER_EPS * (fabs(k * p) + fabs(m))) {
            // Linear
            bb[i] = vplus / kplus;
            cc[i] = dd[i] = 0.0;
          } else {
            // general
            bb[i] = 0.0;
            cc[i] = -D / (k * (p - m)) / kplus;
            dd[i] = -p * m * (k - 1) / D * vscale;
          }
        } else {
          // Non-monotonic
          if (fabs(D2) <= STALKER_EPS * (fabs(k * k * p) + fabs(m))) {
            // Parabola
            bb[i] = vplus / kplus / kplus;
            dd[i] = 0.0;
            cc[i] = NA_REAL; // signals parabola
          } else {
            // general
            bb[i] = p * m * (k - 1) / D2 * vscale / kplus;
            cc[i] = -D2 / (k * D) / kplus;
            dd[i] = -k * p * m * (k - 1) * D / (D2 * D2) * vscale;
          }
        }
      }
    }
    a[index] -= dd[i];
    stride *= dims[i];
  }
}

static double eval_hstalker(const double *x, const int ndim, const int *dims, const double **grid, const double *val,
                            blend_type blend, const double *a, const double *b, const double *c, const double *d) {
  // loop over the corners of the hypercube surrounding x.
  // For each corner, evaluate the function there.
  // Make weighted sum.

  R_xlen_t stride = 1;
  // the blended weights of the lower and upper end of the interval in each dimension
  double wlo[ndim], wup[ndim];
  int gbase[ndim];
  R_xlen_t llcorner = 0;
  // Find the lower left corner for each dimension, store it in gbase
  for (int g = 0; g < ndim; g++) {
    const double *gr = grid[g];
    int gp = find_interval(gr, dims[g], x[g]);
    // gp is below x
    gbase[g] = gp;
    double weight = 1 - (x[g] - gr[gbase[g]]) / (gr[gbase[g] + 1] - gr[gbase[g]]);
    // This ensures extrapolation is linear
    if (weight < 0)
      weight = 0;
    if (weight > 1)
      weight = 1;
    wlo[g] = blend_fun(weight, blend);
    wup[g] = blend_fun(1 - weight, blend);
    llcorner += stride * gp;
    stride *= dims[g];
  }

  double V = 0.0;
  double zx[ndim];
  double sumw = 0.0;
  for (int i = 0; i < (1 << ndim); i++) {
    // i is a corner represented as ndim bits
    // bit j is 1 if above in dimension j, 0 if below
    // Find the position and weight of the corner
    R_xlen_t stride = 1;
    R_xlen_t corner = llcorner;
    double cw = 1.0;
    for (int g = 0; g < ndim; g++) {
      const double *gr = grid[g];
      if (i & (1 << g)) {
        cw *= wup[g];
        corner += stride;
        // zero based in this corner
        zx[g] = x[g] - gr[gbase[g] + 1];
      } else {
        zx[g] = x[g] - gr[gbase[g]];
        cw *= wlo[g];
      }
      stride *= dims[g];
    }
    if (cw == 0.0)
      continue;
    // Find b, c, and d for the corner. These have stride 'ndim'
    const double *bb = &b[ndim * corner], *cc = &c[ndim * corner], *dd = &d[ndim * corner];
    // These are for a zero-based function, i.e. based on the corner being 0,
    // and f(0)=0
    // Compute the hyperbolic stalker in x
    // The zero-based x is in zx
    double fval = val[corner] + a[corner];
    for (int j = 0; j < ndim; j++) {
      // It might be a parabola, but typically a hyperbola
      if (!ISNAN(cc[j])) {
        fval += bb[j] * zx[j] + dd[j] / (1.0 + cc[j] * zx[j]);
      } else {
        fval += bb[j] * zx[j] * zx[j];
      }
    }
    sumw += cw;
    V += cw * fval;
  }
  return V / sumw;
}

static double eval_oldstalker(const double *x, const int ndim, const int *dims, const double **grid, const double *val,
                              blend_type blend, const double *degree, precomp *pcomp) {

  if (ndim == 0)
    return val[0];
  const int last = ndim - 1;
  const double *gr = grid[last];
  const int N = dims[last];
  const double xx = x[last];
  const double *rdet = pcomp[last].det, *rpmin = pcomp[last].pmin, *rpplus = pcomp[last].pplus;
  bool uniform = pcomp[last].uniform;
  R_xlen_t stride = 1;
  for (int r = 0; r < last; r++)
    stride *= dims[r];
  // Find the interval of the last coordinate
  const int imin = find_interval(gr, N, xx);

  // Now, find the function values on the two grid points on each side.
  // i.e. imin, imin-1, and imin+1 and imin+2
  // use these to create the stalker splines across this dimension

  // values of the lower dimensional spline
  double v1 = NA_REAL, v2, v3, v4 = NA_REAL;
  if (last > 0) {
    if (imin > 0)
      v1 = eval_oldstalker(x, last, dims, grid, val + (imin - 1) * stride, blend, degree, pcomp);
    v2 = eval_oldstalker(x, last, dims, grid, val + imin * stride, blend, degree, pcomp);
    v3 = eval_oldstalker(x, last, dims, grid, val + (imin + 1) * stride, blend, degree, pcomp);
    if (imin < N - 2)
      v4 = eval_oldstalker(x, last, dims, grid, val + (imin + 2) * stride, blend, degree, pcomp);
  } else {
    // save a recursion step
    if (imin > 0)
      v1 = val[(imin - 1) * stride];
    v2 = val[imin * stride];
    v3 = val[(imin + 1) * stride];
    if (imin < N - 2)
      v4 = val[(imin + 2) * stride];
  }
  double dmin = NA_REAL;
  if (imin > 0)
    dmin = gr[imin] - gr[imin - 1];
  double dplus = gr[imin + 1] - gr[imin];
  const double low =
      stalk_1(xx - gr[imin], v1, v2, v3, dmin, dplus, degree[last], rdet[imin], rpmin[imin], rpplus[imin], uniform);
  dmin = dplus;
  if (imin < N - 2)
    dplus = gr[imin + 2] - gr[imin + 1];
  else
    dplus = NA_REAL;
  const double high = stalk_1(xx - gr[imin + 1], v2, v3, v4, dmin, dplus, degree[last], rdet[imin + 1], rpmin[imin + 1],
                              rpplus[imin + 1], uniform);

  // combine the basis functions
  // map into normalized coordinates
  const double nx = (xx - gr[imin]) / (gr[imin + 1] - gr[imin]);
  double w = 1 - nx;

  if (w == 0)
    return high;
  if (w == 1)
    return low;
  w = blend_fun(w, blend);
  return w * low + (1 - w) * high;
}

SEXP c_make_stalker(SEXP s_val, SEXP s_grid) {
  const int ndim = LENGTH(s_grid);
  int dims[ndim];
  R_xlen_t len = 1;
  for (int r = 0; r < ndim; r++) {
    dims[r] = LENGTH(VECTOR_ELT(s_grid, r));
    if (dims[r] < 2)
      Rf_error("`grid` dimension %d must have at least two points, not %d.", r + 1, dims[r]);
    len *= dims[r];
  }
  if (len != XLENGTH(s_val)) {
    Rf_error("`val` must have one value per grid point (%lld), not %lld.", (long long)len, (long long)XLENGTH(s_val));
  }
  if (len > INT_MAX) {
    Rf_error("grid size (%lld) is too large", (long long)len);
  }
  const double *val = REAL_RO(s_val);
  const double *grid[ndim];
  for (int i = 0; i < ndim; i++)
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
  int point[ndim];
  SEXP s_b = PROTECT(Rf_allocMatrix(REALSXP, ndim, (int)len));
  SEXP s_c = PROTECT(Rf_allocMatrix(REALSXP, ndim, (int)len));
  SEXP s_r = PROTECT(Rf_allocMatrix(REALSXP, ndim, (int)len));
  const char *names[] = {"b", "c", "r", "val", "grid", ""};
  SEXP s_out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(s_out, 0, s_b);
  SET_VECTOR_ELT(s_out, 1, s_c);
  SET_VECTOR_ELT(s_out, 2, s_r);
  SET_VECTOR_ELT(s_out, 3, s_val);
  SET_VECTOR_ELT(s_out, 4, s_grid);
#ifdef HAVE_GSL
  gsl_error_handler_t *gslh = gsl_set_error_handler_off();
#else
  // Check if grid is uniform
  for (int i = 0; i < ndim; i++) {
    const double *gr = grid[i];
    double df = gr[1] - gr[0];
    for (int j = 2; j < dims[i]; j++)
      if (fabs(gr[j] - gr[j - 1] - df) > STALKER_EPS)
        Rf_warning("Non-uniform grid requires linking chebpol with GSL");
  }
#endif

  make_stalker(ndim, dims, grid, val, REAL(s_b), REAL(s_c), REAL(s_r), 0, point);

#ifdef HAVE_GSL
  gsl_set_error_handler(gslh);
#endif

  UNPROTECT(4);
  return s_out;
}

SEXP c_make_hstalker(SEXP s_val, SEXP s_grid) {
  const int ndim = LENGTH(s_grid);
  int dims[ndim];
  R_xlen_t len = 1;
  for (int r = 0; r < ndim; r++) {
    dims[r] = LENGTH(VECTOR_ELT(s_grid, r));
    if (dims[r] < 2)
      Rf_error("`grid` dimension %d must have at least two points, not %d.", r + 1, dims[r]);
    len *= dims[r];
  }
  if (len != XLENGTH(s_val)) {
    Rf_error("`val` must have one value per grid point (%lld), not %lld.", (long long)len, (long long)XLENGTH(s_val));
  }
  if (len > INT_MAX) {
    Rf_error("grid size (%lld) is too large", (long long)len);
  }
  const double *val = REAL_RO(s_val);
  const double *grid[ndim];
  for (int i = 0; i < ndim; i++)
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
  int point[ndim];
  SEXP s_a = PROTECT(Rf_allocVector(REALSXP, len));
  SEXP s_b = PROTECT(Rf_allocMatrix(REALSXP, ndim, (int)len));
  SEXP s_c = PROTECT(Rf_allocMatrix(REALSXP, ndim, (int)len));
  SEXP s_d = PROTECT(Rf_allocMatrix(REALSXP, ndim, (int)len));
  const char *names[] = {"a", "b", "c", "d", "val", "grid", ""};
  SEXP s_out = PROTECT(Rf_mkNamed(VECSXP, names));
  SET_VECTOR_ELT(s_out, 0, s_a);
  SET_VECTOR_ELT(s_out, 1, s_b);
  SET_VECTOR_ELT(s_out, 2, s_c);
  SET_VECTOR_ELT(s_out, 3, s_d);
  SET_VECTOR_ELT(s_out, 4, s_val);
  SET_VECTOR_ELT(s_out, 5, s_grid);

  make_hstalker(ndim, dims, grid, val, REAL(s_a), REAL(s_b), REAL(s_c), REAL(s_d), 0, point);
  UNPROTECT(5);
  return s_out;
}

SEXP c_eval_stalker(SEXP s_x, SEXP s_stalker, SEXP s_blend, SEXP s_threads) {
  const double *x = REAL_RO(s_x);
  const int n_points = Rf_ncols(s_x);
  blend_type blend = (blend_type)Rf_asInteger(s_blend);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  const double *b = REAL_RO(VECTOR_ELT(s_stalker, 0));
  const double *c = REAL_RO(VECTOR_ELT(s_stalker, 1));
  const double *r = REAL_RO(VECTOR_ELT(s_stalker, 2));
  const double *val = REAL_RO(VECTOR_ELT(s_stalker, 3));
  SEXP s_grid = VECTOR_ELT(s_stalker, 4);
  const int ndim = LENGTH(s_grid);
  if (Rf_nrows(s_x) != ndim)
    Rf_error("`x` must have %d rows, one per dimension, not %d.", ndim, Rf_nrows(s_x));
  const double *grid[ndim];
  for (int i = 0; i < ndim; i++)
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
  int dims[ndim];
  for (int r = 0; r < ndim; r++)
    dims[r] = LENGTH(VECTOR_ELT(s_grid, r));

  SEXP s_out = PROTECT(Rf_allocVector(REALSXP, n_points));
  double *out = REAL(s_out);

#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(guided) if (n_points > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    out[i] = eval_stalker(x + i * ndim, ndim, dims, grid, val, blend, b, c, r);
  }

  UNPROTECT(1);
  return s_out;
}

SEXP c_eval_hstalker(SEXP s_x, SEXP s_stalker, SEXP s_blend, SEXP s_threads) {
  const double *x = REAL_RO(s_x);
  const int n_points = Rf_ncols(s_x);
  blend_type blend = (blend_type)Rf_asInteger(s_blend);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  const double *a = REAL_RO(VECTOR_ELT(s_stalker, 0));
  const double *b = REAL_RO(VECTOR_ELT(s_stalker, 1));
  const double *c = REAL_RO(VECTOR_ELT(s_stalker, 2));
  const double *d = REAL_RO(VECTOR_ELT(s_stalker, 3));
  const double *val = REAL_RO(VECTOR_ELT(s_stalker, 4));
  SEXP s_grid = VECTOR_ELT(s_stalker, 5);
  const int ndim = LENGTH(s_grid);
  if (Rf_nrows(s_x) != ndim)
    Rf_error("`x` must have %d rows, one per dimension, not %d.", ndim, Rf_nrows(s_x));
  const double *grid[ndim];
  for (int i = 0; i < ndim; i++)
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
  int dims[ndim];
  for (int r = 0; r < ndim; r++)
    dims[r] = LENGTH(VECTOR_ELT(s_grid, r));

  SEXP s_out = PROTECT(Rf_allocVector(REALSXP, n_points));
  double *out = REAL(s_out);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(guided) if (n_points > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    out[i] = eval_hstalker(x + i * ndim, ndim, dims, grid, val, blend, a, b, c, d);
  }

  UNPROTECT(1);
  return s_out;
}

SEXP c_eval_oldstalker(SEXP s_x, SEXP s_stalker, SEXP s_degree, SEXP s_blend, SEXP s_threads) {
  const double *x = REAL_RO(s_x);
  const int n_rows = Rf_nrows(s_x);
  const int n_points = Rf_ncols(s_x);
#ifdef _OPENMP
  int threads = Rf_asInteger(s_threads);
#endif
  const double *val = REAL_RO(VECTOR_ELT(s_stalker, 0));
  SEXP s_grid = VECTOR_ELT(s_stalker, 1);
  SEXP s_det = VECTOR_ELT(s_stalker, 2);
  SEXP s_pmin = VECTOR_ELT(s_stalker, 3);
  SEXP s_pplus = VECTOR_ELT(s_stalker, 4);
  SEXP s_uniform = VECTOR_ELT(s_stalker, 5);
  blend_type blend = (blend_type)Rf_asInteger(s_blend);
  const int ndim = LENGTH(s_grid);
  if (n_rows != ndim)
    Rf_error("`x` must have %d rows, one per dimension, not %d.", ndim, n_rows);

  int dims[ndim];
  R_xlen_t len = 1;
  for (int r = 0; r < ndim; r++) {
    dims[r] = LENGTH(VECTOR_ELT(s_grid, r));
    if (dims[r] < 2)
      Rf_error("`grid` dimension %d must have at least two points, not %d.", r + 1, dims[r]);
    len *= dims[r];
  }
  if (len != XLENGTH(VECTOR_ELT(s_stalker, 0))) {
    Rf_error("number of values (%lld) does not match grid size (%lld)", (long long)XLENGTH(VECTOR_ELT(s_stalker, 0)),
             (long long)len);
  }
  if (LENGTH(s_degree) != ndim) {
    Rf_error("degree length must match dimension (%d)", ndim);
  }
  SEXP s_deg = PROTECT(Rf_coerceVector(s_degree, REALSXP));
  const double *degree = REAL_RO(s_deg);
  const double *grid[ndim];
  for (int i = 0; i < ndim; i++)
    grid[i] = REAL_RO(VECTOR_ELT(s_grid, i));
  precomp pcomp[ndim];
  for (int i = 0; i < ndim; i++) {
    pcomp[i].det = REAL_RO(VECTOR_ELT(s_det, i));
    pcomp[i].pmin = REAL_RO(VECTOR_ELT(s_pmin, i));
    pcomp[i].pplus = REAL_RO(VECTOR_ELT(s_pplus, i));
    pcomp[i].uniform = LOGICAL_RO(s_uniform)[i];
#ifndef HAVE_GSL
    if (!pcomp[i].uniform && ISNAN(degree[i]))
      Rf_error("`degree = NA` on a non-uniform grid requires chebpol built with GSL.");
#endif
  }
#ifdef HAVE_GSL
  gsl_error_handler_t *gslh = gsl_set_error_handler_off();
#endif
  SEXP s_out = PROTECT(Rf_allocVector(REALSXP, n_points));
  double *out = REAL(s_out);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(guided) if (n_points > 1 && threads > 1)
#endif
  for (R_xlen_t i = 0; i < n_points; i++) {
    out[i] = eval_oldstalker(x + i * ndim, ndim, dims, grid, val, blend, degree, pcomp);
  }
#ifdef HAVE_GSL
  gsl_set_error_handler(gslh);
#endif
  UNPROTECT(2);
  return s_out;
}

SEXP c_has_gsl(void) {
#ifdef HAVE_GSL
  return Rf_ScalarLogical(TRUE);
#else
  return Rf_ScalarLogical(FALSE);
#endif
}
