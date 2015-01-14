// The Portion of my code use the adaptive multidimensional integration written by Steven G. Johnson.
// The following comments are the copyright description in his C package::

/* Adaptive multidimensional integration of a vector of integrands.
 *
 * Copyright (c) 2005-2013 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Adaptive multidimensional integration on hypercubes (or, really,
   hyper-rectangles) using cubature rules.

   A cubature rule takes a function and a hypercube and evaluates
   the function at a small number of points, returning an estimate
   of the integral as well as an estimate of the error, and also
   a suggested dimension of the hypercube to subdivide.

   Given such a rule, the adaptive integration is simple:

   1) Evaluate the cubature rule on the hypercube(s).
   Stop if converged.

   2) Pick the hypercube with the largest estimated error,
   and divide it in two along the suggested dimension.

   3) Goto (1).

   The basic algorithm is based on the adaptive cubature described in
 
   A. C. Genz and A. A. Malik, "An adaptive algorithm for numeric
   integration over an N-dimensional rectangular region,"
   J. Comput. Appl. Math. 6 (4), 295-302 (1980).

   and subsequently extended to integrating a vector of integrands in

   J. Berntsen, T. O. Espelid, and A. Genz, "An adaptive algorithm
   for the approximate calculation of multiple integrals,"
   ACM Trans. Math. Soft. 17 (4), 437-451 (1991).

   Note, however, that we do not use any of code from the above authors
   (in part because their code is Fortran 77, but mostly because it is
   under the restrictive ACM copyright license).  I did make use of some
   GPL code from Rudolf Schuerer's HIntLib and from the GNU Scientific
   Library as listed in the copyright notice above, on the other hand.

   I am also grateful to Dmitry Turbiner <dturbiner@alum.mit.edu>, who
   implemented an initial prototype of the "vectorized" functionality
   for evaluating multiple points in a single call (as opposed to
   multiple functions in a single call).  (Although Dmitry implemented
   a working version, I ended up re-implementing this feature from
   scratch as part of a larger code-cleanup, and in order to have
   a single code path for the vectorized and non-vectorized APIs.  I
   subsequently implemented the algorithm by Gladwell to extract
   even more parallelism by evalutating many hypercubes at once.)

   TODO:

   * Putting these routines into the GNU GSL library would be nice.

   * A Python interface would be nice.  (Also a Matlab interface,
   a GNU Octave interface, ...)

   * For high-dimensional integrals, it would be nice to implement
   a sparse-grid cubature scheme using Clenshaw-Curtis quadrature.
   Currently, for dimensions > 7 or so, quasi Monte Carlo methods win.

   * Berntsen et. al also describe a "two-level" error estimation scheme
   that they claim makes the algorithm more robust.  It might be
   nice to implement this, at least as an option (although I seem
   to remember trying it once and it made the number of evaluations
   substantially worse for my test integrands).

*/

/* USAGE: Call cubature with your function as described in cubature.h.

   To compile a test program, compile cubature.c with
   -DTEST_INTEGRATOR as described at the end. */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

# include <R.h>
# include <Rinternals.h>
# include <Rmath.h>
# include <R_ext/Linpack.h>
# include "cubature.h"




#include "cubature.h"

/* error return codes */
#define SUCCESS 0
#define FAILURE 1

/***************************************************************************/
/* Basic datatypes */

typedef struct {
  double val, err;
} esterr;

static double errMax(unsigned fdim, const esterr *ee)
{
  double errmax = 0;
  unsigned k;
  for (k = 0; k < fdim; ++k)
    if (ee[k].err > errmax) errmax = ee[k].err;
  return errmax;
}

typedef struct {
  unsigned dim;
  double *data;	/* length 2*dim = center followed by half-widths */
  double vol;	/* cache volume = product of widths */
} hypercube;

static double compute_vol(const hypercube *h)
{
  unsigned i;
  double vol = 1;
  for (i = 0; i < h->dim; ++i)
    vol *= 2 * h->data[i + h->dim];
  return vol;
}

static hypercube make_hypercube(unsigned dim, const double *center, const double *halfwidth)
{
  unsigned i;
  hypercube h;
  h.dim = dim;
  h.data = (double *) malloc(sizeof(double) * dim * 2);
  h.vol = 0;
  if (h.data) {
    for (i = 0; i < dim; ++i) {
      h.data[i] = center[i];
      h.data[i + dim] = halfwidth[i];
    }
    h.vol = compute_vol(&h);
  }
  return h;
}

static hypercube make_hypercube_range(unsigned dim, const double *xmin, const double *xmax)
{
  hypercube h = make_hypercube(dim, xmin, xmax);
  unsigned i;
  if (h.data) {
    for (i = 0; i < dim; ++i) {
      h.data[i] = 0.5 * (xmin[i] + xmax[i]);
      h.data[i + dim] = 0.5 * (xmax[i] - xmin[i]);
    }
    h.vol = compute_vol(&h);
  }
  return h;
}

static void destroy_hypercube(hypercube *h)
{
  free(h->data);
  h->dim = 0;
}

typedef struct {
  hypercube h;
  unsigned splitDim;
  unsigned fdim; /* dimensionality of vector integrand */
  esterr *ee; /* array of length fdim */
  double errmax; /* max ee[k].err */
} region;

static region make_region(const hypercube *h, unsigned fdim)
{
  region R;
  R.h = make_hypercube(h->dim, h->data, h->data + h->dim);
  R.splitDim = 0;
  R.fdim = fdim;
  R.ee = R.h.data ? (esterr *) malloc(sizeof(esterr) * fdim) : NULL;
  return R;
}

static void destroy_region(region *R)
{
  destroy_hypercube(&R->h);
  free(R->ee);
  R->ee = 0;
}

static int cut_region(region *R, region *R2)
{
  unsigned d = R->splitDim, dim = R->h.dim;
  *R2 = *R;
  R->h.data[d + dim] *= 0.5;
  R->h.vol *= 0.5;
  R2->h = make_hypercube(dim, R->h.data, R->h.data + dim);
  if (!R2->h.data) return FAILURE;
  R->h.data[d] -= R->h.data[d + dim];
  R2->h.data[d] += R->h.data[d + dim];
  R2->ee = (esterr *) malloc(sizeof(esterr) * R2->fdim);
  return R2->ee == NULL;
}

struct rule_s; /* forward declaration */

typedef int (*evalError_func)(struct rule_s *r,
			      unsigned fdim, integrand_v f, void *fdata,
			      unsigned nR, region *R);
typedef void (*destroy_func)(struct rule_s *r);


typedef struct rule_s {
  unsigned dim, fdim;         /* the dimensionality & number of functions */
  unsigned num_points;       /* number of evaluation points */
  unsigned num_regions; /* max number of regions evaluated at once */
  double *pts; /* points to eval: num_regions * num_points * dim */
  double *vals; /* num_regions * num_points * fdim */
  evalError_func evalError;
  destroy_func destroy;
} rule;

static void destroy_rule(rule *r)
{
  if (r) {
    if (r->destroy) r->destroy(r);
    free(r->pts);
    free(r);
  }
}

static int alloc_rule_pts(rule *r, unsigned num_regions)
{
  if (num_regions > r->num_regions) {
    free(r->pts);
    r->pts = r->vals = NULL;
    r->num_regions = 0;
    num_regions *= 2; /* allocate extra so that
			 repeatedly calling alloc_rule_pts with
			 growing num_regions only needs
			 a logarithmic number of allocations */
    r->pts = (double *) malloc(sizeof(double) * 
			       (num_regions
				* r->num_points * (r->dim + r->fdim)));
    if (r->fdim + r->dim > 0 && !r->pts) return FAILURE;
    r->vals = r->pts + num_regions * r->num_points * r->dim;
    r->num_regions = num_regions;
  }
  return SUCCESS;
}

static rule *make_rule(size_t sz, /* >= sizeof(rule) */
		       unsigned dim, unsigned fdim, unsigned num_points,
		       evalError_func evalError, destroy_func destroy)
{
  rule *r;

  if (sz < sizeof(rule)) return NULL;
  r = (rule *) malloc(sz);
  if (!r) return NULL;
  r->pts = r->vals = NULL;
  r->num_regions = 0;
  r->dim = dim; r->fdim = fdim; r->num_points = num_points;
  r->evalError = evalError;
  r->destroy = destroy;
  return r;
}

/* note: all regions must have same fdim */
static int eval_regions(unsigned nR, region *R, 
			integrand_v f, void *fdata, rule *r)
{
  unsigned iR;
  if (nR == 0) return SUCCESS; /* nothing to evaluate */
  if (r->evalError(r, R->fdim, f, fdata, nR, R)) return FAILURE;
  for (iR = 0; iR < nR; ++iR)
    R[iR].errmax = errMax(R->fdim, R[iR].ee);
  return SUCCESS;
}

/***************************************************************************/
/* Functions to loop over points in a hypercube. */

/* Based on orbitrule.cpp in HIntLib-0.0.10 */

/* ls0 returns the least-significant 0 bit of n (e.g. it returns
   0 if the LSB is 0, it returns 1 if the 2 LSBs are 01, etcetera). */

static unsigned ls0(unsigned n)
{
#if defined(__GNUC__) &&					\
  ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ > 3)
  return __builtin_ctz(~n); /* gcc builtin for version >= 3.4 */
#else
  const unsigned bits[256] = {
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8,
  };
  unsigned bit = 0;
  while ((n & 0xff) == 0xff) {
    n >>= 8;
    bit += 8;
  }
  return bit + bits[n & 0xff];
#endif
}

/**
 *  Evaluate the integration points for all 2^n points (+/-r,...+/-r)
 *
 *  A Gray-code ordering is used to minimize the number of coordinate updates
 *  in p, although this doesn't matter as much now that we are saving all pts.
 */
static void evalR_Rfs(double *pts, unsigned dim, double *p, const double *c, const double *r)
{
  unsigned i;
  unsigned signs = 0; /* 0/1 bit = +/- for corresponding element of r[] */

  /* We start with the point where r is ADDed in every coordinate
     (this implies signs=0). */
  for (i = 0; i < dim; ++i)
    p[i] = c[i] + r[i];

  /* Loop through the points in Gray-code ordering */
  for (i = 0;; ++i) {
    unsigned mask, d;

    memcpy(pts, p, sizeof(double) * dim); pts += dim;

    d = ls0(i);	/* which coordinate to flip */
    if (d >= dim)
      break;

    /* flip the d-th bit and add/subtract r[d] */
    mask = 1U << d;
    signs ^= mask;
    p[d] = (signs & mask) ? c[d] - r[d] : c[d] + r[d];
  }
}

static void evalRR0_0fs(double *pts, unsigned dim, double *p, const double *c, const double *r)
{
  unsigned i, j;

  for (i = 0; i < dim - 1; ++i) {
    p[i] = c[i] - r[i];
    for (j = i + 1; j < dim; ++j) {
      p[j] = c[j] - r[j];
      memcpy(pts, p, sizeof(double) * dim); pts += dim;
      p[i] = c[i] + r[i];
      memcpy(pts, p, sizeof(double) * dim); pts += dim;
      p[j] = c[j] + r[j];
      memcpy(pts, p, sizeof(double) * dim); pts += dim;
      p[i] = c[i] - r[i];
      memcpy(pts, p, sizeof(double) * dim); pts += dim;

      p[j] = c[j];	/* Done with j -> Restore p[j] */
    }
    p[i] = c[i];		/* Done with i -> Restore p[i] */
  }
}

static void evalR0_0fs4d(double *pts, unsigned dim, double *p, const double *c,
			 const double *r1, const double *r2)
{
  unsigned i;

  memcpy(pts, p, sizeof(double) * dim); pts += dim;

  for (i = 0; i < dim; i++) {
    p[i] = c[i] - r1[i];
    memcpy(pts, p, sizeof(double) * dim); pts += dim;

    p[i] = c[i] + r1[i];
    memcpy(pts, p, sizeof(double) * dim); pts += dim;

    p[i] = c[i] - r2[i];
    memcpy(pts, p, sizeof(double) * dim); pts += dim;

    p[i] = c[i] + r2[i];
    memcpy(pts, p, sizeof(double) * dim); pts += dim;

    p[i] = c[i];
  }
}

#define num0_0(dim) (1U)
#define numR0_0fs(dim) (2 * (dim))
#define numRR0_0fs(dim) (2 * (dim) * (dim-1))
#define numR_Rfs(dim) (1U << (dim))

/***************************************************************************/
/* Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
   cubature rule of degree 7 (embedded rule degree 5) due to A. C. Genz
   and A. A. Malik.  See:

   A. C. Genz and A. A. Malik, "An imbedded [sic] family of fully
   symmetric numerical integration rules," SIAM
   J. Numer. Anal. 20 (3), 580-588 (1983).
*/

typedef struct {
  rule parent;

  /* temporary arrays of length dim */
  double *widthLambda, *widthLambda2, *p;

  /* dimension-dependent constants */
  double weight1, weight3, weight5;
  double weightE1, weightE3;
} rule75genzmalik;

#define real(x) ((double)(x))
#define to_int(n) ((int)(n))

static int isqr(int x)
{
  return x * x;
}

static void destroy_rule75genzmalik(rule *r_)
{
  rule75genzmalik *r = (rule75genzmalik *) r_;
  free(r->p);
}

static int rule75genzmalik_evalError(rule *r_, unsigned fdim, integrand_v f, void *fdata, unsigned nR, region *R)
{
  /* lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19) */
  const double lambda2 = 0.3585685828003180919906451539079374954541;
  const double lambda4 = 0.9486832980505137995996680633298155601160;
  const double lambda5 = 0.6882472016116852977216287342936235251269;
  const double weight2 = 980. / 6561.;
  const double weight4 = 200. / 19683.;
  const double weightE2 = 245. / 486.;
  const double weightE4 = 25. / 729.;
  const double ratio = (lambda2 * lambda2) / (lambda4 * lambda4);

  rule75genzmalik *r = (rule75genzmalik *) r_;
  unsigned i, j, iR, dim = r_->dim;
  size_t npts = 0;
  double *diff, *pts, *vals;

  if (alloc_rule_pts(r_, nR)) return FAILURE;
  pts = r_->pts; vals = r_->vals;

  for (iR = 0; iR < nR; ++iR) {
    const double *center = R[iR].h.data;
    const double *halfwidth = R[iR].h.data + dim;
	  
    for (i = 0; i < dim; ++i)
      r->p[i] = center[i];
	  
    for (i = 0; i < dim; ++i)
      r->widthLambda2[i] = halfwidth[i] * lambda2;
    for (i = 0; i < dim; ++i)
      r->widthLambda[i] = halfwidth[i] * lambda4;

    /* Evaluate points in the center, in (lambda2,0,...,0) and
       (lambda3=lambda4, 0,...,0).  */
    evalR0_0fs4d(pts + npts*dim, dim, r->p, center, 
		 r->widthLambda2, r->widthLambda);
    npts += num0_0(dim) + 2 * numR0_0fs(dim);

    /* Calculate points for (lambda4, lambda4, 0, ...,0) */
    evalRR0_0fs(pts + npts*dim, dim, r->p, center, r->widthLambda);
    npts += numRR0_0fs(dim);

    /* Calculate points for (lambda5, lambda5, ..., lambda5) */
    for (i = 0; i < dim; ++i)
      r->widthLambda[i] = halfwidth[i] * lambda5;
    evalR_Rfs(pts + npts*dim, dim, r->p, center, r->widthLambda);
    npts += numR_Rfs(dim);
  }

  /* Evaluate the integrand function(s) at all the points */
  if (f(dim, npts, pts, fdata, fdim, vals))
    return FAILURE;

  /* we are done with the points, and so we can re-use the pts
     array to store the maximum difference diff[i] in each dimension 
     for each hypercube */
  diff = pts;
  for (i = 0; i < dim * nR; ++i) diff[i] = 0;

  for (j = 0; j < fdim; ++j) {
    const double *v = vals + j;
#         define VALS(i) v[fdim*(i)]	       
    for (iR = 0; iR < nR; ++iR) {
      double result, res5th;
      double val0, sum2=0, sum3=0, sum4=0, sum5=0;
      unsigned k, k0 = 0;
      /* accumulate j-th function values into j-th integrals
	 NOTE: this relies on the ordering of the eval functions
	 above, as well as on the internal structure of
	 the evalR0_0fs4d function */

      val0 = VALS(0); /* central point */
      k0 += 1;

      for (k = 0; k < dim; ++k) {
	double v0 = VALS(k0 + 4*k);
	double v1 = VALS((k0 + 4*k) + 1);
	double v2 = VALS((k0 + 4*k) + 2);
	double v3 = VALS((k0 + 4*k) + 3);
		    
	sum2 += v0 + v1;
	sum3 += v2 + v3;
		    
	diff[iR * dim + k] += 
	  fabs(v0 + v1 - 2*val0 - ratio * (v2 + v3 - 2*val0));
      }
      k0 += 4*k;

      for (k = 0; k < numRR0_0fs(dim); ++k)
	sum4 += VALS(k0 + k);
      k0 += k;
	       
      for (k = 0; k < numR_Rfs(dim); ++k)
	sum5 += VALS(k0 + k);
	       
      /* Calculate fifth and seventh order results */
      result = R[iR].h.vol * (r->weight1 * val0 + weight2 * sum2 + r->weight3 * sum3 + weight4 * sum4 + r->weight5 * sum5);
      res5th = R[iR].h.vol * (r->weightE1 * val0 + weightE2 * sum2 + r->weightE3 * sum3 + weightE4 * sum4);
	       
      R[iR].ee[j].val = result;
      R[iR].ee[j].err = fabs(res5th - result);
	       
      v += r_->num_points * fdim;
    }
#         undef VALS
  }


  /* figure out dimension to split: */
  for (iR = 0; iR < nR; ++iR) {
    double maxdiff = 0;
    unsigned dimDiffMax = 0;

    for (i = 0; i < dim; ++i)
      if (diff[iR*dim + i] > maxdiff) {
	maxdiff = diff[iR*dim + i];
	dimDiffMax = i;
      }
    R[iR].splitDim = dimDiffMax;
  }
  return SUCCESS;
}

static rule *make_rule75genzmalik(unsigned dim, unsigned fdim)
{
  rule75genzmalik *r;

  if (dim < 2) return NULL; /* this rule does not support 1d integrals */

  /* Because of the use of a bit-field in evalR_Rfs, we are limited
     to be < 32 dimensions (or however many bits are in unsigned).
     This is not a practical limitation...long before you reach
     32 dimensions, the Genz-Malik cubature becomes excruciatingly
     slow and is superseded by other methods (e.g. Monte-Carlo). */
  if (dim >= sizeof(unsigned) * 8) return NULL;

  r = (rule75genzmalik *) make_rule(sizeof(rule75genzmalik),
				    dim, fdim,
				    num0_0(dim) + 2 * numR0_0fs(dim)
				    + numRR0_0fs(dim) + numR_Rfs(dim),
				    rule75genzmalik_evalError,
				    destroy_rule75genzmalik);
  if (!r) return NULL;

  r->weight1 = (real(12824 - 9120 * to_int(dim) + 400 * isqr(to_int(dim)))
		/ real(19683));
  r->weight3 = real(1820 - 400 * to_int(dim)) / real(19683);
  r->weight5 = real(6859) / real(19683) / real(1U << dim);
  r->weightE1 = (real(729 - 950 * to_int(dim) + 50 * isqr(to_int(dim)))
		 / real(729));
  r->weightE3 = real(265 - 100 * to_int(dim)) / real(1458);

  r->p = (double *) malloc(sizeof(double) * dim * 3);
  if (!r->p) { destroy_rule((rule *) r); return NULL; }
  r->widthLambda = r->p + dim;
  r->widthLambda2 = r->p + 2 * dim;

  return (rule *) r;
}

/***************************************************************************/
/* 1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in
   GNU GSL (which in turn is based on QUADPACK). */

static int rule15gauss_evalError(rule *r,
				 unsigned fdim, integrand_v f, void *fdata,
				 unsigned nR, region *R)
{
  /* Gauss quadrature weights and kronrod quadrature abscissae and
     weights as evaluated with 80 decimal digit arithmetic by
     L. W. Fullerton, Bell Labs, Nov. 1981. */
  const unsigned n = 8;
  const double xgk[8] = {  /* abscissae of the 15-point kronrod rule */
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.000000000000000000000000000000000
    /* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
       xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule */
  };
  static const double wg[4] = {  /* weights of the 7-point gauss rule */
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327
  };
  static const double wgk[8] = { /* weights of the 15-point kronrod rule */
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714
  };
  unsigned j, k, iR;
  size_t npts = 0;
  double *pts, *vals;

  if (alloc_rule_pts(r, nR)) return FAILURE;
  pts = r->pts; vals = r->vals;

  for (iR = 0; iR < nR; ++iR) {
    const double center = R[iR].h.data[0];
    const double halfwidth = R[iR].h.data[1];

    pts[npts++] = center;

    for (j = 0; j < (n - 1) / 2; ++j) {
      int j2 = 2*j + 1;
      double w = halfwidth * xgk[j2];
      pts[npts++] = center - w;
      pts[npts++] = center + w;
    }
    for (j = 0; j < n/2; ++j) {
      int j2 = 2*j;
      double w = halfwidth * xgk[j2];
      pts[npts++] = center - w;
      pts[npts++] = center + w;
    }

    R[iR].splitDim = 0; /* no choice but to divide 0th dimension */
  }

  if (f(1, npts, pts, fdata, fdim, vals))
    return FAILURE;
     
  for (k = 0; k < fdim; ++k) {
    const double *vk = vals + k;
    for (iR = 0; iR < nR; ++iR) {
      const double halfwidth = R[iR].h.data[1];
      double result_gauss = vk[0] * wg[n/2 - 1];
      double result_kronrod = vk[0] * wgk[n - 1];
      double result_abs = fabs(result_kronrod);
      double result_asc, mean, err;

      /* accumulate integrals */
      npts = 1;
      for (j = 0; j < (n - 1) / 2; ++j) {
	int j2 = 2*j + 1;
	double v = vk[fdim*npts] + vk[fdim*npts+fdim];
	result_gauss += wg[j] * v;
	result_kronrod += wgk[j2] * v;
	result_abs += wgk[j2] * (fabs(vk[fdim*npts]) 
				 + fabs(vk[fdim*npts+fdim]));
	npts += 2;
      }
      for (j = 0; j < n/2; ++j) {
	int j2 = 2*j;
	result_kronrod += wgk[j2] * (vk[fdim*npts] 
				     + vk[fdim*npts+fdim]);
	result_abs += wgk[j2] * (fabs(vk[fdim*npts]) 
				 + fabs(vk[fdim*npts+fdim]));
	npts += 2;
      }
	       
      /* integration result */
      R[iR].ee[k].val = result_kronrod * halfwidth;

      /* error estimate 
	 (from GSL, probably dates back to QUADPACK
	 ... not completely clear to me why we don't just use
	 fabs(result_kronrod - result_gauss) * halfwidth */
      mean = result_kronrod * 0.5;
      result_asc = wgk[n - 1] * fabs(vk[0] - mean);
      npts = 1;
      for (j = 0; j < (n - 1) / 2; ++j) {
	int j2 = 2*j + 1;
	result_asc += wgk[j2] * (fabs(vk[fdim*npts]-mean)
				 + fabs(vk[fdim*npts+fdim]-mean));
	npts += 2;
      }
      for (j = 0; j < n/2; ++j) {
	int j2 = 2*j;
	result_asc += wgk[j2] * (fabs(vk[fdim*npts]-mean)
				 + fabs(vk[fdim*npts+fdim]-mean));
	npts += 2;
      }
      err = fabs(result_kronrod - result_gauss) * halfwidth;
      result_abs *= halfwidth;
      result_asc *= halfwidth;
      if (result_asc != 0 && err != 0) {
	double scale = pow((200 * err / result_asc), 1.5);
	err = (scale < 1) ? result_asc * scale : result_asc;
      }
      if (result_abs > DBL_MIN / (50 * DBL_EPSILON)) {
	double min_err = 50 * DBL_EPSILON * result_abs;
	if (min_err > err) err = min_err;
      }
      R[iR].ee[k].err = err;
	       
      /* increment vk to point to next batch of results */
      vk += 15*fdim;
    }
  }
  return SUCCESS;
}

static rule *make_rule15gauss(unsigned dim, unsigned fdim)
{
  if (dim != 1) return NULL; /* this rule is only for 1d integrals */

  return make_rule(sizeof(rule), dim, fdim, 15,
		   rule15gauss_evalError, 0);
}

/***************************************************************************/
/* binary heap implementation (ala _Introduction to Algorithms_ by
   Cormen, Leiserson, and Rivest), for use as a priority queue of
   regions to integrate. */

typedef region heap_item;
#define KEY(hi) ((hi).errmax)

typedef struct {
  size_t n, nalloc;
  heap_item *items;
  unsigned fdim;
  esterr *ee; /* array of length fdim of the total integrand & error */
} heap;

static void heap_resize(heap *h, size_t nalloc)
{
  h->nalloc = nalloc;
  h->items = (heap_item *) realloc(h->items, sizeof(heap_item) * nalloc);
}

static heap heap_alloc(size_t nalloc, unsigned fdim)
{
  heap h;
  unsigned i;
  h.n = 0;
  h.nalloc = 0;
  h.items = 0;
  h.fdim = fdim;
  h.ee = (esterr *) malloc(sizeof(esterr) * fdim);
  if (h.ee) {
    for (i = 0; i < fdim; ++i) h.ee[i].val = h.ee[i].err = 0;
    heap_resize(&h, nalloc);
  }
  return h;
}

/* note that heap_free does not deallocate anything referenced by the items */
static void heap_free(heap *h)
{
  h->n = 0;
  heap_resize(h, 0);
  h->fdim = 0;
  free(h->ee);
}

static int heap_push(heap *h, heap_item hi)
{
  int insert;
  unsigned i, fdim = h->fdim;

  for (i = 0; i < fdim; ++i) {
    h->ee[i].val += hi.ee[i].val;
    h->ee[i].err += hi.ee[i].err;
  }
  insert = h->n;
  if (++(h->n) > h->nalloc) {
    heap_resize(h, h->n * 2);
    if (!h->items) return FAILURE;
  }

  while (insert) {
    int parent = (insert - 1) / 2;
    if (KEY(hi) <= KEY(h->items[parent]))
      break;
    h->items[insert] = h->items[parent];
    insert = parent;
  }
  h->items[insert] = hi;
  return SUCCESS;
}

static int heap_push_many(heap *h, size_t ni, heap_item *hi)
{
  size_t i;
  for (i = 0; i < ni; ++i)
    if (heap_push(h, hi[i])) return FAILURE;
  return SUCCESS;
}

static heap_item heap_pop(heap *h)
{
  heap_item ret;
  int i, n, child;

  if (!(h->n)) {
    // Yumi Erases the following Dec 6 2014 as R CMD gives err
    // fprintf(stderr, "attempted to pop an empty heap\n");
    // exit(EXIT_FAILURE);
  }

  ret = h->items[0];
  h->items[i = 0] = h->items[n = --(h->n)];
  while ((child = i * 2 + 1) < n) {
    int largest;
    heap_item swap;

    if (KEY(h->items[child]) <= KEY(h->items[i]))
      largest = i;
    else
      largest = child;
    if (++child < n && KEY(h->items[largest]) < KEY(h->items[child]))
      largest = child;
    if (largest == i)
      break;
    swap = h->items[i];
    h->items[i] = h->items[largest];
    h->items[i = largest] = swap;
  }

  {
    unsigned i, fdim = h->fdim;
    for (i = 0; i < fdim; ++i) {
      h->ee[i].val -= ret.ee[i].val;
      h->ee[i].err -= ret.ee[i].err;
    }
  }
  return ret;
}

/***************************************************************************/

static int converged(unsigned fdim, const esterr *ee,
		     double reqAbsError, double reqRelError, error_norm norm)
#define ERR(j) ee[j].err
#define VAL(j) ee[j].val
#include "converged.h"

/***************************************************************************/

/* adaptive integration, analogous to adaptintegrator.cpp in HIntLib */

  static int rulecubature(rule *r, unsigned fdim, 
			  integrand_v f, void *fdata, 
			  const hypercube *h, 
			  size_t maxEval,
			  double reqAbsError, double reqRelError,
			  error_norm norm,
			  double *val, double *err, int parallel)
{
  size_t numEval = 0;
  heap regions;
  unsigned i, j;
  region *R = NULL; /* array of regions to evaluate */
  size_t nR_alloc = 0;
  esterr *ee = NULL;

  if (fdim <= 1) norm = ERROR_INDIVIDUAL; /* norm is irrelevant */
  //if (norm < 0 || norm > ERROR_LINF) return FAILURE; /* invalid norm */

  regions = heap_alloc(1, fdim);
  if (!regions.ee || !regions.items) goto bad;

  ee = (esterr *) malloc(sizeof(esterr) * fdim);
  if (!ee) goto bad;
     
  nR_alloc = 2;
  R = (region *) malloc(sizeof(region) * nR_alloc);
  if (!R) goto bad;
  R[0] = make_region(h, fdim);
  if (!R[0].ee
      || eval_regions(1, R, f, fdata, r)
      || heap_push(&regions, R[0]))
    goto bad;
  numEval += r->num_points;
     
  while (numEval < maxEval || !maxEval) {
    if (converged(fdim, regions.ee, reqAbsError, reqRelError, norm))
      break;

    if (parallel) { /* maximize potential parallelism */
      /* adapted from I. Gladwell, "Vectorization of one
	 dimensional quadrature codes," pp. 230--238 in
	 _Numerical Integration. Recent Developments,
	 Software and Applications_, G. Fairweather and
	 P. M. Keast, eds., NATO ASI Series C203, Dordrecht
	 (1987), as described in J. M. Bull and
	 T. L. Freeman, "Parallel Globally Adaptive
	 Algorithms for Multi-dimensional Integration,"
	 http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638
	 (1994). 

	 Basically, this evaluates in one shot all regions
	 that *must* be evaluated in order to reduce the
	 error to the requested bound: the minimum set of
	 largest-error regions whose errors push the total
	 error over the bound.

	 [Note: Bull and Freeman claim that the Gladwell
	 approach is intrinsically inefficent because it
	 "requires sorting", and propose an alternative
	 algorithm that "only" requires three passes over the
	 entire set of regions.  Apparently, they didn't
	 realize that one could use a heap data structure, in
	 which case the time to pop K biggest-error regions
	 out of N is only O(K log N), much better than the
	 O(N) cost of the Bull and Freeman algorithm if K <<
	 N, and it is also much simpler.] */
      size_t nR = 0;
      for (j = 0; j < fdim; ++j) ee[j] = regions.ee[j];
      do {
	if (nR + 2 > nR_alloc) {
	  nR_alloc = (nR + 2) * 2;
	  R = (region *) realloc(R, nR_alloc * sizeof(region));
	  if (!R) goto bad;
	}
	R[nR] = heap_pop(&regions);
	for (j = 0; j < fdim; ++j) ee[j].err -= R[nR].ee[j].err;
	if (cut_region(R+nR, R+nR+1)) goto bad;
	numEval += r->num_points * 2;
	nR += 2;
	if (converged(fdim, ee, reqAbsError, reqRelError, norm))
	  break; /* other regions have small errs */
      } while (regions.n > 0 && (numEval < maxEval || !maxEval));
      if (eval_regions(nR, R, f, fdata, r)
	  || heap_push_many(&regions, nR, R))
	goto bad;
    }
    else { /* minimize number of function evaluations */
      R[0] = heap_pop(&regions); /* get worst region */
      if (cut_region(R, R+1)
	  || eval_regions(2, R, f, fdata, r)
	  || heap_push_many(&regions, 2, R))
	goto bad;
      numEval += r->num_points * 2;
    }
  }

  /* re-sum integral and errors */
  for (j = 0; j < fdim; ++j) val[j] = err[j] = 0;  
  for (i = 0; i < regions.n; ++i) {
    for (j = 0; j < fdim; ++j) { 
      val[j] += regions.items[i].ee[j].val;
      err[j] += regions.items[i].ee[j].err;
    }
    destroy_region(&regions.items[i]);
  }

  /* printf("regions.nalloc = %d\n", regions.nalloc); */
  free(ee);
  heap_free(&regions);
  free(R);
  return SUCCESS;

 bad:
  free(ee);
  heap_free(&regions);
  free(R);
  return FAILURE;
}

static int cubature(unsigned fdim, integrand_v f, void *fdata, 
		    unsigned dim, const double *xmin, const double *xmax, 
		    size_t maxEval, double reqAbsError, double reqRelError, 
		    error_norm norm,
		    double *val, double *err, int parallel)
{
  rule *r;
  hypercube h;
  int status;
  unsigned i;
     
  if (fdim == 0) /* nothing to do */ return SUCCESS;
  if (dim == 0) { /* trivial integration */
    if (f(0, 1, xmin, fdata, fdim, val)) return FAILURE;
    for (i = 0; i < fdim; ++i) err[i] = 0;
    return SUCCESS;
  }
  r = dim == 1 ? make_rule15gauss(dim, fdim)
    : make_rule75genzmalik(dim, fdim);
  if (!r) { 
    for (i = 0; i < fdim; ++i) {
      val[i] = 0;
      err[i] = HUGE_VAL; 
    }
    return FAILURE;
  }
  h = make_hypercube_range(dim, xmin, xmax);
  status = !h.data ? FAILURE
    : rulecubature(r, fdim, f, fdata, &h,
		   maxEval, reqAbsError, reqRelError, norm,
		   val, err, parallel);
  destroy_hypercube(&h);
  destroy_rule(r);
  return status;
}

int hcubature_v(unsigned fdim, integrand_v f, void *fdata, 
                unsigned dim, const double *xmin, const double *xmax, 
                size_t maxEval, double reqAbsError, double reqRelError, 
                error_norm norm,
                double *val, double *err)
{
  return cubature(fdim, f, fdata, dim, xmin, xmax, 
		  maxEval, reqAbsError, reqRelError, norm, val, err, 1);
}

#include "vwrapper.h"

int hcubature(unsigned fdim, integrand f, void *fdata, 
	      unsigned dim, const double *xmin, const double *xmax, 
	      size_t maxEval, double reqAbsError, double reqRelError, 
	      error_norm norm,
	      double *val, double *err)
{
  int ret;
  fv_data d;

  if (fdim == 0) return SUCCESS; /* nothing to do */     
     
  d.f = f; d.fdata = fdata;
  ret = cubature(fdim, fv, &d, dim, xmin, xmax, 
		 maxEval, reqAbsError, reqRelError, norm, val, err, 0);
  return ret;
}

/***************************************************************************/
// Following codes are written by Yumi Kondo to fit Negative Binomial mixed-effect model and compute
// the conditional probability index.


double getrij(
	      const int ivec, 
	      SEXP X,  // reference even w/o &
	      const int wID[],  // reference even w/o &
	      const double beta[],
	      const int Ntot, // = sum_i ni
	      const int p
	      )// reference even w/o &
{
  // This function modifies the values of sizeij to be sizeij = beta0+Xi1*beta1+Xi2*beta2+...:

  // wID:     A vector of length max(ni), containing the position index to indicate the locations of repeated measures 
  //          of the selected patient in the vector ID in the first ni entries. The rest of max(ni) - ni + 1 entries are junks
  // ivec:    A scalar, indicating the ivec^th repeated measure

  double sizeij = 0.0;
  
  for (int ip = 0 ; ip < p ; ip ++ ) 
    sizeij +=  REAL(X)[ wID[ivec]+Ntot*ip ]*beta[ip];

  return(sizeij);

  // sizeij = beta0+Xi1*beta1+Xi2*beta2+...:
}



void getwID(int *wIDsize, 
	    int wID[], // wID is treated as pointer even w/o 
	    SEXP ID, 
	    const int ipat)
{
  // This function modifies the values of wIDsize and wID:
  
  // ID:      A vector of length sum ni, patient labeling starts from ZERO 
  // ipat:    A scalar the patient index that you are interested in (0 to N-1)
  *wIDsize = 0; 
  for (int ivec=0; ivec < length(ID) ;ivec++)
    {
      if( INTEGER(ID)[ivec]==ipat) 
	{
	  wID[*wIDsize] = ivec; // location at which ipat is recorded
	  (*wIDsize)++; 
	}
    }
  //Rprintf("\v ipat %d wIDsize %d  :",ipat,*wIDsize);
  //for (int ii = 0; ii < *wIDsize; ii++) Rprintf(" %d ",wID[ii]);
}


double intJacob1(const double x){
  //double x2 = R_pow(x,2);
  //return((1+x2)/R_pow(1-x2,2));
  return(1/R_pow(1-x,2));
}
double getProb(double gY, double alpha)
{
  return( 1/(gY*alpha + 1));
}
double intChangeOfVar(const double x)
{
  // Change of variable to integrate from -Inf to Inf
  return(1/(1-x)-1);//x/(1-R_pow(x,2)));//1/(1-x[0])-1
}


double dbetaBinom(const double x, const int N, const double al, const double bet,const int logTF)
{
  double temp = lbeta(x+al, N-x+bet)-lbeta(al,bet) + lchoose(N,x);
  if (logTF) return (temp);
  else return (exp(temp));
}

double distRE(const double g,const double theta,const int idistRE)
{
    if (idistRE == 1)
      {
	return(dgamma(g,1/theta,theta,0));
      }else if (idistRE == 2){
      return(dlnorm(g,-log(theta+1)/2,sqrt(log(theta+1)),0));
    }else{
      return(R_NegInf);
    }
}


// Global variables to define the locations of objects in fd
const int Np = 7;
// parameters for Y CEL
const int l_alpha = 0;
const int l_thY = 1;
const int l_delta = 2;
const int l_maxni = 3;
const int l_distRE = 4;
const int l_wIDsizeY = 5;
const int l_typeSummary = 6;
//If maxni = 8 then: 
//rYij: 6-13
//Yij: 14-21
//diff:22-29
//lnpY: 30-37
// ivec must be 0 <= ivec < maxni
int get_rY(const int ivec, const int maxni){  return(ivec + Np);}
int get_Y(const int ivec, const int maxni){   return(ivec + Np + maxni);}
int get_diff(const int ivec,const int maxni){ return(ivec + Np + maxni*2); }
//============ The following is required only when computing the CPI =============
int get_lnpY(const int ivec, const int maxni){return(ivec + Np + maxni*3); }
// The size of fd to compute likelihood 
int get_fdsize(const int maxni){ return( Np + maxni*4 ); }

const double xmin[1] = {0}, xmax[1] = {1};



double ifelse(const int TF, const double Tout,const double Fout)
{
  if (TF){ return(Tout); }else{ return(Fout); }
}

double dnbinomYK(double x, double size, double prob, int give_log)
{
  const double small = 1E-8;
  if (prob <small){ 
    return(ifelse(give_log,R_NegInf,0));
  }else return(dnbinom(x, size, prob, give_log));
}

double pnbinomYK(double x, double size, double prob, int give_log)
{
  const double small = 1E-8;
  if (prob <small){ 
    return(ifelse(give_log,R_NegInf,0));
  }else return(pnbinom(x, size, prob,1, give_log)); // pnbinom(double x, double size, double prob, int lower_tail, int log_p)
}

double PrYijGivenYij_1AndGYAR(double Yik,double Yik_1,
			      double sizeYik,double sizeYik_1,
			      double delta,double prob,int logTF,int CDF)
{
  // Yik | Yi(k-1), GY if CDF = 1 then compute Pr(Yik<=yik|Yi(k-1), GY ) else compute Pr(Yik = yik|Yi(k-1), GY )
  // under AR(1) model
  int support_max = ifelse(Yik < Yik_1,Yik,Yik_1);//imin2(Yik,Yik_1);
  double fYijGivenYik_1 = 0, temp;
  double sizeNB = sizeYik - delta*sizeYik_1;

  //Rprintf("\n\n support_max %d Yik %f Yik_1 %f",support_max,Yik,Yik_1);
  if (Yik <0 || sizeNB <= 0){
    fYijGivenYik_1 = 0;
  }else{
  
    for (int t = 0; t <= support_max;t++ )
      {
	R_CheckUserInterrupt();
	// (double x, double size, double prob, int give_log)
	if (CDF){
	  temp = pnbinom(Yik - t,sizeNB,prob,1,0);
	}else{
	  temp = dnbinom(Yik - t,sizeNB,prob,0);
	}

        fYijGivenYik_1 += dbetaBinom(t,Yik_1,
				     sizeYik_1*delta,
				     sizeYik_1*(1-delta),0)*temp;
	//Rprintf("\n t %d fYijGivenYik_1 %f Yik_1 %f",
	//	t,fYijGivenYik_1,t,Yik_1);
        
      }
  }
  if (logTF){
    return(log(fYijGivenYik_1));
  }else return(fYijGivenYik_1);
}






double densYijGivenYij_1AndGY(double Yik, double Yik_1, double sizeYik, double sizeYik_1,
			      double delta,double prob, int logTF)
{
  // Compute Pr(Yij=yij|Yij,GY)
  double log_fYijGivenYik_1;
  if ( 0 < delta && delta < 1 && sizeYik_1 >= 0 && Yik_1 >= 0){ 
    // AR(1) 0 < delta < 1
    log_fYijGivenYik_1 = PrYijGivenYij_1AndGYAR(Yik,Yik_1,
						sizeYik,sizeYik_1,
						delta,prob,1,0);
  }else{ 
    // IND <=> delta < 0
    // sizeYik_1 == R_NegInf <=> first scan (could be AR(1))
    log_fYijGivenYik_1 = dnbinomYK(Yik,sizeYik,prob,1);
    //Rprintf("\n log_dnbinomYK %f",log_fYijGivenYik_1);
  }
  if (logTF){ return(log_fYijGivenYik_1); }else{ return(exp(log_fYijGivenYik_1));}
}


int intLik(unsigned ndim, const double *x, void *fdata,
	   unsigned fdim, double *fval)
{
  double *fd = (double*) fdata;
  const double alpha = fd[l_alpha], theta = fd[l_thY], delta = fd[l_delta];
  const int maxni = (int) fd[l_maxni];
  const int RE = (int) fd[l_distRE];
  const int wIDsize = (int) fd[l_wIDsizeY];

  const double g = intChangeOfVar(x[0]), prob = getProb(g,alpha);
  double rij_1, rij, lik, delta_dif,dif, yij_1, yij;
  // Random effect part
  lik = distRE(g, theta, RE);
 

  // 2nd to the final observation of single patient
  for (int ivec = 0 ; ivec < wIDsize ; ivec++)
    {
      R_CheckUserInterrupt();
      rij = fd[get_rY(ivec,maxni)];
      rij_1 = ifelse(ivec==0,R_NegInf,fd[get_rY(ivec,maxni) - 1]);
      yij = fd[get_Y(ivec,maxni)];
      yij_1 = ifelse(ivec==0,R_NegInf,fd[get_Y(ivec,maxni)- 1]);
      dif = fd[get_diff(ivec,maxni)];
      delta_dif = pow(delta,dif); // delta_dif is not used for the computation with ivec = 0
      lik *=densYijGivenYij_1AndGY(yij, yij_1,rij, rij_1,delta_dif,prob, 0);
      //Rprintf("\n ivec %d lik %f delta_dif %f ",ivec,lik,delta_dif);
    }

  fval[0]=lik*intJacob1(x[0]);
  //Rprintf("\nfval[0] %f",fval[0]);
  return 0;
  }



SEXP nllk_C(
	    SEXP Y_, 
	    SEXP X, 
	    
	    SEXP ID,
	    SEXP alpha_, 
	    SEXP theta_,
	    SEXP delta_,
	    SEXP betas_, 
	    SEXP maxni_, 
	    SEXP NpatTot_, // The number of patients
	    SEXP RE,
	    SEXP DIF_, // length(dif) == length(Y)
	    SEXP absTol_
	    )
{
  const double *Y = REAL(Y_);
  const double *dif = REAL(DIF_); // Not used if delta == 0
  const int maxni = INTEGER(maxni_)[0];
  const double absTol = REAL(absTol_)[0];
  double fdata[get_fdsize(maxni)];
  fdata[l_alpha]= REAL(alpha_)[0];
  fdata[l_thY]= REAL(theta_)[0];
  fdata[l_delta]= REAL(delta_)[0];
  fdata[l_maxni]= maxni;
  fdata[l_distRE]= INTEGER(RE)[0]; // RE = 1 if "G" and RE = 2 if "LN"
  //fdata[Np- 1]= wIDsize;

  const double *betas = REAL(betas_);  
  const int NpatTot = INTEGER(NpatTot_)[0];
  const int Ntot = length(ID); // 
  const int p = length(X)/Ntot;

  GetRNGstate();
  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  
  SEXP res = PROTECT(allocVector(VECSXP, 1)); // result is stored in res
  SEXP nllk = allocVector(REALSXP, 1); 
  SET_VECTOR_ELT(res, 0, nllk); 

 
  double val,err;

  int ipat,ivec,temp;
  int wID[maxni],wIDsize = 0,wIDs[NpatTot][maxni], wIDsizes[NpatTot];


  for (ipat = 0 ; ipat < NpatTot ; ipat++)
    {
      getwID( &wIDsize,wID, ID, ipat);
      // record the number of scans for each patient
      wIDsizes[ipat] = wIDsize;
      //Rprintf("\n ipat %d ni %d \n",ipat, wIDsizes[ipat]);
      for (ivec = 0 ; ivec < wIDsize ; ivec++) wIDs[ipat][ivec]=wID[ivec];
      for (ivec = wIDsize ; ivec < maxni; ivec++) wIDs[ipat][ivec]=-1000;

    }

  REAL(nllk)[0] =0;
  for (ipat = 0 ; ipat < NpatTot; ipat++)
    {
      fdata[l_wIDsizeY] = wIDsizes[ipat];
      R_CheckUserInterrupt(); 
      
      for (ivec = 0 ; ivec < wIDsizes[ipat]; ivec++)
	{
	  fdata[get_rY(ivec,maxni)] = exp(getrij(ivec,X, wIDs[ipat],betas, Ntot,p))/REAL(alpha_)[0];
	  fdata[get_Y(ivec,maxni)] = Y[wIDs[ipat][ivec]];
	  fdata[get_diff(ivec,maxni)] = dif[wIDs[ipat][ivec]];
	}
      // Rprintf("\n\n\n ipat %d \n\n",ipat);
      //for (ivec = 0 ; ivec < get_fdsize(maxni); ivec++) 
      //Rprintf("\n fd[%d] %f",ivec,fdata[ivec]);
      
      temp = hcubature(1, intLik, fdata, 1, xmin, xmax, 0, 
		       absTol, // absolute tolerance
		       0, // relative tolerance (ignored)
		       ERROR_INDIVIDUAL, &val, &err);
      //Rprintf("\n  pat %d lk %1.3f converge %1.3f err %f",ipat, val, temp,err);
      
      REAL(nllk)[0] -= log(val);
      
      //Rprintf("\n nllk %1.3f",REAL(nllk)[0]);
    }
  
 
  PutRNGstate();
  UNPROTECT(1);

  return res;
}





int getCombForMax(double Yi[],const int NYi, const double K)
{
  const int NYi2 = NYi-2; // The location of the second last CEL 
  int changePos=0, ivec;
  for (ivec = 0; ivec <= NYi2 ;ivec++)
    {
      changePos=NYi2-ivec;
      if (Yi[changePos] < K){
	Yi[changePos]++;
	break;
      }else if (Yi[changePos] == K){
	Yi[changePos] = 0;
      }else{
	error("Yi should not be greater than %f",K);
      }
    }
  if (ivec==NYi2+1){
    for (ivec = 0 ; ivec <= (NYi2+1); ivec++)Yi[ivec] = R_NegInf;
  }
  return(changePos);
  // The code above returns the combination of all Y1,Y2,..,Yn-1 such that max(Y1)<=k...max(Yn-1)<=k and Yn=k
  // in certain order.
  // For example, If n = 4 and K = 3 then the returned patterns are:  
  // 0 0 0 3
  // 0 0 1 3 
  // 0 0 2 3 
  // 0 0 3 3
  // 0 1 0 3
  // 0 1 1 3
  // 0 1 2 3 
  // 0 1 3 3
  // 0 2 0 3
  // 0 2 1 3
  // 0 2 2 3
  // 0 2 3 3 
  // 0 3 0 3 
  // 0 3 1 3 
  // 0 3 2 3 
  // 0 3 3 3 

  // 1 0 0 3
  // 1 0 1 3 
  // 1 0 2 3 
  // 1 0 3 3
  // 1 1 0 3
  // 1 1 1 3
  // 1 1 2 3 
  // 1 1 3 3
  // 1 2 0 3
  // 1 2 1 3
  // 1 2 2 3
  // 1 2 3 3 
  // 1 3 0 3 
  // 1 3 1 3 
  // 1 3 2 3 
  // 1 3 3 3 

  // 2 0 0 3
  // 2 0 1 3 
  // 2 0 2 3 
  // 2 0 3 3
  // 2 1 0 3
  // 2 1 1 3
  // 2 1 2 3 
  // 2 1 3 3
  // 2 2 0 3
  // 2 2 1 3
  // 2 2 2 3
  // 2 2 3 3 
  // 2 3 0 3 
  // 2 3 1 3 
  // 2 3 2 3 
  // 2 3 3 3 

  // 3 0 0 3
  // 3 0 1 3 
  // 3 0 2 3 
  // 3 0 3 3
  // 3 1 0 3
  // 3 1 1 3
  // 3 1 2 3 
  // 3 1 3 3
  // 3 2 0 3
  // 3 2 1 3
  // 3 2 2 3
  // 3 2 3 3 
  // 3 3 0 3 
  // 3 3 1 3 
  // 3 3 2 3 
  // 3 3 3 3 
  // In total there are (K+1)^(n-1) = (3+1)^(4-1)=64 combinations
 }



int getCombForSum(double Yi[],const int NYi, const double K)
{
  const int NYi1 = NYi-1; // The location of the final CEL 
  // Step 1: Check the value of the final CEL is zero
  int changePos=0;
  if (Yi[0]==K){
    for (int ik = 0; ik < NYi; ik++)  Yi[ik]=R_NegInf;
  }else{
    if (Yi[NYi1]==0)
      {
	for (int ivec = 1; ivec < NYi1; ivec++)
	  {
	    //Rprintf("\n\n NYi1:%d  ivec:%d Yi[%d]:%f",
	    //NYi1,ivec,NYi1-ivec,Yi[NYi1-ivec]);
	    if (Yi[NYi1-ivec]>0){
	      
	      Yi[NYi1-ivec]=0;
	      changePos = NYi1-ivec-1;
	      Yi[changePos]++;
	      Yi[NYi1] = K;
	      for (int ik = 0; ik < NYi1; ik++)  Yi[NYi1] -= Yi[ik];
	      break;
	    }
	  }
      }else{
      // Step 2: If The final CEL is greater than zero, 
      // then add value to the second last value
      changePos=NYi1-1;
      Yi[changePos]++;
      Yi[NYi1]--;
    }
    return(changePos);
  }
  return(R_NegInf); // error
  // The code above returns the combination of all Y1,Y2,..,Yn such that Y1+...+Yn = K in certain order.
  // For example, If n = 4 and K = 3 then the returned patterns are:  
  // 0 0 0 3
  // 0 0 1 2 
  // 0 0 2 1 
  // 0 0 3 0
  // 0 1 0 2
  // 0 1 1 1
  // 0 1 2 0 
  // 0 2 0 1
  // 0 2 1 0
  // 0 3 0 0
  // 1 0 0 2 
  // 1 0 1 1
  // 1 0 2 0
  // 1 1 0 1
  // 1 1 1 0
  // 1 2 0 0
  // 2 0 0 1
  // 2 0 1 0
  // 2 1 0 0
  // 3 0 0 0
  // In total, there is choose(n + K-1, K) combinations.
  // Given precious combination Yi and K, the combinations are returned systematically as follows:
  // Step0: If Yi[0]==K, then this is the final pattern. Return a vector containing -Inf. Else go to the next Step1.
  // Step1: If the final column is > 0 then increase the second final column by one, and decrease the final column by one. 
  //        If the final column is 0, then set ivec=2 and go to Step2.
  // Step2: If the ivec^th last is > 0 then set it zero and increase ivec+1^th last column by one. 
  //        Adjust the final column so that the final column = K - sum of all the previous columns 
  //        If the ivec^th last is zero, then increase ivec by one and repeat Step 2.
  //     
}
SEXP getC(SEXP Yi, SEXP K)
{
  double *yi= REAL(Yi);
  const int Ny = length(Yi); 

  SEXP res = PROTECT(allocVector(VECSXP, 2)); // result is stored in res
  SEXP resYi = allocVector(REALSXP, Ny); 
  SET_VECTOR_ELT(res, 0, resYi);
  SEXP changePos = allocVector(INTSXP, 1); 
  SET_VECTOR_ELT(res, 1, changePos); 
  GetRNGstate();
  INTEGER(changePos)[0] = getCombForMax(yi, Ny,REAL(K)[0]);
  for (int iy=0 ; iy < Ny; iy++)
    {
      REAL(resYi)[iy]=yi[iy];
    }
  PutRNGstate();
  UNPROTECT(1);

  return res;
}




int get_NComb(const int wIDsizeYFol, const int K, const int type) // K = q(Yfol) - 1
{
  if (type==1){ // sum
    return( choose(K+wIDsizeYFol-1,K));
  }else if (type ==2 ){ // max
    return( pow(K+1, // + 1 because Yi can take 0
		wIDsizeYFol-1)); // -1 because the final scan is fixed to K
  }else{
    error("Type must be 1 (sum) or 2 (max)");
  }
  return(0);
}

int getComb(double Yi[],const int NYi, const double K, const int typeSummary)
{
  if (typeSummary == 1) return(getCombForSum(Yi,NYi,K));
  else if  (typeSummary == 2) return(getCombForMax(Yi,NYi,K));
  else error("typeSummary must be 1 (=sum) or 2 (=max)!");
}
int init_getComb(double vec[],const int lvec,const double K)
{
  // **Note** The first combination is the same for both getCombForMax and getCombForSum
  // lvec = length(vec) 
  // Initializing vec={0,...,0,K}
  for (int ivec=0; ivec < (lvec-1);ivec++) vec[ivec] = 0;
  vec[lvec-1] = K;
  return(0);
}
















double cdfqYfolGivenYPreAndG(const double fd[], const double K,
			     const double rYFolSum,const double probY) // K = q(Yfol) -1
{
  // Pr (Yfol+ <= K | Ypre=ypre, G=g) or Pr (max(Yfol) <= max(yfol) | Ypre=ypre, G=g) 
  double denYFol=0,lnpY;
  const double delta = fd[l_delta];
  const int wIDsizeY = fd[l_wIDsizeY];
  const int typeSummary = (int) fd[l_typeSummary],maxni = fd[l_maxni];
  if (delta==0 ){
    // IND model
    switch(typeSummary)
      {
      case 1:// summary
	// Pr (Yfol+ <= K | Ypre=ypre, G=g) =  Pr (Yfol+ <= K | G=g) 
	denYFol = pnbinomYK(K,rYFolSum,probY,0);
	break;
      case 2:
	// Pr (maxYfol <= K | Ypre=ypre, G=g) =  Pr (maxYfol <= K | G=g)  
	// = Pr(Yfol1 <= K|G)Pr(Yfol2 <= K|G)....Pr(Yfoln <= K|G)
	denYFol=1;
	for (int ivec = 0; ivec < wIDsizeY; ivec ++){
	  lnpY = fd[get_lnpY(ivec,maxni)];
	  if (lnpY==1){
	    denYFol *= pnbinomYK(K,fd[get_rY(ivec,maxni)],probY,0);
	  }
	}
	break;
      default:
	error("typeSummary must be 1 (sum) or 2 (max)");
    }
  }else{
    // AR(1) model
    int wIDsizeYPre=0,ivec;
    double Yik, Yik_1, sizeYik, sizeYik_1,YFinalPre;
    for (ivec=0;ivec < wIDsizeY;ivec++)
      {
	lnpY = fd[get_lnpY(ivec,maxni)];
	//Rprintf("\n wIDsizeY %d ivec %d lnpY %f",wIDsizeY, ivec, lnpY);

	if (lnpY==0){
	  wIDsizeYPre++;
	  YFinalPre = fd[get_Y(ivec,maxni)];
	}
      }
    const int wIDsizeYFol = wIDsizeY - wIDsizeYPre;
    //Rprintf("\n wIDsizeY %d wIDsizeYfol %d wIDsizeYPre %d",wIDsizeY, wIDsizeYFol,wIDsizeYPre);
    // The combinations differ depending on the summary stat
    const int Ncomb =  get_NComb(wIDsizeYFol, K, typeSummary);
    double YFoli[wIDsizeYFol],denProdUpTo[wIDsizeYFol],denProd,denEach,denProbUpTo1=1,diff;
    int CDFTF, pos = init_getComb(YFoli,wIDsizeYFol, K); 
    // YFoli={0,0,...,0,K} is updated as an initial combination (for both summary type)
    //Rprintf("\n\n within cdfqYfolGivenYPreAndG compute %d combinations of YFol+=%f probY=%f",
    //Ncomb,K,probY);
    denYFol=0;
    for (int icomb = 0 ; icomb < Ncomb;icomb++)
      {
	R_CheckUserInterrupt();
	/* Rprintf("\n          %d st comb. pos=%d YFinalPre=%f  YFoli=",icomb,pos,YFinalPre); */
	/* for (int iii = 0 ; iii < wIDsizeYFol ; iii++) Rprintf(" %f",YFoli[iii]); */
	/* Rprintf("\n            denProdUpTo"); */
	/* for (int ii = 0 ; ii < wIDsizeYFol ; ii++) Rprintf(" %f",denProdUpTo[ii]); */
	//Rprintf("\n icomb %d ",icomb);
	denProd = 1;
	for (int iFol = pos ; iFol < wIDsizeYFol ; iFol ++ )
	  {
	    Yik = YFoli[iFol];
	    Yik_1=ifelse(iFol == 0,YFinalPre,YFoli[iFol-1] );
	    // Rprintf("\n               Yik_1 %f Yik %f",Yik_1,Yik);
	    sizeYik = fd[get_rY(iFol + wIDsizeYPre,maxni)];
	    sizeYik_1 = fd[get_rY(iFol + wIDsizeYPre,maxni)-1];
	    diff = fd[get_diff(iFol + wIDsizeYPre,maxni)];
	    CDFTF = ifelse(iFol==wIDsizeYFol-1,1,0);
	    denEach = PrYijGivenYij_1AndGYAR(Yik, Yik_1, 
					     sizeYik, sizeYik_1,
					     pow(delta,diff),probY, 0,CDFTF);
	    //Rprintf("  denEach %f", denEach);
	    denProbUpTo1 = ifelse(iFol == 0,1, denProdUpTo[iFol-1]);
	    denProdUpTo[iFol] = denProbUpTo1*denEach;
	    //Rprintf(" YFoli[%d] %1.0f",iFol,YFoli[iFol]);
	  }
	//Rprintf("   dens=%1.9f",denProdUpTo[wIDsizeYFol-1]);
	//Rprintf("\n denProdUpTo[%d] %f",wIDsizeYFol-1,denProdUpTo[wIDsizeYFol-1]);
	denYFol += denProdUpTo[wIDsizeYFol-1];
	pos=getComb(YFoli, wIDsizeYFol, K,typeSummary);
      }
  }
  return(denYFol);
}


double update_qYfol(const int typeSummary, const double Yij, double qYfol)
{
  switch(typeSummary){
  case 1: 
    qYfol += Yij; break;
  case 2:
    if (qYfol < Yij)
      qYfol =Yij; 
    break;
  default:
    error("typeSummary must be 1 (=sum), 2(=max)");
  }	
  return(qYfol);
}
// If wIDsize = fd[0] == 1 then Li_integrand becomes the denominator of the CPI
int CPInum_int2(unsigned ndim, const double *x, void *fdata,unsigned fdim, double *fval){
 
  // Pr(sum(YFol) >= sum(yFol),Ypre=ypre| G) 
  // = 1 - Pr( sum(YFol) <sum(yFol),Ypre=ypre| G)
  // = 1 - Pr( sum(YFol) <=sum(yFol)-1|Ypre=ypre,G) Pr(Ypre=ypre| G)
  // = 1 - sum_{all.comb.of.kFol1+...+ kFoln=sum(yFol)-1} Pr(YFol1 = kFol1, YFol2 = kFol2,...,YFoln <= kFoln|Ypre=ypre,G)Pr(Ypre=ypre|G)

  const double *fd = (double*) fdata;
  //const int wIDsizeY = (int)fd[l_wIDsizeY];
  // NXpre = wIDsizeX-1 because the follow-up is only the final BOD
  const int maxni = (int) fd[l_maxni];
  const double alpha = fd[l_alpha],thY = fd[l_thY],delta = fd[l_delta];
  const double gY = intChangeOfVar(x[0]);//1/(1-x[0])-1; Lesion counts
  double rYij,rYij_1,Yij,Yij_1,lnpY,Li_obs,probY = getProb(gY,alpha);
  const int RE = fd[l_distRE];
  Li_obs = distRE(gY,thY,RE);
  int typeSummary = (int) fd[l_typeSummary];
  /* Rprintf("\n\n RE density %f",Li_obs); */
  /* Rprintf("\n \n \n Enter integration function \n \n \n"); */
  /* for (int iv = 0 ; iv < get_fdsize(maxni); iv++)  Rprintf("\n fdata[%d] %f",iv,fd[iv]); */
  /* Rprintf("\n \n \n"); */
  /* Rprintf("\n gY %f gX %f",gY,gX); */

  R_CheckUserInterrupt(); 

  int ivec;
  double rYFolSum=0, qYfol=0,diff;
  for (ivec = 0 ; ivec < maxni; ivec++)
    {
      lnpY = fd[ get_lnpY(ivec,maxni) ];
      rYij = fd[ get_rY(ivec,maxni)   ];
      Yij = fd[  get_Y(ivec,maxni)    ];
      diff = fd[get_diff(ivec,maxni)];
      rYij_1 = ifelse(ivec==0,R_NegInf,fd[get_rY(ivec,maxni) - 1]);
      Yij_1 = ifelse(ivec==0,R_NegInf,fd[get_Y(ivec,maxni) - 1]);

      if (lnpY==0){ 
	Li_obs *= densYijGivenYij_1AndGY(Yij, Yij_1, 
					 rYij, rYij_1, 
					 pow(delta,diff), 
					 probY, 0);

      }else if (lnpY == 1){
	qYfol = update_qYfol(typeSummary, Yij, qYfol);
	rYFolSum += rYij; //sum
      }
    }

  Li_obs *= cdfqYfolGivenYPreAndG(fd,qYfol-1,rYFolSum,probY); 
  //Rprintf("\n\n K %d Li_obs_sum=%f ",K,Li_obs_sum);
    
  // change of variables
  fval[0] = intJacob1(x[0])*Li_obs;
  return 0;
}








double CPISinglePat(double fd[],
		    const double  wIDsizesY_ipat, 
		    const double Y[],
		    const double diff[],
		    const int lnpY[],
		    SEXP ZY,
		    int wIDsY_ipat[],
		    const double betaY[],
		    const int NtotY,
		    const int pY,
		    const int printing,
		    const double absTol)
{
  double anyScanf = 0;
  R_CheckUserInterrupt(); 
  const int maxni = (int) fd[l_maxni];
  int typeSummary = (int)fd[l_typeSummary];
  fd[l_wIDsizeY] = 0; // The number of pre-scans in CEL counts
  double qYfol = 0, val,err,temp, den;
  for (int ivec = 0 ; ivec < maxni; ivec++)
    {
      if (ivec < wIDsizesY_ipat){
	fd[get_rY(ivec,maxni)] = exp(getrij(ivec,ZY, wIDsY_ipat,
					    betaY, NtotY,pY))/fd[l_alpha];
	fd[get_Y(ivec,maxni)] = Y[wIDsY_ipat[ivec]];
	fd[get_diff(ivec,maxni)] = diff[wIDsY_ipat[ivec]];
	fd[get_lnpY(ivec,maxni)] = lnpY[wIDsY_ipat[ivec]];
	if (lnpY[wIDsY_ipat[ivec]]==0) fd[l_wIDsizeY]++;
	else if (lnpY[wIDsY_ipat[ivec]]==1){
	  anyScanf++;
	  qYfol=update_qYfol(typeSummary, Y[wIDsY_ipat[ivec]], qYfol);

	}
      }else{
	fd[get_rY(ivec,maxni)] = -1000;
	fd[get_Y(ivec,maxni)] = -1000;
	fd[get_diff(ivec,maxni)] = -1000;
	fd[get_lnpY(ivec,maxni)] = -1000;
      }
    }
  
  if (printing) Rprintf(" q(Yfol)=%1.0f",qYfol);
  
  if (qYfol==0){
    temp = 1;
  }else{
    //Rprintf("\n\n== Before computation after inputting using loop ===\n");
    //for (ivec = 0 ; ivec < get_fdsize(maxni) ; ivec++) 
    //Rprintf("\n fd[%d] %f",ivec,fd[ivec]);
    //========= denominator if there is at least one observation in the condition ========
    if (fd[l_wIDsizeY] > 0){
      temp = hcubature(1, intLik, fd, 1, xmin, xmax,
		       0, absTol, 0,  // maxEval, reqAbsError, reqRelError
		       ERROR_INDIVIDUAL, &val, &err);
    }else val = 1;
    den = val;    
    //Rprintf(" den= %f",den);
    //Rprintf("\n\n ======ipat=%d======\nComputed Denominator",ipat);
    //for (ivec = 0 ; ivec < get_fdsize(maxni) ; ivec++) Rprintf("\n fd[%d] %f",ivec,fd[ivec]);
    //Rprintf("\n Pr(Xpre,Ypre)=%f",den);
    //==================== numerator ======================
    fd[l_wIDsizeY] = wIDsizesY_ipat;
   
    temp = hcubature(1, CPInum_int2, fd, 1, xmin, xmax,
		     0, absTol, 0,  // maxEval, reqAbsError, reqRelError
		     ERROR_INDIVIDUAL, &val, &err);
    temp = 1.0 - val/den;

  }
  if (printing) Rprintf(" hat.p=%1.3f",temp );
  return(temp);
}


// Compute the conditional probability index for single patient
SEXP CPI_each(SEXP Y_,  // double 
	       SEXP ZY,  // double 
	       
	       SEXP alpha_,  // double
	       SEXP thY_, // double
	       SEXP delta_, // double
	       
	       SEXP betaY_, // double 
	       SEXP distRE_, // double 
	       SEXP lnpY_, // int
	       SEXP diff_, // double
	      SEXP printing_,
	      SEXP typeSummary,
	      SEXP AbsTol)
{
  const int maxni = length(Y_);
  const double *diff= REAL(diff_);
  const double *Y=REAL(Y_);
  const int *lnpY= INTEGER(lnpY_);
  int wIDsY_ipat[maxni];
  for (int ivec = 0; ivec < maxni;ivec++) wIDsY_ipat[ivec] = ivec;
  const double *betaY = REAL(betaY_);
  const int pY = length(betaY_);
  const int printing = INTEGER(printing_)[0];
  const double absTol = REAL(AbsTol)[0];
  double fd[get_fdsize(maxni)];
  fd[l_alpha] = REAL(alpha_)[0];
  fd[l_thY] = REAL(thY_)[0];
  fd[l_delta] = REAL(delta_)[0];
  fd[l_maxni] = (double) maxni;
  fd[l_distRE] = REAL(distRE_)[0];
  fd[l_wIDsizeY] = maxni;
  fd[l_typeSummary] = REAL(typeSummary)[0];
  //If maxni = 8 then: 
  
  SEXP out;
  PROTECT(out = allocVector(REALSXP, 1));
  REAL(out)[0] = CPISinglePat(fd,maxni,Y,diff,
			      lnpY,ZY,wIDsY_ipat,betaY,maxni,
			      pY,printing,absTol);
  UNPROTECT(1);
  return(out);
}


SEXP CPI_ALL(SEXP Y_,  // double 
	     SEXP ZY,  // double
	     
	     SEXP IDY, // integer
	     
	     SEXP alpha_,  // double
	     SEXP thY_, // double
	     SEXP delta_, // double
	     
	     SEXP betaY_, // double 
	     SEXP distRE_, // double 
	     SEXP maxni_,  // int
	     SEXP Npat_, // int 
	     SEXP lnpY_, // int
	     SEXP diff_, // double
	     SEXP printing_,
	     SEXP typeSummary, // double 1 = sum, 2 = max
	     SEXP AbsTol
	     )
{
  const double *Y = REAL(Y_);
  const int maxni = INTEGER(maxni_)[0];
  const int *lnpY = INTEGER(lnpY_);
  const double *diff = REAL(diff_);
  const int printing = INTEGER(printing_)[0];
  // The number of parameters that are passed to the integrant function to compute the likelihood 
  double fd[get_fdsize(maxni)];
  fd[l_maxni]= (double) maxni;
  fd[l_alpha]= REAL(alpha_)[0];
  fd[l_thY]= REAL(thY_)[0];
  fd[l_delta]= REAL(delta_)[0];
  fd[l_distRE]= REAL(distRE_)[0];
  fd[l_typeSummary] = REAL(typeSummary)[0];
  //Rprintf("\n\n=== The very first fd ===\n");
  //for (int ii = 0 ; ii < get_fdsize(maxni) ; ii++) Rprintf("\n fd[%d] %f",ii,fd[ii]);
  const double *betaY = REAL(betaY_);
  const int NpatTot = INTEGER(Npat_)[0];
  const int NtotY = length(IDY); // 
  const int pY = length(ZY)/NtotY;
  const double absTol = REAL(AbsTol)[0];
  GetRNGstate();
  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  
  SEXP res = PROTECT(allocVector(VECSXP, 1)); // result is stored in res
  SEXP CPI = allocVector(REALSXP, NpatTot); 
  SET_VECTOR_ELT(res, 0, CPI); 
  
  int ipat,ivec;
 
  int wIDY[maxni],wIDsizeY = 0, wIDsY[NpatTot][maxni], wIDsizesY[NpatTot];

  for (ipat = 0 ; ipat < NpatTot ; ipat++)
    {
      getwID( &wIDsizeY,wIDY, IDY, ipat);
      wIDsizesY[ipat] = wIDsizeY;
      for (ivec = 0 ; ivec < wIDsizeY ; ivec++) wIDsY[ipat][ivec]=wIDY[ivec];
      for (ivec = wIDsizeY ; ivec < maxni; ivec++) wIDsY[ipat][ivec]=-1000;
    }
  for (ipat = 0; ipat < NpatTot ; ipat++)
    {
      if (printing) Rprintf("\n patient %d",ipat+1);

      REAL(CPI)[ipat] = CPISinglePat(fd, wIDsizesY[ipat],
				     Y,diff,lnpY,ZY,wIDsY[ipat],
				     betaY,NtotY,pY,printing,absTol);
    }
  // check the number of combinations
 
  PutRNGstate();
  UNPROTECT(1);

  return res;
}
















