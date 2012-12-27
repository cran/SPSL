#include <R.h>
#include <Rinternals.h>
// ====================
// License: GPL-3
// Package: Site Percolation on Square Lattice (SPSL)
// Author: Pavel V. Moskalev <moskalefff@gmail.com>
// ================================================
SEXP ssTNd(SEXP pA, SEXP acA, SEXP bA, SEXP eA, SEXP clA) {
// ========================================================
// Function:
// ssTNd() function performs labeling sites on the n-dimensional
// square lattice with anisotropic Moore (1,d)-neighborhood.
// These sites form the cluster connected with a starting subset.
// ==============================================================
// Arguments:
//   bA - number of initial sites;
//   clA - linear indexes of cluster sites;
//   acA - labeled accessibility matrix;
//   eA, pA - linear indexes and relative fractions of 
//            accessibile sites from Moore (1,d)-neighborhood;
// ===========================================================
// Variables:
//   cls, acc, e, b, p - pointers to clA, acA, eA, bA, pA
//   a, *b - indexes of the current cluster perimeter
//           in cls[] vector: from, to;
//   c - index of the current site in cls[] vector;
//   h - index of the current site neighbor in e[] vector;
//   n - number of sites in Moore (1,d)-neighborhood.
// ==================================================
  acA = coerceVector(acA, REALSXP);
  clA = coerceVector(clA, INTSXP);
  pA = coerceVector(pA, REALSXP);
  bA = coerceVector(bA, INTSXP);
  eA = coerceVector(eA, INTSXP);
  double *p, *acc;
  int *b, *e, *cls,
      n=length(eA), a=0, db, c, h, ch;
  cls = INTEGER(clA); e = INTEGER(eA); b = INTEGER(bA); 
  acc = REAL(acA); p = REAL(pA); db = *b;
  while (db>0) {              // While the perimeter is not empty:
    db = 0;                   // The new perimeter is empty.
    for (c=a; c<*b; c++) {    // For all sites in the perimeter:
      for (h=0; h<n; h++) {   // For all sites in the neighborhood:
        ch = cls[c] + e[h];
        if (acc[ch] < p[h]) { // If the site is accessible
          acc[ch] = 2; db++;  // then label this site
          cls[*b+db-1] = ch;  // and save its index.
        }
      }
    }
    a = *b; *b += db;         // Set indexes to the perimeter.
  }
  return(R_NilValue);
}
