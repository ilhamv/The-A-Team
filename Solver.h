#ifndef _SOLVER_HEADER_
#define _SOLVER_HEADER_

#include <vector>

#include "Const.h" // MAX


// Quad solver for sphere surface
// return smallest positive real root if it exists; if it does not, return very big number
// *Seems redundantly hard-coded, yet it saves up the multiplication operation with a
double sphere_quad( const double b, const double c );

// Quad solver for infinite cylinder and cone surface
// return smallest positive real root if it exists; if it does not, return very big number
double solve_quad( const double a, const double b, const double c );

// Binary search a double location in bin grids
// Note:
// 	value < lowest  grid --> -1
// 	value > highest grid --> vector.size - 1 (or number of bins)
// 	value = grid points  --> location of bin whose upper bound is the value
// 	                         (-1 if value = lowest grid)
int Binary_Search( const double x, const std::vector<double>& vec );

#endif
