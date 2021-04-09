#ifndef INCLUDE_searchGraph_h_
#define INCLUDE_searchGraph_h_

#include "base.h"

extern void bfs(int node, int depth, link neighs, int nnode);
// breadth-first search for neighbours of a given node which has the link distance less than or equal to `depth`.
// the neighbours include the center are stored in `neighs`
extern int bfs_dist(int node, int nnode);
extern void graph2adjmat(int *adjmatrix, int M);

#endif
