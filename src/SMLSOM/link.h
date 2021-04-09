#ifndef INCLUDE_link_h_
#define INCLUDE_link_h_

//list structure for undirect graph

#include "base.h"

/*
typedef struct list {
  int node;
  struct list *next;
} list;
typedef list* link;
*/

link *adjsNodes; // ajdsNodes[k] is a list of adjacents of node k
link availNodes; // nodes which is not dead yet

//int nnode;

extern void pushTo(link que, int k);
extern int getTopFrom(link que);
extern boolean getkFrom(link que, int k);
extern int getDim(link que);
extern boolean isIn(link que, int k);
extern void free_que(link que);
extern link make_dummy();
extern void init_list(int nnode);
extern void free_list(int nnode);
extern void createEdge(int k, int node);
extern boolean deleteEdge(int k, int node);
extern boolean deleteNode(int k);
extern void print_list(int nnode);
extern void print_que(link que);

#endif //INCLUDE_link_h_
