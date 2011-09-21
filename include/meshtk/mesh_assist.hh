#ifndef _MESH_ASSIST_HH_
#define _MESH_ASSIST_HH_

#include "mesh_topo.hh"

extern void localchart(Vector &, Vector &, Vector);
extern void localcoord(Vector, Vector, Vector, double*);
extern double prin_curv(double, double, double, double&, double&);


template <class T> 
struct index_cmp {
  index_cmp(const T arr): arr(arr) {}
  bool operator() (const size_t a, const size_t b) const {
    return arr[a] < arr[b];
  }
  const T arr;
};


#endif /* _MESH_ASSIST_HH_ */
