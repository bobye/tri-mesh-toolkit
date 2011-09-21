#ifndef _MESH_ASSIST_HH_
#define _MESH_ASSIST_HH_

#include "mesh_topo.hh"

extern void localchart(Vector &, Vector &, Vector);
extern void localcoord(Vector, Vector, Vector, double*);
extern double prin_curv(double, double, double, double&, double&);


#endif /* _MESH_ASSIST_HH_ */
