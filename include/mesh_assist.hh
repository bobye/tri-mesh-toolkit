#ifndef _MESH_ASSIST_H_
#define _MESH_ASSIST_H_

#include "mesh_topo.h"

extern void localchart(Vector &, Vector &, Vector);
extern void localcoord(Vector, Vector, Vector, double*);
extern double prin_curv(double, double, double, double&, double&);


#endif /* _MESH_ASSIST_H_ */
