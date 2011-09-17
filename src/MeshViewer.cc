#include "MeshViewer.h"

MeshPainter::MeshPainter(TriMesh *pmesh){

}

MeshPainter::MeshPainter(TriMesh *pmesh, Scalar_Fun *psfun){
}

MeshPainter::draw(){
}

MeshMarker::MeshMarker(TriMesh *pmesh, Bool_Fun* pbfun)
  :MeshPainter(pmesh) {
  
}

MeshMarker::MeshMarker(TriMesh *pmesh, Scalar_Fun *psfun, Bool_Fun *pbfun)
  :MeshPainter(pmesh, psfun) {
}


MeshMarker::draw(){
}





MeshViewer::MeshViewer(int width, int height) {
}


void MeshViewer::add_painter(MeshPainter *painter){
  Painters.push_back(painter);
}

void MeshViewer::view(){
}


