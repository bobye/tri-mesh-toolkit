# tri-mesh-toolkit
A toolkit for mesh (sequence) (pre-)processing, shape analysis and visualization

## Build Dependency
- GLUT
- OpenGL
- CGAL (version 4.3)
- lapack
- PETSc (version 3.4.3)
- SLEPc (version 3.4.3)
- ATLAS

## Limitations
tri-mesh-toolkit only supports manifold meshes. 

## Publication

J1: Jianbo Ye and Yizhou Yu, *A fast modal space transform for robust nonrigid shape retrieval*. TVC 2015

C1: Jianbo Ye, Zhicheng Yan, and Yizhou Yu, *Fast nonrigid 3D retrieval using modal space transform*. ICMR 2013

## Dataset
Additional custom built [dataset](https://drive.google.com/file/d/0BzfQ7xg9jFT3b0NaWXpsNU04ckE/view?usp=sharing) used:
standard approach to test a retrieval method is described in the aforementioned paper. It uses 200 models as query 
set, and the rest 40 models plus 600 [SHREC'11 models](http://www.itl.nist.gov/iad/vug/sharp/contest/2011/NonRigid/) as background models. 

_If you use above dataset for research purposes, please cite (C1) or (J1)_

 
