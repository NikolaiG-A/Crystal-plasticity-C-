#ifndef FUNCTIONS_INIT_H
#define FUNCTIONS_INIT_H

#include <vector>
#include <set>
#include <algorithm>
#include <omp.h>

#include "structures.h"
#include "Eigen/Dense"

using namespace Eigen;

// Basis reduction algorithm
Matrix2d basis_reduction(Matrix2d E);
// initialization of the deformation gradients
def_grad def_grad_init(MatrixXd p, Vector3i T, VectorXd U);
// generation of the geometry
geom geom_gen(MatrixXd p, MatrixXi T);
// initialization of the boundary conditions
BCs bound_cond(geom G, BCs_prop B_pro, double s=0., displ_bc u0_val=displ_bc(), displ_bc v0_val=displ_bc());
// initialization of the energy density
en_density en_density_init(en_parameters en_par,Matrix2d C);
// initialization of the energy density and its derivative for an element
en_density en_density_element(en_parameters en_par,MatrixXd p, Vector3i T, VectorXd U);
// initialization of the total energy
gl_energy global_energy_init(en_parameters en_par,MatrixXd p, MatrixXi T, VectorXd U);
// initialization of the total energy with constrained displacements
gl_energy global_energy_constr_init(en_parameters en_par,MatrixXd p, MatrixXi T, BCs BCs_val, BCs_prop B_pro, VectorXd W);
// initialization the displacements
VectorXd total_displ(MatrixXd p,MatrixXi T, BCs BCs_val, BCs_prop B_pro, VectorXd W);

#endif