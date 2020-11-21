#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

struct def_grad
{
	// F - deformation gradient
	Matrix2d F;
	// dF - derivative of the deformation gradient with respect to nodal displacements displacements
    // dF = [dF/du_1; dF/dv_1; dF/du_2; dF/dv_2; dF/du_3; dF/dv_3;]
    // dF/du_i and dF/dv_i are tensors written as columns:
    // dF/du_i=[dF/du_i(1,1) dF/du_i(2,1) dF/du_i(1,2) dF/du_i(2,2)];
	Matrix<double, 4, 6> dF;
	// A - areas
  double A;
  // F_r=F*m - reduced matrix after the Lagrange reduction
  Matrix2d F_r;
  // m - corresponding transform
  Matrix2d m;
  // C - Cauchy-Green tensor
  Matrix2d C;
  // corresponding reduced metric
  Matrix2d C_r;
};

struct en_density
{
	// W - value of the energy density
    double W;
    // dW - derivative with respect to the metric tensor C
    Matrix2d dW;
    // W_T - derivative with respect to the displacements of nodes of the element u1 v1 u2 v2 ...
    VectorXd dW_T;
};

struct en_parameters
{
	// type of the energy density = 'kirchhoff', 'poly'
    std::string en_type;
    // symmetry type 'sq' or 'hex' for polynomial energy (square and hexagonal symmetry, respectively)
    std::string sym;
    // matrix C (for 'kirchhoff')
    Matrix3d C;
    // bulk modulus
    double K;
};

struct stress
{
	// stresses in elements and total values

  // Cauchy stresses in each element
  Matrix2d sigma;
  // 1st Piola-Kirchhoff stresses in each element
  Matrix2d piola_1;
  // 2nd Piola-Kirchhoff stresses in each element
  Matrix2d piola_2;
};

struct gl_energy
{
	// global energy of a material

  // its value
  double E;
  // derivatives with respect to displacements U
  VectorXd dE;
};

struct displ_bc
{
	displ_bc(): top(0),bottom(0), left(0), right(0) { }
	// displacement on the TOP boundary
	VectorXd top;
	// displacement on the BOTTOM boundary
  VectorXd bottom;
  // displacement on the LEFT boundary
  VectorXd left;
  // displacement on the RIGHT boundary
  VectorXd right;
};

struct BCs
{
	// nodes (indices) of the boundaries
	// TOP boundary
	VectorXi top;
	// BOTTOM boundary
  VectorXi bottom;
  // LEFT boundary
  VectorXi left;
  // RIGHT boundary
  VectorXi right;

  // free indices of displacements of the  boundary (even - horizontal displacements, odd - vertical displacements)
  VectorXi ind_free;
	// loading parameter
	double s;
  // HORIZONTAL displacements of the boundaries
  displ_bc u0;
  // VERTICAL displacements of the boundaries
  displ_bc v0;

};

struct BCs_prop
{
	// type of boundary conditons
	// bc_type is a string: 'F' - prescribed deformation gradient, 'hard' - prescribed displacements
	string bc_type;
	// periodic LEFT and RIGHT
  int u_hor_per;
  // periodic TOP and BOTTOM
  int u_ver_per;
  // periodic LEFT and RIGHT
  int v_hor_per;
  // periodic TOP and BOTTOM
  int v_ver_per;

};

struct geom
{
	// generate the geometry from points p and connectivity list T

	// nodes locations (x,y)
	MatrixXd p;
	// connectivity matrix from delaunay triangulation
	MatrixXi T;

	// boundary nodes
	// TOP boundary
	VectorXi top;
	// BOTTOM boundary
  VectorXi bottom;
  // LEFT boundary
  VectorXi left;
  // RIGHT boundary
  VectorXi right;

  // ALL boundary nodes
  VectorXi bound;

  // bulk nodes
  VectorXi bulk;
};
struct model_param
{
  // parameters of the model in one structure

  // parameters of the energy
  en_parameters en_par;
  // points
  MatrixXd points;
  // connectivity matrix
  MatrixXi T;
  // boundary conditions
  BCs BCs_val;
  // properties of the boundary conditions
  BCs_prop BCs_p;

};
#endif
