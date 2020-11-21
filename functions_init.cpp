#include "functions_init.h"
#include "calculation.h"
// // Basis reduction algorithm
// Matrix2d basis_reduction(Matrix2d E)
// {

//     if (E.col(0).norm()>E.col(1).norm())
// 	   E=E(all,seq(1,0,-1));
// 	if (E.col(0).dot(E.col(1))<0)
// 	{
// 		E(0,1)=-E(0,1);
// 	    E(1,1)=-E(1,1);
// 	}
// 	if ((E.col(1)-E.col(0)).norm()<E.col(1).norm())
// 	{
// 		E(0,1)=E(0,1)-E(0,0);
// 		E(1,1)=E(1,1)-E(1,0);
// 		E = basis_reduction(E);
// 	}
// 	return(E);
// }

// Basis reduction algorithm
void basis_reduction(std::vector<double> &e1,std::vector<double> &e2)
{

    if (e1[0]*e1[0]+e1[1]*e1[1]>e2[0]*e2[0]+e2[1]*e2[1])
       e1.swap(e2);
    if (e1[0]*e2[0]+e1[1]*e2[1]<0)
    {
        e2[0]=-e2[0];
        e2[1]=-e2[1];
    }
    if ((e2[0]-e1[0])*(e2[0]-e1[0])+(e2[1]-e1[1])*(e2[1]-e1[1])<e2[0]*e2[0]+e2[1]*e2[1])
    {
        e2[0]=e2[0]-e1[0];
        e2[1]=e2[1]-e1[1];
        basis_reduction(e1,e2);
    }
}











// initialization of the deformation gradients
def_grad def_grad_init(MatrixXd p, Vector3i T, VectorXd U)
{
	////////////// input //////////////
    // p - points of the mesh
    // U - displacements
    // T - connectivity list (of one element)

    ////////////// output //////////////
    // F - structure of the type 'deform_grad' (see structures.h)


	/// calculate element area
	Matrix3d Pe; // construct a 3x3 matrix
    Pe<< 1., 1., 1., // first row
         p(T(0),0), p(T(1),0), p(T(2),0), // second row
         p(T(0),1), p(T(1),1), p(T(2),1); // thirtd row
    // element area
	double A;
	A = 0.5*Pe.determinant();
	/// inverse transform of a current element to a reference one
	Matrix2d invJ;
    invJ<< p(T(2),1)-p(T(0),1),-(p(T(2),0)-p(T(0),0)),
          -(p(T(1),1)-p(T(0),1)),p(T(1),0)-p(T(0),0);
	invJ*=0.5/A;
	/// matrix of dispacement gradients
    Matrix2d def_u;
    def_u<<U(2*T(1))-U(2*T(0)),U(2*T(2))-U(2*T(0)),
           U(2*T(1)+1)-U(2*T(0)+1),U(2*T(2)+1)-U(2*T(0)+1);
    /// calculate the deformation gradient
    Matrix2d F_temp;
    F_temp=Matrix2d::Identity()+def_u*invJ;
    /// the matrix of derivatives comes from differentiation with respect
    Matrix<double,6,4> dF;
    dF(0,0)=-1.0*(invJ(0,0)+invJ(1,0)); dF(0,1)=0; dF(0,2)=-1.0*(invJ(0,1)+invJ(1,1)); dF(0,3)=0;
    dF(1,0)=0; dF(1,1)=-1.0*(invJ(0,0)+invJ(1,0)); dF(1,2)=0; dF(1,3)=-1.0*(invJ(0,1)+invJ(1,1));
    dF(2,0)=invJ(0,0); dF(2,1)=0; dF(2,2)=invJ(0,1); dF(2,3)=0;
    dF(3,0)=0; dF(3,1)=invJ(0,0); dF(3,2)=0; dF(3,3)=invJ(0,1);
    dF(4,0)=invJ(1,0); dF(4,1)=0; dF(4,2)=invJ(1,1); dF(4,3)=0;
    dF(5,0)=0; dF(5,1)=invJ(1,0); dF(5,2)=0; dF(5,3)=invJ(1,1);
    // perform lagrange reduction
    std::vector<double> e1={F_temp(0,0),F_temp(1,0)},e2={F_temp(0,1),F_temp(1,1)};
    basis_reduction(e1,e2);
    Matrix2d F_r;
    //F_r=basis_reduction(F_temp);
    F_r<<e1[0],e2[0],
         e1[1],e2[1];
    // assign matrix m
    Matrix2d m_temp;
    m_temp=F_temp.inverse()*F_r;
    // Cauchy-Green measures
    Matrix2d C;
    C=F_temp.transpose()*F_temp;

    Matrix2d C_r;
    C_r=F_r.transpose()*F_r;
    // define the deformation gradient structure
    def_grad F;
    F.F=F_temp;
    F.dF=dF.transpose();
    F.A=A;
    F.m=m_temp;
    F.F_r=F_r;
    F.C=C;
    F.C_r=C_r;
    return(F);
}





















// generation of the geometry
geom geom_gen(MatrixXd p, MatrixXi T)
{
    ////////////// input //////////////
    // U - displacements
    // T - connectivity list

    ////////////// output //////////////
    // G - structure of the type 'geom' (see structures.h)

    std::vector<int> ind_bound,ind_top,ind_bottom,ind_left,ind_right,ind_bulk;

    // identify the boundaries of the square grid
    for(int i=0;i<p.rows();i++)
    {
        if(abs(p(i,1)-p(all,1).maxCoeff())<0.01)
            ind_top.push_back(i);
        if(abs(p(i,1)-p(all,1).minCoeff())<0.01)
            ind_bottom.push_back(i);
        if(abs(p(i,0)-p(all,0).maxCoeff())<0.01 && (p(i,1)>p(all,1).minCoeff()) && (p(i,1)<p(all,1).maxCoeff()))
            ind_right.push_back(i);
        if(abs(p(i,0)-p(all,0).minCoeff())<0.01 && (p(i,1)>p(all,1).minCoeff()) && (p(i,1)<p(all,1).maxCoeff()))
            ind_left.push_back(i);
        if((p(i,0)>p(all,0).minCoeff()) && (p(i,0)<p(all,0).maxCoeff()) && (p(i,1)>p(all,1).minCoeff()) && (p(i,1)<p(all,1).maxCoeff()))
            ind_bulk.push_back(i);
    }
    // convert these quanntites to the type of Eigen:VectorXi
    Map<VectorXi> top(&ind_top[0], ind_top.size());
    Map<VectorXi> bottom(&ind_bottom[0], ind_bottom.size());
    Map<VectorXi> right(&ind_right[0], ind_right.size());
    Map<VectorXi> left(&ind_left[0], ind_left.size());
    Map<VectorXi> bulk(&ind_bulk[0], ind_bulk.size());

    // sort these vectors
    std::sort(top.data(),top.data()+top.size());
    std::sort(bottom.data(),bottom.data()+bottom.size());
    std::sort(right.data(),right.data()+right.size());
    std::sort(left.data(),left.data()+left.size());
    std::sort(bulk.data(),bulk.data()+bulk.size());

    //collect all the indices into one vector
    std::vector<int> bound_vec=ind_bottom;
    bound_vec.insert(bound_vec.end(),ind_left.begin(),ind_left.end());
    bound_vec.insert(bound_vec.end(),ind_right.begin(),ind_right.end());
    bound_vec.insert(bound_vec.end(),ind_top.begin(),ind_top.end());


    // extract  the unique indices
    std::set<int> bound_set(bound_vec.begin(),bound_vec.end());
    //convert it to the type of Eigen::VectorXi
    std::vector<int> bound_unique(bound_set.begin(),bound_set.end());
    Map<VectorXi> bound(&bound_unique[0], bound_unique.size());
    // make the structure to return the results
    geom G;
    G.p=p;
    G.T=T;
    G.top=top;
    G.bottom=bottom;
    G.right=right;
    G.left=left;
    G.bound=bound;
    G.bulk=bulk;

    return(G);

}









// initialization of the boundary conditions
//BCs bound_cond(geom G, BCs_prop B_pro, double s=0., displ_bc u0_val=displ_bc(), displ_bc v0_val=displ_bc())
BCs bound_cond(geom G, BCs_prop B_pro, double s, displ_bc u0_val, displ_bc v0_val)

{
    ////////////// input //////////////
    // G - structure with geometry properties
    // B_pro - properties of boundary conditions
    // s - loading parameters


    ////////////// output //////////////
    // BCs_val - structure of the type 'BCs' (see structures.h)

    // define the boundary conditions
    BCs BCs_val;
    if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
    {
      BCs_val.top=G.top(seq(1,last-1));
      BCs_val.bottom=G.bottom(seq(1,last-1));
      BCs_val.left=G.left;
      BCs_val.right=G.right;
    }
    else
    {
      BCs_val.top=G.top;
      BCs_val.bottom=G.bottom;
      BCs_val.left=G.left;
      BCs_val.right=G.right;
    }


    displ_bc u0=u0_val;
    displ_bc v0=v0_val;

    if (B_pro.bc_type == "F")
    {
        // deformation gradient to load
        Matrix2d F_load;
        F_load<<1.,s,
                0.,1.;
        // new points
        MatrixXd p_new;
        p_new = G.p*F_load.transpose();

        // horizontal displacements of the bundaries
        u0.top=p_new(BCs_val.top,0)-G.p(BCs_val.top,0);
        u0.bottom=p_new(BCs_val.bottom,0)-G.p(BCs_val.bottom,0);
        u0.left=p_new(BCs_val.left,0)-G.p(BCs_val.left,0);
        u0.right=p_new(BCs_val.right,0)-G.p(BCs_val.right,0);

        // vertical displacements of the bundaries
        v0.top=p_new(BCs_val.top,1)-G.p(BCs_val.top,1);
        v0.bottom=p_new(BCs_val.bottom,1)-G.p(BCs_val.bottom,1);
        v0.left=p_new(BCs_val.left,1)-G.p(BCs_val.left,1);
        v0.right=p_new(BCs_val.right,1)-G.p(BCs_val.right,1);
    }

    // all indices of degrees of freedom for boundary nodes
    std::vector<int> ind_bound;
    VectorXi ind_bound_full(2*G.bound.size());
    ind_bound_full<< 2*G.bound.array(),2*G.bound.array()+1;



    // if HORIZONTAL displacements on the LEFT and RIGHT boundaries are periodic
    if (B_pro.u_hor_per==1)
    {
        for (int i=0;i<ind_bound_full.size();i++)
            for (int j=0;j<BCs_val.left.size();j++)
                if (2*BCs_val.left(j)==ind_bound_full(i))
                    ind_bound.push_back(i);
    }
    // if VERTICAL displacements on the LEFT and RIGHT boundaries are periodic
    if (B_pro.v_hor_per==1)
    {
        for (int i=0;i<ind_bound_full.size();i++)
            for (int j=0;j<BCs_val.left.size();j++)
                if (2*BCs_val.left(j)+1==ind_bound_full(i))
                    ind_bound.push_back(i);
    }

    // if HORIZONTAL displacements on the TOP and BOTTOM boundaries are periodic
    if (B_pro.u_ver_per==1)
    {
        for (int i=0;i<ind_bound_full.size();i++)
            for (int j=0;j<BCs_val.bottom.size();j++)
                if (2*BCs_val.bottom(j)==ind_bound_full(i))
                    ind_bound.push_back(i);
    }
    // if VERTICAL displacements on the TOP and BOTTOM boundaries are periodic
    if (B_pro.v_ver_per==1)
    {
        for (int i=0;i<ind_bound_full.size();i++)
            for (int j=0;j<BCs_val.bottom.size();j++)
                if (2*BCs_val.bottom(j)+1==ind_bound_full(i))
                    ind_bound.push_back(i);
    }

    if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
    {
      // extract only free nodes of the boundaries
      VectorXi ind_free(2*G.bulk.size()+ind_bound.size()+2);;
      //VectorXi ind_free_bound=ind_bound_full(ind_bound);
      ind_free<< 0,1,2*G.bulk.array(),2*G.bulk.array()+1,ind_bound_full(ind_bound).array();
      std::sort(ind_free.data(),ind_free.data()+ind_free.size());
      BCs_val.ind_free=ind_free;
    }
    else
    {
      // extract only free nodes of the boundaries
      VectorXi ind_free(2*G.bulk.size()+ind_bound.size());;
      //VectorXi ind_free_bound=ind_bound_full(ind_bound);
      ind_free<< 2*G.bulk.array(),2*G.bulk.array()+1,ind_bound_full(ind_bound).array();
      std::sort(ind_free.data(),ind_free.data()+ind_free.size());
      BCs_val.ind_free=ind_free;
    }

    BCs_val.u0=u0;
    BCs_val.v0=v0;

    BCs_val.s=s;
    return(BCs_val);
}









// initialization of the energy density
en_density en_density_init(en_parameters en_par,Matrix2d C)
{
    ////////////// input //////////////
    // en_par - structure with energy density properties
    // C - Cauchy - Green strain


    ////////////// output //////////////
    // en_den - structure of the type 'en_density' (see structures.h)

    // constants for the energy densities
    double W_temp=0.,W_v_temp=0.;
    // constants for their derivatives
    double dW_dC11=0.,dW_dC22=0.,dW_dC12=0.,dW_vol_dC11=0.,dW_vol_dC22=0.,dW_vol_dC12=0.;
    // define the components of C
    double C11=C(0,0),C22=C(1,1),C12=C(0,1);
    // bulk modulus
    double K=en_par.K;

    if(en_par.en_type=="kirchhoff")
    {

        // components of the Hooke's stiffness tensor
        double C_h11,C_h22,C_h33,C_h12,C_h13,C_h23;
        C_h11=en_par.C(0,0);
        C_h22=en_par.C(1,1);
        C_h33=en_par.C(2,2);
        C_h12=en_par.C(0,1);
        C_h13=en_par.C(0,2);
        C_h23=en_par.C(1,2);

        // energy density
        W_temp=C12*1.0/sqrt(-C12*C12+C11*C22)*((C_h13*(C11/2.0-1.0/2.0))/((sqrt(-C12*C12+C11*C22))*2.0)+(C_h23*(C22/2.0-1.0/2.0))/((sqrt(-C12*C12+C11*C22))*2.0)+(C12*C_h33)/((sqrt(-C12*C12+C11*C22))*2.0))+(C11/2.0-1.0/2.0)*1.0/sqrt(-C12*C12+C11*C22)*((C_h11*(C11/2.0-1.0/2.0))/((sqrt(-C12*C12+C11*C22))*2.0)+(C_h12*(C22/2.0-1.0/2.0))/((sqrt(-C12*C12+C11*C22))*2.0)+(C12*C_h13)/((sqrt(-C12*C12+C11*C22))*2.0))+(C22/2.0-1.0/2.0)*1.0/sqrt(-C12*C12+C11*C22)*((C_h12*(C11/2.0-1.0/2.0))/((sqrt(-C12*C12+C11*C22))*2.0)+(C_h22*(C22/2.0-1.0/2.0))/((sqrt(-C12*C12+C11*C22))*2.0)+(C12*C_h23)/((sqrt(-C12*C12+C11*C22))*2.0));
        // volumetric part W_v=K*(detC-log(detC)), K - bulk modulus
        W_v_temp=-K*(log(-C12*C12+C11*C22)+C12*C12-C11*C22);

        // derivatives (from symbolic differentiation)
        dW_dC11=1.0/pow((sqrt(-C12*C12+C11*C22)),3.0)*1.0/pow(-C12*C12+C11*C22,3.0/2.0)*((C22*C22)*C_h12*(C12*C12-C11*C22)*2.0+(C22*C22)*C_h22*(C12*C12-C11*C22)*2.0-(C22*C22*C22)*C_h22*(C12*C12-C11*C22)+C22*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)+C22*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+C22*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)-(C22*C22)*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-(C22*C22)*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+(C22*C22*C22)*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)-C22*C_h11*(C12*C12-C11*C22)-C22*C_h12*(C12*C12-C11*C22)*2.0-C22*C_h22*(C12*C12-C11*C22)-C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0-C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0+C11*C22*C_h11*(C12*C12-C11*C22)*2.0+C11*C22*C_h12*(C12*C12-C11*C22)*2.0+C12*C22*C_h13*(C12*C12-C11*C22)*4.0+C12*C22*C_h23*(C12*C12-C11*C22)*4.0+C11*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0+C12*C_h13*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*8.0+C22*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0-(C11*C11)*C22*C_h11*(C12*C12-C11*C22)-C11*(C22*C22)*C_h12*(C12*C12-C11*C22)*2.0-C12*(C22*C22)*C_h23*(C12*C12-C11*C22)*4.0-(C12*C12)*C22*C_h33*(C12*C12-C11*C22)*4.0-C11*C22*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-C11*C22*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-C12*C22*C_h13*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0-C12*C22*C_h23*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0+(C11*C11)*C22*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)+C11*(C22*C22)*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+C12*(C22*C22)*C_h23*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0+(C12*C12)*C22*C_h33*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0-C11*C12*C22*C_h13*(C12*C12-C11*C22)*4.0+C11*C12*C22*C_h13*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0)*(-1.0/1.6E+1);
        dW_dC22=1.0/pow((sqrt(-C12*C12+C11*C22)),3.0)*1.0/pow(-C12*C12+C11*C22,3.0/2.0)*((C11*C11)*C_h11*(C12*C12-C11*C22)*2.0+(C11*C11)*C_h12*(C12*C12-C11*C22)*2.0-(C11*C11*C11)*C_h11*(C12*C12-C11*C22)+C11*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)+C11*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+C11*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)-(C11*C11)*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-(C11*C11)*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+(C11*C11*C11)*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)-C11*C_h11*(C12*C12-C11*C22)-C11*C_h12*(C12*C12-C11*C22)*2.0-C11*C_h22*(C12*C12-C11*C22)-C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0-C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0+C11*C12*C_h13*(C12*C12-C11*C22)*4.0+C11*C22*C_h12*(C12*C12-C11*C22)*2.0+C11*C12*C_h23*(C12*C12-C11*C22)*4.0+C11*C22*C_h22*(C12*C12-C11*C22)*2.0+C11*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0+C12*C_h23*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*8.0+C22*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0-(C11*C11)*C12*C_h13*(C12*C12-C11*C22)*4.0-(C11*C11)*C22*C_h12*(C12*C12-C11*C22)*2.0-C11*(C22*C22)*C_h22*(C12*C12-C11*C22)-C11*(C12*C12)*C_h33*(C12*C12-C11*C22)*4.0-C11*C12*C_h13*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0-C11*C22*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-C11*C12*C_h23*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0-C11*C22*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+(C11*C11)*C12*C_h13*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0+(C11*C11)*C22*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+C11*(C22*C22)*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)+C11*(C12*C12)*C_h33*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0-C11*C12*C22*C_h23*(C12*C12-C11*C22)*4.0+C11*C12*C22*C_h23*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0)*(-1.0/1.6E+1);
        dW_dC12=1.0/pow((sqrt(-C12*C12+C11*C22)),3.0)*1.0/pow(-C12*C12+C11*C22,3.0/2.0)*(C_h13*pow(C12*C12-C11*C22,2.0)*-4.0-C_h23*pow(C12*C12-C11*C22,2.0)*4.0+C11*C_h13*pow(C12*C12-C11*C22,2.0)*4.0+C12*C_h33*pow(C12*C12-C11*C22,2.0)*4.0+C22*C_h23*pow(C12*C12-C11*C22,2.0)*4.0-C12*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)-C12*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-C12*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)+C12*C_h11*(C12*C12-C11*C22)+C12*C_h12*(C12*C12-C11*C22)*2.0+C12*C_h22*(C12*C12-C11*C22)-C11*C12*C_h11*(C12*C12-C11*C22)*2.0-C11*C12*C_h12*(C12*C12-C11*C22)*2.0-C11*C22*C_h13*(C12*C12-C11*C22)*4.0-C12*C22*C_h12*(C12*C12-C11*C22)*2.0-C11*C22*C_h23*(C12*C12-C11*C22)*4.0-C12*C22*C_h22*(C12*C12-C11*C22)*2.0+C12*C_h33*pow((sqrt(-C12*C12+C11*C22)),2.0)*(C12*C12-C11*C22)*4.0+(C11*C11)*C12*C_h11*(C12*C12-C11*C22)+(C11*C11)*C22*C_h13*(C12*C12-C11*C22)*4.0+C11*(C22*C22)*C_h23*(C12*C12-C11*C22)*4.0+C12*(C22*C22)*C_h22*(C12*C12-C11*C22)+C11*C12*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+C11*C12*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+C11*C22*C_h13*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0+C12*C22*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0+C11*C22*C_h23*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0+C12*C22*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-(C11*C11)*C12*C_h11*pow((sqrt(-C12*C12+C11*C22)),2.0)-(C11*C11)*C22*C_h13*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0-C11*(C22*C22)*C_h23*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0-C12*(C22*C22)*C_h22*pow((sqrt(-C12*C12+C11*C22)),2.0)+C11*C12*C22*C_h12*(C12*C12-C11*C22)*2.0+C11*C12*C22*C_h33*(C12*C12-C11*C22)*4.0-C11*C12*C22*C_h12*pow((sqrt(-C12*C12+C11*C22)),2.0)*2.0-C11*C12*C22*C_h33*pow((sqrt(-C12*C12+C11*C22)),2.0)*4.0)*(-1.0/8.0);

        // derivative of the volumetric part
        dW_vol_dC11 = K*(C22+C22/(C12*C12-C11*C22));
        dW_vol_dC22 = K*(C11+C11/(C12*C12-C11*C22));
        dW_vol_dC12 = K*(-2.*C12-(2.*C12)/(C12*C12-C11*C22));
    }
    else if(en_par.en_type=="poly")
    {
        // beta is a parameter to define the symmetry
        double beta;
        if (en_par.sym == "sq")
            beta=-0.25;
        else if(en_par.sym == "hex")
            beta=4.;
        else
            beta=0.;

        // energy density
        W_temp=(1.0/pow(C12*C12-C11*C22,3.0)*(C11*(C12*C12*C12*C12*C12)*1.568E+3-C11*(C22*C22*C22*C22*C22)*8.0-(C11*C11*C11*C11*C11)*C22*8.0+(C12*C12*C12*C12*C12)*C22*1.568E+3+(C12*C12*C12*C12*C12*C12)*beta*2.48E+2-(C12*C12*C12*C12*C12*C12)*8.96E+2-(C11*C11)*(C12*C12*C12*C12)*6.24E+2-(C11*C11*C11)*(C12*C12*C12)*1.6E+2+(C11*C11*C11*C11)*(C12*C12)*4.0E+1-(C11*C11)*(C22*C22*C22*C22)*1.4E+1+(C11*C11*C11)*(C22*C22*C22)*6.0E+1-(C11*C11*C11*C11)*(C22*C22)*1.4E+1+(C12*C12)*(C22*C22*C22*C22)*4.0E+1-(C12*C12*C12)*(C22*C22*C22)*1.6E+2-(C12*C12*C12*C12)*(C22*C22)*6.24E+2+(C11*C11)*(C12*C12*C12*C12)*beta*3.66E+2-(C11*C11*C11)*(C12*C12*C12)*beta*1.6E+2+(C11*C11*C11*C11)*(C12*C12)*beta*4.0E+1+(C11*C11)*(C22*C22*C22*C22)*beta*1.9E+1-(C11*C11*C11)*(C22*C22*C22)*beta*2.8E+1+(C11*C11*C11*C11)*(C22*C22)*beta*1.9E+1+(C12*C12)*(C22*C22*C22*C22)*beta*4.0E+1-(C12*C12*C12)*(C22*C22*C22)*beta*1.6E+2+(C12*C12*C12*C12)*(C22*C22)*beta*3.66E+2-(C11*C11)*(C12*C12)*(C22*C22)*9.6E+2-C11*(C12*C12*C12*C12)*C22*2.88E+3-C11*(C12*C12*C12*C12*C12)*beta*4.12E+2-C11*(C22*C22*C22*C22*C22)*beta*8.0-(C11*C11*C11*C11*C11)*C22*beta*8.0-(C12*C12*C12*C12*C12)*C22*beta*4.12E+2+C11*(C12*C12)*(C22*C22*C22)*1.36E+2+C11*(C12*C12*C12)*(C22*C22)*1.4E+3+(C11*C11)*(C12*C12*C12)*C22*1.4E+3+(C11*C11*C11)*(C12*C12)*C22*1.36E+2+C11*(C12*C12)*(C22*C22*C22)*beta*4.0-C11*(C12*C12*C12)*(C22*C22)*beta*1.84E+2-(C11*C11)*(C12*C12*C12)*C22*beta*1.84E+2+(C11*C11*C11)*(C12*C12)*C22*beta*4.0+(C11*C11)*(C12*C12)*(C22*C22)*beta*9.6E+1+C11*(C12*C12*C12*C12)*C22*beta*3.54E+2))/1.98E+2;
        // volumetric part W_v=K*(detC-log(detC)), K - bulk modulus
        W_v_temp=-K*(log(-C12*C12+C11*C22)+C12*C12-C11*C22);

        // derivatives (from symbolic differentiation)
        dW_dC11=1.0/pow(C12*C12-C11*C22,4.0)*(C11*(C12*C12*C12*C12*C12*C12)*1.248E+3+C11*(C22*C22*C22*C22*C22*C22)*1.6E+1+(C12*C12*C12*C12*C12*C12)*C22*5.568E+3+(C12*C12*C12*C12*C12*C12*C12)*beta*4.12E+2-(C12*C12*C12*C12*C12*C12*C12)*1.568E+3+(C11*C11)*(C12*C12*C12*C12*C12)*4.8E+2-(C11*C11*C11)*(C12*C12*C12*C12)*1.6E+2+(C11*C11)*(C22*C22*C22*C22*C22)*1.4E+1-(C11*C11*C11*C11)*(C22*C22*C22)*1.4E+1-(C11*C11*C11*C11*C11)*(C22*C22)*1.6E+1-(C12*C12)*(C22*C22*C22*C22*C22)*1.12E+2+(C12*C12*C12)*(C22*C22*C22*C22)*4.8E+2+(C12*C12*C12*C12)*(C22*C22*C22)*1.736E+3-(C12*C12*C12*C12*C12)*(C22*C22)*6.104E+3+(C11*C11)*(C12*C12*C12*C12*C12)*beta*4.8E+2-(C11*C11*C11)*(C12*C12*C12*C12)*beta*1.6E+2-(C11*C11)*(C22*C22*C22*C22*C22)*beta*1.9E+1+(C11*C11*C11*C11)*(C22*C22*C22)*beta*1.9E+1-(C11*C11*C11*C11*C11)*(C22*C22)*beta*1.6E+1-(C12*C12)*(C22*C22*C22*C22*C22)*beta*1.12E+2+(C12*C12*C12)*(C22*C22*C22*C22)*beta*4.8E+2-(C12*C12*C12*C12)*(C22*C22*C22)*beta*1.102E+3+(C12*C12*C12*C12*C12)*(C22*C22)*beta*1.42E+3+(C11*C11)*(C12*C12)*(C22*C22*C22)*7.8E+2-(C11*C11)*(C12*C12*C12)*(C22*C22)*1.4E+3+(C11*C11*C11)*(C12*C12)*(C22*C22)*5.6E+1-C11*(C12*C12*C12*C12*C12)*C22*5.936E+3-C11*(C12*C12*C12*C12*C12*C12)*beta*7.32E+2+C11*(C22*C22*C22*C22*C22*C22)*beta*1.6E+1-(C12*C12*C12*C12*C12*C12)*C22*beta*1.098E+3-C11*(C12*C12)*(C22*C22*C22*C22)*2.44E+2-C11*(C12*C12*C12)*(C22*C22*C22)*2.8E+3+C11*(C12*C12*C12*C12)*(C22*C22)*7.68E+3+(C11*C11)*(C12*C12*C12*C12)*C22*2.16E+2+(C11*C11*C11*C11)*(C12*C12)*C22*8.0E+1-C11*(C12*C12)*(C22*C22*C22*C22)*beta*4.6E+1+C11*(C12*C12*C12)*(C22*C22*C22)*beta*3.68E+2-C11*(C12*C12*C12*C12)*(C22*C22)*beta*9.0E+2-(C11*C11)*(C12*C12*C12*C12)*C22*beta*3.78E+2+(C11*C11*C11*C11)*(C12*C12)*C22*beta*8.0E+1-(C11*C11)*(C12*C12)*(C22*C22*C22)*beta*1.2E+1+(C11*C11)*(C12*C12*C12)*(C22*C22)*beta*1.84E+2-(C11*C11*C11)*(C12*C12)*(C22*C22)*beta*7.6E+1+C11*(C12*C12*C12*C12*C12)*C22*beta*1.192E+3)*(-1.0/1.98E+2);
        dW_dC22=1.0/pow(C12*C12-C11*C22,4.0)*(C11*(C12*C12*C12*C12*C12*C12)*5.568E+3+(C11*C11*C11*C11*C11*C11)*C22*1.6E+1+(C12*C12*C12*C12*C12*C12)*C22*1.248E+3+(C12*C12*C12*C12*C12*C12*C12)*beta*4.12E+2-(C12*C12*C12*C12*C12*C12*C12)*1.568E+3-(C11*C11)*(C12*C12*C12*C12*C12)*6.104E+3+(C11*C11*C11)*(C12*C12*C12*C12)*1.736E+3+(C11*C11*C11*C11)*(C12*C12*C12)*4.8E+2-(C11*C11*C11*C11*C11)*(C12*C12)*1.12E+2-(C11*C11)*(C22*C22*C22*C22*C22)*1.6E+1-(C11*C11*C11)*(C22*C22*C22*C22)*1.4E+1+(C11*C11*C11*C11*C11)*(C22*C22)*1.4E+1-(C12*C12*C12*C12)*(C22*C22*C22)*1.6E+2+(C12*C12*C12*C12*C12)*(C22*C22)*4.8E+2+(C11*C11)*(C12*C12*C12*C12*C12)*beta*1.42E+3-(C11*C11*C11)*(C12*C12*C12*C12)*beta*1.102E+3+(C11*C11*C11*C11)*(C12*C12*C12)*beta*4.8E+2-(C11*C11*C11*C11*C11)*(C12*C12)*beta*1.12E+2-(C11*C11)*(C22*C22*C22*C22*C22)*beta*1.6E+1+(C11*C11*C11)*(C22*C22*C22*C22)*beta*1.9E+1-(C11*C11*C11*C11*C11)*(C22*C22)*beta*1.9E+1-(C12*C12*C12*C12)*(C22*C22*C22)*beta*1.6E+2+(C12*C12*C12*C12*C12)*(C22*C22)*beta*4.8E+2+(C11*C11)*(C12*C12)*(C22*C22*C22)*5.6E+1-(C11*C11)*(C12*C12*C12)*(C22*C22)*1.4E+3+(C11*C11*C11)*(C12*C12)*(C22*C22)*7.8E+2-C11*(C12*C12*C12*C12*C12)*C22*5.936E+3-C11*(C12*C12*C12*C12*C12*C12)*beta*1.098E+3+(C11*C11*C11*C11*C11*C11)*C22*beta*1.6E+1-(C12*C12*C12*C12*C12*C12)*C22*beta*7.32E+2+C11*(C12*C12)*(C22*C22*C22*C22)*8.0E+1+C11*(C12*C12*C12*C12)*(C22*C22)*2.16E+2+(C11*C11)*(C12*C12*C12*C12)*C22*7.68E+3-(C11*C11*C11)*(C12*C12*C12)*C22*2.8E+3-(C11*C11*C11*C11)*(C12*C12)*C22*2.44E+2+C11*(C12*C12)*(C22*C22*C22*C22)*beta*8.0E+1-C11*(C12*C12*C12*C12)*(C22*C22)*beta*3.78E+2-(C11*C11)*(C12*C12*C12*C12)*C22*beta*9.0E+2+(C11*C11*C11)*(C12*C12*C12)*C22*beta*3.68E+2-(C11*C11*C11*C11)*(C12*C12)*C22*beta*4.6E+1-(C11*C11)*(C12*C12)*(C22*C22*C22)*beta*7.6E+1+(C11*C11)*(C12*C12*C12)*(C22*C22)*beta*1.84E+2-(C11*C11*C11)*(C12*C12)*(C22*C22)*beta*1.2E+1+C11*(C12*C12*C12*C12*C12)*C22*beta*1.192E+3)*(-1.0/1.98E+2);
        dW_dC12=C12*1.0/pow(C12*C12-C11*C22,4.0)*(C11*(C12*C12*C12*C12*C12)*7.84E+2+C11*(C22*C22*C22*C22*C22)*1.6E+1+(C11*C11*C11*C11*C11)*C22*1.6E+1+(C12*C12*C12*C12*C12)*C22*7.84E+2-(C11*C11)*(C12*C12*C12*C12)*6.24E+2-(C11*C11*C11)*(C12*C12*C12)*2.4E+2+(C11*C11*C11*C11)*(C12*C12)*8.0E+1+(C11*C11)*(C22*C22*C22*C22)*9.4E+1-(C11*C11*C11)*(C22*C22*C22)*7.8E+2+(C11*C11*C11*C11)*(C22*C22)*9.4E+1+(C12*C12)*(C22*C22*C22*C22)*8.0E+1-(C12*C12*C12)*(C22*C22*C22)*2.4E+2-(C12*C12*C12*C12)*(C22*C22)*6.24E+2+(C11*C11)*(C12*C12*C12*C12)*beta*3.66E+2-(C11*C11*C11)*(C12*C12*C12)*beta*2.4E+2+(C11*C11*C11*C11)*(C12*C12)*beta*8.0E+1+(C11*C11)*(C22*C22*C22*C22)*beta*6.1E+1+(C11*C11*C11)*(C22*C22*C22)*beta*1.2E+1+(C11*C11*C11*C11)*(C22*C22)*beta*6.1E+1+(C12*C12)*(C22*C22*C22*C22)*beta*8.0E+1-(C12*C12*C12)*(C22*C22*C22)*beta*2.4E+2+(C12*C12*C12*C12)*(C22*C22)*beta*3.66E+2-(C11*C11)*(C12*C12)*(C22*C22)*7.68E+3-C11*C12*(C22*C22*C22*C22)*2.4E+2-C11*(C12*C12*C12*C12)*C22*5.568E+3-(C11*C11*C11*C11)*C12*C22*2.4E+2-C11*(C12*C12*C12*C12*C12)*beta*2.06E+2+C11*(C22*C22*C22*C22*C22)*beta*1.6E+1+(C11*C11*C11*C11*C11)*C22*beta*1.6E+1-(C12*C12*C12*C12*C12)*C22*beta*2.06E+2-C11*(C12*C12)*(C22*C22*C22)*9.76E+2+C11*(C12*C12*C12)*(C22*C22)*6.02E+3+(C11*C11)*C12*(C22*C22*C22)*2.1E+3+(C11*C11)*(C12*C12*C12)*C22*6.02E+3+(C11*C11*C11)*C12*(C22*C22)*2.1E+3-(C11*C11*C11)*(C12*C12)*C22*9.76E+2+C11*(C12*C12)*(C22*C22*C22)*beta*7.4E+2-C11*(C12*C12*C12)*(C22*C22)*beta*1.306E+3-(C11*C11)*C12*(C22*C22*C22)*beta*2.76E+2-(C11*C11)*(C12*C12*C12)*C22*beta*1.306E+3-(C11*C11*C11)*C12*(C22*C22)*beta*2.76E+2+(C11*C11*C11)*(C12*C12)*C22*beta*7.4E+2+(C11*C11)*(C12*C12)*(C22*C22)*beta*9.0E+2-C11*C12*(C22*C22*C22*C22)*beta*2.4E+2+C11*(C12*C12*C12*C12)*C22*beta*1.098E+3-(C11*C11*C11*C11)*C12*C22*beta*2.4E+2)*(-1.0/9.9E+1);

        // derivative of the volumetric part
        dW_vol_dC11 = K*(C22+C22/(C12*C12-C11*C22));
        dW_vol_dC22 = K*(C11+C11/(C12*C12-C11*C22));
        dW_vol_dC12 = K*(-2.*C12-(2.*C12)/(C12*C12-C11*C22));
    }

    // define the structure to return the results
    en_density en_den;
    // the total energy density
    en_den.W=W_temp+W_v_temp;
    // matrix of derivatives
    en_den.dW(0,0)=dW_dC11+dW_vol_dC11;
    en_den.dW(0,1)=0.5*(dW_dC12+dW_vol_dC12);
    en_den.dW(1,0)=0.5*(dW_dC12+dW_vol_dC12);
    en_den.dW(1,1)=dW_dC22+dW_vol_dC22;

    return(en_den);

}









// initialization of the energy density and its derivative for an element
en_density en_density_element(en_parameters en_par,MatrixXd p, Vector3i T, VectorXd U)
{
    ////////////// input //////////////
    // en_par - structure with energy density properties
    // p - points of the mesh
    // T - connectivity list
    // U - displacements


    ////////////// output //////////////
    // en_den - structure of the type 'en_density' (see structures.h)

    // define the structure for the deformation gradient
    def_grad F_el;
    // define the structure for the energy density
    en_density E_el;
    // matrix for the 1st 1st Piola-Kirchhof stress
    Matrix2d P;
    // turn P into a vector
    Vector4d P_vec;


    // call the deformation gradient
    F_el=def_grad_init(p, T, U);
    // call the energy within an element
    E_el=en_density_init(en_par,F_el.C_r);



    // the 1st Piola-Kirchhof stress
    P=2.*F_el.F*F_el.m*E_el.dW*F_el.m.transpose();
    // initialise the corresponding vector notation
    P_vec(0)=P(0,0);
    P_vec(1)=P(1,0);
    P_vec(2)=P(0,1);
    P_vec(3)=P(1,1);

    // define the structure to return the results
    en_density en_den;

    en_den.W=F_el.A*E_el.W;
    en_den.dW=E_el.dW;
    en_den.dW_T=F_el.A*P_vec.transpose()*F_el.dF; // derivatives for each elements, every row is [dE/du_1 dE/dv_1 dE/du_2 dE/dv_2 dE/du_3 dE/dv_3]


    return en_den;
}







// initialization of the total energy
gl_energy global_energy_init(en_parameters en_par,MatrixXd p, MatrixXi T, VectorXd U)
{
    ////////////// input //////////////
    // en_par - structure with energy density properties
    // p - points of the mesh
    // T - connectivity list
    // U - displacements


    ////////////// output //////////////
    // E_tot - structure of the type 'gl_energy' (see structures.h)

    double E_val=0.;
    VectorXd dE_val=VectorXd::Zero(U.size());

    // variable for the loop
    long int N=T.rows();
    long int i;
    int  k;

    // set up the parallelisation
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        // define the structure for the energy density of each element
        en_density E_el;
        #pragma omp for private(i,k,E_el) reduction(+:E_val)
        for(i=0;i<N;i++)
        {
            // call the energy within an element
            E_el=en_density_element(en_par,p, T(i,all), U);
            // accumulate the total energy
            E_val+=E_el.W;

            // accumulate values for the gradient
            for(k=0;k<T.cols();k++)
            {
                dE_val(2*T(i,k))+=E_el.dW_T(2*k);
                dE_val(2*T(i,k)+1)+=E_el.dW_T(2*k+1);
            }

        }
    }
    // define the structure to return the results
    gl_energy E_tot;
    // initialise the structure
    E_tot.E=E_val;
    E_tot.dE=dE_val;

    return(E_tot);


}









// initialization of the total energy with constrained displacements
gl_energy global_energy_constr_init(en_parameters en_par,MatrixXd p, MatrixXi T, BCs BCs_val, BCs_prop B_pro, VectorXd W)
{
    ////////////// input //////////////
    // en_par - structure with energy density properties
    // p - points of the mesh
    // T - connectivity list
    // BCs_val - boundary conditions with indices of free and constrained nodes
    // B_pro - properties of the boundary conditions
    // W - displacements (sought during the optimisation)

    ////////////// output //////////////
    // E_tot_constr - structure of the type 'gl_energy' (see structures.h)
    // define the free, uknown nodes
    VectorXd U(2*p.rows());

    // new points
    MatrixXd p_new;
    if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
    {
      // deformation gradient to load
      Matrix2d F_load;
      F_load<<1.,BCs_val.s,
              0.,1.;
      p_new = p*F_load.transpose();
      /*
      for(int i=0;i<BCs_val.ind_free.size();i++)
      {
        if(BCs_val.ind_free(i)%2)
          U(BCs_val.ind_free(i))=p_new(BCs_val.ind_free(i)/2,0)-p(BCs_val.ind_free(i)/2,0)+W(i);
        else
          U(BCs_val.ind_free(i))=p_new((BCs_val.ind_free(i)-1)/2,1)-p((BCs_val.ind_free(i)-1)/2,1)+W(i);
      }

      // deformation gradient to load
      Matrix2d F_load;
      F_load<<1.,BCs_val.s,
              0.,1.;
      p_new = p*F_load.transpose();
      for(int i=0;i<BCs_val.ind_free.size();i++)
      {
        if(BCs_val.ind_free(i)%2)
          U(BCs_val.ind_free(i))=W(i);
        else
          U(BCs_val.ind_free(i))=W(i);
      }
      */
    }
    //else
    //  U(BCs_val.ind_free)=W;

    U(BCs_val.ind_free)=W;
    // define the known displacements

    // TOP or BOTTOM - HORIZONTAL DISPLACEMENTS
    if(B_pro.u_ver_per==0)
    {
        U(2*BCs_val.top.array())=BCs_val.u0.top;
        U(2*BCs_val.bottom.array())=BCs_val.u0.bottom;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.top.array())=U(2*BCs_val.bottom.array())+BCs_val.u0.top-BCs_val.u0.bottom;
			else
				U(2*BCs_val.top.array())=U(2*BCs_val.bottom.array());
		}
    // TOP or BOTTOM - VERTICAL DISPLACEMENTS
    if(B_pro.v_ver_per==0)
    {
        U(2*BCs_val.top.array()+1)=BCs_val.v0.top;
        U(2*BCs_val.bottom.array()+1)=BCs_val.v0.bottom;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.top.array()+1)=U(2*BCs_val.bottom.array()+1)+BCs_val.v0.top-BCs_val.v0.bottom;
			else
				U(2*BCs_val.top.array()+1)=U(2*BCs_val.bottom.array()+1);
		}
    // LEFT or RIGHT - HORIZONTAL DISPLACEMENTS
    if(B_pro.u_hor_per==0)
    {
        U(2*BCs_val.left.array())=BCs_val.u0.left;
        U(2*BCs_val.right.array())=BCs_val.u0.right;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.right.array())=U(2*BCs_val.left.array())+BCs_val.u0.right-BCs_val.u0.left;
			else
				U(2*BCs_val.right.array())=U(2*BCs_val.left.array());
		}
    // LEFT or RIGHT - VERTICAL DISPLACEMENTS
    if(B_pro.v_hor_per==0)
    {
        U(2*BCs_val.left.array()+1)=BCs_val.v0.left;
        U(2*BCs_val.right.array()+1)=BCs_val.v0.right;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.right.array()+1)=U(2*BCs_val.left.array()+1)+BCs_val.v0.right-BCs_val.v0.left;
			else
				U(2*BCs_val.right.array()+1)=U(2*BCs_val.left.array()+1);
		}
    if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
    {
      // displacements of the RIGHT-BOTTOM corner point
      // HORIZONTAL
      U(2*(BCs_val.bottom(last)+1))=U(0)+p_new(BCs_val.bottom(last)+1,0)-p(BCs_val.bottom(last)+1,0)-(p_new(0,0)-p(0,0));
      // VERTICAL
      U(2*(BCs_val.bottom(last)+1)+1)=U(1)+p_new(BCs_val.bottom(last)+1,1)-p(BCs_val.bottom(last)+1,1)-(p_new(0,1)-p(0,1));

      // displacements of the LEFT-TOP corner point
      // HORIZONTAL
      U(2*(BCs_val.right(last)+1))=U(0)+p_new(BCs_val.right(last)+1,0)-p(BCs_val.right(last)+1,0)-(p_new(0,0)-p(0,0));
      // VERTICAL
      U(2*(BCs_val.right(last)+1)+1)=U(1)+p_new(BCs_val.right(last)+1,1)-p(BCs_val.right(last)+1,1)-(p_new(0,1)-p(0,1));

      // displacements of the RIGHT-TOP corner point
      // HORIZONTAL
      U(2*(BCs_val.top(last)+1))=U(0)+p_new(BCs_val.top(last)+1,0)-p(BCs_val.top(last)+1,0)-(p_new(0,0)-p(0,0));
      // VERTICAL
      U(2*(BCs_val.top(last)+1)+1)=U(1)+p_new(BCs_val.top(last)+1,1)-p(BCs_val.top(last)+1,1)-(p_new(0,1)-p(0,1));

    }
    // call the global energy function
    gl_energy E_tot=global_energy_init(en_par,p, T, U);


    /////////////  additional terms due to periodicity  /////////////

    // TOP or BOTTOM - HORIZONTAL DISPLACEMENTS
    if(B_pro.u_ver_per==1)
        E_tot.dE(2*BCs_val.bottom.array())+=E_tot.dE(2*BCs_val.top.array());
    // TOP or BOTTOM - VERTICAL DISPLACEMENTS
    if(B_pro.v_ver_per==1)
        E_tot.dE(2*BCs_val.bottom.array()+1)+=E_tot.dE(2*BCs_val.top.array()+1);
    // LEFT or RIGHT - HORIZONTAL DISPLACEMENTS
    if(B_pro.u_hor_per==1)
        E_tot.dE(2*BCs_val.left.array())+=E_tot.dE(2*BCs_val.right.array());
    // LEFT or RIGHT - VERTICAL DISPLACEMENTS
    if(B_pro.v_hor_per==1);
        E_tot.dE(2*BCs_val.left.array()+1)+=E_tot.dE(2*BCs_val.right.array()+1);
    if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
        {
          E_tot.dE(0)+=E_tot.dE(2*(BCs_val.bottom(last)+1))+E_tot.dE(2*(BCs_val.right(last)+1))+E_tot.dE(2*(BCs_val.top(last)+1));
          E_tot.dE(1)+=E_tot.dE(2*(BCs_val.bottom(last)+1)+1)+E_tot.dE(2*(BCs_val.right(last)+1)+1)+E_tot.dE(2*(BCs_val.top(last)+1)+1);
        }
    // initialise the results
    gl_energy E_tot_constr;
    E_tot_constr.E=E_tot.E;
    E_tot_constr.dE=E_tot.dE(BCs_val.ind_free);
return(E_tot_constr);
}









// initialization the displacements
VectorXd total_displ(MatrixXd p,MatrixXi T, BCs BCs_val, BCs_prop B_pro, VectorXd W)
{
    ////////////// input //////////////
    // en_par - structure with energy density properties
    // p - points of the mesh
    // T - connectivity list
    // BCs_val - boundary conditions with indices of free and constrained nodes
    // B_pro - properties of the boundary conditions
    // W - displacements (sought during the optimisation)

    ////////////// output //////////////
    // U - vector of all displacements of the lattice (see structures.h)

    // define the free, uknown nodes
    VectorXd U(2*p.rows());

    // new points
    MatrixXd p_new;
    if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
    {
      // deformation gradient to load
      Matrix2d F_load;
      F_load<<1.,BCs_val.s,
              0.,1.;
      p_new = p*F_load.transpose();
      /*
      for(int i=0;i<BCs_val.ind_free.size();i++)
      {
        if(BCs_val.ind_free(i)%2)
          U(BCs_val.ind_free(i))=p_new(BCs_val.ind_free(i)/2,0)-p(BCs_val.ind_free(i)/2,0)+W(i);
        else
          U(BCs_val.ind_free(i))=p_new((BCs_val.ind_free(i)-1)/2,1)-p((BCs_val.ind_free(i)-1)/2,1)+W(i);
      }

      // deformation gradient to load
      Matrix2d F_load;
      F_load<<1.,BCs_val.s,
              0.,1.;
      p_new = p*F_load.transpose();
      for(int i=0;i<BCs_val.ind_free.size();i++)
      {
        if(BCs_val.ind_free(i)%2)
          U(BCs_val.ind_free(i))=W(i);
        else
          U(BCs_val.ind_free(i))=W(i);
      }
      */
    }
    //else
    //  U(BCs_val.ind_free)=W;

    U(BCs_val.ind_free)=W;
    // define the known displacements

    // TOP or BOTTOM - HORIZONTAL DISPLACEMENTS
    if(B_pro.u_ver_per==0)
    {
        U(2*BCs_val.top.array())=BCs_val.u0.top;
        U(2*BCs_val.bottom.array())=BCs_val.u0.bottom;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.top.array())=U(2*BCs_val.bottom.array())+BCs_val.u0.top-BCs_val.u0.bottom;
			else
				U(2*BCs_val.top.array())=U(2*BCs_val.bottom.array());
		}
    // TOP or BOTTOM - VERTICAL DISPLACEMENTS
    if(B_pro.v_ver_per==0)
    {
        U(2*BCs_val.top.array()+1)=BCs_val.v0.top;
        U(2*BCs_val.bottom.array()+1)=BCs_val.v0.bottom;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.top.array()+1)=U(2*BCs_val.bottom.array()+1)+BCs_val.v0.top-BCs_val.v0.bottom;
			else
				U(2*BCs_val.top.array()+1)=U(2*BCs_val.bottom.array()+1);
		}
    // LEFT or RIGHT - HORIZONTAL DISPLACEMENTS
    if(B_pro.u_hor_per==0)
    {
        U(2*BCs_val.left.array())=BCs_val.u0.left;
        U(2*BCs_val.right.array())=BCs_val.u0.right;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.right.array())=U(2*BCs_val.left.array())+BCs_val.u0.right-BCs_val.u0.left;
			else
				U(2*BCs_val.right.array())=U(2*BCs_val.left.array());
		}
    // LEFT or RIGHT - VERTICAL DISPLACEMENTS
    if(B_pro.v_hor_per==0)
    {
        U(2*BCs_val.left.array()+1)=BCs_val.v0.left;
        U(2*BCs_val.right.array()+1)=BCs_val.v0.right;
    }
    else
		{
			if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
				U(2*BCs_val.right.array()+1)=U(2*BCs_val.left.array()+1)+BCs_val.v0.right-BCs_val.v0.left;
			else
				U(2*BCs_val.right.array()+1)=U(2*BCs_val.left.array()+1);
		}
    if ((B_pro.u_ver_per==1) && (B_pro.u_hor_per==1) && (B_pro.v_ver_per==1) && (B_pro.v_hor_per==1))
    {
      // displacements of the RIGHT-BOTTOM corner point
      // HORIZONTAL
      U(2*(BCs_val.bottom(last)+1))=U(0)+p_new(BCs_val.bottom(last)+1,0)-p(BCs_val.bottom(last)+1,0)-(p_new(0,0)-p(0,0));
      // VERTICAL
      U(2*(BCs_val.bottom(last)+1)+1)=U(1)+p_new(BCs_val.bottom(last)+1,1)-p(BCs_val.bottom(last)+1,1)-(p_new(0,1)-p(0,1));

      // displacements of the LEFT-TOP corner point
      // HORIZONTAL
      U(2*(BCs_val.right(last)+1))=U(0)+p_new(BCs_val.right(last)+1,0)-p(BCs_val.right(last)+1,0)-(p_new(0,0)-p(0,0));
      // VERTICAL
      U(2*(BCs_val.right(last)+1)+1)=U(1)+p_new(BCs_val.right(last)+1,1)-p(BCs_val.right(last)+1,1)-(p_new(0,1)-p(0,1));

      // displacements of the RIGHT-TOP corner point
      // HORIZONTAL
      U(2*(BCs_val.top(last)+1))=U(0)+p_new(BCs_val.top(last)+1,0)-p(BCs_val.top(last)+1,0)-(p_new(0,0)-p(0,0));
      // VERTICAL
      U(2*(BCs_val.top(last)+1)+1)=U(1)+p_new(BCs_val.top(last)+1,1)-p(BCs_val.top(last)+1,1)-(p_new(0,1)-p(0,1));

    }
    return U;
}
