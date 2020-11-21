#include "optimization.h"
#include <omp.h>

#include <iostream>
#include <iomanip>



#include "calculation.h"
#include "functions_init.h"

#include "Eigen/Dense"

using namespace Eigen;
using namespace alglib;

#include <fstream>
#include <chrono>

// libraries to create new directories for results
#include <stdio.h>
#include <sys/stat.h>
#include <sstream>

// the number of threads
int NUM_THREADS;
// size of the grid in x  and y directions
int n_x;
int n_y;
// separation between the points
double h0;
// properties of the boundary conditions
BCs_prop BC_pro;
// intial value of the loading parameter
double s;
// loading step
double ds;
// maximum load
double s_max;
// energy parameters
en_parameters en_par;
// parameters of optimization
double epsg;
double epsf;
double epsx;
int maxits;

void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
    // define the new vector from x to calculate the energy
    VectorXd W((long) x.length());
    // get the parameters of the model (see structures.h)
    model_param m_p=*(static_cast<model_param*>(ptr));
    // define the other parameter
    for(int i=0;i<W.size();i++)
        W(i)=x[i];
    // get the energy values and its gradient
    gl_energy E_tot=global_energy_constr_init(m_p.en_par,m_p.points,m_p.T,m_p.BCs_val,m_p.BCs_p,W);

    // define the values of function
    func = E_tot.E;
    // its gradient
    for(int i=0;i<E_tot.dE.size();i++)
        grad[i]=E_tot.dE(i);
}

int main(int argc, char **argv)
{
  // load the file with parameters
  std::ifstream file_parameters("./Parameters.txt");
  //read it line by line
  std::string str;
  stringstream os;
  // two additional counters
  int index_holder; // within a line
  int line_count=0; // line number
  while (std::getline(file_parameters, str)) {
    index_holder = 0;
    for(std::string::size_type i = 0; i < str.size(); ++i)
      {
          if (str[i] == '\t'){
              if (line_count==0)
                {os.str(str.substr(index_holder, i - index_holder));os>>NUM_THREADS;os.str("");os.clear();}
              else if (line_count==1)
                {os.str(str.substr(index_holder, i - index_holder));os>>n_x;os.str("");os.clear();}
              else if (line_count==2)
                {os.str(str.substr(index_holder, i - index_holder));os>>n_y;os.str("");os.clear();}
              else if (line_count==3)
                {os.str(str.substr(index_holder, i - index_holder));os>>h0;os.str("");os.clear();}
              else if (line_count==4)
                BC_pro.bc_type=str.substr(index_holder, i - index_holder);
              else if (line_count==5)
                {os.str(str.substr(index_holder, i - index_holder));os>>BC_pro.u_hor_per;os.str("");os.clear();}
              else if (line_count==6)
                {os.str(str.substr(index_holder, i - index_holder));os>>BC_pro.v_hor_per;os.str("");os.clear();}
              else if (line_count==7)
                {os.str(str.substr(index_holder, i - index_holder));os>>BC_pro.u_ver_per;os.str("");os.clear();}
              else if (line_count==8)
                {os.str(str.substr(index_holder, i - index_holder));os>>BC_pro.v_ver_per;os.str("");os.clear();}
              else if (line_count==9)
                {os.str(str.substr(index_holder, i - index_holder));os>>s;os.str("");os.clear();}
              else if (line_count==10)
                {os.str(str.substr(index_holder, i - index_holder));os>>ds;os.str("");os.clear();}
              else if (line_count==11)
                {os.str(str.substr(index_holder, i - index_holder));os>>s_max;os.str("");os.clear();}
              else if (line_count==12)
                en_par.en_type=str.substr(index_holder, i - index_holder);
              else if (line_count==13)
                en_par.sym=str.substr(index_holder, i - index_holder);
              else if (line_count==14)
                {os.str(str.substr(index_holder, i - index_holder));os>>en_par.K;os.str("");os.clear();}
              else if (line_count==15)
                {os.str(str.substr(index_holder, i - index_holder));os>>epsg;os.str("");os.clear();}
              else if (line_count==16)
                {os.str(str.substr(index_holder, i - index_holder));os>>epsf;os.str("");os.clear();}
              else if (line_count==17)
                {os.str(str.substr(index_holder, i - index_holder));os>>epsx;os.str("");os.clear();}
              else if (line_count==18)
                {os.str(str.substr(index_holder, i - index_holder));os>>maxits;os.str("");os.clear();}
              line_count++;
              break;
          }

      }
    }
    // total number of points
    int n_tot=n_x*n_y;
    // x and y coordinates of points
    double p_x[n_x*n_y];
    double p_y[n_x*n_y];
    // points of the grid
    for(int i=0;i<n_y;i++)
        for(int j=0;j<n_x;j++)
        {
          p_x[j+i*n_x]=j*h0;
          p_y[j+i*n_x]=i*h0;
        }
    // put everyting into one matrix
    MatrixXd p(n_tot,2);
    for(int i=0;i<n_tot;i++)
        {
            p(i,0)=p_x[i];
            p(i,1)=p_y[i];
        }
    // make the connectivity matrix
    MatrixXi T(2*(n_x-1)*(n_y-1),3);
	for(int i=0;i<n_y-1;i++)
		for(int j=0;j<n_x-1;j++)
			{
				T(2*(j+i*(n_x-1)),0)=j+i*n_x;
				T(2*(j+i*(n_x-1)),1)=(j+1)+i*n_x;
				T(2*(j+i*(n_x-1)),2)=j+(i+1)*n_x;

				T(2*(j+i*(n_x-1))+1,0)=(j+1)+(i+1)*n_x;
				T(2*(j+i*(n_x-1))+1,1)=j+(i+1)*n_x;
				T(2*(j+i*(n_x-1))+1,2)=(j+1)+i*n_x;
			}
    // generate the geometry structure
    geom G=geom_gen(p,T);

    // deformation gradient for the initial positions
    Matrix2d F_s0;
    F_s0<<1.,s,
          0.,1.;
    // deformation gradient for the increments
    Matrix2d F_ds;
    F_ds<<1.,ds,
          0.,1.;
    // nodes of the current location
    MatrixXd p_cur(n_tot,2);
    p_cur=p*F_s0.transpose();

    // generate the load for boundary conditions
    BCs BCs_load;
    BCs_load=bound_cond(G, BC_pro, s);

    // collect all the parameters into one variable
    model_param mod_par;
    mod_par.en_par=en_par;
    mod_par.points=p;
    mod_par.T=T;
    mod_par.BCs_val=BCs_load;
    mod_par.BCs_p=BC_pro;

    model_param * ptr_mod_par;
    ptr_mod_par=&mod_par;

    // vector of all displacements
    VectorXd U=VectorXd::Zero(2*n_tot);


    for(int i=0;i<n_tot;i++)
    {
      U(2*i)=p_cur(i,0)-p(i,0);
      U(2*i+1)=p_cur(i,1)-p(i,1);
    }

    // vector of free displacements
    VectorXd W(BCs_load.ind_free.size());
    W=U(BCs_load.ind_free);
    // First, we create optimizer object and tune its properties
    real_1d_array x_opt;
    real_1d_array s_opt;
    x_opt.setlength(W.size()); // vector of initial guess
    s_opt.setlength(W.size()); // scaling of components of the function gradient

    for(int i=0;i<W.size();i++)
        s_opt[i]=1.;


    minlbfgsstate state;
    minlbfgsreport rep;

    // set the parameters
    minlbfgscreate(1, x_opt, state);
    minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    minlbfgssetscale(state, s_opt);



	  //// to check the gradient calculation
	  //optguardreport ogrep;
    //minlbfgsoptguardsmoothness(state);
    //minlbfgsoptguardgradient(state, 0.00001);


    // global energy - constrained (which is minimized)
    gl_energy E,E_tot;

    //creating new directories for results
    std::stringstream nx_str,ny_str,set_number_str; // string structures
    nx_str<<n_x; // convert the nx size to a string
    ny_str<<n_y; // convert the nx size to a string
    std::string data_folder="./data"; // main folder with results
    std::string size_folder=data_folder+"/"+nx_str.str()+"x"+ny_str.str(); // folder for each data size
    std::string set_folder=data_folder+"/"+nx_str.str()+"x"+ny_str.str()+"/set "; // folder for each calculation set
    std::string set_folder_new; // structure to create a new folder if necessary

    struct stat sb_0; // structure to check if folders exist
    struct stat sb_1; // structure to check if folders exist
    struct stat sb_2; // structure to check if folders exist

    // check if the main folder exists (create if not)
    if (stat(data_folder.c_str(), &sb_0) == 1 || S_ISDIR(sb_0.st_mode)==0)
        mkdir(data_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // check if the size folder exists (create if not)
    if (stat(size_folder.c_str(), &sb_1) == 1 || S_ISDIR(sb_1.st_mode)==0)
        mkdir(size_folder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // computation set number
    int set_number=1;
    set_number_str<<set_number; // convert to a string
    set_folder_new=set_folder+set_number_str.str(); // find the set number to create a new folder


    while (stat(set_folder_new.c_str(), &sb_2) == 0 && S_ISDIR(sb_2.st_mode)) // while the number exist
	{
		set_number++;
		set_number_str.str("");
    set_number_str.clear();
		set_number_str<<set_number;
		set_folder_new=set_folder+set_number_str.str();
	}
    mkdir(set_folder_new.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // create a directory

    // files with the output results
    std::ofstream file_s(set_folder_new+"/s.dat", std::ios::out);
    std::ofstream file_U(set_folder_new+"/U.dat", std::ios::out);
    // copy the file content of the current calculation.cpp with the settings
    std::string srce_str="./calculation.cpp";
    std::string dest_str=set_folder_new+"/calculation.cpp";
    std::ifstream srce(srce_str.c_str(), std::ios::binary );
    std::ofstream dest(dest_str.c_str(), std::ios::binary ) ;
    dest << srce.rdbuf() ;
    // copy the file content of the current Parameters.txt with the settings
    std::string srce_2_str="./Parameters.txt";
    std::string dest_2_str=set_folder_new+"/Parameters.txt";
    std::ifstream srce_2(srce_2_str.c_str(), std::ios::binary );
    std::ofstream dest_2(dest_2_str.c_str(), std::ios::binary ) ;
    dest_2 << srce_2.rdbuf();
    // count the time
    auto begin_time = std::chrono::system_clock::now();
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds;
    // make the loop
    std::ofstream file_U0(set_folder_new+"/U0.dat", std::ios::out);
    for(int i=0;i<U.size();i++)
        {
            if(i<U.size()-1)
                file_U0 << std::fixed << std::setprecision(7) << U(i) << "\t";
            else
                file_U0 << std::fixed << std::setprecision(7) << U(i) << "\n";
        }
    while(s<=s_max)
    {
        //set  the initial guess

        if ((BC_pro.u_hor_per==1) && (BC_pro.v_hor_per==1) && (BC_pro.u_ver_per==1) && (BC_pro.v_ver_per==1))
  			  {
            for(int i=0;i<n_tot;i++)
            {
              // update the locations
              p_cur(i,0)=p(i,0)+U(2*i);
              p_cur(i,1)=p(i,1)+U(2*i+1);
              // add increment
              p_cur(i,0)=F_ds(0,0)*p_cur(i,0)+F_ds(0,1)*p_cur(i,1);
              p_cur(i,1)=F_ds(1,0)*p_cur(i,0)+F_ds(1,1)*p_cur(i,1);
              // make the initial guess for the displacements
              U(2*i)=p_cur(i,0)-p(i,0);
              U(2*i+1)=p_cur(i,1)-p(i,1);
            }
            // initial guess for the displacements of the free nodes
            W=U(BCs_load.ind_free);
            W+=1.0E-2*ds*VectorXd::Random(BCs_load.ind_free.size());
          }
    		else
    		  W+=1.0E-2*ds*VectorXd::Random(BCs_load.ind_free.size());


        // set up the boundary conditions following the load
        BCs_load=bound_cond(G, BC_pro, s);
        mod_par.BCs_val=BCs_load;

        // set up the initial guess
        for(int i=0;i<W.size();i++)
            x_opt[i]=W(i);

        begin_time = std::chrono::system_clock::now();
        alglib::minlbfgsoptimize(state, function1_grad,NULL,ptr_mod_par);
        minlbfgsresults(state, x_opt, rep);
        end_time = std::chrono::system_clock::now();
        elapsed_seconds=end_time-begin_time;

        for(int i=0;i<W.size();i++)
            W(i)=x_opt[i];
        // calculate the function value
        E=global_energy_constr_init(en_par,p,T,BCs_load,BC_pro,W);

        // get the full vector of displacements
        U=total_displ(p,T,BCs_load, BC_pro, W);
        // get the total energy of the system
        E_tot=global_energy_init(en_par,p,T,U);

        //// to check the gradient calculation
        //minlbfgsoptguardresults(state, ogrep);
	      //printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
	      //printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
	      //printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false


        // output line
        std::cout<<"Loading step "<< s << ", Maximum load " << s_max << ", Norm of the gradient " << E.dE.norm() << ", Energy " << E_tot.E<<std::endl;
        std::cout<<"Time of the optimisation (sec) " << elapsed_seconds.count()<<std::endl;
        std::cout<< (long)rep.terminationtype<<std::endl;
        // save the results into files
        file_s<< std::fixed << std::setprecision(7) << s << "\t" << E_tot.E<<"\n";
        for(int i=0;i<U.size();i++)
            {
                if(i<U.size()-1)
                    file_U << std::fixed << std::setprecision(7) << U(i) << "\t";
                else
                    file_U << std::fixed << std::setprecision(7) << U(i) << "\n";
            }


        // make an increment
        s+=ds;

        minlbfgsrestartfrom(state, x_opt);
    }
    file_s.close();
    file_U.close();
    return 0;
}
