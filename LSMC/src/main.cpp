#include "../include/LSMC.hpp"
#include <iostream>
#include <Eigen/Dense>
#include "gnuplot-iostream.h"
#include <vector>

using namespace std;

int main() {
    int I(10000), M(50), N(5);
    double r(.05), sigma(0.3), S_0(100.), T(1.0), K(90.);

    LSMC lsmc(I, M, N, r, sigma, S_0, T, K);
    Eigen::VectorXd P = lsmc.pricePut();



    Gnuplot gp;

    gp << "set style line 1 linecolor rgb '#00008B' linetype 1 dt 2  linewidth 2\n";
    
    double dt = T / (double)(M - 1);

    vector <pair<double, double> > out;
    
    
    for(int i = 0; i < P.size(); i++) {
        double t_0 = i * dt;
        
        out.push_back(make_pair(t_0, P(i)));
        
    }
    gp << "plot"<< gp.file1d(out) << "w l ls 1 title 'LSMC American Put Price'" << endl;



}