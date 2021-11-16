#include "../include/LSMC.hpp"

using namespace std;

Eigen::VectorXd LSMC::pricePut() const{
        // % Generate trajs
        Eigen::MatrixXd S = trajs();

        double dt = T / (double)(S.cols() - 1);
        double df = exp(-r * dt); // Discount factor
        Eigen::MatrixXd H = Payoffs(S);
        Eigen::MatrixXd V(H.rows(), H.cols());
        V = Eigen::MatrixXd::Zero(H.rows(), H.cols());
        V(Eigen::all, Eigen::last) = H(Eigen::all, Eigen::last);

        for(int t = H.cols() - 1; t >= 0; t--) {
            vector<int> valid = validPaths(H.col(t));
            Eigen::MatrixXd A = buildA(S(valid, t));
            Eigen::VectorXd X = V(valid, t + 1) * df;
            Eigen::VectorXd b = (A.transpose()) * X;
            Eigen::MatrixXd mat = A.transpose() * A;
            Eigen::VectorXd beta = mat.colPivHouseholderQr().solve(b);
            Eigen::VectorXd C = A * beta;

            vector<int> exercice = validPaths(H(valid, t) - C);
            V(exercice, t) = H(exercice, t);

            for(int i = t + 1; i < V.cols(); i++) 
                V(exercice, Eigen::all).col(i) = Eigen::VectorXd::Zero(exercice.size());

            vector<int> discount = discountPaths(V.col(t));
                V(discount, t) = V(discount, t + 1) * df;
        }

        Eigen::VectorXd price(M - 1);
        for(int j = 0; j < M - 1; j++) {
            
            price(j) = df * mean(V.col(j + 1));
        }

        return price;
}

Eigen::MatrixXd LSMC::getTrajs() const {
    Eigen::MatrixXd S = trajs();
    return S;
}


