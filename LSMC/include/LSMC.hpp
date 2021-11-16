#ifndef LSMC_H
#define LSMC_H

#include <iostream>
#include <random>
#include <ctime>
#include <Eigen/Dense>
#include <string>
#include <vector>

class LSMC {
    private:
        int I, M, N;
        double r, sigma, S_0, T, K;

        // % Generate trajectories
        Eigen::MatrixXd trajs() const {
            std::mt19937 G(time(NULL));
            std::normal_distribution <double> N(0,1);

            double dt = T /(double)(M - 1);
            Eigen::MatrixXd out(I, M);

            for(int i=0; i < I; i++) {
                out(i, 0) = S_0;

                for(int j=1; j < M; j++) 
                    out(i, j) = out(i, j - 1) * exp(-(r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * N(G));


            }

            return out;
                        }
        // % Validpaths (payoff > 0)
        std::vector<int> validPaths(const Eigen::VectorXd& H) const {
            std::vector<int> idx;

            for(int i = 0; i < H.size(); i++) {
                if(H(i) > 0) idx.push_back(i);
                else continue;  
            }
            return idx;
        }

        // % Paths where we didn't exercise
        std::vector<int> discountPaths (const Eigen::VectorXd& H) const {
            std::vector<int> idx;

            for(int i = 0; i < H.size(); i++) {
                if(H(i) == 0) idx.push_back(i);
                else continue;  
            }
            return idx;
        }

        // % Payoffs matrix
        Eigen::MatrixXd Payoffs(const Eigen::MatrixXd& trajs) const {
            Eigen::MatrixXd H(trajs.rows(), trajs.cols());
            for(int j = 0; j < H.cols(); j++) {
                for(int i = 0; i < H.rows(); i++ ) 
                    
                    H(i, j) = std::max(K - trajs(i, j), 0.);
                    
                
        }
            return H;
        }

        // % Build A_t = ((S_t^i)^j)
        Eigen::MatrixXd buildA(const Eigen::VectorXd& price) const {
            Eigen::MatrixXd A(price.size(), N + 1);
            for(int i=0; i < A.rows(); i++)
                for(int j=0; j < A.cols(); j++) {
                    A(i, j) = pow(price(i), j);
                }
            
            return A;
        }

        double mean(const Eigen::VectorXd& V) const{
            int N = V.size();

            double sum(0.);

            for(int i = 0; i < N; i++) {
                sum += V(i);
            }

            double m = sum / (double) (N);
            return m;
        }
    
    public:
        LSMC() {}
        LSMC(const int& i, const int& m, const int& n, 
             const double& rate, const double& sig, const double& S0,
             const double& t, const double& k): 
             I(i), M(m), N(n), r(rate), sigma(sig), S_0(S0), T(t), K(k) {}
        
        ~LSMC() {}
        Eigen::VectorXd pricePut() const;
        Eigen::MatrixXd getTrajs() const;
};
#endif