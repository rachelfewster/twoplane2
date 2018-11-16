
#define _GLIBCXX_DEBUG

#include <Rcpp.h>

#include <algorithm>
#include <vector>
//#include <tuple>
#include <cmath>

#ifdef NDEBUG
// Make sure assertions are switched on. 
#undef NDEBUG
#endif

#include <assert.h>
#include <stdio.h>
#include <iostream>
#include "to_string.h"
#include <fenv.h>

using namespace std;

typedef long double num;
typedef vector<vector<num> > matrix;

class dfs_compute_likelihood
{
    int numpairs;
    
    vector<pair<int,int> > closepairs;  // gives the indices of the two observations that are close
    
    public:
    vector<int> varstate; // 0,1, -1
    private:
    int potential_matches_s1, potential_matches_s2;
    
    
    vector<vector<int> >  constraints; //  indexed by [len(closepairs)] 
    // contains vector of indices into varstates, to set to 0.
    
    
    vector<int> btstack;  // indices into varstates to set back to -1.
    
    unsigned long long solcount;
    
    
    // Numerical part.
    num const Pi;
    
    //Parameters.
    
    num theta_1;    //   Log density in animals / plane timestep.
    num theta_3;    //   Log of animal speed in km/timestep.
    
    num p1, p2;
    num k;         //   Separation of observers in timesteps. 
    num tL;         //   Length of flight in plane timesteps
    num planespd;  //   aircraft speed in km/timestep.
    
    
    vector<num> s1;    // Location (in plane seconds) along the line that plane 1 observed animals
    vector<num> s2;    // Same for plane 2. 
    num dmax;          // Maximum separation (in plane seconds) where two observations may be paired.
    
    // Vars calculated at root node
    num p10k, p01k, p11k;
    int nplus10, nplus01;
    
public:
    num Lk_accumulator;  // at the end, this equals Lk.
    
    // Ctor
    dfs_compute_likelihood(vector<num> _s1, vector<num> _s2, num _dmax,
        num _theta_1, num _gamma_12, num _theta_3, num _gamma_21, 
        num _p1, num _p2, num _k, num _tL, num _planespd, num _p10k, num _p01k, num _p11k) : 
    potential_matches_s1(0), potential_matches_s2(0),
     solcount(0),
    
    Pi(4*atan(1)),
    
    theta_1(_theta_1), theta_3(_theta_3),
    p1(_p1), p2(_p2), k(_k), tL(_tL), planespd(_planespd),
    s1(_s1), s2(_s2), dmax(_dmax), p10k(_p10k), p01k(_p01k), p11k(_p11k),
    Lk_accumulator(-std::numeric_limits<num>::infinity())    //  Log likelihood starts at 0. -inf in logs.
    {
        numpairs=0;
        for(unsigned int i=0; i<s1.size(); i++) {
            bool potential_match=false;
            for(unsigned int j=0; j<s2.size(); j++) {
                if( abs(s1[i]-s2[j])<=dmax ) {
                    closepairs.push_back(make_pair(i,j));
                    numpairs++;
                    potential_match=true;
                }
            }
            if(potential_match) potential_matches_s1++;
        }
        
        for(unsigned int i=0; i<s2.size(); i++) {
            bool potential_match=false;
            for(unsigned int j=0; j<s1.size(); j++) {
                if( abs(s1[j]-s2[i])<=dmax ) {
                    potential_match=true;
                }
            }
            if(potential_match) potential_matches_s2++;
        }
        
        varstate.resize(numpairs, -1);
        
        // set up the constraints vector.
        // all constraints are binary: not(x) or not(y)
        
        constraints.resize(numpairs);
        
        for(unsigned int i=0; i<closepairs.size(); i++) {
            for(unsigned int j=0; j<closepairs.size(); j++) {
                if(i!=j && (closepairs[i].first == closepairs[j].first 
                    || closepairs[i].second == closepairs[j].second)) {
                    constraints[i].push_back(j);
                    constraints[j].push_back(i);
                }
            }
        }
        
    }
  
    num dfs() {
        rootnode_calculations();
        
        assert(p1>=0.0 && p1<=1.0);
        assert(p2>=0.0 && p2<=1.0);
        
        dfs_inner(0);
        
        return Lk_accumulator;  // This is the log likelihood
    }
    
    void dfs_inner(int depth) {
        
        if(depth==numpairs) {
            leafnode_calculations();
            
            solcount++;
            return;
        }
        
        if(varstate[depth]==-1) {
            // branch for 0 and 1.
            int btmarker=btstack.size();
            assign_propagate(depth, 0);
            dfs_inner(depth+1);
            backtrack_to(btmarker);
            
            assign_propagate(depth, 1);
            dfs_inner(depth+1);
            backtrack_to(btmarker);
        }
        else {
            dfs_inner(depth+1);
        }
        return;
    }
    
    void assign_propagate(int var, int val) {
        // do the assignment and propagate consequences, bt stack all changes. 
        varstate[var]=val;
        btstack.push_back(var);
        
        if(val==1) {
            for(unsigned int i=0; i<constraints[var].size(); i++) {
                int var2=constraints[var][i];
                assert(varstate[var2] != 1);    // This would mean something else was not propagated properly.
                if(varstate[var2]==-1) {   // if unassigned:
                    varstate[var2]=0;       // set to 0.
                    btstack.push_back(var2);
                }
            }
        }
    }
    
    void backtrack_to(int btmarker) {
        while(btstack.size()> (unsigned int)btmarker) {
            varstate[btstack.back()]=-1;
            btstack.pop_back();
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
    //
    //    The numerical part.
    
    void rootnode_calculations() {
    }
    
    ////////////////////////////////////////////////////////////////////////////
    //   To be called for each matching.
    
    void leafnode_calculations() {
        calculate_Ljk();
    }
    
    // eqn (1) -- PDF of seeing an animal given distance moved d in plane timesteps.
    // Uses theta_3 which is log(sigma), where sigma is in km/timestep.
    num calculate_f(num d) {
      //   t=k+d always.
      num t=k+d;
      num f_firstterm = 1.0/(sqrt(2.0*Pi*t*exp(2.0*(theta_3-log(planespd)))));
      
      num f_secondterm= 2.0*t*exp(2.0*(theta_3-log(planespd)));
      
      return f_firstterm*exp((-(d*d))/( f_secondterm ));
    }
    
    num factorial(int i) {
      num temp=1;
      for(int j=i; j>1; j--) {
        temp=temp*j;
      }
      return temp;
    }
    
    void calculate_Ljk() {
        // Implements formula 19.
        //std::cout << "In calculate_Ljk" <<std::endl;
        // varstate and numpairs...
        
        int n11=0;
        
        for(unsigned int i=0; i<varstate.size();i++) {
            if(varstate[i]==1) n11++;
        }
        int n10=s1.size()-n11;
        int n01=s2.size()-n11;
        
        num prod1=-(p01k+p10k+p11k)*exp(theta_1)*tL;
        
        num prod2=theta_1*(n01+n10+n11);
        
        num prod3=log(p01k)*n01;
        prod3=prod3+log(p10k)*n10;
        prod3=prod3+log(p11k)*n11;
        
        num prod4=0.0;
        for(unsigned int i=0; i<varstate.size();i++) {
            if(varstate[i]==1) {
                pair<int,int> a=closepairs[i];
                num di=s2[a.second]-s1[a.first];  // distance (in plane seconds)
                
                prod4=prod4+log(calculate_f(di));
            }
        }
        
        Lk_accumulator=logadd(Lk_accumulator, prod1+prod2+prod3+prod4);
    }
    
    // matrix multiply
    
    matrix multiply(matrix a, matrix b) {
        
        int cols=b[0].size();
        int rows=a.size();
        
        int acols=a[0].size();
        assert(acols==(int)b.size());
        
        // make the return matrix.
        matrix ret=mkmatrix(rows, cols, 0.0);
        
        for(int i=0; i<rows; i++) {
            for(int j=0; j<cols; j++) {
                num acc=0.0;
                for(int k=0; k<acols; k++) {
                    acc+=(a[i][k]*b[k][j]);
                }
                ret[i][j]=acc;
            }
        }
        
        return ret;
    }
    
    // add scalar multiply and addition.
    
    matrix mkmatrix(int rows, int cols, num value) {
        matrix a;
        a.resize(rows);
        for(int i=0; i<rows; i++) {
            a[i].resize(cols, value);
        }
        return a;
    }
    
    num log_check0(num a) {
        if(a==0.0) {
            return -1e50;    // Huge negative penalty value. 
        }
        else {
            return log(a);
        }
    }
    
    num logadd(num x, num y) {
        if(x == -std::numeric_limits<num>::infinity()) {
            return y;
        }
        if(y == -std::numeric_limits<num>::infinity()) {
            return x;
        }
        return max(x, y) + log1p(exp( -fabs(x - y) ));
    }
    
    void printmatrix(matrix a) {
        std::cout << "[";
        for(unsigned int i=0; i<a.size(); i++) {
            std::cout << "[";
            for(unsigned int j=0; j<a[i].size(); j++) {
                std::cout << a[i][j] << ",";
            }
            std::cout <<"]"<<std::endl;
        }
        std::cout << "]";
    }
};

RcppExport SEXP compute_likelihood(SEXP _s1, SEXP _s2, SEXP _dmax, SEXP _theta_1, SEXP _theta_2, SEXP _theta_3, SEXP _theta_4, 
SEXP _p1, SEXP _p2, SEXP _k, SEXP _tL, SEXP _planespd, SEXP _p10k, SEXP _p01k, SEXP _p11k) {
    
    vector<num> s1=Rcpp::as<vector<num> >(_s1);
    vector<num> s2=Rcpp::as<vector<num> >(_s2);
    
    num dmax=Rcpp::as<num>(_dmax);
    num theta_1=Rcpp::as<num>(_theta_1);
    num theta_2=Rcpp::as<num>(_theta_2);
    num theta_3=Rcpp::as<num>(_theta_3);
    num theta_4=Rcpp::as<num>(_theta_4);
    num p1=Rcpp::as<num>(_p1);
    num p2=Rcpp::as<num>(_p2);
    num k=Rcpp::as<num>(_k);
    num tL=Rcpp::as<num>(_tL);
    num planespd=Rcpp::as<num>(_planespd);
    
    num p10k=Rcpp::as<num>(_p10k);
    num p01k=Rcpp::as<num>(_p01k);
    num p11k=Rcpp::as<num>(_p11k);
    
    dfs_compute_likelihood d(s1, s2, dmax, theta_1, theta_2, theta_3, theta_4, p1, p2, k, tL, planespd, p10k, p01k, p11k);
    
    // Clear any floating point exceptions that may have come from R.
    //feclearexcept(FE_ALL_EXCEPT);
    
    //feenableexcept(FE_DIVBYZERO| FE_INVALID ); //|FE_OVERFLOW);
    
    num loglik=d.dfs();
    
    //fedisableexcept(FE_DIVBYZERO| FE_INVALID ); //|FE_OVERFLOW);
    
    return Rcpp::wrap(loglik);
}

RcppExport SEXP calculate_Ljk(SEXP _s1, SEXP _s2, SEXP _dmax_time, SEXP _theta_1, SEXP _theta_2, SEXP _theta_3, SEXP _theta_4, 
SEXP _p1, SEXP _p2, SEXP _k, SEXP _tL, SEXP _planespd, SEXP _p10k, SEXP _p01k, SEXP _p11k, SEXP _varstate_fixed) {
    
    vector<num> s1=Rcpp::as<vector<num> >(_s1);
    vector<num> s2=Rcpp::as<vector<num> >(_s2);
    
    num dmax=Rcpp::as<num>(_dmax_time);
    num theta_1=Rcpp::as<num>(_theta_1);
    num theta_2=Rcpp::as<num>(_theta_2);
    num theta_3=Rcpp::as<num>(_theta_3);
    num theta_4=Rcpp::as<num>(_theta_4);
    num p1=Rcpp::as<num>(_p1);
    num p2=Rcpp::as<num>(_p2);
    num k=Rcpp::as<num>(_k);
    num tL=Rcpp::as<num>(_tL);
    num planespd=Rcpp::as<num>(_planespd);
    
    num p10k=Rcpp::as<num>(_p10k);
    num p01k=Rcpp::as<num>(_p01k);
    num p11k=Rcpp::as<num>(_p11k);
    
    vector<int> varstate_fixed=Rcpp::as<vector<int> >(_varstate_fixed);
    
    dfs_compute_likelihood d(s1, s2, dmax, theta_1, theta_2, theta_3, theta_4, p1, p2, k, tL, planespd, p10k, p01k, p11k);
    
    d.Lk_accumulator=-std::numeric_limits<num>::infinity();
    d.varstate=varstate_fixed;
    d.calculate_Ljk();
    
    return Rcpp::wrap(d.Lk_accumulator);
}
