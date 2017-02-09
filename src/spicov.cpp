/////////////////////////////////////////////////////////////////
// Spicov.cpp
//
// Version: 1.0
// Author: Didier Chételat (dc623@cornell.edu)
//

// [[Rcpp::depends(RcppArmadillo)]]
#include <sstream> 
#include <RcppArmadillo.h>
#include <string>
#include <stack>
#include <vector>
#include <queue>
#include <algorithm>

using namespace Rcpp;


/****************************************************************
 Declarations
****************************************************************/

/////////////////////////////////////////////////////////////////
// Covariance estimator class (i.e. "CE")
class CE
{
  int n;
  int p;
  int r;
  double* l;
  double* g;
  double** dg;
  double** d2g;
  double s2;
  double* ds2;
  double* d2s2;
  double A;
  double B;
  
  void coef();

public:
  double ure();
  arma::vec eigs();
  
  CE(arma::vec l_t, int n_t, int r_t);
  ~CE();
};

/****************************************************************
 Function Calls
****************************************************************/

/////////////////////////////////////////////////////////////////
// spiCov(S,n,r,verbose)
//
// Computes the robust spiked covariance estimate.
// See section 4 and proposition 1 for the construction.
//
// Arguments
// --------------------------------------
// S       : pxp sample covariance matrix
// n       : sample size (n>=p)
// r       : spiked rank at which to evaluate the estimator (0 to p). 
//           If unspecified or NULL, r is chosen so as to minimize 
//           the unbiased risk estimator.
// verbose : specifies if verbose (TRUE) or not (FALSE). Default is FALSE.
//
// Return value
// --------------------------------------
// hS : the covariance estimator
//

// [[Rcpp::export]]
NumericMatrix spiCov(arma::mat S, int  n, SEXP r=R_NilValue, bool verbose=false)
{
  arma::mat hS, O;
  arma::vec l, l_inc, l_new;
  int p=S.n_rows;
  int R=0;
  double u=0, up=((double)p+1)/n;
  
  if(p>n) // Check that n>=p
  {
    Rcout << "Error: currently, p>n settings are not supported.\n";
    return wrap(arma::zeros<arma::mat>(p,p));
  }
  if(verbose)
    Rcout << "Unbiased risk estimate of the sample covariance matrix: " << up << ".\n";
 
  // -------------------------------------------
  // 1) Compute the spectral decomposition of S.
  // By default, Armadillo's spectral decomposition algorithm outputs
  // eigenvalues in decreasing order. To keep with the notation in the 
  // paper, we need to sort them in decreasing order.
  arma::eig_sym(l_inc,O,S);
  
  l.zeros(p);
  for(int i=0;i<p;++i)
    l[i]=l_inc[p-i-1];
  
  // -------------------------------------------
  // 2) Computes r, if necessary.
  
  if(!Rf_isNull(r))        // If a rank r was specified, use it.
  {
    R=Rcpp::as<int>(r);
    CE e(l,n,R);
    l_new=e.eigs(); // New corrected eigenvalues
    if(verbose)     // Print out info if verbose
    {
      Rcout << "Computing estimate for spiked rank r=" << R << "...\n\n";
      u=e.ure();
      Rcout << "Unbiased risk estimate: " << u << ".\n";
    }
  }
  else                     // Else search for a r that minimizes the ure.
  {
    char out_buf[1+8+3+12+2+1];
    
    if(verbose)
      Rcout << "Searching minimal spiked rank r with improved U.R.E. ...\n\n"
            << "     Rank |       U.R.E.      \n"
            << "------------------------------\n";
            
    // 2.1) Draw up a list of potential ranks
    std::queue<int> R_list;
    int k; double a=0;
    for(k=0;k<p;++k)
    {
      for(int c=k;c<p;++c) 
        a+=l[c];
      a*=(1+sqrt((double)p/n))*(1+sqrt((double)p/n))/(p-k)/l[k];
      if(a>=1)
        R_list.push(k);
    }
    R_list.push(p);
    
    // 2.2) For each potential spiked rank 0<=k<p, compute the ure
    //    until we find one smaller than (p+1)/n.
    while(!R_list.empty())
    {
      /// Obtain next potential rank k
      k=R_list.front();
      
      if(verbose) {
        sprintf(out_buf," %8d | ",k); Rcout << out_buf;
      }
      
      // Compute corresponding ure
      CE e(l,n,k);
      u=e.ure();
      
      if(verbose) {
        sprintf(out_buf," %12f \n",u); Rcout << out_buf;
      }
          
      // Check if minimum
      if(u<=up || k==p)
      {
        R=k;
        l_new=e.eigs();
        break;
      }
      else R_list.pop();
    }
    
    if(verbose)
      Rcout << "\nSelected rank: " << R << " with U.R.E: " << u <<".\n";
  }
  
  // -------------------------------------------
  // 3) Reverse the corrected eigenvalues in increasing order, and output
  // the eigenvalue estimator.
   for(int i=0;i<p;++i)
     l_inc[i]=l_new[p-i-1];
     
   hS=O*arma::diagmat(l_inc)*O.t();
   return wrap(hS);
}

/////////////////////////////////////////////////////////////////
// CE:CE()
//
// Constructor of the covariance estimator class. Computes the 
// corrected eigenvalues at the specified spiked rank.
//
// Arguments
// --------------------------------------
// l_t    : px1 vector of eigenvalues, in decreasing order
// n_t    : sample size (n_t>=p)
// r_t    : spiked rank of the covariance estimator (0<=r<=p)
//
// Careful with r=0!!!!!!!!!!
// l should be imputed in decreasing order

CE::CE(arma::vec l_t, int n_t, int r_t)
{
  double N=n_t;
  double P=l_t.n_elem;
  
  // -------------------------------------------
  // 1) Allocation.
  n=n_t; p=l_t.n_elem; r=r_t;

  l = new double[p];
  s2 = 0;
  ds2 = new double[p];
  d2s2 = new double[p-r];
  
  for(int k=0;k<p;++k)
    l[k]=l_t[k];

  
  if(r) // When r>0.
  {
    g = new double[r];
    dg = new double*[r];
    d2g = new double*[r];
    for(int i=0;i<r;++i)
    {
      dg[i]=new double[p];
      d2g[i]=new double[p];
    }
  }
  else { g=NULL; dg = NULL; d2g = NULL; } // When r=0.
  
  // -------------------------------------------
  // 2) Computation of l, g, s2. The steps follow quite
  // closely the paper, including the definition of A and B.
  if(r<p)
  {    
    // Compute g.
    double a1=0;
    double a2=0;
    for(int k=0;k<r;++k)
    {
      a1=0;a2=0;
      for(int c=r;c<p;++c)
      {
        a1+=l[c];
        a2+=(l[c])/(l[k]-l[c]);
      }
      g[k]=a1/a2;
    }
    
    // Compute A, B and then s2.
    A=0;
    for(int c=0;c<p;++c)
    {
      A+=(N-P-1)/N/P/l[c];
      for(int k=0;k<r;++k)
        A+=(N-P-1)/N/N/P*g[k]/l[k]/l[c];
    }
    for(int k=0;k<r;++k)
    {
      A-=(N-P-1)*(N-P-2)/N/N/P*g[k]/l[k]/l[k];
      for(int c=r;c<p;++c)
      {
        A-=2*(N-P-1)/N/N/P*g[k]/l[c]/(l[k]-l[c]);
        for(int b=0;b<r;++b)
          if(b!=k)
          {
            A+=3/N/N/P*g[b]/(l[k]-l[c])/(l[b]-l[c]);
            A-=3/N/N/P*(g[k]-g[b])/(l[k]-l[b])/(l[k]-l[c]);
          }
        for(int d=r;d<p;++d)
          if(d!=c)
            A-=3/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[d]);
      }
    }
    B=0;
    for(int c=0;c<p;++c)
    {
      B+=(N-P-1)*(N-P-2)/N/N/P/l[c]/l[c];
      for(int d=0;d<p;++d)
        B-=(N-P-1)/N/N/P/l[c]/l[d];
    }
    s2=A/B;
  }
}

/////////////////////////////////////////////////////////////////
// CE:~CE()
//
// Destructor of the covariance estimator class. Nothing exciting.
//

CE::~CE()
{
  delete[] l;
  delete[] ds2;
  delete[] d2s2;
  if(r)
  {
    delete[] g;
    for(int i=0; i<r; ++i)
    {
      delete[] dg[i];
      delete[] d2g[i];
    }
    delete[] dg;
    delete[] d2g;
  }
}

/////////////////////////////////////////////////////////////////
// CE:eigs()
//
// Returns the corrected eigenvalues.
//
// Return value
// --------------------------------------
// hl  : px1 vector with the corrected eigenvalues, in decreasing 
//       order.
//

arma::vec CE::eigs()
{
  arma::vec hl(p);
  if(r<p)
  {
    for(int i=0;i<r;++i)
      hl[i]=g[i];
    for(int i=r;i<p;++i)
      hl[i]=s2;
  }
  else 
  {
    for(int i=0;i<p;++i)
      hl[i]=l[i];
  }
  return hl;
}

/////////////////////////////////////////////////////////////////
// CE:ure()
//
// Computes the unbiased risk estimate for the Haff loss of the
// covariance estimator. See Section 2 in the paper for more details.
//
// Return value
// --------------------------------------
// u  : the unbiased risk estimate for the estimator.
//

double CE::ure(void)
{
  double* psi=new double[p];
  double** dpsi=new double*[p];
  for(int i=0;i<p;++i)
    dpsi[i]=new double[p];
  double* d2psi=new double[p];
  double u=0;
  double N=(double)n;
  double P=(double)p;
  
  // If r=p, save time by returning closed-form expression.
  if(r==p) return (P+1)/N;
  
  // -------------------------------------------
  // 1) Compute dg, ddg, ds2, dds2. The computation is 
  // delegated to another function for cleanliness.
  
  coef();
  
  // -------------------------------------------
  // 2) Construct psi, dpsi, ddpsi.
  
  for(int i=0;i<r;++i)
    psi[i]=g[i];
  for(int i=r;i<p;++i)
    psi[i]=s2;
    
  for(int i=0;i<r;++i)
    for(int j=0;j<p;++j)
      dpsi[i][j]=dg[i][j];
  for(int i=r;i<p;++i)    
    for(int j=0;j<p;++j)
      dpsi[i][j]=ds2[j];
  
  for(int i=0;i<r;++i)
    d2psi[i]=d2g[i][i];
  for(int i=r;i<p;++i)  
    d2psi[i]=d2s2[i-r];
  
  // -------------------------------------------
  // 3) Computation of the ure.
  
  double a1=0;double a2=0;
  u+=1;
  for(int k=0;k<p;++k)
  {
    u+=(N-P-1)*(N-P-2)/N/N/P*psi[k]*psi[k]/l[k]/l[k];
    u+=8/N/N/P*dpsi[k][k]*dpsi[k][k];
    u+=8/N/N/P*psi[k]*d2psi[k];
    u+=8*(N-P-1)/N/N/P*psi[k]/l[k]*dpsi[k][k];
    u-=2*(N-P-1)/N/P*psi[k]/l[k];
    u-=4/N/P*dpsi[k][k];
    
    for(int b=0;b<p;++b)
    {
      u-=(N-P-1)/N/N/P*psi[k]*psi[b]/l[k]/l[b];
      if(b!=k)
      {
        u+=4*(N-P-1)/N/N/P*psi[k]/l[k]*(psi[k]-psi[b])/(l[k]-l[b]);
        u+=8/N/N/P*dpsi[k][k]*(psi[k]-psi[b])/(l[k]-l[b]);
        u+=4/N/N/P*psi[k]*(dpsi[k][k]-dpsi[b][b])/(l[k]-l[b]);
        u+=4/N/N/P*psi[k]*(dpsi[k][k]-dpsi[b][k])/(l[k]-l[b]);
        u-=2/N/P*(psi[k]-psi[b])/(l[k]-l[b]);
        for(int e=0;e<p;++e) 
          if(e!=b && e!=k)
          {
            u+=2/N/N/P*psi[k]/(l[k]-l[b])*
              ((psi[k]-psi[e])/(l[k]-l[e])-(psi[b]-psi[e])/(l[b]-l[e]));
            u+=2/N/N/P*(psi[k]-psi[b])/(l[k]-l[b])*(psi[k]-psi[e])/(l[k]-l[e]);
          }
      }
    }
    
  }
  
  // -------------------------------------------
  // 4) Clean up and output the ure.
  
  delete[] psi;
  for(int i=0;i<p;++i)
    delete[] dpsi[i];
  delete[] dpsi;
  delete[] d2psi;
  
  return u;
}

/////////////////////////////////////////////////////////////////
// CE:coef()
//
// Computes auxiliary coefficients (first and second derivatives)
// that are necessary to compute the unbiased risk estimator. 
// See CE::ure().
//

void CE::coef()
{
  double N=(double)n;
  double P=(double)p;
  double *dA=new double[p];
  double *dB=new double[p];
  double *d2A=new double[p-r];
  double *d2B=new double[p-r];
  
  // -------------------------------------------
  // 1) Compute dg and d2g.
  
  double a1=0;
  for(int c=r;c<p;++c)
    a1+=l[c];
  double a2=0;
  double a3=0;
  
  for(int k=0;k<r;++k)
  {
    // Auxiliary terms
    a2=0;a3=0;
    for(int c=r;c<p;++c)
    {
      a2+=l[c]/(l[k]-l[c])/(l[k]-l[c]);
      a3+=l[c]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c]);
    }
    
    // Main
    dg[k][k]=g[k]*g[k]*a2/a1;
    d2g[k][k]=2*dg[k][k]*dg[k][k]/g[k]-4*g[k]*g[k]*a3/a1;
    for(int b=0;b<r;++b)
      if(b!=k)
      {
        dg[k][b]=0;
        d2g[k][b]=0;
      }
    for(int c=r;c<p;++c)
    {
      dg[k][c]=g[k]/a1*
        (1-g[k]/(l[k]-l[c])-g[k]*l[c]/(l[k]-l[c])/(l[k]-l[c]));
      d2g[k][c]=dg[k][c]/a1*(1-2*g[k]/(l[k]-l[c])-
                                  2*g[k]*l[c]/(l[k]-l[c])/(l[k]-l[c]))-
                g[k]/a1/a1*(1-g[k]/(l[k]-l[c])-
                         g[k]*l[c]/(l[k]-l[c])/(l[k]-l[c]))-
                2*g[k]/a1*(g[k]/(l[k]-l[c])/(l[k]-l[c])+
                         g[k]*l[c]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c]));
    }
  }
  
  
  // -------------------------------------------
  // 2) Compute dA and d2A.
  for(int k=0;k<r;++k)
  {
    dA[k]=0;
    dA[k]+=(N-P-1)/N/P/l[k]/l[k];
    dA[k]+=2*(N-P-1)*(N-P-2)/N/N/P*g[k]/l[k]/l[k]/l[k];
    for(int b=0;b<r;++b)
    {
      dA[k]-=(N-P-1)*(N-P-2)/N/N/P*dg[b][k]/l[b]/l[b];
      for(int c=0;c<p;++c)
      {
        dA[k]+=(N-P-1)/N/N/P*dg[b][k]/l[b]/l[c];
      }
    }
    for(int c=0;c<p;++c)
      dA[k]-=(N-P-1)/N/N/P*g[k]/l[k]/l[k]/l[c];
    for(int c=r;c<p;++c)
    { 
      dA[k]+=2*(N-P-1)/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[c])/l[c];
      for(int b=0;b<r;++b)
      { 
        dA[k]=dA[k]-2*(N-P-1)/N/N/P*dg[b][k]/l[c]/(l[b]-l[c]);
        for(int e=0;e<r;++e)
          if(e!=b)
          {
            dA[k]=dA[k]+3/N/N/P*dg[b][k]/(l[e]-l[c])/(l[b]-l[c]);
            dA[k]=dA[k]+3/N/N/P*dg[b][k]/(l[e]-l[b])/(l[b]-l[c]);
            dA[k]=dA[k]+3/N/N/P*dg[b][k]/(l[e]-l[b])/(l[e]-l[c]);
          }
          if(b!=k)
          { 
            dA[k]=dA[k]-3/N/N/P*g[b]/(l[b]-l[c])/(l[k]-l[c])/(l[k]-l[c]);
            dA[k]=dA[k]-3/N/N/P*g[k]/(l[b]-l[c])/(l[k]-l[c])/(l[k]-l[c]);
            dA[k]=dA[k]+3/N/N/P*(g[k]-g[b])/(l[k]-l[b])/(l[k]-l[b])/(l[k]-l[c]);
            dA[k]=dA[k]+3/N/N/P*(g[k]-g[b])/(l[k]-l[b])/(l[k]-l[b])/(l[b]-l[c]);
            dA[k]=dA[k]+3/N/N/P*(g[k]-g[b])/(l[k]-l[b])/(l[k]-l[c])/(l[k]-l[c]);
          }
        for(int d=r;d<p;++d)
          if(d!=c)
            dA[k]=dA[k]-3/N/N/P*dg[b][k]/(l[b]-l[c])/(l[b]-l[d]);
      }
      for(int d=r;d<p;++d)
          if(d!=c)
            dA[k]=dA[k]+6/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[d]);
    }
  }
  for(int c=r;c<p;++c)
  {
    dA[c]=0;
    dA[c]-=(N-P-1)/N/P/l[c]/l[c];
    for(int k=0;k<r;++k)
    {
      dA[c]-=2*(N-P-1)*(N-P-2)/N/N/P*dg[k][c]/l[k]/l[k];
      dA[c]-=(N-P-1)/N/N/P*g[k]/l[k]/l[c]/l[c];
      dA[c]+=+2*(N-P-1)/N/N/P*g[k]/(l[k]-l[c])/l[c]/l[c];
      dA[c]-=2*(N-P-1)/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[c])/l[c];
      for(int d=0;d<p;++d)
        dA[c]+=(N-P-1)/N/N/P*dg[k][c]/l[k]/l[d];
      for(int b=0;b<r;++b)
        if(b!=k)
        {
          dA[c]+=3/N/N/P*g[b]/(l[k]-l[c])/(l[k]-l[c])/(l[b]-l[c]);
          dA[c]+=3/N/N/P*g[b]/(l[k]-l[c])/(l[b]-l[c])/(l[b]-l[c]);
          dA[c]-=3/N/N/P*(g[k]-g[b])/(l[k]-l[b])/(l[k]-l[c])/(l[k]-l[c]);
          for(int d=r;d<p;++d)
          {
           dA[c]-=3/N/N/P*dg[k][c]/(l[k]-l[b])/(l[k]-l[d]);
           dA[c]+=3/N/N/P*dg[b][c]/(l[k]-l[b])/(l[k]-l[d]);
           dA[c]+=3/N/N/P*dg[k][c]/(l[k]-l[d])/(l[b]-l[d]);
          }
        }
      for(int d=r;d<p;++d)
      {
        dA[c]=dA[c]-2*(N-P-1)/N/N/P*dg[k][c]/(l[k]-l[d])/l[d];
        if(d!=c)
          dA[c]-=6/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[d]);
        for(int e=r;e<p;++e) 
          if(e!=d)
            dA[c]-=3/N/N/P*dg[k][c]/(l[k]-l[e])/(l[k]-l[d]);
      }
    }
  }
  
  for(int c=r;c<p;++c)
  {
    d2A[c-r]=0;
    d2A[c-r]+=2*(N-P-1)/N/P/l[c]/l[c]/l[c];
    for(int k=0;k<r;++k)
    {
      d2A[c-r]-=(N-P-1)*(N-P-2)/N/N/P*d2g[k][c]/l[k]/l[k];
      d2A[c-r]-=2*(N-P-1)/N/N/P*dg[k][c]/l[k]/l[c]/l[c];
      d2A[c-r]+=2*(N-P-1)/N/N/P*g[k]/l[k]/l[c]/l[c]/l[c];
      d2A[c-r]-=4*(N-P-1)/N/N/P*g[k]/(l[k]-l[c])/l[c]/l[c]/l[c];
      d2A[c-r]+=4*(N-P-1)/N/N/P*dg[k][c]/(l[k]-l[c])/l[c]/l[c];
      d2A[c-r]+=4*(N-P-1)/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[c])/l[c]/l[c];
      d2A[c-r]-=4*(N-P-1)/N/N/P*dg[k][c]/(l[k]-l[c])/(l[k]-l[c])/l[c];
      d2A[c-r]+=4*(N-P-1)/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c])/l[c];
      for(int d=0;d<p;++d)
        d2A[c-r]+=(N-P-1)/N/N/P*d2g[k][c]/l[k]/l[d];
      for(int b=0;b<r;++b)
        if(b!=k)
        {
          d2A[c-r]+=6/N/N/P*g[b]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c])/(l[b]-l[c]);
          d2A[c-r]+=6/N/N/P*dg[k][c]/(l[k]-l[c])/(l[k]-l[c])/(l[b]-l[c]);
          d2A[c-r]+=6/N/N/P*g[b]/(l[k]-l[c])/(l[k]-l[c])/(l[b]-l[c])/(l[b]-l[c]);
          d2A[c-r]+=6/N/N/P*dg[k][c]/(l[k]-l[c])/(l[b]-l[c])/(l[b]-l[c]);
          d2A[c-r]+=6/N/N/P*g[b]/(l[k]-l[c])/(l[b]-l[c])/(l[b]-l[c])/(l[b]-l[c]);
          d2A[c-r]-=6/N/N/P*dg[k][c]/(l[k]-l[b])/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c]);
          d2A[c-r]+=6/N/N/P*dg[b][c]/(l[k]-l[b])/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c]);
          d2A[c-r]-=6/N/N/P*(g[k]-g[b])/(l[k]-l[b])/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c]);
          for(int d=r;d<p;++d)
          {
            d2A[c-r]+=3/N/N/P*d2g[k][c]/(l[k]-l[d])/(l[b]-l[d]);
            d2A[c-r]-=3/N/N/P*d2g[k][c]/(l[k]-l[b])/(l[k]-l[d]);
            d2A[c-r]+=3/N/N/P*d2g[b][c]/(l[k]-l[b])/(l[k]-l[d]);
          }
        }
      for(int d=r;d<p;++d)
      {
        d2A[c-r]-=2*(N-P-1)/N/N/P*d2g[k][c]/(l[k]-l[d])/l[d];
        if(d!=c)
        {
          d2A[c-r]-=12/N/N/P*g[k]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[d]);
          d2A[c-r]-=12/N/N/P*dg[k][c]/(l[k]-l[c])/(l[k]-l[c])/(l[k]-l[d]);
        }
        for(int e=r;e<p;++e)
          if(e!=d)
            d2A[c-r]-=3/N/N/P*d2g[k][c]/(l[k]-l[e])/(l[k]-l[d]);
      }
    }
  }
  
  // -------------------------------------------
  // 3) Compute dB and d2B.
  
    for(int c=0;c<p;++c)
    {
      dB[c]=0;
      dB[c]-=2*(N-P-1)*(N-P-2)/N/N/P/l[c]/l[c]/l[c];
      for(int d=0;d<p;++d)
        dB[c]+=2*(N-P-1)/N/N/P/l[c]/l[c]/l[d];
    }
    for(int c=r;c<p;++c)
    {
      d2B[c-r]=0;
      d2B[c-r]+=6*(N-P-1)*(N-P-2)/N/N/P/l[c]/l[c]/l[c]/l[c];
      d2B[c-r]-=2*(N-P-1)/N/N/P/l[c]/l[c]/l[c]/l[c];
      for(int d=0;d<p;++d)
        d2B[c-r]-=4*(N-P-1)/N/N/P/l[c]/l[c]/l[c]/l[d];
    }
    
    
  // -------------------------------------------
  // 4) Compute ds2 and d2s2
  
  for(int c=0;c<p;++c)
    ds2[c]=dA[c]/B-A/B/B*dB[c];
  for(int c=r;c<p;++c)
    d2s2[c-r]=d2A[c-r]/B-2/B/B*dA[c]*dB[c]+2*A/B/B/B*dB[c]*dB[c]-A/B/B*d2B[c-r];
    
  // -------------------------------------------
  // 5) Clean up.
  
  delete[] dA;
  delete[] dB;
  delete[] d2A;
  delete[] d2B;
}

/****************************************************************
 Extra (for simulations)
****************************************************************/

/////////////////////////////////////////////////////////////////
// Stein(S,n)
//
// Computes Stein's isotonized estimator, according to the algorithm described in 
//     Lin, S. and Perlman, M. D. (1985). A Monte Carlo comparison of
//     four estimators of a covariance matrix. Multivariate Analysis, 6:411-429.
//
// Arguments
// --------------------------------------
// S       : pxp sample covariance matrix
// n       : sample size (n>=p)
//
// Return value
// --------------------------------------
// hS : the covariance estimator
//

NumericMatrix Stein(arma::mat S, int  n)
{
  int p=S.n_cols;
  
  arma::mat hS, O;
  arma::vec l_inc(p);
  std::vector<double> l_hat(p), alpha(p);

  // Get the eigenvalues in decreasing order
  arma::eig_sym(l_inc,O,S);
  for(int i=0;i<p;++i)
    l_hat[i]=l_inc[p-i-1];

  // --------------------------------
  // 1) Raw estimates
  for(int i=0;i<p;i++) {
    alpha[i]=n-p-1;
    for(int j=0;j<p;j++) {
      if(j!=i)
        alpha[i]+=2*l_hat[i]/(l_hat[i]-l_hat[j]);
    }
  }

  // --------------------------------
  // 2) Isotonization
  std::stack<int> pooled;
  int i=0;
  
  // Part 1
  // Aggregation of negative alphas
  int neg_alpha=0;
  bool neg_alpha_found=true;
  while(neg_alpha_found) {
    // Get last negative, if there isn't any
    // then break the while loop
    neg_alpha_found=false;
    for(;i<alpha.size();i++) {
      neg_alpha=alpha.size()-1-i;
      if(alpha[neg_alpha]<0) {
        pooled.push(neg_alpha-1);
        
        l_hat[neg_alpha-1]+=l_hat[neg_alpha];
        alpha[neg_alpha-1]+=alpha[neg_alpha];
        l_hat.erase(l_hat.begin()+neg_alpha);
        alpha.erase(alpha.begin()+neg_alpha);
        
        neg_alpha_found=true;
        break;
      }
    }
  }
  
  // Part 2
  // Monotonization  
  int inc_alpha=0;
  bool inc_alpha_found=true;
  while(inc_alpha_found) {
    // Get last negative, if there isn't any
    // then break the while loop
    inc_alpha_found=false;
    for(i=0;i<alpha.size()-1;i++) {
      inc_alpha=alpha.size()-2-i;
      if(l_hat[inc_alpha+1]/alpha[inc_alpha+1]>l_hat[inc_alpha]/alpha[inc_alpha]) {
        pooled.push(inc_alpha);
        
        l_hat[inc_alpha]+=l_hat[inc_alpha+1];
        alpha[inc_alpha]+=alpha[inc_alpha+1];
        l_hat.erase(l_hat.begin()+inc_alpha+1);
        alpha.erase(alpha.begin()+inc_alpha+1);
        
        inc_alpha_found=true;
        break;
      }
    }
  }

  // Part 3
  // Lengthening
  int new_index=0;
  while(!pooled.empty())
  {
    new_index=pooled.top();pooled.pop();
    l_hat.insert(l_hat.begin()+new_index+1,l_hat[new_index]);
    alpha.insert(alpha.begin()+new_index+1,alpha[new_index]);
  }
  
  // --------------------------------

  // Reverse the eigenvalues increasingly
   for(i=0;i<p;++i)
     l_inc[i]=l_hat[p-i-1];
     
   hS=O*arma::diagmat(l_inc)*O.t();
   return wrap(hS);
}

/////////////////////////////////////////////////////////////////
// LedoitWolf(X)
//
// Computes the Ledoit-Wolf estimator, according to the paper
//     Ledoit, O. and Wolf, M. (2004). A well-conditioned estimator for 
//     large-dimensional covariance matrices. Journal of Multivariate 
//     Analysis, 88(2):365-411.
//
// Arguments
// --------------------------------------
// X       : nxp data matrix
//
// Return value
// --------------------------------------
// hS : the covariance estimator
//

NumericMatrix LedoitWolf(arma::mat X)
{
  arma::mat hS, S;
  int n=X.n_rows, p=X.n_cols;
  double N=n, P=p, m, a2, b2, 
          bar_b2, d2, bar_d2, temp=0;
  
  // --------------------------------
  // 1) Constants
  
  S=X.t()*X/N;
  m=trace(S)/P;
  temp=arma::norm(S-m*arma::eye(p,p),"fro");
  d2=temp*temp/P;
  bar_b2=0;
  for(int k=0;k<n;k++)
  {
    temp=arma::norm(X.row(k).t()*X.row(k)-S,"fro");
    bar_b2+=temp*temp/N/N/P;
  }
  b2=std::min(bar_b2,d2);
  a2=d2-b2;
  
  // --------------------------------
  // 2) Computation of estimator
  
  hS=b2/d2*m*arma::eye(p,p)+a2/d2*S;
  return wrap(hS);
}
