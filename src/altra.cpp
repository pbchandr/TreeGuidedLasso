#include <stdio.h>
#include <math.h>
#include <string.h>
#include <Rcpp.h>
using namespace Rcpp;

 /* -------------------------------------------------------------------
 *                       Functions and parameter
 * -------------------------------------------------------------------
 *
 * altra solves the following problem
 *
 * 1/2 \|x-v\|^2 + \sum \lambda_i \|x_{G_i}\|,
 *
 * where x and v are of dimension n,
 *       \lambda_i >=0, and G_i's follow the tree structure
 *
 * It is implemented in Matlab as follows:
 *
 * x=altra(v, n, ind, nodes);
 *
 * ind is a 3 x nodes matrix.
 *       Each column corresponds to a node.
 *
 *       The first element of each column is the starting index,
 *       the second element of each column is the ending index
 *       the third element of each column corrreponds to \lambbda_i.
 *
 * -------------------------------------------------------------------
 *                       Notices:
 * -------------------------------------------------------------------
 *
 * 1. The nodes in the parameter "ind" should be given in the 
 *    either
 *           the postordering of depth-first traversal
 *    or 
 *           the reverse breadth-first traversal.
 *
 * 2. When each elements of x are penalized via the same L1 
 *    (equivalent to the L2 norm) parameter, one can simplify the input
 *    by specifying 
 *           the "first" column of ind as (-1, -1, lambda)
 *
 *    In this case, we treat it as a single "super" node. Thus in the value
 *    nodes, we only count it once.
 *
 * 3. The values in "ind" are in [1,n].
 *
 * 4. The third element of each column should be positive. The program does
 *    not check the validity of the parameter. 
 *
 *    It is still valid to use the zero regularization parameter.
 *    In this case, the program does not change the values of 
 *    correponding indices.
*/

// [[Rcpp::export]]
NumericVector altra(NumericVector v, int n, NumericMatrix ind, int nodes){
    //Rcpp::Rcout << "ind.nrows = " << ind.nrow() << std::endl;
    //Rcpp::Rcout << "ind.ncols = " << ind.ncol() << std::endl;

    NumericVector x(n);
    int i, j;
    double lambda, twoNorm, ratio;
    
    //Rcpp::Rcout << ind(0, 0) << std::endl;
    //Rcpp::Rcout << ind(0, 1) << std::endl;
    //Rcpp::Rcout << ind(0, 2) << std::endl;

    //Rcpp::Rcout << ind;
    
    if (ind(0,0) ==-1){
        //Recheck whether ind[1] equals to zero
        //Rcpp::Rcout << ind(0,0) << std::endl;
        if (ind(1, 0)!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }        
        
        lambda=ind(2,0);
        
        for(j=0;j<n;j++){
            if (v[j]>lambda)
                x[j]=v[j]-lambda;
            else
                if (v[j]<-lambda)
                    x[j]=v[j]+lambda;
                else
                    x[j]=0;
        }
        
        i=1;

    }
    else{
        for(j=0;j<n;j++){
            x[j] = v[j];
        }
        i = 0;
    }

    for(; i < nodes; i ++) {
        // compute L2 norm of this group
        twoNorm = 0;
        for(j = ind(0, i) - 1; j < ind(1, i); j++)
            twoNorm += x[j] * x[j];
        twoNorm = sqrt(twoNorm);

        lambda = ind(2, i);

        if(twoNorm > lambda) {
            ratio = (twoNorm - lambda)/twoNorm;

            /*
             * shrinkage this group by ratio
             */
            for(j= ind(0,i)-1;j< ind(1, i);j++)
                x[j] *= ratio;
        }
        else {
            for(j= ind(0,i)-1;j< ind(1, i);j++)
                x[j] = 0;
        }
    }
    return(x);

}


/*
  * -------------------------------------------------------------------
  *                       Function and parameter
  * -------------------------------------------------------------------
  *
  * findLambdaMax compute
  * 
  * the lambda_{max} that achieves a zero solution for
  *
  *     min  1/2 \|x-v\|^2 +  \lambda_{\max} * \sum  w_i \|x_{G_i}\|,
  *
  * where x is of dimension n,
  *       w_i >=0, and G_i's follow the tree structure
  *
  * The file is implemented in the following in Matlab:
  *
  * lambdaMax=findLambdaMax(v, n, ind,nodes);
*/

// [[Rcpp::export]]         
double findLambdaMax(NumericVector v, int n, NumericMatrix  ind, int nodes){
  
  double lambdaMax;
  int i;
  double lambda=0,squaredWeight=0, lambda1,lambda2;
  NumericVector x(n);
  NumericMatrix ind2(ind.nrow(), ind.ncol());
  
  //Rcpp::Rcout << "ind.nrows = " << ind.nrow() << std::endl;
  //Rcpp::Rcout << "ind.ncols = " << ind.ncol() << std::endl;
  
  //Rcpp::Rcout << "ind2.nrows = " << ind2.nrow() << std::endl;
  //Rcpp::Rcout << "ind2.ncols = " << ind2.ncol() << std::endl;
  
  int num=0;
  
  for(i = 0; i < n; i++){
    lambda += v[i]*v[i];
  }
  
  if(ind(0,0) == -1)
    squaredWeight += n*ind(2,0)*ind(2,0);
  else
    squaredWeight += ind(2,0)*ind(2,0);
  
  for(i = 1; i < nodes; i++) {
    squaredWeight += ind(3,i)*ind(3,i); 
  }
  
  //set lambda to an initial guess
  lambda=sqrt(lambda/squaredWeight);
  
  //Rcpp::Rcout << "Squared Weight = " << squaredWeight << std::endl;
  //Rcpp::Rcout << "lambda = " << lambda << std::endl;
  
  //Copy ind to ind2 and scale the weight
  for(i = 0; i < nodes; i++){
    ind2(0,i) = ind(0,i);
    ind2(1,i) = ind(1,i);
    ind2(2,i) = ind(2,i)*lambda;
  }
  
  //Test Weather Solution is zero or not
  x = altra(v, n, ind2, nodes);
  //Rcpp::Rcout << x << std::endl;
  
  for(i = 0; i< n; i++){
    if(x[i] != 0)
      break;
  }
  if(i >= n){
    //x is a zero vector
    lambda2 = lambda;
    lambda1 = lambda;
    
    num = 0;
    
    while(1) {
      num++;
      lambda2 = lambda;
      lambda1 = lambda1/2;
      
      //update ind2
      for(i = 0; i < nodes; i++){
        ind2(2,i) = ind(2,i)*lambda1;
      }
      
      //compute and test whether x is zero
      x = altra(v, n, ind2, nodes);
      for(i = 0; i< n; i++){
        if(x[i] != 0)
          break;
      }
      
      if (i<n){
        break; //x is not zero vector. we have found lambda1
      }
    }
  }
  
  else {
    //x is a non-zero vector
    lambda2 = lambda;
    lambda1 = lambda;
    
    num = 0;
    
    while(1) {
      num++;
      lambda1 = lambda2;
      lambda2 = lambda2*2;
      
      //update ind2
      for(i = 0; i < nodes; i++){
        ind2(2,i) = ind(2,i)*lambda2;
      }
      
      //compute and test whether x is zero
      x = altra(v, n, ind2, nodes);
      for(i = 0; i< n; i++){
        if(x[i] != 0)
          break;
      }
      
      if (i>=n){
        break; //x is zero vector. we have found lambda2
      }
    }
  }
  
//  Rcpp::Rcout<<"lambda1 = " << lambda1 << std::endl;
//  Rcpp::Rcout<<"lambda2 = " << lambda2 << std::endl;
//  Rcpp::Rcout<<"num = " << num << std::endl;

  while(fabs(lambda2-lambda1) > lambda2 * 1e-10) {
    num++;
    lambda = (lambda1 + lambda2)/2;

    //update ind2
    for (i = 0; i < nodes; i++)
    {
      ind2(2,i) = ind(2,i)*lambda;
    }

    //compute and test whether x is zero
    x = altra(v, n, ind2, nodes);
    for(i = 0; i< n; i++){
      if(x[i] != 0)
        break;
    }

    if (i>=n) {
      lambda2=lambda;
    }
    else {
      lambda1=lambda;
    }
  }

  //Rcpp::Rcout<<"lambda1 = " << lambda1 << std::endl;
  //Rcpp::Rcout<<"lambda2 = " << lambda2 << std::endl;

  lambdaMax = lambda2;

  return(lambdaMax);
}



/*
 * -------------------------------------------------------------------
 *                       Function and parameter
 * -------------------------------------------------------------------
 *
 * treeNorm compute
 *
 *        \sum \lambda_i \|x_{G_i}\|,
 *
 * where x is of dimension n,
 *       \lambda_i >=0, and G_i's follow the tree structure
*/

// [[Rcpp::export]]
double treeNorm(NumericVector x, int n, NumericMatrix ind, int nodes) {
    int i, j;
    double tree_norm = 0, twoNorm, lambda;
    
    // Check if 1st node is special
    if (ind(0,0) ==-1){
        //Recheck whether ind[1] equals to zero
        //Rcpp::Rcout << ind(0,0) << std::endl;
        if (ind(1, 0)!=-1){
            printf("\n Error! \n Check ind");
            exit(1);
        }        
        
        lambda=ind(2,0);
        
        for(j=0;j<n;j++) {
            tree_norm += fabs(x[j]);
        }

        tree_norm = tree_norm * lambda;
        i = 1;
    }
    else {
        i = 0;
    }

    // Sequentially process each node
    for (; i < nodes; i++) {
        // Compute L2 norm of the group
        twoNorm = 0;
        for(j = ind(0, i) - 1; j < ind(1, i); j++)
            twoNorm += x[j] * x[j];
        twoNorm = sqrt(twoNorm);

        lambda = ind(2, i); 
        tree_norm = tree_norm + lambda*twoNorm;    
    }

    return(tree_norm);
}
