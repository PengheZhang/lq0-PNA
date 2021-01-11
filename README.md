# lq,0-PNA (lq,0-proximal Newton algorithm)
Matalab code for the paper: "Solving Multinomial Logistic Regression via lq,0-Proximal Newton Algorithm"  
Author: Penghe Zhang, Rui Wang, Naihua Xiu  

====================================================================  
CONTENTS:  
* `grad_Pw.m`: Matlab function for computing gradient and probability matrix.  
* `hessian_w.m`: Matlab function for product of the Hessian and a vector.  
* `main_simulate.m`: Matlab script file. Perform the numerical test on simulated data.  
* `Matvecn.m`: Matlab function for computing the gradient in the linear CG.  
* `multi_logistic_fun.m`: Matlab function for computing the average multinomial logistic loss.  
* `partial_hessian.m`: Matlab function for calculating submatrix of the Hesssian.  
* `pcgn.m`: Matlab function of conjugate gradient method. 
* `PNA.m`: Matlab function for lq,0-proximal gradient algorithm.
* `prop2lam.m`: Matlab function file. Given the proportion of the significant features, estimate the corresponding regularization parameter lambda.
* `Random_sam.m`: Matlab function file for generating simulated data.
