#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

/**************************************************************************************************
*
* 						MLE Initial Estimate 
*
****************************************************************************************************/
// [[Rcpp::export]]
double llk_h2_rho(double rho,double h2_1,double h2_2,double sigma_1,double sigma_2,double n_1,double n_2,double num_SNP,arma::vec z_1,arma::vec z_2, arma::mat R_1,arma::mat R_2){
	arma::mat V_mat(2,2);
	V_mat(0,0) = h2_1;
    V_mat(1,1) = h2_2;
    V_mat(0,1) = rho*sqrt(h2_1*h2_2);
    V_mat(1,0) = rho*sqrt(h2_1*h2_2);
	V_mat = V_mat/num_SNP;
	
	arma::mat XtX(2*num_SNP,2*num_SNP,fill::zeros);
	XtX.submat(0,0,num_SNP-1,num_SNP-1) = R_1*n_1/sigma_1;
	XtX.submat(num_SNP,num_SNP,2*num_SNP-1,2*num_SNP-1) = R_2*n_2/sigma_2;
	
	
	arma::vec Xty(2*num_SNP);
	Xty.subvec(0,num_SNP-1) = z_1/sigma_1*sqrt(n_1);
	Xty.subvec(num_SNP,2*num_SNP-1) = z_2/sigma_2*sqrt(n_2);
	
    arma::mat V_inv = inv_sympd(V_mat);
	arma::mat sigma_beta_inv = kron(V_inv,arma::eye(num_SNP,num_SNP))+XtX;
    arma::mat sigma_beta = inv_sympd(sigma_beta_inv);
	arma::vec mu_beta = sigma_beta*Xty;
	double llk = -n_1/2*log(sigma_1)-n_2/2*log(sigma_2)-num_SNP/2*log_det_sympd(V_mat)-0.5*log_det_sympd(sigma_beta_inv)-0.5*n_1/sigma_1-0.5*n_2/sigma_2+0.5*as_scalar(mu_beta.t()*sigma_beta_inv*mu_beta);
	return(-llk);
};

// [[Rcpp::export]]
double llk_h2_est(double h2,double sigma,arma::vec z, arma::mat R,double n, double num_SNP){
    
  double sigma_update,h2_update;
  sigma_update = sigma;
 
  arma::mat XtX = R*n;
 
  // Define Xty, UtXty, and mu_beta, lambda
  arma::vec Xty(num_SNP);
 
  
  Xty = z*sqrt(n);
  arma::mat posterior_var(num_SNP,num_SNP);
   arma::vec h2_vec(num_SNP);
   h2_update =  h2/num_SNP;
   h2_vec.fill(1.0/h2_update);
	posterior_var = XtX/sigma_update + diagmat(h2_vec);
	 // cout<<"Post var is "<<posterior_var(0,0)<<endl;
   //   llk_2(i) = -n/2*log(sigma_update)-num_SNP/2*log(h2_update)- log_det_sympd(posterior_var)/2-n/sigma_update/2+ as_scalar(Xty.t()* inv_sympd(posterior_var)*Xty/2/sigma_update/sigma_update);
	 double llk = 0.5*(-n*log(sigma_update)-num_SNP*log(h2_update)- log_det_sympd(posterior_var)-n/sigma_update+ as_scalar(Xty.t()* inv_sympd(posterior_var)*Xty/sigma_update/sigma_update));
	return(-llk);
}

// [[Rcpp::export]]
double llk_sigma_e_est(double sigma,arma::vec z, arma::mat R,double n, double num_SNP){
    
  double sigma_update,h2_update;
  sigma_update = sigma;
 
  arma::mat XtX = R*n;
 
  // Define Xty, UtXty, and mu_beta, lambda
  arma::vec Xty(num_SNP);
  Xty = z*sqrt(n);
  arma::mat posterior_var(num_SNP,num_SNP);
  posterior_var = XtX/sigma_update ;
	 // cout<<"Post var is "<<posterior_var(0,0)<<endl;
   //   llk_2(i) = -n/2*log(sigma_update)-num_SNP/2*log(h2_update)- log_det_sympd(posterior_var)/2-n/sigma_update/2+ as_scalar(Xty.t()* inv_sympd(posterior_var)*Xty/2/sigma_update/sigma_update);
	 double llk = 0.5*(-n*log(sigma_update)- log_det_sympd(posterior_var)-n/sigma_update+ as_scalar(Xty.t()* inv_sympd(posterior_var)*Xty/sigma_update/sigma_update));
	return(-llk);
}

// [[Rcpp::export]]
double llk_h2_trans_est(double h2,double sigma,arma::vec Utz, arma::vec eig_val,double n, double num_SNP){
    
  double sigma_update,h2_update;
  arma::vec mu(Utz.size(),fill::zeros);
  arma::vec posterior_var = h2*n/num_SNP*pow(eig_val,2)+sigma*eig_val;
  double llk = accu(log_normpdf(Utz,  mu, posterior_var));
  return(-llk);
}

// [[Rcpp::export]]
double llk_sigma_e_trans_est(double sigma,arma::vec Utz, arma::vec eig_val){  
  arma::vec mu(Utz.size(),fill::zeros);
  arma::vec posterior_var = sigma*eig_val;
  double llk = accu(log_normpdf(Utz,  mu, posterior_var));
  return(-llk);
}
/**************************************************************************************************
*
* 						Null Estimate for testing of h2
*
****************************************************************************************************/
// [[Rcpp::export]]
SEXP est_h2_null_score(double h2_1,double h2_2,double sigma_1,double sigma_2,arma::vec z_1,arma::vec z_2, arma::mat R_1,arma::mat R_2,double n_1, double n_2,double num_SNP,int max_iter = 100, bool fix_intercept = false){
  
    double h2_1_update,h2_2_update,sigma_1_update,sigma_2_update,alpha_1,alpha_2;
  
    h2_1_update = h2_1/num_SNP;
    h2_2_update = h2_2/num_SNP;
    sigma_1_update = sigma_1;
    sigma_2_update = sigma_2;
    // XtX and its eigenvalue and eigen decomposition
    
    arma::mat XtX_1 = R_1*n_1;
    arma::mat XtX_2 = R_2*n_2;
    
    arma::vec eigval_1,eigval_2;
    arma::mat eigvec_1,eigvec_2;
    
    eig_sym(eigval_1,eigvec_1,XtX_1);
    eig_sym(eigval_2,eigvec_2,XtX_2);
  
    arma::mat eigvec_1_trans = trans(eigvec_1);
    arma::mat eigvec_2_trans = trans(eigvec_2);
    // Define Xty, UtXty, and mu_beta, lambda
    arma::vec Xty_1(num_SNP),UtXty_1(num_SNP),UtXty_sigma_1(num_SNP),mu_beta_1(num_SNP),lambda_1(num_SNP);
    arma::vec Xty_2(num_SNP),UtXty_2(num_SNP),UtXty_sigma_2(num_SNP),mu_beta_2(num_SNP),lambda_2(num_SNP);
    
    
    Xty_1 = z_1*sqrt(n_1);
    UtXty_1 = eigvec_1.t()*Xty_1;
    
    Xty_2 = z_2*sqrt(n_2);
    UtXty_2 = eigvec_2.t()*Xty_2;
    
    arma::vec llk_1(max_iter+1);
    arma::vec llk_2(max_iter+1);
    
    llk_1(0) = -datum::inf;
    llk_2(0) = -datum::inf;
    
     for(int i=0;i<max_iter;i++){
        lambda_1 = 1.0/(eigval_1/sigma_1_update + 1.0/h2_1_update);
        UtXty_sigma_1 = UtXty_1/sigma_1_update;
        
        llk_1(i+1) = -n_1/2*log(sigma_1_update)-num_SNP/2*log(h2_1_update)+accu(log(lambda_1))/2-n_1/2/sigma_1_update+accu(pow(UtXty_sigma_1,2)%lambda_1)/2;
        if(abs(llk_1(i+1)-llk_1(i))<pow(10,-5)){
            llk_1.resize(i+2);
          break;
      }
        
        h2_1_update = (accu(pow(UtXty_sigma_1%lambda_1,2)) + accu(lambda_1))/num_SNP;
      //  alpha_1 = accu(UtXty_sigma_1%lambda_1%UtXty_1)/(accu(pow(UtXty_sigma_1%lambda_1,2)%eigval_1)+accu(lambda_1%eigval_1));
        alpha_1 = accu(pow(UtXty_sigma_1,2)%lambda_1)/(accu(pow(UtXty_sigma_1%lambda_1,2)%eigval_1)/sigma_1_update+accu(lambda_1%eigval_1)/sigma_1_update);//may lack a sigma_update here
      //  sigma_1_update = 1 - alpha_1*accu(UtXty_sigma_1%lambda_1%UtXty_1)/n_1;
      if(fix_intercept==false){
          sigma_1_update = 1 + (pow(alpha_1,2)*(accu(pow(UtXty_sigma_1%lambda_1,2)%eigval_1)+accu(lambda_1%eigval_1))-2.0*alpha_1*sigma_1_update*accu(pow(UtXty_sigma_1,2)%lambda_1))/n_1;
      }
        h2_1_update = pow(alpha_1,2)*h2_1_update;
        
     }
     
     
    for(int i=0;i<max_iter;i++){ 
        lambda_2 = 1.0/(eigval_2/sigma_2_update + 1.0/h2_2_update);
        UtXty_sigma_2 = UtXty_2/sigma_2_update;
        llk_2(i+1) = -n_2/2*log(sigma_2_update)-num_SNP/2*log(h2_2_update)+accu(log(lambda_2))/2-n_2/2/sigma_2_update+accu(pow(UtXty_sigma_2,2)%lambda_2)/2;
        if(abs(llk_2(i+1)-llk_2(i))<pow(10,-5)){
          llk_2.resize(i+2);
          break;
        }
        h2_2_update = (accu(pow(UtXty_sigma_2%lambda_2,2)) + accu(lambda_2))/num_SNP;
       // alpha_2 = accu(UtXty_sigma_2%lambda_2%UtXty_2)/(accu(pow(UtXty_sigma_2%lambda_2,2)%eigval_2)+accu(lambda_2%eigval_2));
      //  sigma_2_update = 1 - alpha_2*accu(UtXty_sigma_2%lambda_2%UtXty_2)/n_2;
        
        alpha_2 = accu(pow(UtXty_sigma_2,2)%lambda_2)/(accu(pow(UtXty_sigma_2%lambda_2,2)%eigval_2)/sigma_2_update+accu(lambda_2%eigval_2)/sigma_2_update);
        if(fix_intercept==false){
        sigma_2_update = 1 + (pow(alpha_2,2)*(accu(pow(UtXty_sigma_2%lambda_2,2)%eigval_2)+accu(lambda_2%eigval_2))-2.0*alpha_2*sigma_2_update*accu(pow(UtXty_sigma_2,2)%lambda_2))/n_2;
        }
        h2_2_update = pow(alpha_2,2)*h2_2_update;    
    }
    
    
      //sanity check if h2 is greater than 0.1 warning and then set to the initial
      bool fix_h2;
      if(h2_1_update>0.1||h2_2_update>0.1){
          fix_h2 = true;
      }else{
          fix_h2 = false;
      }
      if(h2_1_update>0.1){
          h2_1_update = h2_1/num_SNP;
          sigma_1_update = sigma_1;
      }
      if(h2_2_update>0.1){
          h2_2_update = h2_2/num_SNP;
          sigma_2_update = sigma_2;
      }
     
     double llk;
      if(fix_h2){
           lambda_1 = 1.0/(eigval_1/sigma_1_update + 1.0/h2_1_update);
           UtXty_sigma_1 = UtXty_1/sigma_1_update;
           lambda_2 = 1.0/(eigval_2/sigma_2_update + 1.0/h2_2_update);
           UtXty_sigma_2 = UtXty_2/sigma_2_update;
           llk = -n_1/2*log(sigma_1_update)-num_SNP/2*log(h2_1_update)+accu(log(lambda_1))/2-n_1/2/sigma_1_update+accu(pow(UtXty_sigma_1,2)%lambda_1)/2 -n_2/2*log(sigma_2_update)-num_SNP/2*log(h2_2_update)+accu(log(lambda_2))/2-n_2/2/sigma_2_update+accu(pow(UtXty_sigma_2,2)%lambda_2)/2;
      }else{
          llk= max(llk_1) + max(llk_2);
      }
    
    
    arma::vec est(4);
      est(0) = h2_1_update*num_SNP;
      est(1) = h2_2_update*num_SNP;
      est(2) = sigma_1_update;
      est(3) = sigma_2_update;
  
      
     if(fix_intercept==false){
          Xty_1 = Xty_1/sigma_1_update;
          Xty_2 = Xty_2/sigma_2_update;
      }else{
          Xty_1 = Xty_1;
          Xty_2 = Xty_2;
      }
    arma::mat XtX_1_quad_1_XtX_1(num_SNP,num_SNP),XtX_2_quad_2_XtX_2(num_SNP,num_SNP),B_mat(num_SNP,num_SNP);
    double h2_1_score,h2_2_score,rho_score;
     h2_1_score = accu(pow(Xty_1,2))/2;
     h2_2_score = accu(pow(Xty_2,2))/2;
    
    
   if(fix_intercept==false){
     XtX_1_quad_1_XtX_1 = eigvec_1*diagmat(eigval_1/sigma_1_update - pow(eigval_1/sigma_1_update,2)/(eigval_1/sigma_1_update + 1.0/h2_1_update))*eigvec_1.t();
     XtX_2_quad_2_XtX_2 = eigvec_2*diagmat(eigval_2/sigma_2_update - pow(eigval_2/sigma_2_update,2)/(eigval_2/sigma_2_update + 1.0/h2_2_update))*eigvec_2.t();
     rho_score = dot(Xty_1.t()-UtXty_1.t()/sigma_1_update*(eigvec_1_trans.each_col()%(eigval_1/sigma_1_update/(eigval_1/sigma_1_update + 1.0/h2_1_update))),Xty_2.t()-UtXty_2.t()/sigma_2_update*(eigvec_2_trans.each_col()%(eigval_2/sigma_2_update/(eigval_2/sigma_2_update + 1.0/h2_2_update))));
    
   }else{
     XtX_1_quad_1_XtX_1 = eigvec_1*diagmat(eigval_1 - pow(eigval_1,2)/(eigval_1 + 1.0/h2_1_update))*eigvec_1.t();
     XtX_2_quad_2_XtX_2 = eigvec_2*diagmat(eigval_2 - pow(eigval_2,2)/(eigval_2 + 1.0/h2_2_update))*eigvec_2.t(); 
     rho_score = dot(Xty_1.t()-UtXty_1.t()*(eigvec_1_trans.each_col()%(eigval_1/(eigval_1 + 1.0/h2_1_update))),Xty_2.t()-UtXty_2.t()*(eigvec_2_trans.each_col()%(eigval_2/(eigval_2 + 1.0/h2_2_update))));
   }
  //  double score_stat_raw = dot(Xty_1.t()-h2_1_update*Xty_1.t()*inv_sympd(arma::eye(num_SNP, num_SNP)+ h2_1_update/sigma_1_update*XtX_1)*XtX_1/sigma_1_update,Xty_2.t()-h2_2_update*Xty_2.t()*inv_sympd(arma::eye(num_SNP, num_SNP)+ h2_2_update/sigma_2_update*XtX_2)*XtX_2/sigma_2_update);
    B_mat = XtX_1_quad_1_XtX_1 * XtX_2_quad_2_XtX_2;
    
    cx_vec eigval =eig_gen(B_mat); 
    vec evals = real(eigval);
    vec pos_evals = sqrt(evals(find(evals>0)))/2.0;
    pos_evals = join_cols(-pos_evals, pos_evals);
  
    arma::vec beta_1,beta_2;
    beta_1 = eigvec_1*(lambda_1%UtXty_1);
    beta_2 = eigvec_2*(lambda_2%UtXty_2);
  
    List output = List::create(
      _["est"] = est,
      _["h2_1_score"] = h2_1_score,
      _["h2_2_score"] = h2_2_score,
      _["rho_score"] = rho_score,
      _["h2_1_eigvals"] = eigval_1/2,
      _["h2_2_eigvals"] = eigval_2/2,
      _["rho_eigvals"] = pos_evals,
      _["llk_1"] = llk_1,
      _["llk_2"] = llk_2,
      _["llk_null"] = llk,
      _["fix_h2"] = fix_h2,
    _["beta_1"] = beta_1,
    _["beta_2"] = beta_2
    );
    return(output);
  }
/**************************************************************************************************
*
* 						Parameter Expansion EM algorithm
*
****************************************************************************************************/
double alt_rho_cpp(double rho,double h2_1,double h2_2,double sigma_1,double sigma_2,double n_1,double n_2,double num_SNP,arma::mat XtX,arma::vec Xty){
	arma::mat V_mat(2,2);
	V_mat(0,0) = h2_1;
    V_mat(1,1) = h2_2;
    V_mat(0,1) = rho;
    V_mat(1,0) = rho;
    arma::mat V_inv = inv_sympd(V_mat);
	arma::mat sigma_beta_inv = kron(V_inv,arma::eye(num_SNP,num_SNP))+XtX;
    arma::mat sigma_beta = inv_sympd(sigma_beta_inv);
	arma::vec mu_beta = sigma_beta*Xty;
	double llk = -n_1/2*log(sigma_1)-n_2/2*log(sigma_2)-num_SNP/2*log_det_sympd(V_mat)-0.5*log_det_sympd(sigma_beta_inv)-0.5*n_1/sigma_1-0.5*n_2/sigma_2+0.5*as_scalar(mu_beta.t()*sigma_beta_inv*mu_beta);
	return(-llk);
};
// [[Rcpp::export]]
SEXP PX_EM_alt(double h2_1,double h2_2,double rho,double sigma_1,double sigma_2,arma::vec z_1,arma::vec z_2, arma::mat R_1,arma::mat R_2,double n_1, double n_2,double num_SNP,int n_iter =1000,bool fix_intercept = false,bool fix_h2 = false){
  
    double h2_1_update,h2_2_update,rho_update,sigma_1_update,sigma_2_update,alpha_1,alpha_2;
    h2_1_update = h2_1/num_SNP;
    h2_2_update = h2_2/num_SNP;
    rho_update = rho*sqrt(h2_1)*sqrt(h2_2)/num_SNP;
    sigma_1_update = sigma_1;
    sigma_2_update = sigma_2;
    
    arma::mat V_mat(2,2),V_inv(2,2);
    arma::mat XtX(2*num_SNP,2*num_SNP,fill::zeros),sigma_beta_inv(2*num_SNP,2*num_SNP),sigma_beta(2*num_SNP,2*num_SNP);
    arma::mat I_p = arma::eye(num_SNP,num_SNP);
    
    arma::vec Xty(2*num_SNP),mu_beta(2*num_SNP);
    arma::vec llk(n_iter+1);
    double max_llk;
    llk(0) = -datum::inf;
    double trace_R_1_Sigma_1,trace_R_2_Sigma_2,mu_1_Sigma_mu_1,mu_2_Sigma_mu_2,mu_1_z_1,mu_2_z_2;
      if(fix_h2==true){
           XtX.submat(0,0,num_SNP-1,num_SNP-1) = R_1*n_1/sigma_1_update;
           XtX.submat(num_SNP,num_SNP,2*num_SNP-1,2*num_SNP-1) = R_2*n_2/sigma_2_update;
           Xty.subvec(0,num_SNP-1) = z_1/sigma_1_update*sqrt(n_1);
           Xty.subvec(num_SNP,2*num_SNP-1) = z_2/sigma_2_update*sqrt(n_2);
    
                // Extract R's optim function
              Rcpp::Environment stats("package:stats"); 
              Rcpp::Function optim = stats["optim"];
               Rcpp::List opt_results = optim(Rcpp::_["par"]    = rho_update,
                                   // Make sure this function is not exported!
                                   Rcpp::_["fn"]     = Rcpp::InternalFunction(&alt_rho_cpp),
                                   Rcpp::_["method"] = "Brent",
                                   // Pass in the other parameters as everything
                                   // is scoped environmentally
                                   Rcpp::_["h2_1"] = h2_1_update,
                                   Rcpp::_["h2_2"] = h2_2_update,
                                   Rcpp::_["sigma_1"] = sigma_1_update,
                                   Rcpp::_["sigma_2"] = sigma_2_update,                           
                                   Rcpp::_["n_1"] = n_1,
                                   Rcpp::_["n_2"] = n_2,
                                   Rcpp::_["num_SNP"] = num_SNP,
                                   Rcpp::_["XtX"] = XtX,
                                   Rcpp::_["Xty"] = Xty,
                                   Rcpp::_["lower"] = -1*sqrt(h2_1_update*h2_2_update),
                                   Rcpp::_["upper"] = 1*sqrt(h2_1_update*h2_2_update));
    
              // Extract out the estimated parameter values
          rho_update = as_scalar(Rcpp::as<arma::vec>(opt_results[0]));
          max_llk= Rcpp::as<double>(opt_results["value"]);
          llk(1) = max_llk;
          llk.resize(2);
          cout<<"rho update is "<<rho_update;
          cout<<"max_llk is "<<max_llk;
      }
  if(fix_h2==false){
    for(int i=0;i<n_iter;i++){
      
      V_mat(0,0) = h2_1_update;
      V_mat(1,1) = h2_2_update;
      V_mat(0,1) = rho_update;
      V_mat(1,0) = rho_update;
      V_inv = inv_sympd(V_mat);
      
      
      XtX.submat(0,0,num_SNP-1,num_SNP-1) = R_1*n_1/sigma_1_update;
      XtX.submat(num_SNP,num_SNP,2*num_SNP-1,2*num_SNP-1) = R_2*n_2/sigma_2_update;
      sigma_beta_inv = kron(V_inv,I_p)+XtX;
      sigma_beta = inv_sympd(sigma_beta_inv);
      
      Xty.subvec(0,num_SNP-1) = z_1/sigma_1_update*sqrt(n_1);
      Xty.subvec(num_SNP,2*num_SNP-1) = z_2/sigma_2_update*sqrt(n_2);
      mu_beta = sigma_beta*Xty;
      
  
      llk(i+1) = -n_1/2*log(sigma_1_update)-n_2/2*log(sigma_2_update)-num_SNP/2*log_det_sympd(V_mat)-0.5*log_det_sympd(sigma_beta_inv)-0.5*n_1/sigma_1_update-0.5*n_2/sigma_2_update+0.5*as_scalar(mu_beta.t()*sigma_beta_inv*mu_beta);
      if(abs(llk(i+1)-llk(i))<pow(10,-5)){
          llk.resize(i+1);
          break;	
      } 
      if(fix_h2==false){
          h2_1_update = (sum(square(mu_beta.subvec(0,num_SNP-1))) + trace(sigma_beta.submat(0,0,num_SNP-1,num_SNP-1)))/num_SNP;
          h2_2_update = (sum(square(mu_beta.subvec(num_SNP,2*num_SNP-1))) + trace(sigma_beta.submat(num_SNP,num_SNP,2*num_SNP-1,2*num_SNP-1)))/num_SNP;
      }
  
      rho_update = (dot(mu_beta.subvec(0,num_SNP-1),mu_beta.subvec(num_SNP,2*num_SNP-1)) + trace(sigma_beta.submat(0,num_SNP,num_SNP-1,2*num_SNP-1)))/num_SNP;
      
      
      mu_1_z_1 = dot(mu_beta.subvec(0,num_SNP-1),z_1)/sqrt(n_1);
      mu_2_z_2 = dot(mu_beta.subvec(num_SNP,2*num_SNP-1),z_2)/sqrt(n_2);  
      trace_R_1_Sigma_1 = trace(R_1*sigma_beta.submat(0,0,num_SNP-1,num_SNP-1));
      trace_R_2_Sigma_2 = trace(R_2*sigma_beta.submat(num_SNP,num_SNP,2*num_SNP-1,2*num_SNP-1));
      mu_1_Sigma_mu_1 = as_scalar(trans(mu_beta.subvec(0,num_SNP-1))*R_1*mu_beta.subvec(0,num_SNP-1));
      mu_2_Sigma_mu_2 = as_scalar(trans(mu_beta.subvec(num_SNP,2*num_SNP-1))*R_2*mu_beta.subvec(num_SNP,2*num_SNP-1));
      
      alpha_1 = mu_1_z_1/(mu_1_Sigma_mu_1 + trace_R_1_Sigma_1);
      alpha_2 = mu_2_z_2/(mu_2_Sigma_mu_2 + trace_R_2_Sigma_2);
      cout<<"alpha_update "<<alpha_1<<" "<<alpha_2<<endl;
      
      if(fix_intercept==false){
      sigma_1_update = 1 + pow(alpha_1,2)*mu_1_Sigma_mu_1-alpha_1*mu_1_z_1*2+pow(alpha_1,2)*trace_R_1_Sigma_1;
      sigma_2_update = 1 + pow(alpha_2,2)*mu_2_Sigma_mu_2-alpha_2*mu_2_z_2*2+pow(alpha_2,2)*trace_R_2_Sigma_2;
      }
      
   cout<<"sigma_1_update "<<sigma_1_update<<"sigma_2_updat "<<sigma_2_update<<endl;
      if(fix_h2==false){
      h2_1_update = pow(alpha_1,2)*h2_1_update;
      h2_2_update = pow(alpha_2,2)*h2_2_update;
      rho_update = alpha_1*alpha_2*rho_update;
      }
      cout<<"h2_1_update "<<h2_1_update<<"h2_2_update "<<h2_2_update<<endl;
    }
    max_llk = max(llk);
  }
    arma::vec est(5);
    est(0) = h2_1_update*num_SNP;
    est(1) = h2_2_update*num_SNP;
    est(2) = rho_update*num_SNP;
    est(3) = sigma_1_update;
    est(4) = sigma_2_update;
    
   // cout<<"llk "<<llk<<endl;
    List output = List::create(
      _["llk"] = llk,
      _["llk_alt"] = max_llk,	
      _["est"] = est,
      _["beta_1"] = mu_beta.subvec(0,num_SNP-1),
      _["beta_2"] = mu_beta.subvec(num_SNP,2*num_SNP-1)
    );
    return(output);
    
  }
  