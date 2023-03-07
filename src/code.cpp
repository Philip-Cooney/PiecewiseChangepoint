// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]



#include<RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
NumericMatrix surv_nochange(NumericVector time, int num_zero, NumericVector lambda_df) {

    int time_length = time.length();

    NumericMatrix St( time_length , num_zero);


    for (int i=0; i < num_zero; ++i){
        for(int j = 0; j < time_length; ++j ){
            St(j,i) = exp(-lambda_df(i)*time(j));
        }
    }
    return St;
}



// [[Rcpp::export]]
NumericMatrix surv_change(NumericVector time, int nrow_curr, NumericMatrix lambda_df, NumericMatrix changepoint_df, NumericMatrix surv_df) {

    int time_length = time.length();

    NumericMatrix St(time_length, nrow_curr);


    for (int i=0; i < nrow_curr; ++i){
        for(int j = 0; j < time_length; ++j ){
            LogicalVector less_chng1;
            less_chng1 = time(j) < changepoint_df(i,_);
            if(all(less_chng1).is_true()){
                St(j,i) =  exp(-lambda_df(i,1)*time(j));
            }else{

                int index;
                LogicalVector less_chng = (time(j) >= changepoint_df(i,_));
                index = sum(less_chng)-1;

                double time_diff = time(j) - changepoint_df(i,index);
                St(j,i) = surv_df(i,index)*exp(-lambda_df(i,index)*time_diff);

            }
        }
    }
    return St;
}


// [[Rcpp::export]]
NumericMatrix exposure_death_alt(NumericMatrix total_time_diff, IntegerVector change_nums) {

    int n_event; int num_breaks; NumericVector total_time_diff1;
    NumericVector total_time_diff2;NumericVector total_time_diff3;

    NumericVector n_event_vec = total_time_diff(_,1);
    NumericVector total_time_diff_vec = total_time_diff(_,0);
    n_event = sum(n_event_vec);

    if(IntegerVector::is_na(change_nums(0))){
        num_breaks = 0;
    }else{
        num_breaks =  change_nums.length();
    }

    NumericMatrix res_matrix(num_breaks+1,2);

    if(num_breaks == 0){
        res_matrix(0,0) = n_event;
        res_matrix(0,1) = sum(total_time_diff_vec);
    }else if(num_breaks == 1){

        res_matrix(0,0) = sum(n_event_vec[Range(0,change_nums(0)-1)]);
        res_matrix(1,0) = n_event - res_matrix(0,0);

        total_time_diff1 = total_time_diff_vec[Range(0,change_nums(0)-1)];
        res_matrix(0,1) = sum(total_time_diff1);


        total_time_diff2 = total_time_diff_vec[Range(change_nums(0), total_time_diff_vec.length()-1)];
        res_matrix(1,1) = sum(total_time_diff2);

    }else{
        for (int i=0; i <= num_breaks; ++i){
            //Rcout << "The value of i : " << i << "\n";

            if(i == 0){
                res_matrix(i,0) = sum(n_event_vec[Range(0,change_nums(0)-1)]);
                total_time_diff1 = total_time_diff_vec[Range(0,change_nums(0)-1)];
                res_matrix(i,1) = sum(total_time_diff1);

            }else if (i == num_breaks){
                res_matrix(i,0) = sum(n_event_vec[Range(change_nums(i-1),n_event_vec.length()-1)]) ;
                res_matrix(i,1) = sum(total_time_diff_vec[Range(change_nums(i-1),total_time_diff_vec.length()-1)]);

            }else{
                res_matrix(i,0) = sum(n_event_vec[Range(change_nums(i-1),change_nums(i)-1)]) ;
                res_matrix(i,1) = sum(total_time_diff_vec[Range(change_nums(i-1),change_nums(i)-1)]) ;
            }
        }
    }

    return res_matrix;

}


//
//double marg_lik_eval_log( NumericVector x, int alpha, int beta) {

//  long double part1 = lgamma(x[0] + 1) ;
//  long double part2 = (alpha + x[0])*log(beta + x[1]);

//  return part1-part2;
//}


double marg_lik_eval_log( NumericVector x, double alpha, double beta) {

    long double part1 = lgamma(x[0] + alpha) - lgamma(alpha);
    long double part2 = alpha*log(beta) - (alpha + x[0])*log(beta + x[1]);

    return part1+part2;
}

// [[Rcpp::export]]
double margin_lik_calc_log(NumericMatrix time_diffs,
                           IntegerVector k_indexes, double alpha, double beta) {

    NumericMatrix res_matrix = exposure_death_alt(time_diffs, k_indexes);

    double marg_lik_tot = 0;
    for(int i = 0; i < res_matrix.nrow(); ++i){
        NumericVector res_slice = NumericVector::create(res_matrix(i,0),res_matrix(i,1));
        double marg_lik_slice = marg_lik_eval_log(res_slice,alpha,beta);
        marg_lik_tot += marg_lik_slice;
    }
    return  marg_lik_tot;

}



IntegerVector drop_val(IntegerVector x, IntegerVector val){

    IntegerVector vec = match(x,val);
    IntegerVector vec2(vec.length()-1);
    int index = 0;
    for(int i = 0; i < vec.length(); ++i){
        if(vec(i) != 1){
            vec2(index) = x(i);
            index++ ;
        }
    }

    return vec2;
}


// [[Rcpp::export]]
double pos_prob(IntegerVector k_indexes, int n_events, LogicalVector MLE){

    double length_k = k_indexes.length();

    if(all(MLE).is_true()){
        //return 1;

        //IntegerVector  k_indexes2 (length_k+1,0);


        IntegerVector k_indexes2 =  k_indexes;

        k_indexes2.push_front(0);
        NumericVector log_prob_eval(length_k, 0);
        for(int i  =0; i < length_k; ++i){
            log_prob_eval(i) = log(1/(n_events -k_indexes2(i) -(length_k-i)));
        }

        return(exp(sum(log_prob_eval)));

    }else{



        IntegerVector part1 = k_indexes;
        IntegerVector part2 (length_k+1,0);

        for(int i = 0; i < part1.length(); ++i){
            part2(i) = part1(i);
        }
        part2(part1.length()) = n_events;

        IntegerVector part3 = diff(part2);

        IntegerVector k_diff (length_k+1,0);
        k_diff(0) = k_indexes(0)-1;
        for(int i = 1; i <= part3.length(); ++i){
            k_diff(i) =  part3(i-1)-1;
        }

        int neg_val = 0;

        for(int i = 0; i < k_diff.length(); ++i){
            if(k_diff(i)<0){
                neg_val++ ;
            }
        }

        if(neg_val >0 ){
            return(0);
        }else{
            int vn = n_events-1;
            int vk = 2*length_k +1;

            double prod_val =1;
            for(int i = 0; i < k_diff.length(); ++i){
                prod_val= prod_val*k_diff(i);
            }

            return((1/Rf_choose(vn,vk)*prod_val));
        }

    }

}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//



IntegerVector return_pos_df(NumericVector LL, NumericVector pos_prob_vec, IntegerMatrix eval_df){
    if(eval_df.nrow() == 1){
        return  eval_df(0,_);
    }else{
        IntegerVector new_loc;
        NumericVector LL_ratio (LL.length());
        NumericVector LL_final (LL.length());
        double LL_max = max(LL);

        double sum = 0;

        for(int i = 0; i < LL_ratio.length(); ++i){

            //if(NumericVector::is_na(LL[i]))
            //Rprintf("v[%i] is NA.\n", i);
            //if(Rcpp::traits::is_nan<REALSXP>(LL[i]))
            //Rprintf("v[%i] is NaN.\n", i);
            //if(Rcpp::traits::is_infinite<REALSXP>(LL[i]))
            //Rprintf("v[%i] is Inf or -Inf.\n", i);
            if(NumericVector::is_na(LL[i])){
                LL_ratio(i) =0;
            }else{
                LL_ratio(i) = exp(LL(i)- LL_max)*pos_prob_vec(i);
            }

            sum = LL_ratio(i) + sum  ;
        }

        for(int i = 0; i < LL_ratio.length(); ++i){
            LL_final(i) = LL_ratio(i)/sum ;
        }
        //Rcout << "The value of v : " << eval_df << "\n";
        //Rcout << "The value of v : " << LL_ratio << "\n";

        // Rcout << "The value of v : " << LL_ratio << "\n";
        IntegerVector eval_index = seq(0,LL_ratio.length()-1);
        IntegerVector chosen_index = sample(eval_index, 1,false,LL_final);

        return eval_df(chosen_index(0),_);
    }

}



IntegerVector return_pos_vec(NumericVector LL, NumericVector pos_prob_vec, IntegerVector evals){

    if(evals.length() == 1){
        return  IntegerVector::create(evals(0));
    }else{
        IntegerVector new_loc;
        NumericVector LL_ratio (LL.length());
        NumericVector LL_final (LL.length());
        double LL_max = max(LL);

        double sum = 0;

        for(int i = 0; i < LL_ratio.length(); ++i){
            if(NumericVector::is_na(LL[i]))
                Rprintf("v[%i] is NA.\n", i);
            if(Rcpp::traits::is_nan<REALSXP>(LL[i]))
                Rprintf("v[%i] is NaN.\n", i);
            if(Rcpp::traits::is_infinite<REALSXP>(LL[i]))
                Rprintf("v[%i] is Inf or -Inf.\n", i);

            LL_ratio(i) = exp(LL(i)- LL_max)*pos_prob_vec(i);
            sum = LL_ratio(i) + sum  ;
        }

        for(int i = 0; i < LL_ratio.length(); ++i){
            LL_final(i) = LL_ratio(i)/sum ;
        }

        IntegerVector eval_index = seq(0,LL_ratio.length()-1);

        IntegerVector chosen_index = sample(eval_index, 1,false,LL_final);

        return IntegerVector::create(evals(chosen_index(0)));

    }


}


double sumprod(arma::vec x, arma::vec y){

    arma::vec res = x%y;

    return sum(res);
}




double piecewise_loglik_alt_cpp(NumericMatrix time_diffs, IntegerVector change_nums, NumericVector lambda) {

    NumericMatrix res_matrix = exposure_death_alt(time_diffs, change_nums);

    //Rcout << "The value of res_matrix : " << res_matrix << "\n";
    //Rcout << "The value of n_events : " << l1 << "\n";

    if(NumericVector::is_na(lambda(0))){


        NumericVector lambda_MLE (res_matrix.nrow());


        int iters = res_matrix.nrow()-1 ;
        for(int j = 0; j <= iters; ++j){

            lambda_MLE(j) = res_matrix(j,0)/res_matrix(j,1);
        }

        //Rcout << "The value of lambda_MLE : " << lambda_MLE << "\n";

        double log_lik = sumprod(res_matrix(_,0),log(lambda_MLE)) - sumprod(res_matrix(_,1),lambda_MLE);

        return(log_lik);

    }else{

        double log_lik = sumprod(res_matrix(_,0),log(lambda)) - sumprod(res_matrix(_,1),lambda);
        return(log_lik);
    }



}


double beta_samp(NumericMatrix time_diffs,
                 IntegerVector k_indexes, double alpha, double curr_beta,double beta1, double beta2) {
    double beta_new;
    //Rcout << "k_indexes? : " << k_indexes << "\n";
    //Rcout << "isna? : " <<NumericVector::is_na(k_indexes(0))<< "\n";
    if(IntegerVector::is_na(k_indexes(0))){
        NumericMatrix res_matrix = exposure_death_alt(time_diffs, k_indexes);

        double shape_par = alpha + res_matrix(0,0);
        double rate_par =  1/(curr_beta + res_matrix(0,1));// 1/rate
        double lambda_curr;

        lambda_curr= R::rgamma(shape_par, rate_par);

        double shape_beta = alpha + beta1;
        double rate_beta =  1/(lambda_curr + beta2);// 1/rate
        beta_new = R::rgamma(shape_beta, rate_beta);

        //Rcout << "The value of lambda_curr : " << lambda_curr << "\n";
        //Rcout << "The value of shape_par : " << shape_par << "\n";
        //Rcout << "The value of rate_par : " << rate_par << "\n";
        //Rcout << "The value of res_matrix(1,0) : " << res_matrix(1,0)<< "\n";


    }else{
        IntegerVector k_indexes_omit = k_indexes[!is_na(k_indexes)];
        int num_breaks = k_indexes_omit.length();
        NumericVector lambda_curr(num_breaks+1);

        //Rcout << "The value of time_diffs : " << kt_curr << "\n";
        NumericMatrix res_matrix = exposure_death_alt(time_diffs, k_indexes_omit);


        for(int q = 0; q <= num_breaks; ++q ){

            double shape_par = alpha + res_matrix(q,0);
            double rate_par =  1/(curr_beta + res_matrix(q,1));// 1/rate
            lambda_curr(q) =  R::rgamma(shape_par, rate_par);
            //Rcout << "The value of lambda_curr_vec : " << lambda_curr << "\n";
            //Rcout << "The value of shape_par : " << shape_par << "\n";
            //Rcout << "The value of rate_par : " << rate_par << "\n";

        }

        double lambda_curr_sum = sum(lambda_curr);
        //Rcout << "The value of sum lambda_curr_sum : " << lambda_curr_sum << "\n";
        double shape_beta = alpha*(num_breaks+1) + beta1;
        //Rcout << "The value of sum shape_beta : " << shape_beta << "\n";

        double rate_beta =  1/(lambda_curr_sum + beta2);// 1/rate
        beta_new = R::rgamma(shape_beta,rate_beta);

    }

    return beta_new;

}





IntegerVector gibbs_update(NumericMatrix time_diffs, IntegerVector kt, NumericVector lambda, int j_pos, LogicalVector MLE, double alpha, double beta){

    NumericVector LL;
    NumericVector n_event_vec = time_diffs(_,1);
    int n_events = sum(n_event_vec);
    int events_length = n_event_vec.length();
    //Rcout << "The value of n_events : " << n_events << "\n";
    IntegerVector pos_change = seq(2,events_length-1);
    IntegerVector evals;
    IntegerMatrix evals_df;
    IntegerVector drop_value;
    IntegerVector remain_vec;
    IntegerVector mat;
    IntegerMatrix mat_true2;
    int curr_num_breaks = kt.length();
    int j = j_pos;

    if(curr_num_breaks == 1){ // One Interval
        evals = pos_change;
        NumericVector LL (evals.length(),-10000000);
        NumericVector pos_prob_vec(evals.length());

        for(int i = 0; i < evals.length(); ++i){
            IntegerVector evals_i = IntegerVector::create(evals(i));

            if(NumericVector::is_na(lambda(0))){
                LL(i) = margin_lik_calc_log(time_diffs,evals_i,alpha, beta);
            }else{
                LL(i) = piecewise_loglik_alt_cpp(time_diffs,evals_i,lambda);
                //Rcout << "The value of evals_i : " << evals_i << "\n";

            }
            pos_prob_vec(i) = pos_prob(evals_i, events_length,MLE);
        }
        //Rcout << "The value of LL : " << LL << "\n";
        //Rcout << "The value of pos_prob_vec : " << pos_prob_vec << "\n";

        return return_pos_vec(LL, pos_prob_vec, evals);

    }else if (curr_num_breaks != 1 && j == curr_num_breaks){ //Final Interval

        evals = seq(kt(j-2)+1,events_length-1);
        //Rcout << "The value of evals : " << evals << "\n";
        drop_value = {kt(j-1)};
        remain_vec = drop_val(kt, drop_value);

        mat = rep(remain_vec, evals.length());
        mat.attr("dim") = Dimension( remain_vec.length(),evals.length());

        IntegerMatrix mat_true = as<IntegerMatrix>(mat);
        mat_true2 = transpose(mat_true);
        evals_df = cbind(mat_true2,evals);

        NumericVector LL (evals.length(),-10000000);
        NumericVector pos_prob_vec(evals.length());
        for(int i = 0; i < evals.length(); ++i){
            IntegerVector evals_i = evals_df(i,_);
            //Rcout << "The value of evals_i : " << evals_i << "\n";


            if(NumericVector::is_na(lambda(0))){
                LL(i) = margin_lik_calc_log(time_diffs,evals_i, alpha, beta);
            }else{
                LL(i) = piecewise_loglik_alt_cpp(time_diffs,evals_i,lambda);
            }

            pos_prob_vec(i) = pos_prob(evals_i, events_length, MLE);
        }

        // Rcout << "The value of LL : " << LL << "\n";
        //Rcout << "The value of pos_prob_vec : " << pos_prob_vec << "\n";
        return return_pos_df(LL, pos_prob_vec, evals_df);

        //return evals_df;
        //return pos_prob_vec;
        //return LL;
        //evals_df
        //drop the value you don't need. Rep the row a number of times and create a matrix.
        // cbind the result and apply the sort function
        // rep()

    }else if (curr_num_breaks != 1 && j == 1){ //First Interval
        evals = seq(2,kt(j)-1);
        drop_value = {kt(j-1)};
        remain_vec = drop_val(kt, drop_value);

        mat = rep(remain_vec, evals.length());
        mat.attr("dim") = Dimension( remain_vec.length(),evals.length());

        IntegerMatrix mat_true = as<IntegerMatrix>(mat);
        mat_true2 = transpose(mat_true);
        IntegerMatrix mat_true3 = cbind(mat_true2,evals);
        IntegerMatrix evals_df( evals.length(), kt.length());

        for(int i = 0; i < evals.length(); ++i){

            IntegerVector order_row = sort_unique(mat_true3(i,_));
            evals_df(i,_) = order_row;

        }

        NumericVector LL (evals.length(),-10000000);
        NumericVector pos_prob_vec(evals.length());

        for(int i = 0; i < evals.length(); ++i){
            IntegerVector evals_i = evals_df(i,_);

            if(NumericVector::is_na(lambda(0))){
                //Rcout << "The value of evals_i : " << evals_i << "\n";
                LL(i) = margin_lik_calc_log(time_diffs,evals_i, alpha, beta);
            }else{
                //Rcout << "The value of evals_i : " << evals_i << "\n";
                LL(i) = piecewise_loglik_alt_cpp(time_diffs,evals_i,lambda);
                //if(NumericVector::is_na(LL[i])){
                //}
            }

            pos_prob_vec(i) = pos_prob(evals_i, events_length, MLE);

        }

        //Rcout << "The value of LL : " << LL << "\n";
        //Rcout << "The value of pos_prob_vec : " << pos_prob_vec << "\n";

        return return_pos_df(LL, pos_prob_vec, evals_df);


    }else{ //Middle interval

        evals = seq(kt(j-2)+1,kt(j)-1);
        //Rcout << "The value of evals : " << evals << "\n";
        drop_value = {kt(j-1)};
        remain_vec = drop_val(kt, drop_value);

        mat = rep(remain_vec, evals.length());
        mat.attr("dim") = Dimension( remain_vec.length(),evals.length());

        IntegerMatrix mat_true = as<IntegerMatrix>(mat);
        mat_true2 = transpose(mat_true);
        IntegerMatrix mat_true3 = cbind(mat_true2,evals);
        IntegerMatrix evals_df( evals.length(), kt.length());

        for(int i = 0; i < evals.length(); ++i){

            IntegerVector order_row = sort_unique(mat_true3(i,_));
            evals_df(i,_) = order_row;

        }

        NumericVector LL (evals.length(),-10000000);
        NumericVector pos_prob_vec(evals.length());

        for(int i = 0; i < evals.length(); ++i){
            IntegerVector evals_i = evals_df(i,_);
            //Rcout << "The value of evals_i : " << evals_i << "\n";

            if(NumericVector::is_na(lambda(0))){
                LL(i) = margin_lik_calc_log(time_diffs,evals_i, alpha, beta);
            }else{
                LL(i) = piecewise_loglik_alt_cpp(time_diffs,evals_i,lambda);
            }

            pos_prob_vec(i) = pos_prob(evals_i, events_length,MLE);
        }

        return return_pos_df(LL, pos_prob_vec, evals_df);
        //return evals_df;
        //return pos_prob_vec;
        //return LL;
    }
    //return LL;
}


// [[Rcpp::export]]
List RJMCM_core(IntegerMatrix k,int max_breaks,
                NumericMatrix time_diffs, LogicalVector MLE, double alpha, double beta1, double beta2, int lambda =1){

    NumericVector n_event_vec = time_diffs(_,1);
    int n_events = sum(n_event_vec);
    int events_length = n_event_vec.length();
    IntegerVector pos_change = seq(2,events_length-1);

    NumericVector beta_vec (k.nrow());

    //double shape_par = alpha + res_matrix(q,0);
    //double rate_par =  1/(curr_beta + res_matrix(q,1));// 1/rate
    //Initialize beta
    beta_vec(0) =  R::rgamma(beta1, (1/beta2));

    //Rcout << "The value of  beta_vec(0) : " <<  beta_vec(0) << "\n";

    //int i = 1;
    for(int i = 1; i < k.nrow(); ++i){
        //Rcout << "The value of  i : " <<  i << "\n";

        double a_k; double d_k; double a_k_1; double d_k_1; double a_k1;
        double d_k1; IntegerVector k_add; IntegerVector k_propose; LogicalVector add_del_ind;
        IntegerVector index_del; double log_marg_lik_current; double log_marg_lik_propose;
        double pos_prob_propose; double pos_prob_current; double pos_prob_factor; double prob_factor;
        IntegerVector k_current; int curr_num_breaks; double change_prior;
        double P_k_to_k1;
        double P_k1_to_k;
        double A;
        double A_calc;
        NumericVector A_vec;
        double P_k_1_to_k;
        double P_k_to_k_1;
        NumericVector lambda_na = NumericVector::create(NA_REAL);
        double curr_beta = beta_vec(i-1);
        //Rcout << "The value of  curr_beta : " <<  curr_beta << "\n";

        IntegerVector curr_breaksna = k(i-1, _ );
        LogicalVector which_na = is_na(curr_breaksna);
        if(all(which_na).is_true()){
            k_current = {NA_INTEGER};
            curr_num_breaks = 0;
        }else{
            k_current = curr_breaksna[!is_na(curr_breaksna)];
            curr_num_breaks = k_current.length();
        }

        //Choose to add or delete a changepoint
        //a_k (add), d_k (delete)

        if(curr_num_breaks == 0){
            a_k = 1;
        }else if(curr_num_breaks == max_breaks){
            a_k = 0;
        }else{
            a_k = 0.5 ;
        }

        d_k = 1 - a_k;

        //Now for the K-1 changepoint we need to know what
        //the add and deletion proabiblities are

        if(curr_num_breaks - 1 == 0){
            a_k_1 = 1 ;
        }else{
            a_k_1 = 0.5 ;
        }

        d_k_1 = 1- a_k_1 ;

        //Now for the K+1 changepoint.
        if(curr_num_breaks +1 == max_breaks){
            a_k1 = 0;
        }else{
            a_k1 = 0.5;
        }

        d_k1 =  1- a_k1;//deletion probability of k+1 changepoints

        double u = runif(1)(0);

        if(u < a_k){
            //Sample a changepoint to add

            k_add = sample(pos_change,1);

            if(curr_num_breaks == 0){
                k_propose = k_add;
            }else{
                IntegerVector k_propose_unsort = union_(k_current,k_add);
                k_propose = sort_unique(k_propose_unsort);
            }
            add_del_ind = true; //Used later (indicator)

        }else{
            if(curr_num_breaks == 1){
                k_propose = NA_INTEGER;
            }else{
                index_del = sample(k_current, 1);

                k_propose = drop_val(k_current, index_del);

            }
            add_del_ind = false;

        }
        log_marg_lik_current = margin_lik_calc_log(time_diffs, k_current, alpha, curr_beta);
        log_marg_lik_propose = margin_lik_calc_log(time_diffs, k_propose, alpha, curr_beta);
        //Rcout << "The value of  log_marg_lik_current : " <<  log_marg_lik_current << "\n";
        //Rcout << "The value of  log_marg_lik_propose : " <<  log_marg_lik_propose << "\n";


        if(IntegerVector::is_na(k_propose(0))){
            change_prior = R::dpois(0,lambda,false);
            pos_prob_propose = 1*change_prior;
        }else{
            change_prior = R::dpois(k_propose.length(),lambda,false);
            pos_prob_propose = pos_prob(k_propose,events_length,MLE)*change_prior;
        }


        if(IntegerVector::is_na(k_current(0))){
            change_prior = R::dpois(0,lambda,false);
            pos_prob_current = 1;
        }else{
            change_prior = R::dpois(k_current.length(),lambda,false);
            pos_prob_current = pos_prob(k_current,events_length,MLE)*change_prior;
        }


        pos_prob_factor = pos_prob_propose/pos_prob_current ;

        prob_factor = pos_prob_factor*(exp(log_marg_lik_propose-log_marg_lik_current));


        if(all(add_del_ind).is_true() ){
            P_k_to_k1 = a_k/(events_length-curr_num_breaks-1);

            P_k1_to_k = d_k1/(curr_num_breaks+1);

            A_calc = prob_factor*(P_k1_to_k/P_k_to_k1);
            A_vec = {A_calc, 1};
            A = min(A_vec);
            //Rcout << "The value of  A : " <<  A << "\n";

        }else{
            P_k_1_to_k = a_k_1/(events_length-curr_num_breaks);
            P_k_to_k_1 = d_k/curr_num_breaks;

            A_calc = prob_factor*(P_k_1_to_k/P_k_to_k_1);
            A_vec = {A_calc, 1};
            A = min(A_vec);
            //Rcout << "The value of  A : " <<  A << "\n";

        }

        double u2 = runif(1)(0);

        IntegerVector kt = rep(NA_INTEGER, max_breaks);

        if(u2 < A){

            LogicalVector which_na2 = is_na(k_propose);
            if(all(which_na2).is_true()){

            }else{

                IntegerVector samp = seq(1, k_propose.length());
                IntegerVector j_vec = sample(samp,1);
                IntegerVector k_propose_new = gibbs_update(time_diffs,k_propose,lambda_na, j_vec(0), MLE, alpha, curr_beta);
                //Rcout << "The value of  k_propose_new : " <<  k_propose_new << "\n";


                for(int j = 0; j < k_propose_new.length(); ++j){
                    kt(j) = k_propose_new(j);
                }
            }

        }else{

            LogicalVector which_na2 = is_na(k_current);
            if(all(which_na2).is_true()){

            }else{

                IntegerVector samp = seq(1, k_current.length());
                IntegerVector j_vec = sample(samp,1);
                IntegerVector k_current_new = gibbs_update(time_diffs,k_current,lambda_na,j_vec(0),MLE,alpha, curr_beta);
                //Rcout << "The value of  k_current_new : " <<  k_current_new << "\n";


                for(int j = 0; j < k_current_new.length(); ++j){

                    //Rcout << "The value of  k_current_new.length()" << k_current_new.length() << "\n";
                    kt(j) = k_current_new(j);
                    //Rcout << "The value of  kt(i)" << kt(j) << "\n";
                }
            }
        }


        beta_vec(i) = beta_samp(time_diffs,kt,alpha,curr_beta,beta1,beta2);
        k(i,_) =  kt;


    }
    List res;
    res["k"] = k;
    res["beta"] = beta_vec;
    return res;
}


// [[Rcpp::export]]
List Gibbs_core(IntegerMatrix k,int num_breaks,
                NumericMatrix time_diffs,NumericMatrix alpha_array, NumericMatrix beta_array, LogicalVector MLE, double alpha, double beta){

    IntegerVector kt;
    NumericMatrix res_array;
    NumericMatrix lambda(k.nrow(),k.ncol()+1);
    NumericVector max_lik(k.nrow());
    NumericVector marg_lik_log(k.nrow());
    NumericVector lambda_na  = NumericVector::create(NA_REAL);
    NumericVector n_event_vec = time_diffs(_,1);
    //int n_events = sum(n_event_vec);

    for(int i = 1; i < k.nrow(); ++i){

        kt = k(i-1,_);

        IntegerVector kt_curr; IntegerVector kt_new;

        for(int j = 0; j < num_breaks; ++j){

            if(j == 0){
                kt_curr = kt;
            }else{
                kt_curr = kt_new;
            }
            //Rcout << "The value of time_diffs : " << kt_curr << "\n";
            NumericMatrix res_array = exposure_death_alt(time_diffs, kt_curr);
            NumericVector lambda_curr(num_breaks+1);

            for(int q = 0; q <= num_breaks; ++q ){

                double shape_par = alpha_array(i,q) + res_array(q,0);
                double rate_par =  1/(beta_array(i,q) + res_array(q,1));// 1/rate
                lambda_curr(q) =  R::rgamma(shape_par, rate_par);

                if(Rcpp::traits::is_infinite<REALSXP>(lambda_curr(q))){
                    //Rcout << "The value of res_array : " << res_array << "\n";
                    //Rcout << "The value of beta_array(i,q) : " << beta_array(i,q) << "\n";
                    //Rcout << "The value of alpha_array(i,q) : " << alpha_array(i,q) << "\n";

                }

            }
            //Rcout << "The value of time_diffs : " << kt_curr << "\n";
            kt_new = gibbs_update(time_diffs, kt_curr,lambda_curr, j + 1, MLE, alpha, beta);
            //Rcout << "The value of time_diffs : " << time_diffs << "\n";
        }

        k(i,_) = kt_new;

        max_lik(i) = piecewise_loglik_alt_cpp(time_diffs,kt_new,lambda_na);

        // Set pos prob to be equal to 1
        marg_lik_log(i) = margin_lik_calc_log(time_diffs,kt_new, alpha, beta)+log(pos_prob(kt_new, n_event_vec.length(), true));

        //Final lambda

        NumericMatrix res_array_final = exposure_death_alt(time_diffs, kt_curr);
        NumericVector lambda_curr_final(num_breaks+1);

        for(int q = 0; q <= num_breaks; ++q ){

            double shape_par = alpha_array(i,q) + res_array_final(q,0);
            double rate_par =  1/(beta_array(i,q) + res_array_final(q,1));// 1/rate
            lambda_curr_final(q) =  R::rgamma(shape_par, rate_par);

        }

        lambda(i,_) = lambda_curr_final;


    }

    List res;
    res["k"] = k;
    res["max_lik"] = max_lik;
    res["marg_lik_log"] = marg_lik_log;
    res["lambda"] = lambda;

    return res;


}





