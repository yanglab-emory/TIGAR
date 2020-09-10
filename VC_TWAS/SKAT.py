# %%
#!/usr/bin/env python
# -*- coding: utf-8 -*-
##########################################################################################
# import packages needed
import pandas as pd
import numpy as np
from numpy import mean, sqrt, std, where, diag,sum
import statsmodels.api as sm
import scipy.stats 
import qfc_checked
#######################################################################################

# %%
############# get s^2,residual,and cov matrix for continuous phenotype #############
def get_linear_residual(pheno, cov, outtype):       
    if len(cov) != 0:
        cov_x = sm.add_constant(cov)
    else:
        cov_x = np.ones(len(pheno))
    lm = sm.OLS(pheno,cov_x).fit()
    s2 = lm.mse_resid
    res = lm.resid
    return s2, res, cov_x       

# %%
############# get cov matrix , mu, pi_1, res, res_out(resampling if sample size <2000) dichotomous phenotype #############
def get_logistic_residual(pheno, cov, n_resampling):       
    if len(cov) != 0:
        cov_x = sm.add_constant(cov)
    else:
        cov_x = np.ones(len(pheno))
    ####from logistic regression#####
    lm = sm.Logit(pheno,cov_x).fit()
    mu = lm.predict()
    pi_1 = mu*(1-mu)
    res = pheno.values.flatten() - mu 
    n1 = len(res)
    n_case = sum(pheno)
    res_out = 0
    ####resampling####
    if n_resampling > 0:
        res_out = np.zeros((n1,n_resampling))
        for i in range(n_resampling):
            res_out1 = np.random.binomial(1,mu,n1)
            res_out2 = np.random.binomial(1,mu,n1)
            id_case1 = where(res_out1==1)
            id_case2 = where(res_out2==1)
            id_c1 = set(id_case1[0]).intersection(id_case2[0]) # intersection 
            id_c1 = np.array(list(id_c1))
            set1 = set(id_case1[0]).difference(set(id_case2[0]))
            set2 = set(id_case2[0]).difference(set(id_case1[0]))
            id_c2 = set1.union(set2)
            id_c2 = np.array(list(id_c2))
            if n_case[0] <= len(id_c1):
                id_case = np.random.choice(id_c1, size = n_case[0])
            elif (n_case[0] > len(id_c1)) and n_case[0] <= (len(id_c1) + len(id_c2)):
                id_c3 = np.random.choice(id_c2, size = n_case[0] - len(id_c1),replace = False)
                id_case = np.concatenate((id_c1,id_c3))
            else:
                id_case3 = set(id_c1).union(set(id_c2))
                id_c4 = set(range(n1)).difference(id_case3)
                n_needed = n_case[0] - len(np.array(list(id_case3)))
                id_c4 = np.array(list(id_c4))
                id_c5 = np.random.choice(id_c4, size = n_needed, p = mu[id_c4]/sum(mu[id_c4]))
                id_case = id_case3.union(set(id_c5))
                id_case = np.array(list(id_case))
            res_out[id_case, i] = 1
            res_out [:,i] = res_out [:,i] - mu                 
    return cov_x, mu, pi_1, res, res_out

# %%
############# get residual (resampling & without resampling) #############
def SKAT_Null_Model_MomentAdjust(pheno, cov, n_resampling):
    re1 = get_logistic_residual(pheno, cov, 0)
    re2 = get_logistic_residual(pheno, cov, n_resampling)
    return (re1, re2)

# %%
############# SKAT Null Model #############
def SKAT_Null_Model(pheno, cov, outtype):
    regdata = pd.concat([pheno, cov], axis=1)
    samplesize = len(regdata.dropna())
    n_resampling = 0
    if outtype == "D":
        if samplesize < 2000:
            n_resampling=10000
            re = SKAT_Null_Model_MomentAdjust(pheno, cov, n_resampling)
        else:
            re = get_logistic_residual(pheno, cov, 0)
    if outtype =="C":
        re = get_linear_residual(pheno, cov, outtype)
    return re, n_resampling

# %%
############# logistic model with adjustment (sample size < 2000) #############
def SKAT_With_NullModel_ADJ(re, weight, geno):
    obj_res = re[0]
    res = obj_res[3]
    res2 = re[1][4]
    X1 = obj_res[0]
    pi_1 = obj_res[2]
    mu = obj_res[1]
    res_moment = res2
    Z1 = weight * geno
    Q1 = (res).dot(Z1)
    Q = Q1.dot(Q1.T) / 2
    Z2 = (sqrt(pi_1) * Z1.T).T - ((sqrt(pi_1) * X1.T).T).dot(np.linalg.inv((X1.T).dot((pi_1 * X1.T).T))).values.dot((X1.T).dot((pi_1 * Z1.T).T).values)
    Q_Temp_res1 = (res_moment.T).dot(Z1)
    Q_sim = pow(Q_Temp_res1, 2).sum(axis = 1) / 2
    p_value = SKAT_PValue_Logistic_VarMatching(Q, Z2 /sqrt(2), mu, Q_sim)
    return p_value    

# %%
############# Get_Lambda_U_From_Z #############
def Get_Lambda_U_From_Z(Z):
    out_svd = np.linalg.svd(Z)
    lambda_org = pow(out_svd[1],2)
    IDX = where(lambda_org > mean(lambda_org)/100000)
    lambda_org = pow(out_svd[1],2)
    Lambda = lambda_org[IDX]
    u = out_svd[0]
    U = u[:,IDX[0]]
    return Lambda, U

# %%
############# SKAT_Get_Cov_Param #############
def SKAT_Get_Cov_Param(Lambda,p_all,U):
    p_m = len(Lambda)
    m4 = p_all * (1 - p_all) * (3 * pow(p_all, 2) - 3 * p_all + 1) / pow(p_all * (1 - p_all), 2)
    zeta = np.zeros(len(Lambda))
    var_i = np.zeros(len(Lambda))
    varQ = 0
    for i in range(p_m):
        temp_M1 = pow(sum(pow(U[:,i], 2)), 2) - sum(pow(U[:,i], 4))
        zeta[i] = sum(m4 * pow(U[:,i], 4)) + 3 * temp_M1
        var_i[i] = zeta[i] - 1 
    if p_m == 1:
        cov_mat = zeta * pow(Lambda, 2)
    elif p_m > 1:
        cov_mat = np.diag(zeta * pow(Lambda, 2))
        for i in range(p_m - 1):
            for j in range(i+1, p_m):
                cov_mat[i,j] = SKAT_Get_Var_Elements(m4, p_all, U[:,i], U[:,j])
                cov_mat[i,j] = cov_mat[i,j] * Lambda[i] * Lambda[j]
    
    cov_mat = cov_mat + cov_mat.T
    cov_mat2 = cov_mat
    cov_mat2[np.diag_indices_from(cov_mat2)] = np.diag(cov_mat)/2
    varQ = sum(cov_mat2) - pow(sum(Lambda), 2)
    muQ = sum(Lambda)
    Lambda_new = Lambda * sqrt(var_i) / sqrt(2)
    return zeta, var_i, varQ, muQ, Lambda_new, p_m

# %%
############# SKAT_GET_kurtosis #############
def SKAT_GET_kurtosis(x):
    if std(x) == 0:
        return -100
    m4 = mean(pow(x - mean(x), 4))
    kurt = m4 / pow(std(x), 4) - 3
    return kurt

# %%
############# SKAT_GET_skewness #############
def SKAT_GET_skewness(x):
    m3 = mean(pow(x - mean(x), 3))
    skew = m3 / pow(std(x), 3)
    return skew

# %%
############# SKAT_Get_DF_Sim #############
def SKAT_Get_DF_Sim(Q_sim):
    s2_sim = SKAT_GET_kurtosis(Q_sim)
    df_sim = 12 / s2_sim
    if s2_sim <= 0 :
        df_sim = 100000
    elif df_sim < 0.01 :
        s1_sim = SKAT_GET_skewness(Q_sim)
        df_sim = 8 / pow(s1_sim, 2)
    return df_sim

# %%
############# Get_Liu_Params_Mod #############
def Get_Liu_Params_Mod(c1):
    muQ = c1[0]
    sigmaQ = sqrt(2 * c1[1])
    s1 = c1[2] / pow(c1[1], 3/2)
    s2 = c1[3] / pow(c1[1], 2)
    if pow(s1, 2) > s2:
        a = 1 / (s1 - sqrt(pow(s1, 2) - s2))
        d = s1 *a^3 - a^2
        l = a^2 - 2*d
    else: 
        l = 1/s2
        a = sqrt(l)
        d = 0
    muX = l + d
    sigmaX = sqrt(2) * a
    return l, d, muQ, muX, sigmaQ, sigmaX

# %%
############# SKAT_Logistic_VarMatching_GetParam #############
def SKAT_Logistic_VarMatching_GetParam(Lambda, U, p_all, Q_sim):
    re = SKAT_Get_Cov_Param(Lambda, p_all, U)
    Lambda_new = re[4]
    muQ = re[3]
    varQ = re[2]
    df = SKAT_Get_DF_Sim(Q_sim)
    c1 = np.zeros(4)
    for i in range(4):
        c1[i] = sum(pow(Lambda, i + 1))
    para = Get_Liu_Params_Mod(c1)
    n_Lambda = len(Lambda_new)
    return muQ, varQ, df, Lambda_new, para, n_Lambda

# %%
############# SKAT_Logistic_VarMatching_GetParam1 #############
def SKAT_Logistic_VarMatching_GetParam1(Z, p_all, Q_sim):
    out_svd = Get_Lambda_U_From_Z(Z)
    Lambda = out_svd[0]
    U = out_svd[1]
    para = SKAT_Logistic_VarMatching_GetParam(Lambda, U, p_all, Q_sim)
    return para

# %%
############# SKAT_Get_Var_Elements #############
def SKAT_Get_Var_Elements(m4,p_all,u1,u2):
    temp1 = pow(u1, 2) * pow(u2, 2)
    a1 = sum(m4 * temp1)
    a2 = sum(pow(u1, 2))*sum(pow(u2, 2)) - sum(temp1)
    a3 = pow(sum(u1 * u2), 2) - sum(temp1)
    a3 = a3 * 2
    return (a1 + a2 + a3)

# %%
############# SKAT_PValue_Logistic_VarMatching (get p_value and p_value without adjustment) logistic #############
def SKAT_PValue_Logistic_VarMatching(Q, Z, p_all, Q_sim):
    para = SKAT_Logistic_VarMatching_GetParam1(Z, p_all, Q_sim)
    if para[1] == 0:
        p_value = np.ones(1, len(Q))
    param = para[4]
    Q_norm = (Q - para[0]) / sqrt(para[1])
    Q_norm1 = Q_norm * sqrt(2 * para[2]) + para[2]
    p_value =  scipy.stats.distributions.chi2.sf(x = Q_norm1, df = para[2])
    if len(param) > 1:
        if param[4] == 0:
            p_value_noadj = np.ones(1, len(Q))
        else:
            Q_norm = (Q - param[2]) / param[4]
            Q_norm1 = Q_norm * param[5] + param[3]
            p_value_noadj =  scipy.stats.distributions.chi2.sf(x = Q_norm1, df = param[0])                
    return p_value, p_value_noadj

# %%
############# SKAT logistic Model(get Q,W) ################
def SKAT_Logistic_Model(geno, weights, re):
    X = re[0]
    mu = re[1]
    pi_1 = re[2]
    res = re[3]
    Z1 = weights * geno
    Q1 = (res.T).dot(Z1)
    Q = Q1.dot(Q1.T) / 2
    W = (Z1.T).dot((pi_1 * Z1.T).T.values) - (pi_1 * Z1.T).dot(X).dot(np.linalg.inv((X.T).dot((pi_1 * X.T).T.values))).dot((X.T).values.dot((pi_1 * Z1.T).T.values)) 
    return Q,W 

# %%
############# SKAT Linear Model(get Q,W) ################
def SKAT_Linear_Model(geno, weights, re):
    s2 = re[0]
    res = re[1]
    X = re[2]
    Z1 = weights * geno
    Q1 = (res.T).dot(Z1)
    Q = Q1.dot(Q1.T) / s2 / 2
    W = (Z1.T).dot(Z1)  - ((Z1.T).dot(X)).dot(np.linalg.inv((X.T).dot(X))).dot(((X.T).dot(Z1)))     
    return Q,W 

# %%
############# Get_Lambda #############
def Get_Lambda(K):
    lambda1 = np.linalg.eigvalsh(K, UPLO = 'L')
    a = mean(lambda1[lambda1 >= 0] / 100000)
    lambda2 = lambda1[lambda1 > a]
    #len(lambda2)==0, stop?
    return lambda2

# %%
############# Get_Liu_Params_Mod_Lambda (function from Get_Liu_PVal_MOD_Lambda) #############
def Get_Liu_Params_Mod_Lambda(Lambda):
    c1 = np.zeros(shape = (4))
    for i in range (0,4):
        c1[i] = sum(pow(Lambda, i + 1))
    muQ = c1[0]
    sigmaQ = sqrt(2 * c1[1])
    s1 = c1[2] / (pow(c1[1], 3 / 2))
    s2 = c1[3] / (pow(c1[1], 2))
    if pow(s1, 2) > s2:
        a = 1 / (s1 - sqrt(s1^2 - s2))
        d = s1 * pow(a, 3) - pow(a, 2)
        l = pow(a, 2) - 2 * d
    else:
        l = 1 / s2
        a = sqrt(l)
        d = 0
    muX = l + d
    sigmaX = sqrt(2) * a 
    return l, d, muQ, muX, sigmaQ, sigmaX

# %%
############# Get_Liu_PVal_MOD_Lambda (get_pvalue) #############
def Get_Liu_PVal_MOD_Lambda(Q,Lambda):
    l, d, muQ, muX, sigmaQ, sigmaX = Get_Liu_Params_Mod_Lambda(Lambda)
    Q_Norm = (Q - muQ) / sigmaQ
    Q_Norm1 = Q_Norm * sigmaX + muX
    p_value = scipy.stats.distributions.chi2.sf(x = Q_Norm1, df = l)
    return p_value    

# %%
############# SKAT_davies method to get p_value #############
def SKAT_davies(Q, Lambda):
    nc = np.zeros(len(Lambda))
    n = np.ones(len(Lambda))
    r = len(Lambda) # int r1
    sigma = 0 # double sigma
    lim = 10000 # int lim1
    acc = pow(10,-6) # double acc
    trace = np.zeros(7) # double trace
    ifault = 0 # int ifault
    res = 0 # double
    trace, res, ifault = qfc_checked.qfc(Lambda, nc, n, r, sigma, Q, lim, acc, trace, ifault, res)
    return trace, res, ifault

# %%
############# caculate p_value #############
def Get_PValue_Lambda(Lambda,Q):
    p_val_liu = Get_Liu_PVal_MOD_Lambda(Q, Lambda)
    trace, res, ifault = SKAT_davies(Q, Lambda)
    p_val = 1 - res
    if len(Lambda) == 1:
        p_val = p_val_liu
    if p_val > 1 or p_val <= 0:
        p_val = p_val_liu 
    return p_val

# %%
############# combine to get p_value #############
def Get_pvalue(W, Q):
    K = W / 2
    Lambda = Get_Lambda(K)
    p_val = Get_PValue_Lambda(Lambda,Q)
    return p_val

# %%
#############  phenotype, cov, genotype with samples without missing  #############
def missingcheck(geno, pheno, cov):
    regdata = pd.concat([pheno, cov], axis=1)
    regdata["id"] = range(len(regdata))
    regdata_nona = regdata.dropna()
    if len(regdata) != len(regdata_nona):
        keep_col = regdata_nona["id"]
        geno_mat_nona = geno[keep_col,:]
        pheno_nona = regdata_nona[pheno.columns]
        cov_nona = regdata_nona[cov.columns]
        return geno_mat_nona, pheno_nona, cov_nona
    else:
        return geno, pheno, cov

# %%
############# SKAT method main function #############
def SKAT(geno, pheno, cov, weights, outtype):
    geno, pheno, cov = missingcheck(geno, pheno, cov)
    if outtype == "C":
        re, n_resampling = SKAT_Null_Model(pheno, cov, outtype)
        Q,W = SKAT_Linear_Model(geno, weights, re)
        p_val = Get_pvalue(W, Q)
    if outtype =="D":
        re, n_resampling = SKAT_Null_Model(pheno, cov, outtype)
        if n_resampling > 0:
            p_val = SKAT_With_NullModel_ADJ(re, weights, geno)
        else:
            Q,W = SKAT_Logistic_Model(geno, weights, re)
            p_val = Get_pvalue(W, Q)
    return p_val

# %%
############# SKAT method for summary statistics #############
def SKAT_summary(beta_var, beta_estimate, weight, sample_size, COV, D):
    weight_square = pow(weight,2)
    y_estimate = np.zeros(len(weight))
###estimate y'y 
    for i in range(len(weight)):
        y_estimate[i] = (sample_size-1) *D[i]*beta_var[i]*(sample_size-1)+ D[i]*pow(beta_estimate[i],2)*(sample_size-1)
###get median as y'y and variance of y
    y_square_estimate =np.median(y_estimate)
    y_estimate_var =y_square_estimate/(sample_size-1)
###get W
    W = (1/y_estimate_var)*(sample_size-1)*diag(weight).dot(COV).dot(diag(weight))
###get G'y
    res_geno_temp = diag(D*sample_size).dot(beta_estimate)
    res_geno = res_geno_temp/y_estimate_var
    square_res_geno = pow(res_geno,2)
    Q = square_res_geno.dot(weight_square)
    p_val = Get_pvalue(2 * W, Q)
    return p_val

# %%
############# SKAT method for summary statistics #############
def SKAT_summary(beta_var, beta_estimate, weight, sample_size, COV, D):
    weight_square = pow(weight,2)
    y_estimate = np.zeros(len(weight))
###estimate y'y 
    for i in range(len(weight)):
        y_estimate[i] = (sample_size-1) *D[i]*beta_var[i]*(sample_size-1)+ D[i]*pow(beta_estimate[i],2)*(sample_size-1)
###get median as y'y and variance of y
    y_square_estimate =np.median(y_estimate)
    y_estimate_var =y_square_estimate/(sample_size-1)
###get W
    W = (1/y_estimate_var)*(sample_size-1)*diag(weight).dot(COV).dot(diag(weight))
###get G'y
    res_geno_temp = diag(D*sample_size).dot(beta_estimate)
    res_geno = res_geno_temp/y_estimate_var
    square_res_geno = pow(res_geno,2)
    Q = square_res_geno.dot(weight_square)
    p_val = Get_pvalue(2 * W, Q)
    return p_val
# %%
############# SKAT method for summary statistics #############
def SKAT_summary_gety(beta_var,beta_estimate, sample_size, COV, D):
    y_estimate = np.zeros(len(beta_var))
###estimate y'y 
    for i in range(len(beta_var)):
        y_estimate[i] = (sample_size-1) *D[i]*beta_var[i]*(sample_size-1)+ D[i]*pow(beta_estimate[i],2)*(sample_size-1)
###get median as y'y and variance of y
    y_square_estimate =np.median(y_estimate)
    return y_square_estimate


# %%
############# SKAT method for summary statistics #############
def SKAT_summary_withy(y_est,beta_estimate, weight, sample_size, COV, D):
    weight_square = pow(weight,2)
###get median as y'y and variance of y
    y_estimate_var =y_est/(sample_size-1)
###get W
    W = (1/y_estimate_var)*(sample_size-1)*diag(weight).dot(COV).dot(diag(weight))
###get G'y
    res_geno_temp = diag(D*sample_size).dot(beta_estimate)
    res_geno = res_geno_temp/y_estimate_var
    square_res_geno = pow(res_geno,2)
    Q = square_res_geno.dot(weight_square)
    p_val = Get_pvalue(2 * W, Q)
    return p_val
