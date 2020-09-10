#!/usr/bin/env python

####################################################
from math import atan, pow

import pandas as pd
import numpy as np
from numpy import exp, floor, log, sin, sqrt

####################################################
# %%
pi = 3.14159265358979
log28 = 0.0866 

# %%
############# exp #############
def exp1(x):
    x1 = 0 if x <- 50 else exp(x)
    return x1

# %%
############# log #############
def log1(x, first):
    if abs(x) > 0.1:
        s = log(1 + x) if first == True else (log(1 + x) - x)
    else:
        y = x / (2 + x)
        term = 2 * pow(y, 3)
        k = 3
        s = 2 * y if first == True else -x * y
        y = pow(y, 2)
        s1 = s + term / k
        while s1 != s:
            k = k + 2
            term = term * y
            s = s1
            s1 = s + term / k
    return s

# %%
#############errbd#############
#############stop when count>lim#############
def errbd(u, sigsq, lb, nc, n, r):
    #counter count>lim break
    xconst = u * sigsq
    sum1 = u * xconst
    u = 2 * u
    j = r - 1
    while j >= 0:
        x = u * lb[j]
        y = 1 - x
        xconst = xconst + lb[j] * (nc[j] / y + n[j]) / y
        sum1 = sum1 + nc[j] * pow(x / y, 2) + n[j] * (pow(x, 2) / y + log1(-x, False))
        j = j - 1
    cx = xconst
    return_sum = exp1(-0.5*sum1)
    return cx, return_sum

# %%
def ctff(accx, upn, sigsq, lb, nc, n, r, lmin, lmax, mean):
    u2 = upn; u1 = 0; c1 = mean;
    rb = 2* lmax if u2 > 0 else 2*lmin
    #first statement
    u = u2 / (1.0 + u2 * rb)
    #condition
    c2, return_sum = errbd(u, sigsq, lb, nc, n, r)
    while return_sum > accx:
        u1 = u2;  c1 = c2;  u2 = 2.0 * u2
        u = u2 / (1.0 + u2 * rb)
        c2, return_sum = errbd(u, sigsq, lb, nc, n, r)
    ###############
    #first statement
    u = (c1 - mean) / (c2 - mean)
    #condition
    while u < 0.9:
        u = (u1 + u2) / 2.0
        xconst, return_sum = errbd(u / (1.0 + u * rb),sigsq, lb, nc, n, r)
        if return_sum > accx:
            u1 = u; c1 = xconst;
        else:
            u2 = u; c2 = xconst;
        u = (c1 - mean) / (c2 - mean)
    upn = u2; 
    return upn, c2

# %%
#############truncation#############checked
def truncation(u, tausq, sigsq, r, lb, nc, n):
    #counter count>lim break
    sum1 = 0
    prod2 = 0
    prod3 = 0
    s = 0
    sum2 = (sigsq + tausq) * pow(u, 2)
    prod1 = 2 * sum2
    u = 2 * u
    #first statement
    j = 0
    #condition
    while j < r:
        x = pow(u * lb[j], 2)
        sum1 = sum1 + nc[j] * x / (1 + x)
        if x > 1:
            prod2 = prod2 + n[j] * log(x)
            prod3 = prod3 + n[j] * log1(x, True)
            s = s + n[j]
        else:
            prod1 = prod1 + n[j] * log1(x, True)
        j = j + 1
    sum1 = 0.5 * sum1
    prod2 = prod1 + prod2
    prod3 = prod1 + prod3
    x = exp1(-sum1 - 0.25 * prod2) / pi
    y = exp1(-sum1-0.25*prod3) / pi
    err1 = 1 if s == 0 else x * 2 / s
    err2 = 2.5 * y if prod3 > 1 else 1
    if err2 < err1:
        x = 0.5 * sum2
        err2 = 1 if x <= y else y / x
    err = err1 if err1 < err2 else err2
    return err        

# %%
#############findu#############checked
def findu(utx, accx, sigsq, r, lb, nc, n):
    divis = [2, 1.4, 1.2, 1.1]
    ut = utx #why point in c++????
    u = ut / 4
    if truncation(u, 0, sigsq, r, lb, nc, n) > accx:
        #first statement
        u = ut
        #conidtion
        while truncation(u, 0, sigsq, r, lb, nc, n) > accx:
            ut = ut * 4
            u = ut
    else:
        ut = u
        #first statement
        u = u / 4
        #condition
        while truncation(u, 0, sigsq, r, lb, nc, n) <= accx:
            ut = u
            u = u / 4
    #first statement
    i = 0 
    #condition
    while i < 4:
        u = ut / divis[i]
        if truncation(u, 0, sigsq, r, lb, nc, n) <= accx:
            ut = u
        i = i + 1
    utx = ut
    return utx  

# %%
#############integrate#############
def integrate(nterm, interv, tausq, mainx, sigsq, c, r, n, lb, nc, intl, ersm):
    inpi = interv / pi
    #first statement
    k = nterm
    #condition
    while k >= 0:
        u = (k + 0.5) * interv
        sum1 = -2 * u * c
        sum2 = abs(sum1)
        sum3 = -0.5 * sigsq * pow(u, 2)
        #first statement
        j = r-1
        #condition
        while j >= 0:
            x = 2 * lb[j] * u
            y = pow(x, 2)
            sum3 = sum3 - 0.25 * n[j] * log1(y, True)
            y = nc[j] * x / (1 + y)
            z = n[j] * atan(x) + y
            sum1 = sum1 + z
            sum2 = sum2 + abs(z)
            sum3 = sum3 - 0.5 * x * y
            j = j - 1
        x = inpi * exp1(sum3) / u
   
        if mainx == False:
            x = x * (1 - exp1(-0.5 * tausq * pow(u, 2)))
        sum1 = sin(0.5 * sum1) * x
        sum2 = 0.5 * sum2 * x
        intl = intl + sum1
        ersm = ersm + sum2
        k = k - 1
    return intl, ersm

# %%
#############order num absolute value of lb, increasing order#############
def order(lb):
    sortnum = np.argsort(-np.abs(lb))
    return sortnum

# %%
#############cfe#############
def cfe(x, fail,n, lb, nc, r):
    #counter count>lim break
    sortnum= order(lb)
    axl = abs(x)
    sxl= 1 if x > 0 else -1
    sum1 = 0
    j = r - 1
    while j >= 0:
        t = sortnum[j]
        if lb[t] * sxl > 0:
            lj = abs(lb[t])
            axl1 = axl - lj * (n[t] + nc[t])
            axl2 = lj / log28
            if (axl1 > axl2):
                axl = axl1
            else:
                if axl > axl2:
                    axl = axl2
                sum1 = (axl - axl1) / lj
                    #print(sum1)
                k = j - 1
                while k >= 0:
                    sum1 = sum1 + n[sortnum[k]] + nc[sortnum[k]]
                    k = k - 1
                if sum1 > 100:
                    fail = True
                    return_value = 1
                else:
                    return_value = pow(2.0, (sum1 / 4)) / (pi * pow(axl, 2))
                return return_value, fail
        j = j - 1
    if sum1 > 100:
        fail = True
        return_value = 1
    else:
        return_value = pow(2.0, (sum1 / 4)) / (pi * pow(axl, 2))
    return return_value, fail                 

# %%
def qfc(lb, nc, n, r, sigma, c, lim, acc, trace, ifault, res):
    acc1 = acc
    qfval = -1.0
    rats = [1,2,4,8]
    ifault = 0; count = 0; intl = 0; ersm = 0
    fail = False
    sigsq = pow(sigma, 2)
    sd = sigsq
    xlim = lim
    lmax = 0; lmin = 0; mean = 0
    #first statement
    j = 0
    #condition
    while j < r:
        if n[j] < 0 or nc[j] < 0:
            trace[6] = count
            res = qfval
            ifault = 3
            return trace, res, ifault
        sd = sd + pow(lb[j], 2)*(2*n[j] + 4*nc[j])
        mean = mean + lb[j]*(n[j] + nc[j])
        if lmax < lb[j]:
            lmax = lb[j]
        elif lmin > lb[j]:
                lmin = lb[j]
        j = j + 1
    #######################
    if sd == 0:
        qfval = 1 if c > 0 else 0
        trace[6] = count
        res = qfval
        return trace, res, ifault
    if lmin == 0 and lmax == 0 and sigma == 0:
        trace[6] = count
        res = qfval
        ifault = 3
        return trace, res, ifault
    sd = sqrt(sd)
    almx = -lmin if lmax < -lmin else lmax
    utx1 = 16 / sd; up = 4.5 / sd; un = -up
    accx = 0.5 * acc1
    utx = findu(utx1, accx,sigsq, r, lb, nc, n)
    if c != 0 and almx > 0.07 * sd:
        return_value, fail = cfe(c, fail, n, lb, nc, r)
        tausq = 0.25 * acc / return_value
        if fail:
            fail = False
        err = truncation(utx, tausq, sigsq, r, lb, nc, n)
        if err < 0.2 * acc1:
            sigsq = sigsq + tausq
            utx=findu(utx, 0.25 * acc1, sigsq, r, lb, nc, n)
            trace[5] = sqrt(tausq)
    trace[4] = utx;  acc1 = 0.5 * acc1
    #######################
    flag = True
    while flag == True:
        up, c2 = ctff(acc1, up, sigsq, lb, nc, n, r, lmin, lmax, mean)
        d1 = c2 - c
        if d1 < 0:
            qfval = 1
            trace[6] = count
            res = qfval
            return trace, res, ifault
        un, c3 = ctff(acc1, un, sigsq, lb, nc, n, r, lmin, lmax, mean)
        d2 = c - c3
        if d2 < 0:
            qfval = 1
            trace[6] = count
            res = qfval
            return trace, res, ifault
        intv = 2 * pi / d1 if d1 > d2 else 2 * pi / d2
        xnt = utx / intv
        xntm = 3 / sqrt(acc1)
        if xnt > xntm * 1.5:
            if xntm > xlim:
                ifault = 1
                qfval = 1
                trace[6] = count
                res = qfval
                return trace, res, ifault
            ntm = floor(xntm + 0.5)
            intv1= utx / ntm; x = 2.0 * pi / intv1
            if x <= abs(c):
                flag=False
            else:
                x1, fail = cfe(c-x, fail, n, lb, nc, r)
                x2, fail = cfe(c+x, fail, n, lb, nc, r)
                tausq = 0.33 *acc1 / (1.1 * (x1 + x2))
            if fail:
                flag=False
            else:
                acc1 = 0.67 * acc1
                intl, ersm = integrate(ntm, intv1, tausq, False, sigsq, c, r, n, lb, nc, intl, ersm)
                xlim = xlim - xntm; sigsq = sigsq + tausq
                trace[2] = trace[2] + 1; trace[1] = trace[1] + ntm + 1
                utx = findu(utx, 0.25 * acc1, sigsq, r, lb, nc, n)
                acc1 = 0.75*acc1
                flag = True
        else:
            flag = False      
    #######################
    trace[3] = intv
    if xnt > xlim:
        ifault = 1
        qfval = 1
        trace[6] = count
        res = qfval
        return trace, res, ifault
    nt = floor(xnt + 0.5)
    intl, ersm = integrate(nt, intv, 0, True, sigsq, c, r, n, lb, nc, intl, ersm)
    trace[2] = trace[2] + 1; trace[1] = trace[1] + nt + 1
    qfval = 0.5 - intl
    trace[0] = ersm
    up = ersm; x = up + acc1 / 10.0
    j = 0 
    while j < 4:
        if rats[j] * x == rats[j] * up:
            ifault = 2
        j = j + 1
    #######################
    trace[6] = count
    res = qfval
    return trace, res, ifault   




