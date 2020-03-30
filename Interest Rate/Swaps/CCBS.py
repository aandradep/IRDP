import numpy as np
import pandas as pd
import datetime

import sys
sys.path.append('D:/Documents/IRDP/Interest Rate')

from Bonds.bonds import InterpolationMethods

term = np.arange(1,11)
coupons = np.array([0.05, 0.051, 0.052, 0.053, 0.054, 0.055, 0.056, 0.057, 
                    0.058, 0.059])
basis = np.array([-0.001, -0.0012, -0.0014, -0.0016, -0.0018, -0.002, -0.0022,
                  -0.0024, -0.0026, -0.0028])


def DFBasicSwap(term, coupons):
    """Calculates the discount factors of a basic swap based on term and coupons
    # Arguments:
        term {np.array} -- Numpy array with the terms of the swaps in years. 
        coupons {np.array} -- Numpy array with the coupons of the swaps in number.    
    # Return: The result is a Numpy array with the discount factors
    # Description: The function asums that the swap pays yearly and the day
    count is 30/360 so the payments are done each 1 year. This simplifies greatly the 
    calculation of the numerator because it does not have to be multiply by its 
    appropiet df.
    """
    
    number_instruments = len(term)
    discount_factors = np.zeros(number_instruments)
    years_btw_terms = np.diff(np.append([0],term))
    
    for instrument in range(number_instruments):

        if instrument == 0:
            numerator = 1
        else:
            numerator = 1 - coupons[instrument]*sum(discount_factors[0:instrument])
            
        denominator = 1+coupons[instrument]*years_btw_terms[instrument]
        discount_factors[instrument] = numerator/denominator
    
    return discount_factors

def DFBasisSwapWrong(term, coupons, basis):
    """Calculates the discount factors of a Basis swap based as PRACTITIONERS DO (WRONG !!!)
    # Arguments:
        term {np.array} -- Numpy array with the terms of the swaps in years. 
        coupons {np.array} -- Numpy array with the coupons of the swaps in number.   
        basis {np.array} -- Numpy array with the basis of the swaps in number.
    # Return: The result is a Numpy array with the discount factors
    # Description: The function asums that the swap pays yearly and the day
    count is 30/360 so the payments are done each 1 year. This simplifies greatly the 
    calculation of the numerator because it does not have to be multiply by its 
    appropiet df.
    """
    number_instruments = len(term)
    discount_factors = np.zeros(number_instruments)
    years_btw_terms = np.diff(np.append([0],term))
    basic_df = DFBasicSwap(term, coupons)
    fwd_rate = np.append([1],basic_df)[:-1]/basic_df - 1
    
    for instrument in range(number_instruments):
        
        if instrument == 0:
            numerator = 1
        else:
            sumatory = sum([(fwd_rate[j]+basis[instrument])*discount_factors[j] 
            for j in range(term[instrument]-1)])
            numerator = 1-sumatory
        
        denominator = 1+(fwd_rate[instrument]+basis[instrument])*years_btw_terms[instrument]
        discount_factors[instrument] = numerator/denominator
        
    return discount_factors

def DFBasisSwap(term, coupons, basis):
    """Calculates the discount factors of a Basis swap based as HfB paper
    # Arguments:
        term {np.array} -- Numpy array with the terms of the swaps in years. 
        coupons {np.array} -- Numpy array with the coupons of the swaps in number.   
        basis {np.array} -- Numpy array with the basis of the swaps in number.
    # Return: The result is a Numpy array with the discount factors 
              and an other with forward rates.
    # Description: The function asums that the swap pays yearly and the day
    count is 30/360 so the payments are done each 1 year. This simplifies greatly the 
    calculation of the numerator because it does not have to be multiply by its 
    appropiet df.
    """
    
    number_instruments = len(term)
    basic_df = np.zeros(number_instruments)
    years_btw_terms = np.diff(np.append([0],term))
    
    for i in range(number_instruments):
        
        if i == 0:
            numerator = 1
        else:
            numerator = 1-(coupons[i]+basis[i])*sum(basic_df[0:i])*years_btw_terms[i]
        
        denominator = 1+(coupons[i]+basis[i])*years_btw_terms[i]
        basic_df[i] = numerator/denominator
    
    basis_df = np.cumsum(basic_df*years_btw_terms)*basis+basic_df
    F_t = basic_df - basis_df
    df_part = (np.append([1],basic_df[:-1])/basic_df-1)/years_btw_terms
    margin_func_part = (F_t-np.append([0],F_t[:-1]))/(years_btw_terms*basic_df)
    fwd_rates = df_part + margin_func_part
    
    return basis_df, fwd_rates    
    
    