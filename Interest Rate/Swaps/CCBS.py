import numpy as np
import pandas as pd

import sys
sys.path.append('D:/Documents/IRDP/Interest Rate')
from Bonds.bonds import InterpolationMethods

class BasisSwapPractitioners(InterpolationMethods):
    """Calibrates a basis swap curve as practitioners do."""
    
    def __init__(self,
                 terms: np.array,
                 coupons: np.array,
                 periodicities: np.array,
                 basis: np.array):
        self._terms = terms
        self._coupons = coupons
        self._periodicities = periodicities
        self._basis = basis
        
    
    def BasisSwapCurve(self) -> pd.DataFrame:
        """Calculates the discount factors of a Basis swap as PRACTITIONERS DO (WRONG !!!)
        # Arguments:
            term {np.array} -- Numpy array with the terms of the swaps in years. 
            coupons {np.array} -- Numpy array with the coupons of the swaps in number. 
            periodicities {np.array} -- Numpy array with the coupons periodicities. 1 for year, 2 for semesets 4 for quarters
            basis {np.array} -- Numpy array with the basis of the swaps in number.
        # Return: The result is a pandas dataframe with the calibrated curve
        # Description: The function asums that the day count is 30/360.
        """
        
        number_instruments = len(self._terms)
        basic_df = np.zeros(number_instruments)
        ccs_curve = np.array(pd.DataFrame({'term': self._terms, 
                                           'rate':np.mean(self._coupons)}))
        fwd_rates = np.zeros(number_instruments)
        discount_factors = np.zeros(number_instruments)
        
        for j in range(30):
            for i in range(number_instruments):
                temp_rate = self._SwapRate(term=self._terms[i], coupon=self._coupons[i],
                                          periodicity=self._periodicities[i], 
                                          ccs_curve=ccs_curve)
                ccs_curve[i,1] = temp_rate
                basic_df[i] = np.exp(-temp_rate*self._terms[i])
                delta = 1/self._periodicities[i]
                
                if i == 0:
                    fwd_rates[i] = (1/basic_df[i]-1)/delta
                    numerator = 1
                else:
                    fwd_rates[i] = (basic_df[i-1]/basic_df[i]-1)/delta
                    numerator = 1- sum((fwd_rates[0:i]+self._basis[i])*discount_factors[0:i])
                    
                denominator = 1+(fwd_rates[i]+self._basis[i])*delta
                discount_factors[i] = numerator/denominator
        
        rates = -np.log(discount_factors)/self._terms
        basis_curve = pd.DataFrame({'term': self._terms, 'rate': rates,
                                    'basis': self._basis, 'forward_rates': fwd_rates})
        
        return basis_curve
            
    @staticmethod
    def _SwapRate(term: float, coupon: float, periodicity: int, 
                  ccs_curve: np.array) -> float:
        """Calculates the rate of a basic swap based on term, coupon and coupon periodicity
        # Arguments:
            term {float} -- term in years of the swap, It can be 0.5 for 6 months for example 
            coupons {float} -- coupon of the swap for a swap with price 1, e.j 5 is 0.05. 
            periodicity {int} -- coupon periodicity. Usually 1 for yearly, 2 for semester 4 for quarterly
            ccs_curve {np.array} -- numpy array with a curve
        # Return: The result is float with the rate of the swap
        # Description: The function asums that the day count is 30/360 .
        """
        
        coupons_vec = np.repeat(coupon,term*periodicity)
        deltas = np.repeat(1/periodicity,len(coupons_vec))
        rates_curve = [InterpolationMethods(ccs_curve[:,0], ccs_curve[:,1], delta).RawInterpolation()
                       for delta in np.cumsum(deltas)]
        disc_factors = np.exp(-np.cumsum(deltas)*rates_curve)
        numerator = 1 - coupons_vec[::-1][0]*sum(disc_factors[:-1]*deltas[:-1])
        denominator = 1+coupon*deltas[::-1][0]
        rate = -np.log(numerator/denominator)/term
        
        return rate
        
class BasisSwap(InterpolationMethods):
    """Calibrates a basis swap curve."""
    
    def __init__(self,
                 terms: np.array,
                 coupons: np.array,
                 periodicities: np.array,
                 basis: np.array):
        self._terms = terms
        self._coupons = coupons
        self._periodicities = periodicities
        self._basis = basis
    
    def BasisSwapCurve(self) -> pd.DataFrame:
    
        number_instruments = len(self._terms)
        basic_df = np.zeros(number_instruments)
        basis_df = np.zeros(number_instruments)
        deltas = np.zeros(number_instruments)
        F_t = np.zeros(number_instruments)
        fwd_rates = np.zeros(number_instruments)
        basis_curve = np.array(pd.DataFrame({'term': self._terms, 
                                             'rate': np.mean(self._coupons+self._basis)}))
        
        for j in range(30):
            for i in range(number_instruments):
                temp_rate = self._SwapRate(term=self._terms[i], coupon=self._coupons[i], 
                                          periodicity=self._periodicities[i], basis=self._basis[i],
                                          basis_curve=basis_curve)
                basis_curve[i,1] = temp_rate
                basic_df[i] = np.exp(-temp_rate*self._terms[i])
                deltas[i] = 1/self._periodicities[i]
                
                if i == 0:
                    numerator = 1
                else:
                    numerator = 1-(self._coupons[i]+self._basis[i])*deltas[i]*sum(basic_df[0:i])
            
                denominator = 1+(self._coupons[i]+self._basis[i])*deltas[i]
                basic_df[i] = numerator/denominator
                basis_df[i] = basic_df[i] + self._basis[i]*sum(deltas[0:(i+1)]*basic_df[0:(i+1)])
                F_t[i] = basic_df[i] - basis_df[i]
                
                if i == 0:
                    basic_fwd = (1/basic_df[i]-1)/deltas[i]
                    margin_func = (F_t[i] - 0)/(deltas[i]*basic_df[i])
                else:
                    basic_fwd = (basic_df[i-1]/basic_df[i]-1)/deltas[i]
                    margin_func = (F_t[i] - F_t[i-1])/(deltas[i]*basic_df[i])
                
                fwd_rates[i] = basic_fwd + margin_func
        
        rates = -np.log(basis_df)/self._terms
        basis_curve = pd.DataFrame({'term': self._terms, 'rate': rates, 'basis': self._basis, 
                                    'forward_rates': fwd_rates})
        
        return basis_curve
        
    @staticmethod
    def _SwapRate(term: float, coupon: float, periodicity: int, basis: float ,
                  basis_curve: np.array) -> float:
        """Calculates the rate of a basic swap based on term, coupon and coupon periodicity
        # Arguments:
            term {float} -- term in years of the swap, It can be 0.5 for 6 months for example 
            coupons {float} -- coupon of the swap for a swap with price 1, e.j 5 is 0.05. 
            periodicity {int} -- coupon periodicity. Usually 1 for yearly, 2 for semester 4 for quarterly
            basis {float} -- basis of the swap.
            basis_curve {np.array} -- numpy array with a curve
        # Return: The result is float with the rate of the swap
        # Description: The function asums that the day count is 30/360 .
        """
        coupons_vec = np.repeat(coupon, term*periodicity)
        deltas = np.repeat(1/periodicity,len(coupons_vec))
        rates_curve = [InterpolationMethods(basis_curve[:,0], basis_curve[:,1], delta).RawInterpolation()
                       for delta in np.cumsum(deltas)]
        disc_factors = np.exp(-np.cumsum(deltas)*rates_curve)
        numerator = 1 - sum((coupons_vec[:-1]*deltas[:-1]+basis*1/periodicity)*disc_factors[:-1])
        denominator = 1+coupon*deltas[::-1][0]+basis*1/periodicity
        rate = -np.log(numerator/denominator)/term
        
        return(rate)    
   