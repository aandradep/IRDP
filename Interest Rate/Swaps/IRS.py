import numpy as np
import pandas as pd
import datetime

import sys
sys.path.append('D:/Documents/IRDP/Interest Rate')

from Bonds.bonds import InterpolationMethods

class IRSCurveCalibration:
    def __init__(self,
                 date_val: datetime.date = datetime.date(2002, 12, 31),
                 swaps_rates: list = [],
                 swaps_term: list = [],
                 curve_term: list = []
                 ):
        self.date_val = date_val
        self.swaps_rates = swaps_rates
        self.swaps_term = swaps_term
        self.curve_term = curve_term
        self._curve_rate = (np.repeat(np.mean(swaps_rates),
                                      len(swaps_rates))/100).tolist()
    def _IRSRate(self, term, market_rate) -> float:
        """Calculates the rate of an interest rate bullet swap.
        # Arguments:
            term {flaot} -- Term of  the swap in years, ej 6 months is 0.5.
            market_rate {float} -- interest rate that the swap is trading, 
            a.k.a cupon of the swap.
        # Return: the result is the rate used to calculate DF of the swap.
        """
        rate = -np.log(100/(100+market_rate*term))/term
        
        return rate

    def _OISRate(self, term, market_rate, curve_term, curve_rates) -> float:
        """Calculates the rate of a OIS swap that pays every 3 months
        # Arguments:
            term {float} -- Term of the swap in years.
            market_rate {float} -- interest rate that the swap is trading, 
            a.k.a cupon of the swap.
            date_val {datetime.date} -- date of pricing
            curve_term {list} -- Array with the terms of the curve
            curve_rates {list} -- Rates for each termn of the curve. 
        # Return: the result is the rate used to calculate DF of the swap.
        """
        
        date_end =  self.date_val + pd.DateOffset(months = term*12)
        coupon_dates = pd.date_range(self.date_val, date_end, freq='3M')
        coupon_dates = [date.replace(day = self.date_val.day) for date in coupon_dates]
        coupon_dates.append(pd.Timestamp(date_end))
        term_btw_coupons = (pd.Series(
                np.diff(coupon_dates)).values.view('<i8')/10**9 / 
            (24 * 60 * 60)) / 360        
        coupon_dates.pop(0)
        coupon_terms = ((pd.Series(coupon_dates)-
                         pd.Timestamp(self.date_val)).values.view('<i8')/10**9 / 
        (24 * 60 * 60)) /360  
        N = len(coupon_dates)
        discount_factors = [np.exp(-term*InterpolationMethods(
                curve_term, curve_rates, float(term)).RawInterpolation()) 
        for term in coupon_terms]
        numerator = 1 - np.dot(
                market_rate/100 *  term_btw_coupons[0:(N-1)],
                np.array(discount_factors)[0:(N-1)])
        denominator = 1 + market_rate/100 * term_btw_coupons[N-1]
        rate = -np.log(numerator/denominator)/term
        
        return rate
        
    def CurveCalibration(self) -> np.array:
        """Calibrates the term structure of interest rate swaps
        # Arguments:
            swaps_rates {list} -- List with the rates of the swaps from the market, a.k.a their coupons
            swaps_term {list} -- List with the term of the swaps in years
            date_val {datetime.date} -- Date of pricing
            curve_rates {list} -- List with interest rate of the curve
            curve_term {list} -- List of the terms of the 
        # Return: The result is a np.array with the terms and the rates of the curve
        calibrated.
        """
        
        number_instruments = len(self.swaps_rates)
        new_curve_term = self.curve_term.copy()
        new_curve_rate = self._curve_rate.copy()
        
        for j in range(30):
            # Loops 30 time to be sure but in the paper it says that in 5 it converges
            for instrument in range(number_instruments):
                if self.swaps_term[instrument] < 2:
                    new_curve_rate[instrument] = self._IRSRate(
                            self.swaps_term[instrument], self.swaps_rates[instrument])
                else:
                    new_curve_rate[instrument] = self._OISRate(
                            self.swaps_term[instrument], self.swaps_rates[instrument], 
                            new_curve_term, new_curve_rate)
                                   
        curve = np.c_[new_curve_term, new_curve_rate]
        
        return curve
    
    
class IRSwap(IRSCurveCalibration):
    def __init__(self,
                 curve: np.array,
                 market_rate: float = 0.05,
                 term: float = 5,
                 date_val: datetime.date = datetime.date(2019, 3, 1),
                 date_end: datetime.date = datetime.date(2024, 3, 21)
                 ):
        self._market_rate = market_rate
        self._term = term
        self._date_val = date_val
        self._date_end = date_end
        self._curve = curve

    def IRSPricing(self):
        """Price a interest rate swap with a curve.
        # Arguments:
            curve {np.array} -- Swap curve calibraded with a RawInterpolation method
            market_rate {float} -- interest rate that the swap is trading, 
            a.k.a cupon of the swap so it is not a percentage
            term {float} -- Term of  the swap in years, ej 6 months is 0.5.
            date_val {datetime.date} -- Date of pricing. 
            date_end {datetime.date} -- Date when the swap expirces. 
        # Return: The price of the swap.
        """   
        
        if self._term > 2:
            # If it is OIS 
            coupon_dates = pd.date_range(self._date_val, self._date_end, freq='3M')
            coupon_dates = [date.replace(day = self._date_end.day) for date in coupon_dates]
            coupon_dates.append(pd.Timestamp(self._date_end))
            coupon_dates[0] = pd.Timestamp(self._date_val)
            term_btw_coupons = (pd.Series(
                    np.diff(coupon_dates)).values.view('<i8')/10**9 / 
                        (24 * 60 * 60)) / 360        
            coupon_dates.pop(0)
        else:
            # If it is a bullet swap
            coupon_dates = self._date_end
            term_btw_coupons = self._term
            
        term_df = ((pd.Series(coupon_dates)-
                    pd.Timestamp(self._date_val)).values.view('<i8')/10**9 / 
                        (24 * 60 * 60)) /360  
        discount_factors = [np.exp(-self._term*InterpolationMethods(
                self._curve[:,0], self._curve[:,1], float(self._term)).RawInterpolation()) 
            for term in term_df]
        price = (np.dot(self._market_rate/100 * term_btw_coupons,np.array(discount_factors)) + 
                 discount_factors[len(discount_factors)-1])
        
        return price
   