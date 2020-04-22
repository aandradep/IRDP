import numpy as np
import scipy.stats as si

class Option:
    
    def __init__(self,
                 strike: float,
                 term: float,
                 option_type: str = 'call',
                 option_style: str = 'european'
                 ):
        self._strike = strike
        self._term = term
        self._type = option_type
        self._style = option_style
    
    
    def price(self, spot: float, risk_free: float, volatility: float,
              price_method: str = 'bsm', steps: int = None)  -> float:
       
        if (price_method == 'binomial') and (self._style == 'european' or self._type == 'call'):
            dt = self._term/steps
            A = 0.5*(np.exp(-risk_free*dt)+np.exp((risk_free+volatility**2)*dt))
            Up = A + np.sqrt(A**2-1)
            Down = 1 / Up
            price = self.EuropeanBinomialTree(spot=spot,steps=steps,risk_free=risk_free,
                                              Up=Up,Down=Down)
            
        elif(price_method == 'bsm') and (self._style == 'european'):
            price = self.BlackScholesMerton(spot=spot, risk_free=risk_free,
                                            volatility=volatility)
            
        else:
            dt = self._term/steps
            A = 0.5*(np.exp(-risk_free*dt)+np.exp((risk_free+volatility**2)*dt))
            Up = A + np.sqrt(A**2-1)
            Down = 1 / Up
            price = self.AmericanBinomialTree(spot=spot,steps=steps,risk_free=risk_free,
                                              Up=Up,Down=Down)
        
        return price
            
    def EuropeanBinomialTree(self, spot: float, steps: int, risk_free: float, 
                             Up: float, Down: float) -> float:
        """Pricing of an european option, call or put based on binomial tree.
        # Arguments:
            spot {float} -- spot of the underlying 
            steps {int} -- number of steps to simulate the price, e.j 12 for monthly
            risk_free {float} -- risk free rate in number, e.j for 5% then 0.05
            Up {float} -- 1 + percentage the stock goes up when it increases in price,
                          e.j if it increases 5% each time it goes up then 1.05
            Down {float} -- 1 - percentage the stock goes down when its price decreases, 
                            e.j if it decreases 10% each time it goes down then 0.9
        
        # Return: the result is a float with the price of the option
                            
        # The pricing is done following the code 8 describe in:
          Nine Ways to Implement the Binomial Method for Option Valuation in MATLAB
        """
        
        strike = self._strike
        option_type = self._style
        dt = self._term/steps
        p = (np.exp(risk_free*dt)-Down)/(Up-Down)
        q = 1-p
        # steps + 1 because at the end we have steps + 1 prices
        W = spot * (Down **np.arange(steps,-1,-1)) * (Up **np.arange(0,steps+1))
        # Options values at  T 
        if option_type == 'call':
            W = np.maximum(W - strike, 0)
        else:
            W = np.maximum(strike - W, 0)
            
        # price based on log pascal triangles
        lsc = np.cumsum(np.log(np.append(1,np.arange(1,steps+1))))
        temp = (lsc[steps] - lsc - lsc[::-1] + np.log(p) * \
                np.arange(0,steps+1) + np.log(q)*np.arange(steps,-1,-1))
        price = np.exp(-risk_free*self._term) * sum(np.exp(temp)*W)
                
        return price
    
    
    def AmericanBinomialTree(self, spot: float, steps: int, risk_free: float,
                             Up: float, Down: float):
        """Pricing of an american option, call or put based on binomial tree.
        # Arguments:
             spot {float} -- spot of the underlying 
             steps {int} -- number of steps to simulate the price, e.j 12 for monthly
             risk_free {float} -- risk free rate in number, e.j for 5% then 0.05
             Up {float} -- 1 + percentage the stock goes up when it increases in price,
                           e.j if it increases 5% each time it goes up then 1.05
            Down {float} -- 1 - percentage the stock goes down when its price decreases, 
                            e.j if it decreases 10% each time it goes down then 0.9
        
        # Return: the result is a float with the price of the option
        # The pricing is done following the code 5 describe in 
        Nine Ways to Implement the Binomial Method for Option Valuation in MATLAB
        """
        
        strike = self._strike
        option_type = self._style
        dt = self._term/steps
        p = (np.exp(risk_free*dt)-Down)/(Up-Down)
        q = 1-p
        
        dpowers = Down ** np.arange(steps,-1,-1)
        upowers = Up ** np.arange(0,steps+1)
        scale1 = p * np.exp(-risk_free*dt)
        scale2 = q * np.exp(-risk_free*dt)
        
        # steps + 1 because at the end we have steps + 1 prices
        W = spot*dpowers*upowers
        # Options values at  T 
        if option_type == 'call':
            W = np.maximum(W - strike, 0)
        else:
            W = np.maximum(strike - W, 0)
            
        # backward valuation
        for i in np.arange(steps, 0,-1):
            Si = spot*dpowers[(steps-i+1):steps+1]*upowers[0:i]
            if option_type == 'call':
                W = np.maximum(np.maximum(Si - strike ,0), scale1*W[1:(i+1)] + scale2*W[1:i])
            else:
                W = np.maximum(np.maximum(strike - Si, 0), scale1*W[1:(i+1)] + scale2*W[0:i])
        
        return W[0]
    
    def BlackScholesMerton(self, spot: float, risk_free: float,
                           volatility: float) -> float:
        """Pricing of an european option, call or put based on BSM formula.
        # Arguments:
            spot {float} -- spot of the underlying 
            risk_free {float} -- risk free rate in number, e.j for 5% then 0.05
        
        # Return: the result is a float with the price of the option
                            
        # The pricing is done following the code 8 describe in:
          Nine Ways to Implement the Binomial Method for Option Valuation in MATLAB
        """
    
        d1 = (np.log(spot/self._strike) + (risk_free + 0.5 * volatility**2)/self._term)/(volatility * np.sqrt(self._term))
        d2 = d1 - volatility * np.sqrt(self._term)
        
        if self._style == 'call':
            price = spot * si.norm.cdf(d1, 0, 1) - self._strike * np.exp(-risk_free * self._term) * si.norm.cdf(d2,0,1)
        else:
            price = self._strike * np.exp(-risk_free * self._term) * si.norm.cdf(-d2,0,1) - spot * si.norm.cdf(-d1, 0, 1)
        
        return price
        

put_option = Option(strike = 10, term = 1, option_type = 'put', option_style = 'european')
european_price = put_option.price(spot = 5, risk_free = 0.06, volatility = 0.3,
                                  price_method = 'binomial', steps = 256)
bs_price = put_option.price(spot = 5 , risk_free = 0.06, volatility = 0.3, 
                            price_method = 'bsm')

put_american = Option(strike = 10, term = 1, option_type = 'put', option_style = 'american')
put_american.price(spot = 9, risk_free = 0.06, volatility = 0.3, price_method = 'binomial', 
                 steps = 256)





