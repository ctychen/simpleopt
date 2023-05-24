class OptModel_Box:
    def __init__(self, threshold = 0.5, gprev = 10000, gnow = 1, ):
        #all placeholder values for now I have no idea what they should be yet - maybe these shouldn't be here and instead in opt? 
        self.threshold_err = 0.5
        self.g_prev = 10000
        self.g_curr = 1
        self.delstep = 0.1
        self.del_e = 1
        return
    
    def calculateDelE(self): 
        self.del_e = self.g_prev - self.g_curr
        return self.del_e
    
    def findTransform(self):
        theta_new = theta - delstep*del_e

    

