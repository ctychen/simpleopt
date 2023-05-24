import time
import sys

import Solid
import ForwardModel_Box
import OptModel_Box

box = Solid()
fwd = ForwardModel_Box()
opt = OptModel_Box()

if __name__ == '__main__':

    #set up scenario
    #read config and get initial values/models
    #run fwd for first time (or maybe opt handles this? or maybe all part of loop?)
    #run opt loop after that 

    #rev:
    #mesh the solid
    #calculate q on solid 
    #calculate g with mesh q's
    #calculate dele between this and prev g
    #check if dele > threshold
    #if not, manipulate the solid and repeat
    #log time for vibes

    t0 = time.time()

    while (del_e > threshold_err): 
        q = fwd.calculateQ()
        g = fwd.calculateG()
        del_e = opt.calculateDelE()

        if del_e > threshold_err: 
            break
        else: 
            box = box.manipulate() #modify it, should probably specify how eventually

        


    
