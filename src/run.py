import time
import sys
import numpy as np

import Solid
import ForwardModel_Box
import OptModel_Box


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

    #objective function
    g_obj = lambda qvals: np.maximum(qvals)

    stpPath = "test.stp"
    qDirIn = [0, 1, 0] #[m]
    qMagIn = 10 #[W/m^2]

    box = Solid(stpPath)
    fwd = ForwardModel_Box(g_obj, box, qMagIn, qDirIn) 
    opt = OptModel_Box()

    t0 = time.time()

    box.loadSTEP()

    # while (del_e > threshold_err): 
    #     q = fwd.calculateQ()
    #     g = fwd.calculateG()
    #     del_e = opt.calculateDelE()

    #     if del_e > threshold_err: 
    #         break
    #     else: 
    #         box = box.manipulate() #modify it, should probably specify how eventually

    while (opt.del_e > opt.threshold_err):
        fwd.processCADModel()
        q_mesh_all = fwd.calcQMesh()
        g_now = fwd.calcObjective(q_mesh_all)
        opt.updategValues(g_now)
        del_e = opt.calculateDelE()

        if (del_e > opt.threshold_err):
            print(f"found opt, time elapsed: {time.time() - t0}, err: {opt.del_e}")
            break
        else: 
            del_theta = opt.doTransform(box) #this prob doesn't match yet, gonna fix
            print(f"transformed: [err]: {opt.del_e} [g_x now]: {opt.g_curr} [g_x prev]: {opt.g_prev} [theta]: {del_theta}")





        


    
