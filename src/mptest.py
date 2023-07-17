from multiprocessing import Process
import numpy as np

def printvals(i):
    print(f"value: {i}") 
    return    

if __name__ == '__main__':
    #test using multiprocessing
    valuesToPrint = np.arange(10) 

    processes = []
    for i in valuesToPrint:
        p = Process(target=printvals, args=(i,))
        processes.append(p)
        p.start()   

    # p = Process(target=printvals, args=(valuesToPrint,))
    # p.start()
    # p.join()
    
