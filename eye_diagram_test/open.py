import pickle
import matplotlib.pyplot as plt
import os
def show(filename):
    f = open(filename,'rb')
    pickle.load(f)
    plt.show()
    f.close()
    return 

if __name__=='__main__':
    # os.chdir("D:\Master_degree\paper\微分方程\couple-mode-theory-in-time\eye_diagram_test")
    show('voltage')
    