import pickle
import matplotlib.pyplot as plt
def show(filename):
    f = open(filename,'rb')
    pickle.load(f)
    plt.show()
    f.close()
    return 

if __name__=='__main__':
    show('b_test')
    show('Free_carrier_density')
    