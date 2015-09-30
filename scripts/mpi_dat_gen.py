import numpy as np

data = np.array([0,0,0,0],dtype=float)
text_file = open('mpi_dat.dat', 'w')
for i in np.arange(1,9):
    counter = 0
    for j in np.array([10,50,100,500]):
        data[counter] = np.loadtxt(str(i)+'_core'+str(j)+'.dat')
        counter += 1
    text_file.write(str(i)+' '+str(data[0])+' '+str(data[1])+' '+str(data[2])+' '+str(data[3])+'\n')
#        for item in np.vstack(data[counter]):
#            print(item)
text_file.close()