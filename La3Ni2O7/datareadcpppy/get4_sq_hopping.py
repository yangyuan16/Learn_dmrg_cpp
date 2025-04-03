import numpy as np
import pandas as pd
import matplotlib.pylab as plt

if __name__ == "__main__":
    print()
    Lz = 2
    Ly = 3
    Lx = 48
    dop = 0.25
    t = 3
    J = 1
    Jz = 1
    dim = 6000 # dim cutoff
    workpath = "E:\\WORK\\Work\\Project\\La3Ni2O7"
    filepath1 = "\\data_dmrgcpp\\Lz%d_Ly%d_Lx%d\\dop%g" % (Lz, Ly, Lx,dop)
    filepath2 = "\\t%d_J%d_Jz%g_dim%d" % (t, J, Jz, dim)
    filepath3 = "\\measurement_sq_hopping.dat"
    filename = workpath + filepath1 + filepath2 + filepath3
    print(filename)
    # load the data
    #df = pd.read_csv(filename, header=None, sep='\t',encoding='utf-8')
    df = pd.read_csv(filename, header=None)
    #df.rename(columns={0: "site1", 1: "site2", 2: "corre"},inplace=True)
    #df = df[(df["site2"]-df["site1"]) % (Lz * Ly) == 0] 
    #df.sort_values(["site1","site2"],inplace=True)
    #print(df.head())
    #print(df.tail())
    print(len(df))