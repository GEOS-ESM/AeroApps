from numpy import load, mean

if __name__ == "__main__":

    sigF = 0.45  # constant sigF assumed in PSAS' solve.x

    for inst in ( 'modo', 'mydo', 'misr' ):
        sigO = ()
        for ch in ( 870, 660, 550, 470 ):
            f = load('a0005.%s-%d.200806.npz'%(inst,ch))
            sigO_ = mean(f['Alpha'],0)[0] 
            sigF_ = mean(f['Alpha'],0)[1] 
            sigO = sigO + (round(sigO_*sigF/sigF_,2),)
        print("%s %6.2f %6.2f %6.2f %6.2f"%(inst,sigO[0],sigO[1],sigO[2],sigO[3]))

