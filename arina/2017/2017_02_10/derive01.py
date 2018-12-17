def derive(X,Y,step=40,cutoff=200, cutoff_s=10):
    for j in range(cutoff_s+step, len(Y)-cutoff):
        yield ((Y[j] -Y[j-step])/(X[j] -X[j-step]), X[j])
