#def fivepoint(x,y):
def derive(x,y):
    lenx = len(x)
    for _j in range(lenx):
        if _j<2:
            yield (-25.*y[_j] +48.*y[_j+1] -36.*y[_j+2] +16.*y[_j+3] -3.*y[_j+4])/((x[_j+4]-x[_j])*3.)
        elif _j>(lenx-3):
            yield (25.*y[_j] -48.*y[_j-1] +36.*y[_j-2] -16.*y[_j-3] +3.*y[_j-4])/((x[_j]-x[_j-4])*3.)
        else:
            yield (y[_j-2] -8.*y[_j-1] +8.*y[_j+1] -y[_j+2])/((x[_j+2]-x[_j-2])*3.)