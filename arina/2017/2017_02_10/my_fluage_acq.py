# -*- coding: utf-8 -*-

from os import getenv, path, listdir
import numpy
from math import factorial, log10
import matplotlib
from cStringIO import StringIO
import matplotlib.pyplot as plt

import madnex
from madnex.data.evolution import Evolution
from madnex.api.api_local import Pivot as Pivot
from madnex.madnex_constant import datatype, TYPE_CDX2HDF, MADNEX_PATH
import madnex.IO as IO


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = numpy.linspace(-4, 4, 500)
    y = numpy.exp( -t**2 ) + numpy.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, numpy.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    try:
        window_size = numpy.abs(numpy.int(window_size))
        order = numpy.abs(numpy.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = numpy.mat([[k**i for i in order_range]
                      for k in range(-half_window, half_window+1)])
    m = numpy.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))
    return numpy.convolve( m[::-1], y, mode='valid')

def reader(fichier_acqui, detect_rupt, dump, lissage_v = 0):
    """ read the acquisition raw file
        and build four madnex evolutions """
    def __prepare(_lines):
        nbl = _lines.__len__()
        ev_MT = None
        ev_MC = None
        ev_F = None
        ev_F2 = None
        irupt = 0
        _il = 0
        while _il < nbl :
            if "MISE EN TEMPERATURE" in _lines[_il] :
                _il += 3
                a_MT_t = numpy.array([])
                a_MT_T1 = numpy.array([])
                a_MT_T2 = numpy.array([])
                a_MT_T3 = numpy.array([])
                while _lines[_il] not in ["\r\n", "\r", "\n"] and _il!=(nbl-1) :
                    _line = _lines[_il].strip('\n').split('\t')
                    a_MT_t = numpy.append(a_MT_t, float(_line[1]) * 60.)
                    a_MT_T1 = numpy.append(a_MT_T1, float(_line[2]))
                    a_MT_T2 = numpy.append(a_MT_T2, float(_line[3]))
                    a_MT_T3 = numpy.append(a_MT_T3, float(_line[4]))
                    _il += 1
                a_MT_T = (a_MT_T1 + a_MT_T2 + a_MT_T3) / 3.
                a_MT = numpy.vstack((a_MT_t, a_MT_T1, a_MT_T2, a_MT_T3, a_MT_T))
                titre_MT = ['Temps', 'Temperature 1', 'Temperature 2',
                            'Temperature 3', 'Temperature moyenne']
                unite_MT = ['s','°C','°C','°C','°C']
                descrip_MT = ['Temps ecoule', 'Thermocouple 1', 'Thermocouple 2'
                              'Thermocouple 3', 'Moyenne thermocouples']
                ev_MT = Evolution(a_MT, titre_MT, unite_MT, descrip_MT)
            elif "MISE EN CHARGE" in _lines[_il] :
                _il += 3
                a_MC_C = numpy.array([])
                a_MC_t = numpy.array([])
                a_MC_A1 = numpy.array([])
                a_MC_A2 = numpy.array([])
                a_MC_T1 = numpy.array([])
                a_MC_T2 = numpy.array([])
                a_MC_T3 = numpy.array([])
                while _lines[_il] not in ["\r\n", "\r", "\n"] and _il!=(nbl-1) :
                    _line = _lines[_il].strip('\n').split('\t')
                    a_MC_C = numpy.append(a_MC_C, float(_line[2]))
                    a_MC_t = numpy.append(a_MC_t, float(_line[3]))
                    a_MC_A1 = numpy.append(a_MC_A1, float(_line[4]) / 100.)
                    a_MC_A2 = numpy.append(a_MC_A2, float(_line[5]) / 100.)
                    a_MC_T1 = numpy.append(a_MC_T1, float(_line[6]))
                    a_MC_T2 = numpy.append(a_MC_T2, float(_line[7]))
                    a_MC_T3 = numpy.append(a_MC_T3, float(_line[8]))
                    _il += 1
                a_MC_T = (a_MC_T1 + a_MC_T2 + a_MC_T3) / 3.
                a_MC_A = (a_MC_A1 + a_MC_A2) / 2.
                a_MC = numpy.vstack((a_MC_t, a_MC_C, a_MC_T1, a_MC_T2, a_MC_T3, a_MC_T,
                                     a_MC_A1, a_MC_A2, a_MC_A))
                titre_MC = ['Temps', 'Contrainte', 'Temperature 1', 'Temperature 2',
                            'Temperature 3', 'Temperature moyenne', 'Allongement 1',
                            'Allongement 2', 'Allongement moyen']
                unite_MC = ['s','MPa','°C','°C','°C','°C','mm/mm','mm/mm','mm/mm']
                descrip_MC = ['Temps ecoule', 'Chargement applique', 'Thermocouple 1',
                              'Thermocouple 2', 'Thermocouple 3', 'Moyenne thermocouples',
                              'LVDT 1', 'LVDT 2', 'LVDT moyenne']
                ev_MC = Evolution(a_MC, titre_MC, unite_MC, descrip_MC)
            elif "FLUAGE" in _lines[_il] :
                _il += 3
                a_F_C = numpy.array([])
                a_F_t = numpy.array([])
                a_F_t2 = numpy.array([])
                a_F_A1 = numpy.array([])
                a_F_A2 = numpy.array([])
                a_F_T1 = numpy.array([])
                a_F_T2 = numpy.array([])
                a_F_T3 = numpy.array([])
                a_F_rupt = numpy.array([])
                while _lines[_il] not in ["\r\n", "\r", "\n"] and _il!=(nbl-1) :
                    _line = _lines[_il].strip('\n').split('\t')
                    a_F_C = numpy.append(a_F_C, float(_line[2]))
                    a_F_t = numpy.append(a_F_t, float(_line[3]) * 3600.)
                    a_F_t2 = numpy.append(a_F_t2, float(_line[3]))
                    a_F_A1 = numpy.append(a_F_A1, float(_line[4]) / 100.)
                    a_F_A2 = numpy.append(a_F_A2, float(_line[5]) / 100.)
                    a_F_T1 = numpy.append(a_F_T1, float(_line[6]))
                    a_F_T2 = numpy.append(a_F_T2, float(_line[7]))
                    a_F_T3 = numpy.append(a_F_T3, float(_line[8]))
                    a_F_rupt = numpy.append(a_F_rupt, float(_line[9]))
                    _il += 1
                a_F_T = (a_F_T1 + a_F_T2 + a_F_T3) / 3.
                a_F_A = (a_F_A1 + a_F_A2) / 2.
                a_F_At = a_MC_A[-1] + a_F_A
                a_F = numpy.vstack((a_F_t, a_F_C, a_F_T1, a_F_T2, a_F_T3, a_F_T,
                                    a_F_A1, a_F_A2, a_F_A, a_F_rupt))
                titre_F = ['Temps', 'Contrainte', 'Temperature 1', 'Temperature 2',
                           'Temperature 3', 'Temperature moyenne', 'Allongement 1',
                           'Allongement 2', 'Allongement moyen']
                unite_F = ['s','MPa','°C','°C','°C','°C','mm/mm','mm/mm','mm/mm']
                descrip_F = ['Temps ecoule', 'Chargement applique', 'Thermocouple 1',
                             'Thermocouple 2', 'Thermocouple 3', 'Moyenne thermocouples',
                             'LVDT 1', 'LVDT 2', 'LVDT moyenne', 'Rupture detectee']
                ev_F = Evolution(a_F, titre_F, unite_F, descrip_F)

                if detect_rupt :
                    irupt = find_rupture(a_F_rupt, a_F_A)
                else :
                    irupt = find_rupture(a_F_rupt)

                if irupt != 0 :
                    lid = numpy.array(xrange(irupt))
                else :
                    lid = numpy.array(xrange(a_F_t.size+1))
                if lissage_v == 0 :
                    autoliss = 2*int(0.01*lid.size)+1
                else :
                    autoliss = lissage_v
                all_s = savitzky_golay(a_F_A[lid], autoliss ,1, 0, 1)
                csrate = []
                step = 10
                cutoff = 0
                for _ii in range(all_s.size-step-cutoff):
                    csrate.append((all_s[_ii+step+cutoff] - all_s[_ii+cutoff])\
                                  /(a_F_t[_ii+step+cutoff] - a_F_t[_ii+cutoff]))
                a_F_v = savitzky_golay(csrate, autoliss, 1, 0, 1)
                a_F_v = numpy.append([numpy.nan]*step, a_F_v)
                a_F2 = numpy.vstack((a_F_t[lid], a_F_C[lid], a_F_T[lid],
                                     a_F_A[lid], a_F_t2[lid], a_F_v,
                                     a_F_At[lid]))
                titre_F2 = ['Time', 'Stress', 'Temperature',
                            'Creep strain', 'Time (h)', 'Strain rate',
                            'Total strain']
                unite_F2 = ['s','MPa','°C','mm/mm', 'h','s-1','mm/mm']
                descrip_F2 = ['Temps ecoule','Chargement applique',
                              'Moyenne thermocouples', 'LVDT moyenne',
                              'Temps ecoule (h)', 'Vitesse moyenne LVDT'
                              'LVDT total moyenne']
                ev_F2 = Evolution(a_F2, titre_F2, unite_F2, descrip_F2)
            _il += 1
        return [ev_MT, ev_MC, ev_F, ev_F2], irupt

    if fichier_acqui is None:
        from cStringIO import StringIO
        from IPython.display import display
        import fileupload
        def _cbk(change,):
            dump[0],dump[1] = __prepare(StringIO(change['new']).readlines())
            print "File '%s' uploaded successfully!" % change['owner'].filename
        #
        _upload_widget = fileupload.FileUploadWidget()
        _upload_widget.observe(_cbk, names='data')
        display(_upload_widget)
        return
    else:
        with open(fichier_acqui, "r") as _file:
            dump[0],dump[1] = __prepare(_file.readlines())
        return dump[1]

def find_rupture(isruprt, allmoy = []):
    idr = isruprt.size - numpy.argmax(isruprt[::-1] == 0)
    if len(allmoy) != 0 and idr != 0 :
        marge = 40
        llid = xrange(max(idr-marge,0),min(idr+marge,isruprt.size))
        look = allmoy[numpy.array(llid)]\
               - (allmoy[numpy.array(numpy.add(llid,-1))] +\
                  allmoy[numpy.array(numpy.add(llid,-2))] +\
                  allmoy[numpy.array(numpy.add(llid,-3))])/2.
        idr2 = idr - marge + numpy.argmax(look>0.)
        if ((idr - marge) < idr2) and (idr2 < (idr + marge)) :
            return idr2
    return idr

def compute_mt_images(workpath, evol):
    """ cree une image png de l evolution de la temperature
        pendant la mise en temperature de l essai """
    arr = evol.data
    tit = evol.title
    uni = evol.unit
    title_font = {'fontname':'Arial', 'size':'16',
                  'color':'black', 'weight':'normal'}
    axis_font = {'fontname':'Arial', 'fontsize':'12'}
    dheight = 0.05 * arr[1:5].max()
    plt.figure(figsize = (18, 7), dpi = 300)
    lit = xrange(0, arr[0].size, 4)
    plt.subplot(121)
    plt.plot(arr[0,lit] / 60., arr[1,lit], "b+", label=tit[1])
    plt.plot(arr[0,lit] / 60., arr[2,lit], "g+", label=tit[2])
    plt.plot(arr[0,lit] / 60., arr[3,lit], "r+", label=tit[3])
    plt.plot(arr[0] / 60., arr[4], "k-", label=tit[4])
    plt.grid(False)
    plt.title('Mise en temperature', **title_font)
    plt.xlabel('%s (%s)' %(tit[0],'min'), labelpad=20, **title_font)
    plt.ylabel('Temperature (deg)', **title_font)
    plt.xticks(xrange(0, int(max(arr[0])), 100), **axis_font)
    plt.yticks(xrange(0, int(arr[1:5].max() + dheight), 100), **axis_font)
    plt.axis([0, max(arr[0] / 60.),
              arr[1:5].min(), arr[1:5].max() + dheight])
    plt.legend(loc=5, prop={'size':title_font['size'],
                            'family':title_font['fontname']})
    plt.subplot(122)
    plt.plot(arr[0] / 60., arr[1], "b+", label=tit[1])
    plt.plot(arr[0] / 60., arr[2], "g+", label=tit[2])
    plt.plot(arr[0] / 60., arr[3], "r+", label=tit[3])
    plt.plot(arr[0] / 60., arr[4], "k-", label=tit[4])
    plt.grid(False)
    plt.title('Mise en temperature [ZOOM]', **title_font)
    plt.xlabel('%s (%s)' %(tit[0],'min'), labelpad=20, **title_font)
    plt.ylabel('Temperature (deg)', **title_font)

    plt.xticks(xrange(0, int(max(arr[0]) / 60.), 20), **axis_font)
    plt.yticks(xrange(0, int(arr[1:5].max() + dheight), 1), **axis_font)
    plt.axis([arr[0][-30] / 60., max(arr[0]) / 60.,
              arr[1:5,-30:].min() - 0.15*dheight,
              arr[1:5,-30:].max() + 0.10*dheight])
    plt.legend(loc=4, prop={'size':title_font['size'],
                            'family':title_font['fontname']})
#    plt.tight_layout(pad=1., w_pad=1., h_pad=5.)
    png_filename = "0_chauffage_temperature.png"
    #plt.savefig(path.join(workpath, png_filename))
    plt.show()
    #return [png_filename]

def compute_mc_images(workpath, evol):
    """ cree trois images png des evolutions de la temperature,
        de la contrainte et de l allongement pendant
        la mise en charge de l essai """
    arr = evol.data
    tit = evol.title
    uni = evol.unit
    title_font = {'fontname':'Arial', 'size':'16',
                  'color':'black', 'weight':'normal'}
    axis_font = {'fontname':'Arial', 'fontsize':'12'}
    dheight = 0.05 * arr[2:6].max()
    plt.figure(figsize = (8, 7), dpi = 300)
    plt.plot(arr[0], arr[2], "b+", label=tit[2])
    plt.plot(arr[0], arr[3], "g+", label=tit[3])
    plt.plot(arr[0], arr[4], "r+", label=tit[4])
    plt.plot(arr[0], arr[5], "k-", label=tit[5])
    plt.grid(False)
    plt.title('Mise en charge', **title_font)
    plt.xlabel('%s (%s)' %(tit[0],uni[0]), labelpad=20, **title_font)
    plt.ylabel('Temperature (deg)', **title_font)
    plt.xticks(xrange(int(arr[0,0]), int(arr[0,-1]), 1), **axis_font)
    plt.yticks(xrange(int(arr[2:6].min() - dheight),
                      int(arr[2:6].max() + dheight), 1), **axis_font)
    plt.axis([min(arr[0]), max(arr[0]),
              arr[2:6].min() - 0.15*dheight,
              arr[2:6].max() + 0.10*dheight])
    plt.legend(loc=4, prop={'size':title_font['size'],
                            'family':title_font['fontname']})
    png_filename1 = "1_charge_temperature.png"
    #plt.savefig(path.join(workpath, png_filename1))
    plt.show()
    yield
    ###
    dheight = 0.05 * max(arr[1])
    plt.figure(figsize = (8, 7), dpi = 300)
    plt.plot(arr[0], arr[1], "k-", label=tit[1])
    plt.grid(False)
    plt.title('Mise en charge', **title_font)
    plt.xlabel('%s (%s)' %(tit[0],uni[0]), labelpad=20, **title_font)
    plt.ylabel('%s (%s)' %(tit[1],uni[1]), **title_font)
    plt.xticks(xrange(int(arr[0,0]), int(arr[0,-1]), 1), **axis_font)
    plt.yticks(xrange(0, int(max(arr[1]) + dheight), 10), **axis_font)
    plt.axis([min(arr[0]), max(arr[0]), 0, max(arr[1]) + dheight])
    png_filename2 = "1_charge_contrainte.png"
    #plt.savefig(path.join(workpath, png_filename2))
    plt.show()
    yield
    ###
    dheight = 0.05 * arr[6:9].max()
    plt.figure(figsize = (8, 7), dpi = 300)
    plt.plot(arr[0], arr[6], "b+", label=tit[6])
    plt.plot(arr[0], arr[7], "r+", label=tit[7])
    plt.plot(arr[0], arr[8], "k-", label=tit[8])
    plt.grid(False)
    plt.title('Mise en charge', **title_font)
    plt.xlabel('%s (%s)' %(tit[0],uni[0]), labelpad=20, **title_font)
    plt.ylabel('Allongement (%s)' %uni[6], **title_font)
    plt.xticks(xrange(int(arr[0,0]), int(arr[0,-1]), 1), **axis_font)
    plt.yticks(**axis_font)
    plt.axis([min(arr[0]), max(arr[0]),
              arr[6:9].min() - dheight,
              arr[6:9].max() + dheight])
    plt.legend(loc=2, prop={'size':title_font['size'],
                            'family':title_font['fontname']})
    png_filename3 = "1_charge_allongement.png"
    #plt.savefig(path.join(workpath, png_filename3))
    plt.show()
    
    #return [png_filename1, png_filename2, png_filename3]

def compute_fl_images(workpath, evol, evol2, irupt = 0):
    """ cree trois images png des evolutions de la temperature,
        de l allongement et de l allongement moyen pendant
        l essai de fluage """
    arr = evol.data
    tit = evol.title
    uni = evol.unit
    arr2 = evol2.data
    tit2 = evol2.title
    uni2 = evol2.unit
    title_font = {'fontname':'Arial', 'size':'16',
                  'color':'black', 'weight':'normal'}
    axis_font = {'fontname':'Arial', 'fontsize':'12'}
    if irupt != 0 :
        lid = xrange(irupt)
    else :
        lid = xrange(arr[0].size+1)
    dheight = 0.05 * arr[2:6,lid].max()
    plt.figure(figsize = (8, 7), dpi = 300)
    plt.plot(arr[0,lid] / 3600., arr[2,lid], "b+", label=tit[2])
    plt.plot(arr[0,lid] / 3600., arr[3,lid], "g+", label=tit[3])
    plt.plot(arr[0,lid] / 3600., arr[4,lid], "r+", label=tit[4])
    plt.plot(arr[0,lid] / 3600., arr[5,lid], "k-", label=tit[5])
    plt.grid(False)
    plt.title('Fluage', **title_font)
    plt.xlabel('%s (h)' %tit[0], labelpad=20, **title_font)
    plt.ylabel('Temperature (deg)', **title_font)
    plt.xticks(**axis_font)
    plt.yticks(xrange(int(arr[2:6,0].mean() - dheight),
                      int(arr[2:6,0].mean() + dheight), 1), **axis_font)
    plt.axis([arr[0,0] / 3600., arr[0,lid[-1]]*1.05 / 3600.,
              arr[2:6,0].mean() - 0.15*dheight,
              arr[2:6,0].mean() + 0.10*dheight])
    plt.legend(loc=4, prop={'size':title_font['size'],
                            'family':title_font['fontname']})
    png_filename1 = "2_fluage_temperature.png"
    #plt.savefig(path.join(workpath, png_filename1))
    plt.show()
    yield
    ###
    dheight = 0.05 * arr[6:9,lid].max()
    plt.figure(figsize = (8, 7), dpi = 300)
    plt.plot(arr[0,lid] / 3600., arr[6,lid], "b+", label=tit[6])
    plt.plot(arr[0,lid] / 3600., arr[7,lid], "r+", label=tit[7])
    plt.plot(arr[0,lid] / 3600., arr[8,lid], "k-", label=tit[8])
    plt.grid(False)
    plt.title('Fluage', **title_font)
    plt.xlabel('%s (h)' %tit[0], labelpad=20, **title_font)
    plt.ylabel('Allongement (%s)' %uni[6], **title_font)
    plt.xticks(**axis_font)
    plt.yticks(**axis_font)
    plt.axis([arr[0,0] / 3600., arr[0,lid[-1]]*1.05 / 3600.,
              arr[6:9,lid].min() - dheight,
              arr[6:9,lid].max() + dheight])
    plt.legend(loc=2, prop={'size':title_font['size'],
                            'family':title_font['fontname']})
    png_filename2 = "2_fluage_allongement.png"
    #plt.savefig(path.join(workpath, png_filename2))
    plt.show()
    yield
    ###
    plt.figure(figsize = (8, 7), dpi = 300)
    plt.plot(arr2[4], arr2[3], "k-", label=tit[8])
    plt.grid(False)
    plt.title('Fluage', **title_font)
    plt.xlabel('Temps (%s)' %uni2[4], labelpad=20, **title_font)
    plt.ylabel('Allongement (%s)' %uni2[3], **title_font)
    plt.xticks(**axis_font)
    plt.yticks(**axis_font)
    plt.axis([arr2[4,0], arr2[4,-1]*1.05,
              arr[6:9,lid].min() - dheight,
              arr[6:9,lid].max() + dheight])
    png_filename3 = "2_fluage_allongement_moy.png"
    #plt.savefig(path.join(workpath, png_filename3))
    plt.show()
    yield
    ###
    plt.figure(figsize = (8, 7), dpi = 300)     
    plt.semilogy(arr2[4], arr2[5], "k-", label=tit2[4])
    plt.grid(False)
    plt.title('Fluage', **title_font)
    plt.xlabel('Temps (%s)' %uni2[4], labelpad=20, **title_font)
    plt.ylabel('Vitesse de fluage (%s)' %uni2[5], **title_font)
    plt.xticks(**axis_font)
    plt.yticks(**axis_font)
    plt.axis([arr2[4,0], arr2[4,-1] * 1.05,
              10**(-10), 10**(-4)])
    png_filename4 = "2_fluage_vitesse_moy.png"
    #plt.savefig(path.join(workpath, png_filename4))
    plt.show()
    yield
    #return [png_filename1, png_filename2, png_filename3, png_filename4]

