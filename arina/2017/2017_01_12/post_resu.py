#!/usr/bin/python
import numpy as np
import scipy.signal.spline as spline
import pylab as plb
import pdb

def lire_exp (_nom_essai,_nom_rep) :
	_f = open("./" + _nom_rep + "/" + _nom_essai +'.prn', "r")
	_L1 = _f.readlines()
	_f.close()

	_entete = _L1[0].split()
	del _L1[0]

	_donnees = {}
	_donnees['nom'] = _nom_essai
	for x in _entete :
		_donnees[x] = []

	for x in _L1 :
		_L2 = x.split()
		for y in range(len(_L2)) :
			_donnees[_entete[y]].append(_L2[y])

	return _donnees

def lire_exp_fluage (_nom_essai,_nom_rep) :
    _f = open("./essais_fluage/" + _nom_rep + "/" + _nom_essai+'.acquisition', "r")
    _L1 = _f.readlines()
    _f.close()
    index = -1
    for x in range(len(_L1)):
        if _L1[x].find('FLUAGE')==0:
            index=x
    for j in range(index):
        del _L1[0]
    _type=(_L1[0].split('\t'))[0].strip('\r\n')
    del _L1[0]
    _entete_1=_L1[0].split('\t')
    del _L1[0]
    _entete_2=_L1[0].split('\t')
    del _L1[0]
    _entete_1[-1]=_entete_1[-1].strip('\r\n')
    _entete_2[-1]=_entete_2[-1].strip('\r\n')
    _entete=_entete_2
    for i in range(len(_entete_2)):
        _entete[i]=''.join(_entete_2[i].split())+'_'+str(i)


    #print _entete
    _donnees = {}
    _donnees['nom'] = _nom_essai
    for x in _entete :
        _donnees[x] = []
    for x in _L1 :
        _L2 = x.split('\t')
        for y in range(len(_L2)) :
            #print _donnees.keys()
            _donnees[_entete[y]].append(_L2[y])
            #_donnees[y].append(_L2[y])

    return _donnees

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
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
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
    #import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

liste_num_essai_1 = [2163, 2165, 2169, 2171, 2172, 2173, 2195]
liste_num_essai_1 = [2266, ]
#liste_label_1 = ['130MPa-F','130MPa-C','130MPa-C/E','130MPa-F/E','100MPa-C','100MPa-C','130MPa-C2']
liste_label_1 = ['130-C','130-WM','130-WM-N','130-C-N','100-WM-1','100-WM-2','140-H-1']

liste_num_essai_2 = [2208,2209,2197,2187,2188,2198, 2196, 2202]
liste_num_essai_2 = []
#liste_label_2 = ['120MPa-C2-6','160MPa-C2-E3','120MPa-C2-4','130MPa-F-E12','130MPa-F-E13','160MPa-C2-E1','130MPa-F-E5']
liste_label_2 = ['120-H-6','160-H-N-3','120-H-4','130-C-N-12','130-C-N-13','160-H-N-1', '140-H-2', '140-H-3']

#liste_num_essai_2 = [2197,2187,2188,2198,2171]
#liste_label_2 = ['120MPa-C2-4','130MPa-F-12','160MPa-C2-E1','130MPa-F-5']

liste_num_essai=liste_num_essai_1+liste_num_essai_2
liste_label = liste_label_1+liste_label_2


#don = lire_tab_aster(liste_nom_essai[0],liste_nom_rep[0])

#print don.keys()

fig_width_pt = 1024.  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
lw = 1.5 # epaisseur des courbes
params = {'backend': 'ps',
    'axes.labelsize': 24,
    'text.fontsize': 24,
    'legend.fontsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'text.usetex': True,
    'figure.figsize': fig_size}
plb.rcParams.update(params)

for i in range(np.size(liste_num_essai)) :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=np.array(don[cle_t],dtype=float)
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    flu_s=savitzky_golay(flu,11,1,deriv=0,rate=1)

    plot = plb.plot(tps,flu,'-',linewidth=2,label=(don['nom'].split('.'))[0]);
    plot = plb.plot(tps,flu_s,'.k',linewidth=0.5);

plb.xlabel('Crrep time (h)')
plb.ylabel('Creep Strain \%')
plb.title('Creep curves');
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='upper left');
plb.savefig("./figs/all_creep.png")
plb.clf()



for i in range(np.size(liste_num_essai)) :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=don[cle_t]
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    if i==2 or i==3:
        flutot=0.5*(flu1+flu2)
        l_ref=4*np.sqrt(1.2**2-0.2**2)
        #pour la deformation de fluage on doit corriger car elle est calculee pour la longueur entre collerettes de 36mm, les entailles font 1.2mm de rayon et de profondeur 1mm,
        flu=flutot*36/l_ref
    else:
        flu=0.5*(flu1+flu2)
    plot = plb.plot(tps,flu,'-',linewidth=2,label=liste_label[i]);


plb.xlabel('Creep time (h)')
plb.ylabel('Creep Strain \%')
plb.title('Creep curves');
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='upper left');
plb.savefig("./figs/all_creep_labelized.png")
plb.savefig("./figs/all_creep_labelized.pdf")
plb.clf()

for i in [0, 1, 4, 5,6,13,14,7,9] :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=don[cle_t]
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    plot = plb.plot(tps,flu,'-',linewidth=2,label=liste_label[i]);

plb.xlabel('Creep time (h)')
plb.ylabel('Creep Strain \%')
plb.title('Creep curves');
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='upper left');
plb.savefig("./figs/all_creep_labelized_smooth.png")
plb.savefig("./figs/all_creep_labelized_smooth.pdf")
plb.clf()


for i in [7,9,0,6,13,14] :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=don[cle_t]
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    plot = plb.plot(tps,flu,'-',linewidth=2,label=liste_label[i]);

plb.xlabel('Creep time (h)')
plb.ylabel('Creep Strain \%')
plb.title('Creep curves');
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='upper left');
plb.savefig("./figs/all_creep_labelized_smooth_select.png")
plb.savefig("./figs/all_creep_labelized_smooth_select.pdf")
plb.clf()


for i in [4,5,1] :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=don[cle_t]
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    plot = plb.plot(tps,flu,'-',linewidth=2,label=liste_label[i]);

plb.xlabel('Creep time (h)')
plb.ylabel('Creep Strain \%')
plb.title('Creep curves');
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='upper left');
plb.savefig("./figs/all_creep_labelized_smooth_WM.png")
plb.savefig("./figs/all_creep_labelized_smooth_WM.pdf")
plb.clf()




for i in [3,10,11,12,8] :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=don[cle_t]
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    l_ref=4*np.sqrt(1.2**2-0.2**2)
    #pour la deformation de fluage on doit corriger car elle est calculee pour la longueur entre collerettes de 36mm, les entailles font 1.2mm de rayon et de profondeur 1mm,
    flu_corr=flu*36/l_ref
    plot = plb.plot(tps,flu_corr,'-',linewidth=2,label=liste_label[i]);


plb.xlabel('Creep time (h)')
plb.ylabel('Creep Strain \%')
plb.ylim(0,10)
plb.title('Creep curves');
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='upper left');
plb.savefig("./figs/all_creep_labelized_notched.png")
plb.savefig("./figs/all_creep_labelized_notched.pdf")
plb.clf()


for i in [0, 1, 4, 5,6,13,14,7,9] :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=np.array(don[cle_t],dtype=float)
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    flu_s=savitzky_golay(flu,21,1,deriv=0,rate=1)
    csrate=[]
    tpsrate=[]
    step=40
    cutoff=200
    cutoff_s=10
    for j in range(len(flu_s)-step-cutoff-cutoff_s):
        csrate.append(0.01*(flu_s[j+step+cutoff_s]-flu_s[j+cutoff_s])/(tps[j+step+cutoff_s]-tps[j+cutoff_s]))
        tpsrate.append(tps[j+step+cutoff_s])
    csrate_s=savitzky_golay(csrate,11,1,deriv=0,rate=1)
    plot = plb.semilogy(tpsrate,csrate_s,'-',linewidth=2,label=liste_label[i]);

plb.xlabel('Creep time (h)')
plb.ylabel('Creep rate (h$^{-1}$)')
plb.title('Creep curves');
plb.ylim(1e-7,1e-3);
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='lower right');
plb.savefig("./figs/creep_rate_smooth.png")
plb.savefig("./figs/creep_rate_smooth.pdf")
plb.clf()

for i in [7,9,0,6,13,14] :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=np.array(don[cle_t],dtype=float)
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    flu_s=savitzky_golay(flu,21,1,deriv=0,rate=1)
    csrate=[]
    tpsrate=[]
    step=40;
    cutoff=200
    cutoff_s=10
    for j in range(len(flu_s)-step-cutoff-cutoff_s):
        csrate.append(0.01*(flu_s[j+step+cutoff_s]-flu_s[j+cutoff_s])/(tps[j+step+cutoff_s]-tps[j+cutoff_s]))
        tpsrate.append(tps[j+step+cutoff_s])
    csrate_s=savitzky_golay(csrate,11,1,deriv=0,rate=1)
    plot = plb.semilogy(tpsrate,csrate_s,'-',linewidth=2,label=liste_label[i]);

plb.xlabel('Creep time (h)')
plb.ylabel('Creep rate (h$^{-1}$)')
plb.title('Creep curves');
plb.ylim(1e-7,1e-3)
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='lower right');
plb.savefig("./figs/creep_rate_smooth_select.png")
plb.savefig("./figs/creep_rate_smooth_select.pdf")
plb.clf()

for i in [4,5,1] :
    nom_essai=str(liste_num_essai[i])
    don=lire_exp_fluage(nom_essai,'.')
    cles=don.keys()
    #print cles
    cles_l=[]
    cles_t=[]
    for x in cles:
        if (x.find('Allongement')>=0) and (x.find('%')>=0):
            cles_l.append(x)
        if x.find('Tempsheure')>=0:
            cle_t=x

    #print cles_l,cle_t
    tps=np.array(don[cle_t],dtype=float)
    #print 'label defo = ',cles[4]
    #flu=0.5*(don[cles[4]]+don[cles[6]])
    flu1=np.array(don[cles_l[0]],dtype=float)
    flu2=np.array(don[cles_l[1]],dtype=float)
    flu=0.5*(flu1+flu2)
    flu_s=savitzky_golay(flu,21,1,deriv=0,rate=1)
    csrate=[]
    tpsrate=[]
    step=40;
    cutoff=200
    cutoff_s=10
    for j in range(len(flu_s)-step-cutoff-cutoff_s):
        csrate.append(0.01*(flu_s[j+step+cutoff_s]-flu_s[j+cutoff_s])/(tps[j+step+cutoff_s]-tps[j+cutoff_s]))
        tpsrate.append(tps[j+step+cutoff_s])
    csrate_s=savitzky_golay(csrate,11,1,deriv=0,rate=1)
    plot = plb.semilogy(tpsrate,csrate_s,'-',linewidth=2,label=liste_label[i]);

plb.xlabel('Creep time (h)')
plb.ylabel('Creep rate (h$^{-1}$)')
plb.title('Creep curves');
plb.ylim(1e-7,1e-3)
#plb.legend(liste_legende,loc='upper right');
plb.legend(loc='lower right');
plb.savefig("./figs/creep_rate_smooth_WM.png")
plb.savefig("./figs/creep_rate_smooth_WM.pdf")
plb.clf()



list_vmin=[6E-6,7E-6,6E-6,2E-5,2E-5,7E-5,5E-5,2E-5,5E-5]
list_tr=[6500,5800,4000,2416,1872,794,760,1650,840]
#liste_label2 = ['120MPa-C2','120MPa-C2','130MPa-F','130MPa-C','100MPa-C','100MPa-C','130MPa-C2']
liste_label2 = ['120-H-4','120-H-6','130-C','130-WM','100-WM-1','100-WM-2','140-H-1','140-H-2','140-H-3']
list_sig=[120,120,130,130,100,100,140,140,140]


list_sig_C_S=[130]
list_tr_C_S=[4000]
list_vmin_C_S=[6E-6]

list_sig_H_S=[120,120,140,140,140]
list_tr_H_S=[6500,5800,760,1650,840]
list_vmin_H_S=[6E-6,7E-6,5E-5,2E-5,5E-5]

list_sig_WM_S=[100,100,130]
list_tr_WM_S=[794,1872,2416]
list_vmin_WM_S=[7E-5,2E-5,2E-5]







list_tr_calc=np.linspace(800,7000,10)
list_vmin_calc_EM=[]
list_vmin_calc_HP40=[]
for i in range(len(list_tr_calc)):
    tr=list_tr_calc[i]
    list_vmin_calc_EM.append(0.103/(np.power(tr,1.08)))
    list_vmin_calc_HP40.append(0.225/(np.power(tr,1.19)))

plb.loglog(list_tr_calc,list_vmin_calc_EM,'-k',linewidth=2,label='MG E. Molinie')
plb.loglog(list_tr_calc,list_vmin_calc_HP40,'--k',linewidth=2,label='MG H-P40')

#plb.loglog(list_tr,list_vmin,'sk',linewidth=5,markersize=12,label='ex-service 15CDV4-10')
plb.loglog(list_tr_C_S,list_vmin_C_S,'sb',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (C)')
plb.loglog(list_tr_WM_S,list_vmin_WM_S,'dm',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (WM)')
plb.loglog(list_tr_H_S,list_vmin_H_S,'or',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (H)')

plb.xlim(500,10000)
plb.ylim(4E-6,1E-4)
plb.xlabel('Time to failure (h)')
plb.ylabel('Min. Creep Strain Rate (h-1)')
plb.legend(loc='upper right');
plb.savefig("./figs/Monkman_Grant.png")
plb.savefig("./figs/Monkman_Grant.pdf")
plb.clf()


expe=lire_exp('German1986','.')

list_sig_ref=np.array(expe['SIG'],dtype='float')
list_tr_ref=np.array(expe['tr'],dtype='float')

list_vmin_ref=0.103/np.power(list_tr_ref,1.08)



list_vmin_F=6E-6
list_sig_F=130.

list_vmin_C=[7E-5,2E-5,2E-5]
list_sig_C=[100,100,130]

list_vmin_C2=[6E-6,6.5E-6,6E-5,5E-5,2E-5,5E-5]
list_sig_C2=[120,120,130,140,140,140]

plb.loglog(list_sig_ref,list_vmin_ref,'^k',linewidth=5,markersize=12,label='German 1986 560C')

plb.loglog(list_sig_F,list_vmin_F,'sb',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (C)')

plb.loglog(list_sig_C,list_vmin_C,'dm',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (WM)')

plb.loglog(list_sig_C2,list_vmin_C2,'or',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (H)')

#plb.xlim(80,250)
plb.xlim(50,1000)
plb.ylim(1E-7,1E-2)
plb.xlabel('Creep Stress (MPa)')
plb.ylabel('Min. Creep Strain Rate (h-1)')
plb.legend(loc='lower right');
plb.savefig("./figs/Norton.png")
plb.savefig("./figs/Norton.pdf")
plb.clf()



plb.loglog(list_tr_ref,list_sig_ref,'^k',linewidth=5,markersize=12,label='German 1986 560C')
plb.loglog(list_tr_C_S,list_sig_C_S,'sb',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (C)')
plb.loglog(list_tr_WM_S,list_sig_WM_S,'dm',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (WM)')
plb.loglog(list_tr_H_S,list_sig_H_S,'or',linewidth=5,markersize=12,label='ex-service 15CDV4-10 (H)')

plb.xlim(10,100000)
plb.ylim(40,1000)
plb.xlabel('Creep Time (h)')
plb.ylabel('Creep Stress (MPa)')
plb.legend(loc='lower left');
plb.savefig("./figs/sig_tr_all.png")
plb.savefig("./figs/sig_tr_all.pdf")
plb.clf()

# Comparison with the "best model" : Hayhurst Sp2, with data imported from matlab

list_sig_Hsp2=[80,100,120,140,160,180,200,220,240,260,280,300]
list_tr_Hsp2=[1.0833742e+05,4.3487565e+04,1.7395040e+04,6.8359608e+03,2.7160911e+03,1.0742487e+03,4.1507240e+02,1.6361089e+02,6.3516562e+01,2.4457396e+01,9.4216667e+00,3.6632812e+00]

list_sig_C_S_av=[130]
list_tr_C_S_av=[4000]


list_sig_H_S_av=[120,140]
list_tr_H_S_av=[6150,1157]


list_sig_WM_S_av=[100,130]
list_tr_WM_S_av=[1333,2416]




plb.loglog(list_tr_Hsp2,list_sig_Hsp2,'-k',linewidth=2,markersize=12,label=r'Hayhurst $H_{sp2}$')
plb.loglog(list_tr_C_S_av,list_sig_C_S_av,'sb',linewidth=5,markersize=12,label='average (C)')
plb.loglog(list_tr_WM_S_av,list_sig_WM_S_av,'dm',linewidth=5,markersize=12,label='average (WM)')
plb.loglog(list_tr_H_S_av,list_sig_H_S_av,'or',linewidth=5,markersize=12,label='average (H)')

plb.xlim(100,100000)
plb.yticks([80,90,100,110,120,130,140,150,160,170,180,190,200],('', '', '100', '', '120','','140','','160','','','','') )
plb.ylim(80,200)
plb.xlabel('Creep Time (h)')
plb.ylabel('Creep Stress (MPa)')
plb.legend(loc='lower left');
plb.savefig("./figs/sig_tr_compar.png")
plb.savefig("./figs/sig_tr_compar.pdf")
plb.clf()



