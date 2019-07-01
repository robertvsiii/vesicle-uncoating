from SloppyCell.ReactionNetworks import *
import csv
import numpy as np
import pickle

def ImportData(PROBE, folder, CLTAmax, probechoice, control, singletrack):
    
    #Truncate initial steady-state region, and increase size of error bars for control data.
    trunc = 12
    ControlErr = 1
    
    if singletrack:
        t, probe, probeerr, CLTA, CLTAerr = LoadData(PROBE,folder,tracktype='singletrack')
    else:
        t, probe, probeerr, CLTA, CLTAerr = LoadData(PROBE,folder,tracktype='normal')
        #Compensate for wrong error calculation in export file
        probeerr = probeerr/2.
        CLTAerr = CLTAerr/2. 
        
    
    #CLTA, CLTAerr = Rescale(CLTA,CLTAerr,CLTAmax)
    CLTAdict = dict((t[k], (CLTA[k],CLTAerr[k])) for k in range(trunc,len(t)))
    CLTAdict = RemoveNaNs(CLTAdict)
    
    expt = Experiment(PROBE)
    
    if probechoice == None:
        expt.set_data(dict([('net' + PROBE,dict([('CLTA',CLTAdict)]))]))
        #expt.set_fixed_sf({'CLTA' : 1})
    else:
        probedict = dict((t[k], (probe[k],probeerr[k])) for k in range(trunc,len(t)))
        if control:
            tNC, probeNC, probeerrNC, CLTANC, CLTAerrNC = LoadData(PROBE,folder,tracktype='control')
            NoClathrindict = dict((tNC[k], (probeNC[k],ControlErr)) for k in range(trunc,len(tNC)))
            NoClathrinCLTAdict = dict((tNC[k], (CLTANC[k],CLTAerrNC[k])) for k in range(trunc,len(tNC)))
            NoClathrinCLTAdict = RemoveNaNs(NoClathrinCLTAdict)
     
            tNL, probeNL, probeerrNL, CLTANL, CLTAerrNL = LoadData('Aux1TRUNC',folder,tracktype='control')
            NoLipiddict = dict((tNL[k], (probeNL[k],ControlErr)) for k in range(trunc,len(tNL)))
            NoLipidCLTAdict = dict((tNL[k], (CLTANL[k],CLTAerrNL[k])) for k in range(trunc,len(tNL)))
            NoLipidCLTAdict = RemoveNaNs(NoLipidCLTAdict)
            
            if PROBE in ['FYVE','TAPP1']:
                expt.set_data(dict([('net' + PROBE,dict([(PROBE,probedict),('CLTA',CLTAdict)])),
                                    ('net' + PROBE + 'NoClathrin',dict([(PROBE,NoClathrindict),
                                                                        ('CLTA',NoClathrinCLTAdict)])),
                                    ('net' + PROBE + 'NoLipid',dict([(PROBE,NoLipiddict),
                                                                     ('CLTA',NoLipidCLTAdict)]))]))
            else:
                expt.set_data(dict([('net' + PROBE,dict([(PROBE,probedict),('CLTA',CLTAdict)])),
                                    ('net' + PROBE + 'NoClathrin',dict([(PROBE,NoClathrindict),
                                                                        ('CLTA',NoClathrinCLTAdict)]))]))
        else:
            expt.set_data(dict([('net' + PROBE,dict([(PROBE,probedict),('CLTA',CLTAdict)]))]))
    return expt, [t, probe, probeerr], [t, CLTA, CLTAerr]

def RemoveNaNs(datadict):
    for item in datadict.keys():
        if np.isnan(datadict[item][0]):
            del datadict[item]
    return datadict

def MakeProbeModel(netdict, exptdict, priors):
    for id in netdict.keys():
        if id in priors.keys():
            val = priors[id]
            exptdict[id].set_sf_prior(id, prior_type = 'gaussian in log sf', prior_params = (np.log(val[0]),val[1]))    
        else:
            exptdict[id].set_fixed_sf({id : 1})
            
    net = []
    expt = []
    for item in netdict:
        net.append(netdict[item].copy('net' + item + 'NoClathrin'))
        net[-1].remove_component('bindpc')
        net[-1].remove_component('bindpcl')
        if net[-1].probe in ['FYVE','TAPP1']:
            net[-1].remove_component('bindpcl2')
            net.append(netdict[item].copy('net' + item + 'NoLipid'))
            net[-1].remove_component('bindpl')
            net[-1].remove_component('bindplc')
            net[-1].remove_component('bindpl')
            net[-1].remove_component('bindplc')
        net.append(netdict[item])
        expt.append(exptdict[item])
    
    model = Model(expt,net)
    for id, val in priors.items():
        if id not in netdict.keys():
            model.set_var_optimizable(id, True)
            model.AddResidual(Residuals.PriorInLog('prior_on_%s' % id, id, np.log(val[0]), val[1]))

    params = model.get_params()
    cost = model.cost(params)
    
    return model, params, cost

def LoadData(PROBE,folder,tracktype='normal'):
    if tracktype == 'control':
        groupname = PROBE + '_control'
    elif tracktype == 'singletrack':
        groupname = 'SingleTrack_' + PROBE
    else:
        groupname = PROBE
            
    f = open(folder + groupname + '_t.csv')
    reader0 = csv.reader(f)
    for row in reader0:
        t = np.asarray(row,dtype=float)
    f.close()

    f = open(folder + groupname + '_' + PROBE + '.csv')
    reader0 = csv.reader(f)
    for row in reader0:
        probe = np.asarray(row,dtype=float)
    f.close()
    
    f = open(folder + groupname + '_' + PROBE + 'err.csv')
    reader0 = csv.reader(f)
    for row in reader0:
        probeerr = np.asarray(row,dtype=float)
    f.close()
    
    f = open(folder + groupname + '_CLTA.csv')
    reader0 = csv.reader(f)
    for row in reader0:
        CLTA = np.asarray(row,dtype=float)
    f.close()
    
    f = open(folder + groupname + '_CLTAerr.csv')
    reader0 = csv.reader(f)
    for row in reader0:
        CLTAerr = np.asarray(row,dtype=float)
    f.close()
    
    return t, probe, probeerr, CLTA, CLTAerr

def Rescale(data,err,datamax):
    sf = datamax*1./max(data)
    scaled_data = sf * data
    scaled_err = sf * err
    
    return scaled_data, scaled_err

def StitchFiles(filenameout, filenamesin):
    content = None
    for filename in filenamesin:
        if content == None:
            content = Utility.load(filename)
        else:
            temp = Utility.load(filename)
            for j in range(len(content)):
                content[j] = content[j] + temp[j]
    Utility.save(content, filenameout)
    
def CleanParamList(filename, extra_params = []):
    """
    Delete derivative parameters and dynamic variables from full parameter list.
    """
    paramlist = ['pit','N3k','N4k','N5p','N3k4','kc','kcon','klon','CLTA','CLTAfree','switch0','switch1','PI','PI3P','PI4P','PI45P2','PI34P2','CLTAmax','klonb','kconb','a','turnover','probe','probeCLTA','probePI4P','probePI3P','probePI45P2','probePI34P2','probeCLTAPI4P','probeCLTAPI3P','probeCLTAPI45P2','probeCLTAPI34P2','circ','ku'] + extra_params
    
    var_vals = Utility.load(filename)
    for param in paramlist:
        try:
            var_vals.del_by_key(param)
        except:
            missing = param
    Utility.save(var_vals, filename)
