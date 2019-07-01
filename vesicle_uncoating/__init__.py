from SloppyCell.ReactionNetworks import *
from NetworkDef import net, MakeProbeNetwork, AddLipidAffinity
from DataImport import ImportData, MakeProbeModel, CleanParamList, StitchFiles, LoadData
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
from functools import partial
from matplotlib.backends import backend_pdf as bpdf
import copy

net_base = net.copy()
defaultfolder = '/Users/Bobby/Dropbox (MIT)/MIT/EnglandGroup/Membranes/Kirchhausen/Data for Bobby/ProbeData02012016/'
#net.set_var_vals({'CLTAmax_true' : 70})
#CLTA_scale = 70*(70./66.7)
CLTA_scale = 70.

class CoatedVesicle:
    """
    Creates network, experiment and model for given lipid probe and optimizable variables.
    
    priors = dictionary specifying prior probabilities. For example, {'ktest' : (3, 1)} would 
    make ktest an optimizable parameter, set 3 as the expected value, and allow the search
    to explore deviations from this value up to a factor of exp(1).
    """

    def __init__(self, priors = {}, folder = defaultfolder, 
                 filename_params = 'ProbeParams' + str(date.today()) + '.bp',
                 filename_ens = 'ProbeEns' + str(date.today()) + '.en', colors = ['red','blue','green'],
                 probechoice = 'All', control = {'FYVE':True,'DrrA':True,'TAPP1':True,'PH':False}, saved = False):
        
        if probechoice in ['All', None]:
            self.probelist = ['FYVE', 'DrrA', 'PH', 'TAPP1']
        elif probechoice == 'NoPH':
            self.probelist = ['FYVE', 'DrrA', 'TAPP1']
        elif probechoice in ['FYVE', 'DrrA', 'PH', 'TAPP1']:
            self.probelist = [probechoice]
        else:
            print 'Invalid probe choice.'
        self.folder = folder
        self.filename_params = folder + filename_params + '.bp'
        self.filename_ens = folder + filename_ens + '.en'
        self.filename_ens_traj = folder + '_traj_' + filename_ens
        self.cpunum = cpu_count()
        self.colors = colors
        
        self.priors = KeyedList(priors.copy())
        self.net_base = net_base.copy()
        self.net_probe = {}
        self.expt = {}
        self.probedata = {}
        self.CLTAdata = {}
        for probe in self.probelist:
            self.net_probe[probe] = MakeProbeNetwork(probe, self.net_base)
            self.expt[probe], self.probedata[probe], self.CLTAdata[probe] = ImportData(probe, self.folder, CLTA_scale, probechoice, control[probe])
        self.model, self.params, self.cost = MakeProbeModel(self.net_probe, self.expt, self.priors)

        if saved:
            self.params = Utility.load(self.filename_params)
            
        self.model.params.update(self.params)
        self.params = self.model.params
        self.cost = self.model.cost(self.params)
        
        self.jtj_log = None
        self.ens = None
        self.tmax = 110
        
 
            
    def AddLipidAffinity(self, new_probed_lipid):
        AddLipidAffinity(new_probed_lipid,self.net_probe)
        
    def GetTrajectory(self, probe):
        self.net_probe[probe].update_optimizable_vars(self.params)
        #self.net_probe[probe].set_var_vals(self.params)
        traj = Dynamics.integrate(self.net_probe[probe], [0, self.tmax])
        return traj.get_times(),traj.get_var_traj('CLTA'),traj.get_var_traj(self.net_probe[probe].probed_lipid),traj.get_var_traj(probe)
        
    def GetBestFit(self, method = 'Nelder-Mead-log', tol = 1e-3):
        if method == 'Powell-log':
            self.params = Optimization.fmin_powell_log_params(self.model, self.params)
        elif method == 'Nelder-Mead-log':
            self.params = Optimization.fmin_log_params(self.model, self.params, xtol = tol)
        #elif method == 'Nelder-Mead':
        #    self.params = Optimization.fmin(self.model, self.params)
        elif method == 'Levenberg-Marquadt-log':
            self.params = Optimization.fmin_lm_log_params(self.model, self.params)
        elif method == 'Least-Squares-log':
            self.params = Optimization.leastsq_log_params(self.model, self.params)
        else:
            print 'Invalid fit method.'
        Utility.save(self.params, self.filename_params)
        
        
    def GetEnsemble(self, maxtime, skip = 0, temperature = 1, plotting = False, load_saved = False, sort = False):
        if load_saved:
            self.ens, self.F, self.ratio = Utility.load(self.filename_ens)
            
        else:
            ensemble_partial = partial(Ensembles.ensemble_log_params, params = self.params, max_run_hours = maxtime,
                                      skip_elems = skip, temperature = temperature)
            models = []
            for n in range(self.cpunum):
                models.append(self.model.copy())
            
            pool = Pool(self.cpunum)
            results = pool.map(ensemble_partial, models)
            pool.close()
            
            self.ens = []
            self.F = []
            self.ratio = []
            for res in results:
                self.ens = self.ens + res[0]
                self.F = self.F + res[1]
                self.ratio.append(res[2])
            Utility.save([self.ens, self.F, self.ratio], self.filename_ens)
        
        if sort:
            temp = dict((self.F[j],self.ens[j]) for j in range(len(self.F)))
            self.ens = [value for (key,value) in sorted(temp.items())]
            self.F = [key for (key,value) in sorted(temp.items())]
        
        if plotting:
            self.F_autocorr = Ensembles.autocorrelation(self.F)[:round(len(self.F)/(2.*self.cpunum))]
            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax1.plot(self.F)
            ax1.set_xlabel('Steps')
            ax1.set_ylabel('Cost')
            
            cutoff = 0.1
            ax2 = fig.add_subplot(212)
            ax2.plot(self.F_autocorr)
            ax2.set_xlabel('Lag')
            ax2.set_ylabel('Cost Autocorrelation')
    
    def PlotModel(self, datachoice = None):
        if datachoice == None:
            datachoice = [self.probelist[0]]
        self.model.CalculateForAllDataPoints(self.params)
        Plotting.figure(1)
        lines, labels = Plotting.plot_model_results(self.model, loc='upper left', data_to_plot = datachoice)
        return lines, labels
        
    def PlotSensitivity(self, variables, n = [100, 100], var_lims = [[0.1,10],[0.1,10]], levels = np.arange(2,50,2),
                       Hessian = True, log = True, vlim = [None, None], flip = False):
        fig1 = plt.figure(figsize = (5,5))
        ax1 = fig1.add_subplot(111)
        
        variable0 = variables[0]
        variable1 = variables[1]
        idx0 = self.model.params.index_by_key(variable0)
        idx1 = self.model.params.index_by_key(variable1)
        if idx0>idx1:
            variables[0] = variable1
            variables[1] = variable0
        
        jtj2_log = self.IntegrateHessian(variables)
        self.sensitivity_variables = variables
        
        
        
        if Hessian:
            if log:
                var_vals = [np.log(self.params.get(variables[0])), np.log(self.params.get(variables[1]))]
                var_ranges = [np.log(var_lims[0][1]) - np.log(var_lims[0][0]), np.log(var_lims[1][1]) - np.log(var_lims[1][0])]
                dtheta = np.asarray(var_ranges)*1./np.asarray(n)
                grd = np.mgrid[0:n[0],0:n[1]]
                theta = [grd[0]*dtheta[0] + np.log(var_lims[0][0]), grd[1]*dtheta[1] + np.log(var_lims[1][0])]
                Deltalogtheta = [theta[0] - var_vals[0], theta[1] - var_vals[1]]
                self.jtj2_log = jtj2_log
                C = 0.5*jtj2_log[0,0]*Deltalogtheta[0]**2 + jtj2_log[0,1]*Deltalogtheta[0]*Deltalogtheta[1] + 0.5*jtj2_log[1,1]*Deltalogtheta[1]**2
                for j in range(2):
                    theta[j] = theta[j]*np.log10(np.exp(1))
                    var_lims[j] = [np.log10(var_lims[j][0]),np.log10(var_lims[j][1])]
                if self.ens != None:
                    theta_ens = [[],[]]
                    cost = []
                    for params in self.ens:
                        theta_ens[0].append(np.log10(params.get(variables[0])))
                        theta_ens[1].append(np.log10(params.get(variables[1])))
                        scat = ax1.scatter(theta_ens[0],theta_ens[1],c=self.F,cmap='hot', 
                                           vmin = vlim[0], vmax = vlim[1])
                    cbar = fig1.colorbar(scat)
                    cbar.ax.set_ylabel('Cost')

            else:
                var_vals = [self.params.get(variables[0]), self.params.get(variables[1])]
                var_ranges = [var_lims[0][1] - var_lims[0][0], var_lims[1][1] - var_lims[1][0]]
                dtheta = np.asarray(var_ranges)*1./np.asarray(n)
                grd = np.mgrid[0:n[0],0:n[1]]
                theta = [grd[0]*dtheta[0] + var_lims[0][0], grd[1]*dtheta[1] + var_lims[1][0]]
                Deltalogtheta = [np.log(theta[0]) - np.log(var_vals[0]), np.log(theta[1]) - np.log(var_vals[1])]
                self.jtj2_log = jtj2_log
                C = 0.5*jtj2_log[0,0]*Deltalogtheta[0]**2 + jtj2_log[0,1]*Deltalogtheta[0]*Deltalogtheta[1] + 0.5*jtj2_log[1,1]*Deltalogtheta[1]**2
                
            
                if self.ens != None:
                    theta_ens = [[],[]]
                    cost = []
                    for params in self.ens:
                        theta_ens[0].append(params.get(variables[0]))
                        theta_ens[1].append(params.get(variables[1]))
                        scat = ax1.scatter(theta_ens[0],theta_ens[1],c=self.F,cmap='hot',vmin = vlim[0], vmax = vlim[1])
                    cbar = fig1.colorbar(scat)
                    cbar.ax.set_ylabel('Cost')
                    
        self.contours = ax1.contour(theta[0], theta[1], np.exp(-C), levels, cmap = 'hot')
        cbar0 = fig1.colorbar(self.contours)
        cbar0.ax.set_ylabel(r'$p(\log$'+variables[0] + ', $\log$ '+variables[1] + ')')
        cbar0.set_ticks(())
        ax1.set_xlim(var_lims[0])
        ax1.set_ylim(var_lims[1])
        
        if flip:
            ax1.set_ylim(ax1.get_ylim()[::-1])
        if log:
            ax1.set_xlabel(r'$\log$ '+variables[0])
            ax1.set_ylabel(r'$\log$ '+variables[1])
        else:
            ax1.set_xlabel(variables[0])
            ax1.set_ylabel(variables[1])
        
            #SAVE FIGURE
        pdf = bpdf.PdfPages(self.filename_ens + '_distribution.pdf')
        pdf.savefig(fig1)
        pdf.close()
        
        

    def IntegrateHessian(self, variables):
        #NOTE: a simpler way to do this would be to just invert the matrix, remove the relevant rows/columns,
        #then invert again. I tested this, and it gives the same answer.
        if self.jtj_log == None:
            self.j_log, self.jtj_log = self.model.GetJandJtJInLogParameters(np.log(self.params))
        n_deleted = 0
        variables_idx = []
        jtj2_log = self.jtj_log.copy()
        
        for item in variables:
            variables_idx.append(self.model.params.index_by_key(item))

        for var0 in range(np.size(self.jtj_log,0)):
            if var0 not in variables_idx:
                var = var0 - n_deleted
                a_log = jtj2_log[var,var]
                b_log = jtj2_log[:,var]
                b2_log = np.meshgrid(b_log,b_log)
                b2_log = b2_log[0]*b2_log[0].T
                jtj2_log = jtj2_log - b2_log/a_log
                jtj2_log = np.delete(jtj2_log, var, axis = 0)
                jtj2_log = np.delete(jtj2_log, var, axis = 1)
    
                n_deleted += 1
        
        return jtj2_log
    
    def SensitivityTable(self, ensemble = False, finite_diff = False, integrate = True):
        if self.jtj_log == None:
            self.j_log, self.jtj_log = self.model.GetJandJtJInLogParameters(np.log(self.model.params))
        self.cov = np.linalg.inv(self.jtj_log)
        if finite_diff:
            eps = []
            for k in range(len(self.model.params)):
                eps.append(0.1*np.sqrt(self.cov[k,k]))
            self.hess_log = self.model.hessian_log_params(self.model.params,eps)
            self.cov = np.linalg.inv(self.hess_log)
        
        self.ParamTable = dict()
        if ensemble:
            for variable in self.params.keys():
                log_values = []
                value = self.params.get(variable)
                for element in self.ens:
                    log_values.append(np.log(element.get(variable)))
                log_values = np.asarray(log_values)
                log_mean = np.mean(log_values)
                log_std = np.std(log_values)
                value_min = np.exp(log_mean-2*log_std)
                value_max = np.exp(log_mean+2*log_std)
                
                if variable in ['NPI3P', 'NPI34P2']:
                    value_min = 0
                    value_max = np.exp(np.max(log_values))
                    
                self.ParamTable[variable] = [value,(value_min, value_max)]
                if variable == 'ku0':
                    value_kex = self.Compute_kex(value)
                    value_min_kex = self.Compute_kex(value_min)
                    value_max_kex = self.Compute_kex(value_max)
                    self.ParamTable['kex'] = [value_kex,(value_min_kex, value_max_kex)]

                
        else:
            if integrate:
                for variable in self.model.params.keys():
                    value = self.model.params.get(variable)
                    idx = self.model.params.index_by_key(variable)
                    log_std = np.sqrt(self.cov[idx,idx])
                    value_min = np.exp(np.log(value)-2*log_std)
                    value_max = np.exp(np.log(value)+2*log_std)
                    self.ParamTable[variable] = [value,(value_min, value_max)]
                    if variable == 'ku0':
                        value_kex = self.Compute_kex(value)
                        value_min_kex = self.Compute_kex(value_min)
                        value_max_kex = self.Compute_kex(value_max)
                        self.ParamTable['kex'] = [value_kex,(value_min_kex, value_max_kex)]
                    if variable in ['N3k0','N4k0','N3k40','N5p0']:
                        kcat = self.model.get_calcs()[0].get_var_val('kcat'+variable[1:-1])
                        value = value*kcat
                        value_min = value_min*kcat
                        value_max = value_max*kcat
                        self.ParamTable[variable+'kcat'] = [value,(value_min,value_max)]
            else:
                for variable in self.model.params.keys():
                    value = self.model.params.get(variable)
                    idx = self.model.params.index_by_key(variable)
                    log_std = np.sqrt(1/self.jtj_log[idx,idx])
                    value_min = np.exp(np.log(value)-2*log_std)
                    value_max = np.exp(np.log(value)+2*log_std)
                    self.ParamTable[variable] = [value,(value_min, value_max)]
                    if variable == 'ku0':
                        value_kex = self.Compute_kex(value)
                        value_min_kex = self.Compute_kex(value_min)
                        value_max_kex = self.Compute_kex(value_max)
                        self.ParamTable['kex'] = [value_kex,(value_min_kex, value_max_kex)]
                    if variable in ['N3k0','N4k0','N3k40','N5p0']:
                        kcat = self.model.get_calcs()[0].get_var_val('kcat'+variable[1:-1])
                        value = value*kcat
                        value_min = value_min*kcat
                        value_max = value_max*kcat
                        self.ParamTable[variable+'kcat'] = [value,(value_min,value_max)]
                    
        return self.ParamTable
    
    def Compute_kex(self, ku0):
        knet = self.model.get_calcs()[0].get_var_val('knet')
        CLTAmax_true = self.model.get_calcs()[0].get_var_val('CLTAmax_true')
        fleak = self.model.get_calcs()[0].get_var_val('fleak')
        kc = ku0 + knet
        a = (fleak*CLTAmax_true*ku0/knet)**2
        CLTAmax = CLTAmax_true*(1 - np.sqrt(np.sqrt((1-4*a)**2)))/(2*a)
        kex = kc * np.sqrt(np.sqrt(((CLTAmax_true/CLTAmax) - (CLTAmax_true/CLTAmax)**2)**2)) / CLTAmax_true
        return kex
    
def IntegrateEnsemble(vesicle, range_dict = {}):
    variables = range_dict.keys()
    
    if range_dict != {}:
        ens_dict = {}
        for name in variables:
            ens_dict[name] = np.linspace(range_dict[name][0],range_dict[name][1],range_dict[name][2])
        vesicle.ens = []
        for j in range(len(ens_dict[variables[0]])):
            vesicle.ens.append(vesicle.params.copy())
            for name in variables:
                vesicle.ens[-1].set(name,ens_dict[name][j])
        MakeTrajPartial = partial(MakeTrajFromParams, variables = variables, model = vesicle.model.copy(), 
                                  tmax = vesicle.tmax)
        pool = Pool()
        results = pool.map(MakeTrajPartial,vesicle.ens)
        pool.close()
    else:
        results = MakeTrajFromParams(vesicle.params.copy(), model = vesicle.model.copy(), tmax = vesicle.tmax)
        results = [results]
        
    Utility.save(results,vesicle.filename_ens_traj)
    
    return results
    
def MakeTrajFromParams(params, variables = None, model = None, tmax = None):
    if model != None and tmax != None:
        model.params.update(params)
        params = model.params
        if variables != None:       
            var_vals = {}
            for variable in variables:
                model.set_var_optimizable(variable, False)
                var_vals[variable]=params.get(variable)
                params.del_by_key(variable)
                model.residuals.del_by_key('prior_on_'+variable)
                for net_probe in model.get_calcs().values():
                    net_probe.set_var_val(variable,var_vals[variable])
            model.params.update(params)
            params = model.params
            params = Optimization.fmin_log_params(model, params, xtol = 1e-3)
            for variable in variables:
                params.set(variable,var_vals[variable])
                
        results = {}
        #Load new parameter values into all nets
        for net_probe in model.get_calcs().values():
            for name in params.keys():
                try:
                    net_probe.set_var_val(name,params.get(name))
                except:
                    e = 1

        scale_factors = model.GetScaleFactors()
        
        #Integrate trajectories and save results
        for net_probe in model.get_calcs().values():
            probed_lipid = net_probe.probed_lipid
            probename = net_probe.probe
            netname = net_probe.id
            results[netname] = {}
            
            #Make network with no probe binding to the relevant lipid, to model wild-type cell
            probes = ['DrrA','FYVE','TAPP1']
            if probename == 'FYVE':
                net_lipid = MakeProbeNetwork('DrrA',net_base) 
            else:
                net_lipid = MakeProbeNetwork('FYVE',net_base) 
                
                
            for name in params.keys():
                try:
                    net_lipid.set_var_val(name,params.get(name))
                except:
                    e = 1
            
            traj_lipid = Dynamics.integrate(net_lipid, [0, tmax+10])
            traj = Dynamics.integrate(net_probe, [0, tmax+10])
            
            results[netname]['t_probe'] = traj.get_times()
            results[netname]['t_lipid'] = traj_lipid.get_times()
            results[netname]['lipid'] = traj_lipid.get_var_traj(probed_lipid)
            results[netname]['probe'] = traj.get_var_traj(probename)*scale_factors[probename][probename]
        
        
        results['params'] = params
        results['scale_factors'] = scale_factors
        return results
    else:
        print 'Missing keyword argument in MakeTrajFromParams'
        return None
    
def ConvertToKD(params,Nf = 1):
    KD_dict = {
        'Cytosolic concentration (uM)' : Nf,
        'Clathrin' : params.get('kcoff')*210*Nf/params.get('kcon0'),
    }

    for lipid in ['PI45P2','PI4P']:
        KD_dict[lipid] = params.get('kloff'+lipid)*Nf/params.get('klon')
    for lipid in ['PI3P','PI34P2']:
        KD_dict[lipid] = params.get('kloff'+lipid)**2 *Nf/(params.get('sl')*params.get('klon')**2)/(1 + params.get('kloff'+lipid)/(params.get('sl')*params.get('klon')))
    return KD_dict