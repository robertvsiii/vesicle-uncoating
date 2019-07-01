from SloppyCell.ReactionNetworks import *

def SetPriors(model,priors):
    """Assigns priors to model, specified as python dict, with (mean, lognormal std)"""
    
    priors = KeyedList(priors)
    
    for id, val in priors:
        model.AddResidual(Residuals.PriorInLog('prior_on_%s' % id, id, scipy.log(val[0]), val[1]))