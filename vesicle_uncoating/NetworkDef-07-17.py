from SloppyCell.ReactionNetworks import *

#INITIATE NETWORK
net = Network('net')
net.add_compartment('pit',initial_size = 1.0)

#DEFINE PARAMETERS

#Initial Numbers
net.add_parameter('NC', 5, is_optimizable = False)
net.add_parameter('NPI4P', 10, is_optimizable = False)
net.add_parameter('NPI3P', 0.01, is_optimizable = False)
net.add_parameter('NPI34P2', 0.01, is_optimizable = False)
net.add_parameter('NPI45P2', 100, is_optimizable = False)
net.add_parameter('NPI', 100, is_optimizable = False)

#Enzyme parameters

net.add_parameter('N3k0', 0.000001, is_optimizable = False)
net.add_parameter('kcat3k',50, is_optimizable = False)
net.add_parameter('KM3k',4, is_optimizable = False)
net.add_parameter('N4k0', 0.000001, is_optimizable = False)
net.add_parameter('kcat4k',50, is_optimizable = False)
net.add_parameter('KM4k',4, is_optimizable = False)
net.add_parameter('N5p0', 0.000001, is_optimizable = False)
net.add_parameter('kcat5p',50, is_optimizable = False)
net.add_parameter('KM5p',4, is_optimizable = False)
net.add_parameter('N3k40', 0.000001, is_optimizable = False)
net.add_parameter('kcat3k4',50, is_optimizable = False)
net.add_parameter('KM3k4',4, is_optimizable = False)

net.add_species('N3k','pit', initial_conc = 'N3k0')
net.add_species('N4k','pit', initial_conc = 'N4k0')
net.add_species('N5p','pit', initial_conc = 'N5p0')
net.add_species('N3k4','pit', initial_conc = 'N3k40')



#Clathrin parameters
net.add_parameter('CLTAmax_true', 70, is_optimizable = False)
net.add_species('CLTAmax', 70)
net.add_parameter('knet',7, is_optimizable = False) #net forward rate parameter before scission
net.add_parameter('kc',107)
net.add_parameter('kcfinal',20, is_optimizable = False) #on-rate after uncoating
net.add_parameter('ku0',100, is_optimizable = False) #off- rate before scission
net.add_assignment_rule('kc', 'ku0 + knet') #on-rate before scission
net.add_parameter('fleak',0.0003, is_optimizable = False) #rate of losing first trimer of complete coat
net.add_parameter('ku1',0, is_optimizable = False) #off-rate after scission
net.add_parameter('ku2', 0.5, is_optimizable = False) #off-rate with Hsc70
#Scission time (replace with dynamin binding later)
net.add_parameter('tsciss', 60)
net.add_parameter('delay', 10, is_optimizable = False)
net.add_parameter('tdst', 70, is_optimizable = False)
net.add_assignment_rule('tsciss', 'tdst-delay')
#Reduce dimensionality by constraining clathrin parameters
#See "Dimensional Reduction" section in Lipid Probe Plots notebook
net.add_parameter('a',0)
net.add_parameter('turnover',0)
net.add_assignment_rule('a', '(fleak*CLTAmax_true*ku0/knet)**2')
net.add_assignment_rule('CLTAmax', 'CLTAmax_true*(1 - sqrt(sqrt((1-4*a)**2)))/(2*a)')
net.add_assignment_rule('turnover', 'kc * sqrt(sqrt(((CLTAmax_true/CLTAmax) - (CLTAmax_true/CLTAmax)**2)**2)) / CLTAmax_true')




#Binding of clathrin-binding domain to clathrin
#Free protein binding
net.add_parameter('kcon0', 0.2, is_optimizable = False)
net.add_parameter('kcon', 0)
net.add_assignment_rule('kcon', 'kcon0/210')
net.add_parameter('kcoff', 0.005, is_optimizable = False)

net.add_parameter('s', 100, is_optimizable = False)


#on-rate for all lipids
net.add_parameter('klon', 0.1, is_optimizable = False)


#ADD SPECIES
net.add_species('CLTA','pit',initial_conc = 0)
net.add_species('CLTAfree','pit',initial_conc = 'NC')
net.add_species('switch0','pit',initial_conc = 0)
net.add_species('switch1','pit',initial_conc = 0)

net.add_species('PI','pit',initial_conc = 'NPI')
net.add_species('PI3P','pit',initial_conc = 'NPI3P')
net.add_species('PI4P','pit',initial_conc = 'NPI4P')
net.add_species('PI45P2','pit',initial_conc = 'NPI45P2')
net.add_species('PI34P2','pit',initial_conc = 'NPI34P2')


#DEFINE REACTIONS FOR LIPID SWITCH
net.addReaction(Reactions.MichaelisMentenReaction, 'K3', S = 'PI', P = 'PI3P', E = 'N3k', k = 'kcat3k * switch0',  
                Km = 'KM3k')
net.addReaction(Reactions.MichaelisMentenReaction, 'K4', S = 'PI3P', P = 'PI34P2', E = 'N4k', k = 'kcat4k * switch0',  
                Km = 'KM4k')
net.addReaction(Reactions.MichaelisMentenReaction, 'P5', S = 'PI45P2', P = 'PI4P', E = 'N5p', k = 'kcat5p * switch0',  
                Km = 'KM5p')
net.addReaction(Reactions.MichaelisMentenReaction, 'K34', S = 'PI4P', P = 'PI34P2', E = 'N3k4', k = 'kcat3k4 * switch0',  
                Km = 'KM3k4')

#DEFINE SCISSION AND UNCOATING EVENTS
net.add_event('scission', 'gt(time, tsciss)', {'switch0': '1'})
net.add_event('destabilization', 'gt(time, tdst)', {'switch1': '1'})


def MakeProbeNetwork(probe, net_base):
           

    if probe == 'FYVE':
        probed_lipid = 'PI3P'
    elif probe == 'DrrA':
        probed_lipid = 'PI4P'
    elif probe == 'PH':
        probed_lipid = 'PI45P2'
    elif probe == 'TAPP1':
        probed_lipid = 'PI34P2'
    elif probe in ['Aux1', 'GAK']:
        probed_lipid = 'PI3P'
            
    net_probe = net_base.copy('net' + probe)
    
    net_probe.probed_lipid = probed_lipid
    net_probe.probe = probe
    net_probe.new_probed_lipids = []
    #Number of free probe molecules, starting concentration in coat, and true max height of CLTA trace
    net_probe.add_parameter('N' + probe, 1, is_optimizable = False)
    
    #Binding of lipid-binding domain
    net_probe.add_parameter('kloff' + probed_lipid, 2, is_optimizable = False)
    
    #Cooperativity factor
    net_probe.add_parameter('s'+probe, 100)
    net_probe.add_assignment_rule('s'+probe,'s')
    net_probe.add_parameter('kconb'+probe,0)
    net_probe.add_assignment_rule('kconb'+probe, 's'+probe+'*kcon')
    net_probe.add_parameter('klonb'+probe, 0)
    net_probe.add_assignment_rule('klonb'+probe, 's'+probe+'*klon')
    
    #Add probe species
    net_probe.add_species(probe,'pit',initial_conc = 0)
    net_probe.add_species(probe + 'CLTA','pit',initial_conc = 0)
    net_probe.add_species(probe + probed_lipid,'pit',initial_conc = 0)
    net_probe.add_species(probe +'CLTA' + probed_lipid,'pit',initial_conc = 0)
    
    #Update assignment rule
    net_probe.add_assignment_rule('CLTA','(CLTAfree + ' + probe + 'CLTA + ' + probe + 'CLTA' + probed_lipid + ')/3')
    net_probe.add_assignment_rule(probe, probe + 'CLTA + ' + probe + probed_lipid + ' + ' + probe + 'CLTA' + probed_lipid)
    
    #DEFINE REACTIONS FOR ADDING AND REMOVING CLATHRIN
    #('switch' parameters say whether vesicle has pinched off and whether uncoating threshold has been reached)
    net_probe.add_species('circ'+probe, 0)
    net_probe.add_species('ku'+probe, 10)
    net_probe.add_assignment_rule('circ'+probe,'sqrt(sqrt(((CLTA/CLTAmax) - (CLTA/CLTAmax)**2)**2))')
    net_probe.add_assignment_rule('ku'+probe,'((1 - switch1) * ((1 - switch0) * ku0 + switch0 * ku1) * (circ'+probe+' + fleak * CLTA) / (CLTA + 0.01)) + switch1 * ku2')

    net_probe.addReaction(Reactions.ProductionReaction, 'add', product = 'CLTAfree',
                          rate = '3 * (kc * (1 - switch0) ) * circ' + probe + '+ 3 * kcfinal * switch1')
    net_probe.addReaction(Reactions.ExponentialDecayReaction, 'removefree', species = 'CLTAfree',
                          rate = 'ku'+probe)
    net_probe.addReaction(Reactions.ExponentialDecayReaction, 'removeprobeCLTA', species = probe + 'CLTA',
                          rate = 'ku'+probe)
    net_probe.addReaction(Reactions.TransformationReaction, 'removeprobeCLTAlipid', 
                          old = probe + 'CLTA' + probed_lipid, new = probe + probed_lipid, rate = 'ku'+probe)

    #DEFINE REACTIONS FOR p BINDING
    net_probe.addReaction(Reactions.TransformationReaction,'bindpc',
                   new = probe + 'CLTA', old = 'CLTAfree', rate = 'kcon * N' + probe)
    net_probe.addReaction(Reactions.TransformationReaction,'bindpl',
                   new = probe + probed_lipid, old = probed_lipid, rate = 'klon * N' + probe)
    net_probe.addReaction(Reactions.HeterodimerizationReaction,'bindplc',
                   dimer = probe + 'CLTA' + probed_lipid, A = probe + 'CLTA', B = probed_lipid, rate = 'klonb'+probe)
    net_probe.addReaction(Reactions.HeterodimerizationReaction,'bindpcl',
                   dimer = probe + 'CLTA' + probed_lipid, A = probe + probed_lipid, B = 'CLTAfree', rate = 'kconb'+probe)
    net_probe.addReaction(Reactions.TransformationReaction,'dispc',
                   old = probe + 'CLTA', new = 'CLTAfree', rate = 'kcoff')
    net_probe.addReaction(Reactions.TransformationReaction,'displ',
                   old = probe + probed_lipid, new = probed_lipid, rate = 'kloff' + probed_lipid)
    net_probe.addReaction(Reactions.HeterodimerDissociationReaction,'displc',
                   dimer = probe + 'CLTA' + probed_lipid, A = probe + 'CLTA', B = probed_lipid, rate = 'kloff' + probed_lipid)
    net_probe.addReaction(Reactions.HeterodimerDissociationReaction,'dispcl',
                   dimer = probe + 'CLTA' + probed_lipid, A = probe + probed_lipid, B = 'CLTAfree', rate = 'kcoff')
    
    if probe in ['Aux1', 'GAK']:
        AddLipidAffinity('PI4P', net_probe)
        AddLipidAffinity('PI45P2', net_probe)
    
    return net_probe

def AddLipidAffinity(new_probed_lipid, net_probe):
    probe = net_probe.probe
    
    net_probe.new_probed_lipids.append(new_probed_lipid)
    
    net_probe.add_parameter('kloff' + new_probed_lipid, 20, is_optimizable = False)
    
    net_probe.add_species(probe + new_probed_lipid,'pit',initial_conc = 0)
    net_probe.add_species(probe + 'CLTA' + new_probed_lipid,'pit',initial_conc = 0)
    
    #Update assignment rule
    CLTAstring = '(CLTAfree + ' + probe + 'CLTA' + ' + ' + probe + 'CLTA' + net_probe.probed_lipid
    probestring = probe + 'CLTA + ' + probe + net_probe.probed_lipid + ' + ' + probe + 'CLTA' + net_probe.probed_lipid
    for lipid in net_probe.new_probed_lipids:
        CLTAstring = CLTAstring + ' + ' + probe + 'CLTA' + lipid
        probestring = probestring + ' + ' + probe + lipid + ' + ' + probe + 'CLTA' + lipid
    CLTAstring = CLTAstring + ')/3'
    
    net_probe.add_assignment_rule('CLTA',CLTAstring)
    net_probe.add_assignment_rule(probe,probestring)
    
    net_probe.addReaction(Reactions.TransformationReaction, 'removeprobeCLTAnewlipid' + new_probed_lipid, 
                    old = probe + 'CLTA' + new_probed_lipid, new = probe + new_probed_lipid, rate = 'ku'+probe)

    #DEFINE REACTIONS FOR p BINDING
    net_probe.addReaction(Reactions.TransformationReaction,'bindplnew' + new_probed_lipid,
                   new = probe + new_probed_lipid, old = new_probed_lipid, rate = 'klon * N' + probe)
    net_probe.addReaction(Reactions.HeterodimerizationReaction,'bindplcnew' + new_probed_lipid,
                   dimer = probe + 'CLTA' + new_probed_lipid, A = probe + 'CLTA', B = new_probed_lipid, rate = 'klonb'+probe)
    net_probe.addReaction(Reactions.HeterodimerizationReaction,'bindpclnew' + new_probed_lipid,
                   dimer = probe + 'CLTA' + new_probed_lipid, A = probe + new_probed_lipid, B = 'CLTAfree', rate = 'kconb'+probe)
    net_probe.addReaction(Reactions.TransformationReaction,'displnew' + new_probed_lipid,
                   old = probe + new_probed_lipid, new = new_probed_lipid, rate = 'kloff' + new_probed_lipid)
    net_probe.addReaction(Reactions.HeterodimerDissociationReaction,'displcnew' + new_probed_lipid,
                   dimer = probe + 'CLTA' + new_probed_lipid, A = probe + 'CLTA', B = new_probed_lipid, rate = 'kloff' + new_probed_lipid)
    net_probe.addReaction(Reactions.HeterodimerDissociationReaction,'dispclnew' + new_probed_lipid,
                   dimer = probe + 'CLTA' + new_probed_lipid, A = probe + new_probed_lipid, B = 'CLTAfree', rate = 'kcoff')