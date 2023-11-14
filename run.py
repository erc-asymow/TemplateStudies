import argparse
import os
import sys
import copy
import math

parser = argparse.ArgumentParser(description='run')

parser.add_argument('--none', action='store_true'  , help = 'none')
parser.add_argument('--dryrun', action='store_true'  , help = 'dry run')
parser.add_argument('--smear', action='store_true'  , help = 'smear')
parser.add_argument('--tag',  default='test'       , help = 'tag')
parser.add_argument('--post_tag',  default='test'  , help = 'post tag')
parser.add_argument('--algo',   default='all'        , help = 'algo')
parser.add_argument('--run',   default='full'        , help = 'run')
parser.add_argument('--nevents', dest = 'nevents'  , type = int,  default=10000000, help='number of events')
parser.add_argument('--rebinX', dest = 'rebinX'  , type = int,  default=-1, help='')
parser.add_argument('--rebinY', dest = 'rebinY'  , type = int,  default=-1, help='')
parser.add_argument('--scale0', action='store_true'  , help = 'scale A0(x,y)')
parser.add_argument('--scale1', action='store_true'  , help = 'scale A1(x,y)')
parser.add_argument('--scale2', action='store_true'  , help = 'scale A2(x,y)')
parser.add_argument('--scale3', action='store_true'  , help = 'scale A3(x,y)')
parser.add_argument('--scale4', action='store_true'  , help = 'scale A4(x,y)')
parser.add_argument('--jacmass', dest = 'jacmass'  , type = int,  default=-1, help='')

args = parser.parse_args()

sample_sizes = {
    '10M':     10000000,
    '100M':   100000000,
    '1G':    1000000000,
    #'4G':    4000000000,
    #'10G':  10000000000,
    #'40G':  40000000000,
    ##'200G': 200000000000,
}

fit_opts = {
    'ULA0A1A2A3A4' :       ' --jUL --j0 --j1 --j2 --j3 --j4',
    'ADDMC_ULA0A1A2A3A4' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert',
    #'FBB_ULA0A1A2A3A4'         : ' --jUL --j0 --j1 --j2 --j3 --j4 --compute_deltachi2 --with_offset',
    #'ADDMCFBB_ULA0A1A2A3A4'    : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --compute_deltachi2 --with_offset',
    #'ADDMCJACEXT_ULA0A1A2A3A4' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --jacobians_from_external',
    #'ADDMCHMCEXT_ULA0A1A2A3A4' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --hMC_from_external',
    ##'ADDMCHMCEXTSCALED_ULA0A1A2A3A4' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --hMC_from_external',
    #'ADDMCMASSUNCERT_ULA0A1A2A3A4' : ' --jUL --j0 --j1 --j2 --j3 --j4 --hMC_from_external --jacobians_from_external --add_mass_uncert=1.0',
    #'ADDMC_ULA0A1A2A4':    ' --jUL --j0 --j1 --j2 --j4 --add_MC_uncert',
    #'ULA0A1A2A4' :         ' --jUL --j0 --j1 --j2 --j4',
    #'ADDMC_ULA0A1A2A3':    ' --jUL --j0 --j1 --j2 --j3 --add_MC_uncert',
    #'ULA0A1A2A3' :         ' --jUL --j0 --j1 --j2 --j3',
    #'ADDMC_ULA0A1A2A3A4_Ymin25' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=25',
    #'ADDMC_ULA0A1A2A3A4_Ymin26' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=26',
    #'ADDMC_ULA0A1A2A3A4_Ymin27' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=27',
    #'ADDMC_ULA0A1A2A3A4_Ymin28' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=28',
    #'ADDMC_ULA0A1A2A3A4_Ymin29' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=29',
    #'ADDMC_ULA0A1A2A3A4_Ymin30' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=30',
    #'ADDMC_ULA0A1A2A3A4_Ymin31' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=31',
    #'ADDMC_ULA0A1A2A3A4_Ymin32' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=32',
    #'ADDMC_ULA0A1A2A3A4_Ymin33' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=33',
    #'ADDMC_ULA0A1A2A3A4_Ymin34' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=34',
    #'ADDMC_ULA0A1A2A3A4_Ymin35' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_min=35',
    #'ADDMC_ULA0A1A2A3A4_Ymax72' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=72',
    #'ADDMC_ULA0A1A2A3A4_Ymax70' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=70',
    #'ADDMC_ULA0A1A2A3A4_Ymax68' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=68',
    #'ADDMC_ULA0A1A2A3A4_Ymax66' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=66',
    #'ADDMC_ULA0A1A2A3A4_Ymax64' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=64',
    #'ADDMC_ULA0A1A2A3A4_Ymax62' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=62',
    #'ADDMC_ULA0A1A2A3A4_Ymax60' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=60',
    #'ADDMC_ULA0A1A2A3A4_Ymax58' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=58',
    #'ADDMC_ULA0A1A2A3A4_Ymax56' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=56',
    #'ADDMC_ULA0A1A2A3A4_Ymax54' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=54',
    #'ADDMC_ULA0A1A2A3A4_Ymax52' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --Y_max=52',
    #'ADDMC_ULA0A1A2A3A4_Xmax2p6' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=2.6',
    #'ADDMC_ULA0A1A2A3A4_Xmax2p4' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=2.4',
    #'ADDMC_ULA0A1A2A3A4_Xmax2p2' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=2.2',
    #'ADDMC_ULA0A1A2A3A4_Xmax2p0' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=2.0',
    #'ADDMC_ULA0A1A2A3A4_Xmax1p8' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=1.8',
    #'ADDMC_ULA0A1A2A3A4_Xmax1p5' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=1.5',
    #'ADDMC_ULA0A1A2A3A4_Xmax1p0' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=1.0',
    #'ADDMC_ULA0A1A2A3A4_Xmax0p8' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=0.8',
    #'ADDMC_ULA0A1A2A3A4_Xmax0p5' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=0.5',
    #'ADDMC_ULA0A1A2A3A4_Xmax0p3' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=0.3',
    #'ADDMC_ULA0A1A2A3A4_Xmax0p2' : ' --jUL --j0 --j1 --j2 --j3 --j4 --add_MC_uncert --X_max=0.2',
}

pol_default = {
    'run'    : "full",
    'corr'   : [10,4],
    'A0'     : [3,2],
    'A1'     : [2,3],
    'A2'     : [3,2],
    'A3'     : [3,4],
    'A4'     : [3,3]
}

pol_systs = []

#pol_systs.append( pol_default )

pol_syst = copy.deepcopy(pol_default)
pol_syst['A0'] = [2,2]
pol_syst['A2'] = [2,2]
#pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['A3'] = [2,4]
#pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['A3'] = [3,2]
#pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['A4'] = [1,3]
#pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['corr'] = [14,4]
#pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['corr'] = [12,4]
#pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['corr'] = [8,4]
#pol_systs.append( pol_syst )

pol_syst= copy.deepcopy(pol_default)
pol_syst['corr'] = [10,6]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [8,6]
pol_syst['A0']   = [3,2]
pol_syst['A1']   = [3,2]
pol_syst['A2']   = [3,2]
pol_syst['A3']   = [3,4]
pol_syst['A4']   = [4,8]
pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [5,4] #[8,6]
pol_syst['A0']   = [2,2] #[3,2]
pol_syst['A1']   = [2,2] #[3,2]
pol_syst['A2']   = [2,2] #[3,2]
pol_syst['A3']   = [2,4] #[3,4]
pol_syst['A4']   = [2,6] #[4,8]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [5,4] #[5,4] #[8,6]
pol_syst['A0']   = [1,2] #[2,2] #[3,2]
pol_syst['A1']   = [1,2] #[2,2] #[3,2]
pol_syst['A2']   = [1,2] #[2,2] #[3,2]
pol_syst['A3']   = [1,4] #[2,4] #[3,4]
pol_syst['A4']   = [2,6] #[2,6] #[4,8]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [5,4]
pol_syst['A0']   = [1,2]
pol_syst['A1']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A3']   = [1,4]
pol_syst['A4']   = [1,6]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [5,4]
pol_syst['A0']   = [1,2]
pol_syst['A1']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A3']   = [1,4]
pol_syst['A4']   = [1,4]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [5,4]
pol_syst['A0']   = [1,2]
pol_syst['A1']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A3']   = [1,2]
pol_syst['A4']   = [1,2]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,4]
pol_syst['A0']   = [1,2]
pol_syst['A1']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A3']   = [1,2]
pol_syst['A4']   = [1,2]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,4]
pol_syst['A0']   = [2,2]
pol_syst['A1']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A3']   = [2,2]
pol_syst['A4']   = [2,2]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [2,2]
pol_syst['A1']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A3']   = [1,2]
pol_syst['A4']   = [2,2]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [2,2]
pol_syst['A1']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A2']   = [2,2]
pol_syst['A3']   = [2,2]
pol_syst['A4']   = [1,2]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [1,2]
pol_syst['A1']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A2']   = [1,2]
pol_syst['A3']   = [1,2]
pol_syst['A4']   = [1,2]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "corr"
pol_syst['corr'] = [3,2]
pol_syst['A0']   = [1,1]
pol_syst['A1']   = [1,1]
pol_syst['A2']   = [1,1]
pol_syst['A2']   = [1,1]
pol_syst['A3']   = [1,1]
pol_syst['A4']   = [1,1]
#pol_systs.append( pol_syst )

pol_syst = copy.deepcopy(pol_default)
pol_syst['run'] = "grid"
pol_syst['corr'] = [8,6]
pol_syst['A0']   = [8,6]
pol_syst['A1']   = [8,6]
pol_syst['A2']   = [8,6]
pol_syst['A2']   = [8,6]
pol_syst['A3']   = [8,6]
pol_syst['A4']   = [8,6]
#pol_systs.append( pol_syst )

pol_default = pol_syst

if args.algo=='jac2':
    command  = './jac2 --nevents='+str(args.nevents) +' --tag='+args.tag+' --run='+pol_default['run']
    for k in pol_default.keys():
        if k=='run': continue
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if args.smear:
        command += ' --smear'
    print(command)
    if not args.dryrun:
        os.system(command)   

elif args.algo=='jac2_systs':
    count = 0
    for syst in pol_systs:
        command  = './jac2 --nevents='+str(args.nevents) +' --tag='+args.tag
        for k in syst.keys():
            if k=='run':
                command += ' --run='+syst[k]
            else:
                command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
        if args.smear:
            command += ' --smear'
        print(command)
        if not args.dryrun:
            os.system(command)
        count += 1
    print('%.0f jobs will be submitted' % count)

elif args.algo=='jac2_systs_vsN':
    count = 0
    for syst in pol_systs:
        for size in sample_sizes.keys():
            command  = './jac2 --nevents='+str(sample_sizes[size]) +' --tag='+args.tag+'_'+size
            for k in syst.keys():
                if k=='run':
                    command += ' --run='+syst[k]
                else:
                    command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
            if args.smear:
                command += ' --smear'
            command += ' --max_y=3.0 --max_x=0.3'
            command += ' --relativistic'
            #command += ' --unweighted'
            command += ' --fullphasespace'
            print(command)
            if not args.dryrun:
                os.system(command)   
            count += 1            
    print('%.0f jobs will be submitted' % count)

elif args.algo=='fit':
    count = 0
    for opt in fit_opts.keys():
        command  = './fit --nevents='+str(args.nevents)+' --tag='+args.tag
        for k in pol_default.keys():
            if k=='run':
                command += ' --run='+pol_default[k]
            else:
                command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
        if (args.rebinX>0 or args.rebinY>0):
            command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
        command += fit_opts[opt]+' --post_tag='+args.post_tag+'_'+opt            
        print(command)
        if not args.dryrun:
            os.system(command)   
        count += 1            
    print('%.0f jobs will be submitted' % count)

elif args.algo=='fit_systs':
    count = 0
    for syst in pol_systs:
        for opt in fit_opts.keys():
            command  = './fit --nevents='+str(args.nevents)+' --tag='+args.tag
            for k in syst.keys():
                if k=='run':
                    command += ' --run='+syst[k]
                else:
                    command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
            if (args.rebinX>0 or args.rebinY>0):
                command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
            command += fit_opts[opt]+' --post_tag='+args.post_tag+'_'+opt            
            print(command)
            if not args.dryrun:
                os.system(command)   
            count += 1            
    print('%.0f jobs will be submitted' % count)

elif args.algo=='fit_systs_vsN':
    count = 0
    for syst in pol_systs:
        for size in sample_sizes.keys():
            for opt in fit_opts.keys():
                command  = './fit --nevents='+str(args.nevents)+' --tag='+args.tag+'_'+size
                for k in syst.keys():
                    if k=='run':
                        command += ' --run='+syst[k]
                    else:
                        command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
                if (args.rebinX>0 or args.rebinY>0):
                    command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
                command += fit_opts[opt]+' --post_tag='+args.post_tag+'_'+opt            
                if 'hMC_from_external' in command and 'jacobians_from_external' not in command:
                    command += ' --scale_MC_uncert='+str(math.sqrt(1.0e+10/sample_sizes[size]))
                if 'add_mass_uncert' in command:
                    command += ' --scale_mass_uncert='+str(math.sqrt(1.0e+10/sample_sizes[size]))
                #command += ' --verbose'
                print(command)
                if not args.dryrun:
                    os.system(command)   
                count += 1            
    print('%.0f jobs will be submitted' % count)

    
elif args.algo=='fit_systs':
    for syst in pol_systs:
        command  = './fit --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure --post_tag='+args.post_tag
        for k in syst.keys():
            command += ' --degs_'+k+'_x='+str(syst[k][0])+' --degs_'+k+'_y='+str(syst[k][1])
        if (args.rebinX>0 or args.rebinY>0):
            command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
        for hel in ['0','1','2','3','4']:
            if ( hasattr(args,"scale"+hel) and getattr(args,"scale"+hel) ):
                command += (' --scale'+hel)
        if args.jacmass>-1:
            command += (' --jacmass='+str(args.jacmass))

        fit_opts = [' --jUL', 
                    ' --jUL --j0', 
                    ' --jUL --j0 --j1', 
                    ' --jUL --j0 --j1 --j2',
                    ' --jUL --j0 --j1 --j2 --j3',
                    ' --jUL --j0 --j1 --j2 --j4',
                    ' --jUL --j0 --j1 --j2 --j3 --j4',
                ]
        for i in fit_opts:
            new_command = command+i
            print(new_command)
            os.system(new_command) 


            
elif args.algo=='fit_fast':
    command  = './fit --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=closure --post_tag='+args.post_tag
    for k in pol_default.keys():
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    fit_opts = ['jUL', 'j0', 'j1', 'j2', 'j3', 'j4']
    #fit_opts = ['jUL', 'j0', 'j1']
    for i in range(len(fit_opts)):
        command += (' --'+fit_opts[i])
    print(command)
    os.system(command) 

elif args.algo=='fit_grid_fast':
    command  = './fit_grid --nevents='+str(args.nevents) +' --prior=0.5 --extraprior=10 --tag='+args.tag+' --run=grid --post_tag='+args.post_tag
    for k in pol_default.keys():
        if 'corr' not in k: continue
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    fit_opts = ['jUL', 'j0', 'j1', 'j2', 'j3', 'j4']
    #fit_opts = ['jUL']
    for i in range(len(fit_opts)):
        command += (' --'+fit_opts[i])
    print(command)
    os.system(command) 

elif args.algo=='fit_grid':
    command  = './fit_grid --nevents='+str(args.nevents) +' --tag='+args.tag+' --run=grid --post_tag='+args.post_tag
    for k in pol_default.keys():
        if 'corr' not in k: continue
        command += ' --degs_'+k+'_x='+str(pol_default[k][0])+' --degs_'+k+'_y='+str(pol_default[k][1])
    if (args.rebinX>0 or args.rebinY>0):
        command += ' --rebinX='+str(args.rebinX)+' --rebinY='+str(args.rebinY)
    fit_opts = [' --jUL', 
                ' --jUL --j0', 
                ' --jUL --j0 --j1', 
                ' --jUL --j0 --j1 --j2',
                ' --jUL --j0 --j1 --j2 --j3',
                ' --jUL --j0 --j1 --j2 --j4',
                ' --jUL --j0 --j1 --j2 --j3 --j4',
            ]
    for i in fit_opts:
        new_command = command+i
        print(new_command)
        os.system(new_command) 
