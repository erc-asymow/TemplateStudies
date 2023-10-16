import os

do_fits = False
do_syst = True

phase_spaces = ["--xf_max=0.4 --yf_max=3.5",
                #"--xf_max=0.3 --yf_max=3.0"
                ]

systs = ["scet",
         "scale",
         "pdf",
         "altpdf"
         ]

procs = [#"UL",
         "A0",
         "A1",
         "A2",
         #"A3",
         "A4"
         ]


if __name__ == '__main__':

    if do_fits:
        for ips in phase_spaces:
            for isys in systs:
                for iproc in procs:
                    if isys=="scet" and iproc!="UL":
                        continue
                    posttag = ' --posttag=add'
                    cmd_jac  = 'python run_getparam.py --algo=run --mt --jac '+ips+(' --do'+iproc if iproc!="UL" else '' )+posttag #+' --dryrun'
                    cmd_fit  = 'python run_getparam.py --algo=run --mt --fit '+ips+(' --do'+iproc if iproc!="UL" else '' )+' --systname='+isys+posttag
                    cmd_plot = 'python run_getparam.py --algo=plot_fitres '+ips+(' --do'+iproc if iproc!="UL" else '' )+' --systname='+isys+' --posttag=add'+' --batch'+posttag
                    print(cmd_jac)
                    print(cmd_fit)
                    print(cmd_plot)
                    os.system(cmd_jac)
                    os.system(cmd_fit)
                    os.system(cmd_plot)
    elif do_syst:
        #posttag = ' --posttag=V0'
        posttag = ' --posttag=add'
        for ips in phase_spaces:
            cms_syst = 'python run_getparam.py --algo=run --mt --syst '+ips+' '+posttag #+' --dryrun'
            print(cms_syst)
            os.system(cms_syst)
