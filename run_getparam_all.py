import os

do_fits = False
do_syst = True

phase_spaces = [#"--xf_max=0.4 --yf_max=3.5",
                "--xf_max=0.3 --yf_max=3.0"
                ]

systs = [
    "scet",
    "scale",
    "pdf",
    "altpdf",
    "scale_31point"
]

procs = [#"UL",
         #"A0",
         #"A1",
         #"A2",
         "A3",
         #"A4"
         ]


if __name__ == '__main__':

    if do_fits:
        for ips in phase_spaces:
            for iproc in procs:
                do_str = (' --do'+iproc if iproc!="UL" else '' )
                posttag = ' --posttag=add '
                #posttag = ''
                cmd_jac  = 'python run_getparam.py --algo=run --mt --jac ' + posttag + ips + do_str #+' --dryrun'
                print(cmd_jac)
                os.system(cmd_jac)
                for isys in systs:                
                    if isys=="scet" and iproc!="UL":
                        continue
                    if isys=="scale_31point" and iproc=="UL":
                        continue
                    cmd_fit  = 'python run_getparam.py --algo=run --mt --fit ' + posttag + ips + do_str + ' --systname='+isys
                    cmd_plot = 'python run_getparam.py --algo=plot_fitres '    + posttag + ips + do_str + ' --systname='+isys + ' --batch'
                    print(cmd_fit)
                    print(cmd_plot)
                    os.system(cmd_fit)
                    os.system(cmd_plot)
    elif do_syst:
        posttag = ' --posttag=V10 '
        for ips in phase_spaces:
            cms_syst = 'python run_getparam.py --algo=run --mt --syst ' + ips + ' ' + posttag #+' --dryrun'
            print(cms_syst)
            os.system(cms_syst)

'''
Log-book:
V0  = start
add = change jac to additive for A3
V1  = increase UL:syst_deg_y 4 -> 6
V2  = decrease all. Smaller ranges
V3  = clip and revert to larger degrees
V4  = switched to V3/ samples. Clip. Decrease all but A4 and UL
V5  = switched to V3/ samples. Clip. Larger degrees
V6  = switched to V3/ samples. Clip. Larger degrees
V7  = switched to V4/ samples (CT18Z). Clip. Minimal degrees: (5,4), (1,2), ..., (2,2), (4,2)
V8  = switched to V4/ samples (CT18Z). Clip. Larger degrees for UL, same for others: (8,6), (1,2), ..., (2,2), (4,2)
V9  = switched to V4/ samples (CT18Z). Fix to A4: (4,2) -> (2,4)
'''
