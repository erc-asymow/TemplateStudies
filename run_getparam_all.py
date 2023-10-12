import os

phase_spaces = ["--xf_max=0.4 --yf_max=3.5",
                #"--xf_max=0.3 --yf_max=3.0"
                ]

systs = [#"scet",
         #"scale",
         #"pdf",
         "altpdf"
         ]

procs = ["UL",
         #"A0",
         #"A1",
         #"A2",
         #"A3",
         #"A4"
         ]


if __name__ == '__main__':
    for ips in phase_spaces:
        for isys in systs:
            for iproc in procs:
                if isys=="scet" and iproc!="UL":
                    continue
                cmd_fit  = 'python run_getparam.py --algo=run --mt --fit '+ips+(' --do'+iproc if iproc!="UL" else '' )+' --systname='+isys
                cmd_plot = 'python run_getparam.py --algo=plot_fitres '+ips+(' --do'+iproc if iproc!="UL" else '' )+' --systname='+isys+' --batch'
                print(cmd_fit)
                print(cmd_plot)
                os.system(cmd_fit)
                #os.system(cmd_plot)
