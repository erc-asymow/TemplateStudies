import os

lumi = 1000000

nevents = [
    100000,
    500000,
    1000000,
    5000000,
    10000000,
    50000000,
    100000000,
    500000000,
]

degs = [
    2,4,6,8
]

if __name__ == '__main__':    
    for deg in degs:
        count = 0
        for nevent in nevents: 
            tag = 'mctstat'+str(count)
            cmd  = './toy --nevents='+str(nevent)+'  --lumi='+str(lumi)+' --degs_x='+str(deg)+' --do_fit --scalejac=1 --tag='+tag
            print(cmd)
            os.system(cmd)
            count += 1

