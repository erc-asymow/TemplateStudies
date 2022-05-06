# TemplateStudies
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.02/x86_64-centos7-gcc48-opt/bin/thisroot.sh
g++ -O2 main.cpp -o main.exe `root-config --cflags --glibs`
./main.exe