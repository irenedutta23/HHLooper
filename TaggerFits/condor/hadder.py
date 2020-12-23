import os,sys,re,fileinput,string,shutil
from datetime import date
##             Dataset        Name   
datasets = ["Data_2018",
            "ttbar_2018",
            "top_2018",
            "QCD_2018",
            "wqq_2018",
            "zqq_2018"
]

pwd = os.getcwd()

outDir=pwd +"/"+sys.argv[1]
print(outDir)
os.chdir(outDir)
for data in datasets:
    haddComm = "hadd -k -f " + data+".root " + data+"_job*.root "
    os.system(haddComm)
