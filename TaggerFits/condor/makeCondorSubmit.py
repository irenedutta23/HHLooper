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

NSections = 10
readFiles = ""

pwd = os.getcwd()

outDir=pwd +"/"+ str(date.today())
print(outDir)
if not os.path.exists(outDir):
    mkdirComm="mkdir " +outDir
    print(mkdirComm)
    os.system(mkdirComm)
    

for data in datasets:
    jobidx = 0
    inputfname = data+".txt"
    with open(inputfname) as inputfile:
        readFiles = inputfile.readlines()
        print "len(readFiles)", len(readFiles)
    NSections = 1


    NFilesTotal = len(readFiles)
    TotalFiles = NFilesTotal

    print "Dataset ",  data, " NFilesTotal ", NFilesTotal
    NFilesDone  = 0

    while( NFilesDone < NFilesTotal ) :
        thisList = readFiles[NFilesDone : NFilesDone+NSections]
        print "NFilesDone ", NFilesDone, "len(thisList)", len(thisList)

        ##you may have to give full path i.e. CurrentDIR/condor_submit/runlist_...
        inputRunListName = pwd+"/condor_submit/"+data+"_"+str(jobidx)+".txt"
        inputRunList = open(inputRunListName, "w")
        for line in thisList:
            inputRunList.write(line)

        condorSubmit = "condor_submit/submitCondor_"+data+"_"+str(jobidx)
        jobName      = str(date.today())+"_"+data+"_job"+str(jobidx)
        outHistFile = data+"_job"+str(jobidx)+".root"
        if("Data" in data):
            isData       =  "1"
        else:
            isData       =  "0"

        if("2018" in data):
            year = "2018"

        shutil.copyfile("proto_condor_submit",condorSubmit)
        for line in fileinput.FileInput(condorSubmit, inplace=1):
            line=line.replace("JOBNAME", jobName)
            line=line.replace("FILELIST",inputRunListName)
            line=line.replace("ROOTOUT",outHistFile)
            line=line.replace("DATANAME",data)
            line=line.replace("ISDATA",isData)
            line=line.replace("OUTDIR",outDir)
            line=line.replace("YEAR",year)
            print line.rstrip()
                        
        submitCommand = "condor_submit "+condorSubmit
        print submitCommand
        os.system(submitCommand)     
        jobidx = jobidx+1
        NFilesDone = NFilesDone + len(thisList)

    print "Final NFilesDone ", NFilesDone
