import os
import stat

def main():

    debug = False

    #File name: 
    #Start and end of run from RCT [STF - 120000:ETF + 120000] 2 minutes buffer, from run_IR_Bfield.txt files

    #General path to add flexibility to the code + period name
    #period = "LHC23_pass4_skimmed_QC1" #pp skimmed QC data of 2023 pass 4
    #period = "LHC23_PbPb_pass3_I-A11" #Pb-Pb dataset - one of the two used for the analyses of Nazar
    #period = "LHC23_PbPb_pass3_fullTPC" #Pb-Pb dataset - other used for the analyses of Nazar
    #period = "LHC22o_pass7_minBias"
    #period = "LHC22_pass7_skimmed"
    #period = "LHC23_pass4_skimmed"
    #period = "LHC23_PbPb_pass4"
    #period = "LHC24_pass1_skimmed"
    #period = "LHC25ae_pass2"
    period = "LHC24_ppref_pass1"

    globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/" #Where txt file with star/end of run is stored

    ccdbPath = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/" #Where CCDB objects are stored

    periodRuns = []
    startTimes = []
    endTimes = []

    #Open list of runs
    with open(globalPath+"run_IR_Bfield.txt") as config:
        runs = config.readlines()
        #print(run)
        for run in runs:
            asList = run.split("\t")
            if int(asList[0]) == 1:
                periodRun = int(asList[1])
                time1 = int(asList[4])
                time2 = int(asList[5])
                periodRuns.append(periodRun)
                startTimes.append(time1)
                endTimes.append(time2)
                
                """
                if period == "LHC22_pass7_skimmed":
                    startTimes.append(time1)
                    endTimes.append(time2)

                else:
                    startTimes.append(time1 - 120000)
                    endTimes.append(time2 + 120000)
                """
    if debug:
        print(periodRuns)
        print(startTimes)
        print(endTimes)

    print("\n\n",endTimes[-1],"\n\n")

    #Create upload.sh script to upload the data to ccdb
    if os.path.exists(ccdbPath+"upload.sh"):
        os.remove(ccdbPath+"upload.sh")
        print("Removing old file")
    
    else:
        print("The file does not exist!")

    outFile = open(ccdbPath+"upload.sh", "a") 

    outFile.write("#!/bin/bash\n")
    outFile.write("\n")
    outFile.write("ccdbhost=http://alice-ccdb.cern.ch\n")
    outFile.write("\n")

    #Open folder where all root files are saved and get them into an array
    files = os.listdir(ccdbPath)
    print(files)

    for file in files:
    # Check if item is a file, not a directory
        if not os.path.isdir(os.path.join(ccdbPath, file)) and not (os.path.join(ccdbPath, file)).endswith('.sh'):
            if debug:
                print(file)
            
            file = file.removesuffix('.root')
            file = file.split("_")
            if debug:
                print(file)
                print(file[1],file[2])
            
            #first and lust run of the merging
            run1 = int(file[1])
            run2 = int(file[2]) 
            
            #placeholder variables
            s = 0
            e = 0

            #Get start time of first run and end time of last run in the merge
            for r1, t1 in zip(periodRuns,startTimes):
                if r1 == run1:
                    s = t1
                    break
            
            for r2, t2 in zip(periodRuns, endTimes):
                if r2 == run2:
                    e = t2
                    break
            
            #print("For runs",file[1],file[2]," the start time is",s," and the time is",e)
            outFile.write("o2-ccdb-upload --host \"$ccdbhost\" -p MID/Calib/ChamberEfficiency -f o2-mid-ChEffCounter_" + str(run1) + "_" + str(run2) + 
                        ".root -k ccdb-object --starttimestamp " + str(s-120000) + " --endtimestamp " + str(e+120000) + " -m \"runNumber=" + str(run1) + "-" + str(run2) +
                        ";JIRA=O2-6423;adjustableEOV=true;Created=" + str(endTimes[-1] + 1000) + "\"\n")

    #Flush content to file and close it
    outFile.flush()
    outFile.close()

    #Make the file executable
    #os.chmod(outFile, os.stat(outFile).st_mode | stat.S_IEXEC)

if __name__ == "__main__":
    main()