import requests
import os

def main():

    #File name: 
    #Start and end of run from RCT [STF - 120000:ETF + 120000] 2 minutes buffer, from "run_dates.txt"

    #General path to add flexibility to the code + period name
    #period = "LHC23_pass4_skimmed_QC1" #pp skimmed QC data of 2023 pass 4
    #period = "LHC23_PbPb_pass3_I-A11" #Pb-Pb dataset - one of the two used for the analyses of Nazar
    #period = "LHC23_PbPb_pass3_fullTPC" #Pb-Pb dataset - other used for the analyses of Nazar
    #period = "LHC22o_pass7_minBias"
    #period = "LHC22_pass7_skimmed"
    #period = "LHC23_pass4_skimmed"
    #period = "LHC23_PbPb_pass4"
    period = "LHC24_pass1_skimmed"
    #period = "LHC25ae_pass2"

    globalPath = "/media/luca/Extreme SSD/MIDefficieincy/"+period+"/" #Where txt file with star/end of run is stored

    ccdbPath = "/home/luca/cernbox/assegnoTorino/MIDefficiency/AO2D/"+period+"/ccdb/" #Where CCDB objects are stored

    periodRuns = []
    startTimes = []

    with open(globalPath+"run_IR_Bfield.txt") as config:
        runs = config.readlines()
        #print(run)
        for run in runs:
            asList = run.split("\t")
            if int(asList[0]) == 1:
                periodRun = int(asList[1])
                time = int(asList[4])
                periodRuns.append(periodRun)
                if period == "LHC22_pass7_skimmed":
                    startTimes.append(time)
                else:
                    startTimes.append(time - 120000)

    print(periodRuns)
    print(startTimes)

    #with open(ccdbPath+"modifyMetaData.sh","w")

    if os.path.exists(ccdbPath+"modifyMetaData.sh"):
        os.remove(ccdbPath+"modifyMetaData.sh")
        print("Removing old file")
    else:
        print("The file does not exist!")

    outFile = open(ccdbPath+"modifyMetaData.sh", "a") 

    outFile.write("#!/bin/bash\n")
    outFile.write("\n")
    outFile.write("ccdbhost=http://alice-ccdb.cern.ch\n")
    outFile.write("\n")

    for startrun,periodRun in zip(startTimes,periodRuns):
        url = "http://alice-ccdb.cern.ch/MID/Calib/ChamberEfficiency/" + str(startrun)
        print(url, periodRun)
        response = requests.get(url, allow_redirects=False)  # don't follow 303 redirect

        etag = response.headers.get("ETag")  # directly extract header
        etag = etag.strip('"')
        print("ETag:", etag)

        outFile.write("curl -i -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem -H -L -X PUT \'http://alice-ccdb.cern.ch/MID/Calib/ChamberEfficiency/" + str(startrun) + "/" + etag + "?Creation=1640991600000\'\n")
        if period == "LHC24_pass1_skimmed":
            outFile.write("curl -i -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem -H -L -X PUT \'http://alice-ccdb.cern.ch/MID/Calib/ChamberEfficiency/" + str(startrun) + "/" + etag + "?JIRA=O2-5780\'\n")

        print("curl -i -k --cert ~/.globus/usercert.pem --key ~/.globus/userkey.pem -H -L -X PUT \'http://alice-ccdb.cern.ch/MID/Calib/ChamberEfficiency/" + str(startrun) + "/" + etag + "?Creation=1640991600000\'")    

    outFile.flush()
    outFile.close()

if __name__ == "__main__":
    main()