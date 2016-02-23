import time, json

folder =  time.strftime("%B_%d").upper()
json_data = open("run_config.json", "r").read()
thisjson  = json.loads(json_data)
path  = "/afs/cern.ch/work/h/hardenbr/2015/DIJET/JETPROB/"

output_folder = thisjson["outputDir"]
split = output_folder.split("/")
datedir = split[-2]
sampledir = split[-1]
print "mkdir " +path + "/" + datedir
print "mkdir " + path + "/" + datedir + "/" + sampledir

