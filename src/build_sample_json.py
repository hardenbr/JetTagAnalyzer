import json, sys

if len(sys.argv) != 2:
    print "build_sample_json.py [list_of_signal files]"
    

file = open(sys.argv[1], "r")
lines = map(lambda x: x.rstrip("\n"), file.readlines()) 

for line in lines:
    #    print line
    label = line.split("_")[0]
    #signal
    x_label = line.split("M-")[1].split("_")[0]
    y_label = line.split("CTau-")[1].split("_")[0].split("mm")[0]
    # qcd
    # x_label = line.split("Pt_")[1].split("to")[0]
    # y_label = line.split("Pt_")[1].split("to")[1].split("_")[0]
    # x_label = line.split("HT")[1].split("to")[0]
    # y_label = line.split("HT")[1].split("to")[1].split("_")[0]

    # neutrino
    #x_label = int(line.split("_M")[1].split("_")[0])
    #y_label = int(line.split("tau")[1].split("um")[0])

    json_val = """{  
          \"runSample\": true,
          \"label\": \"%s_%s_%smm\",
          \"stack\": \"%s_%s_%s\",
          \"path\":  \"%s\",
          \"isMC\":  true,
          \"isSig\": true,
          \"xsec\": 1,
          \"x_limit_label\": %s,
          \"y_limit_label\": %s
       },""" % (label, x_label, y_label, label, x_label, y_label, line, x_label, y_label)
       
    print json_val

       
       

