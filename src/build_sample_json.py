import json, sys

if len(sys.argv) != 2:
    print "build_sample_json.py [list_of_signal files]"
    

file = open(sys.argv[1], "r")
lines = map(lambda x: x.rstrip("\n"), file.readlines()) 

real_samps = {}

ctau0_val_num = [1, 2, 3, 5, 7, 10,
                 20, 30, 50, 70, 100, 
                 200, 300, 400, 500, 700, 1000,
                 2000]
ctau0_vals = []

for ii in ctau0_val_num: ctau0_vals.append(str(ii))


json_vals = {}
#all samples that need no re-weighting
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

    #mh125
    #x_label = line.split("MS")[1].split("_")[0]
    #y_label = line.split("ctauS")[1].split("_")[0]

    # neutrino
    #x_label = int(line.split("_M")[1].split("_")[0])
    #y_label = int(line.split("tau")[1].split("um")[0])

#    print "REAL SAMPS", (x_label,y_label)
    real_samps[(x_label,y_label)] = line

    json_val = """{  
          \"runSample\": true,
          \"label\": \"%s_%s_%smm\",
          \"stack\": \"%s_%s_%s\",
          \"path\":  \"%s\",
          \"isMC\":  true,
          \"isSig\": true,
          \"xsec\": 1,
          \"orignal_lifetime\": %s,
          \"x_limit_label\": %s,
          \"y_limit_label\": %s
       },""" % (label, x_label, y_label, label, x_label, y_label, line, y_label, x_label, y_label)
       
    json_vals[(x_label,y_label)] = json_val

def find_best_reweight(mx, ctau0):

    best_y = -999
    minDeltaY = 9999

    # the highest ctau0 value 
    max_y = -999
    for pair in real_samps.keys():
        if pair[0] == mx: 
            max_y = max(float(pair[1]), float(max_y))
    

    for pair in real_samps.keys():
        x,y = pair
        
        # they must have the same mass value
        if x != mx: continue
        
        # find the closest point larger than the value
        deltaY = float(y) - float(ctau0)
        if (deltaY  > 0 and deltaY < minDeltaY) and (float(ctau0) < max_y):
            best_y    = y
            minDeltaY = deltaY

        if (float(ctau0) > max_y) and abs(deltaY) < minDeltaY:
            best_y    = y
            minDeltaY = abs(deltaY)
            
            
    
    return best_y

for pair in json_vals.keys():
    x,y = pair 

    for ctau0 in ctau0_vals:
        if (x, str(ctau0)) not in json_vals.keys():
            x_label = x
            y_label = ctau0
            
            original_lifetime = find_best_reweight(x_label, y_label)
            line = real_samps[(x_label, original_lifetime)] 

            json_val = """{  
          \"runSample\": true,
          \"label\": \"%s_%s_%smm\",
          \"stack\": \"%s_%s_%s\",
          \"path\":  \"%s\",
          \"isMC\":  true,
          \"isSig\": true,
          \"xsec\": 1,
          \"original_lifetime\": %s,
          \"x_limit_label\": %s,
          \"y_limit_label\": %s
},""" % (label, x_label, y_label, label, x_label, y_label, line, original_lifetime, x_label, y_label)
        
            json_vals[(x_label,ctau0)] = json_val

sorted_keys = sorted(json_vals)

print "{\"samples\": ["
for key in sorted_keys:
    print json_vals[key]
print "]}"            
