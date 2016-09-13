import sys, os, time
nDiv        = int(sys.argv[1])
run_command = False
only_probs  = False
minIndex    = 0
maxIndex    = 9999

if len(sys.argv) > 3:
    minIndex    = int(sys.argv[2])
    maxIndex    = int(sys.argv[3])

if len(sys.argv) > 4:
    run_command = True

print "usage: execComposite.py [nDivisions] [beginning probability index] [end probability index]"

for ii in range(nDiv):
    if ii < minIndex or ii > maxIndex: continue
    for jj in range(nDiv):
        #if only_probs and ((jj != ii+1) or (ii==nDiv-1)): continue

        if ii == jj:
            runString = "./jetTagAnalyzer %i %i %i > /dev/null &" % (nDiv, ii ,jj)
            print runString

            if run_command: 
                print "RUNNING IT.....", runString                
                os.system(runString)
                print "SLEEPING 45 Seconds....."
                time.sleep(20)
            
