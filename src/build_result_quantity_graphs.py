from  optparse  import OptionParser
import ROOT as rt
import os, sys, json
import array, math
parser = OptionParser()

parser.add_option("-s", "--siglist", dest="sig_list",
                  help="list of signal result jsons to be analyzed",default="xx4j.list",
                  action="store",type="string")


parser.add_option("-q", "--quantity", dest="quantity",
                  help="quantity",default="nTagHLTSYSRel",
                  action="store",type="string")


parser.add_option("-n", "--tagbin", dest="tagbin",
                  help="which bin in the result json to check",default=2,
                  action="store",type="string")


parser.add_option("--parse", dest="result_list",
                  help="use this option",  default="xx4j_result.list",
                  action="store",type="string")



# contains the systematic value
class result:    
    def __init__(self, x_label, y_label, 
                 val, valUp, valDn):
        # labeling for slices
        self.x_label = x_label
        self.y_label = y_label

        #value varied up and down
        self.val     = val
        self.valUp   = valUp
        self.valDn   = valDn

#container for dictionaries for each limit containing tuples whose
#first index is the dependent variable. Example: a slice of mass
# gives dictionaries index by mass containing tuples of (ctau, limit) 
class slice_container:
    def __init__(self,name):
        self.name    = name  
        self.nPoints = 0
        self.val     = {}
        self.valUp   = {}
        self.valDn   = {}

    def add_result2slice(self, result, fixed_label, float_label):
        # add the list if its not already in the dictionary
        if fixed_label not in self.val.keys():
            # val Up and Dn
            self.val[fixed_label] = []
            self.valUp[fixed_label] = []      
            self.valDn[fixed_label] = []    
        else:
            print "val", result.val, "valUp", result.valUp, "valDn", result.valDn
            self.val[fixed_label].append((float(float_label), result.val))
            self.valUp[fixed_label].append((float(float_label), result.valUp))
            self.valDn[fixed_label].append((float(float_label), result.valDn))        
        
        self.nPoints += 1

    # sort each list by the first tuple element (the floating label)
    def sort(self):
        for key in self.val.keys():
            self.val[key].sort()
            self.valUp[key].sort()
            self.valDn[key].sort()

# container for a full interpretation -
# from an interpretation you can derive limit bands in slices of parameters
class interpretation:    
    def __init__(self):
        self.results = []

    def add_result(self, json_path, x_label, y_label):        

        sig_json_data = open(json_path).read()
        sig_json      = json.loads(sig_json_data)

        val = 0

        if "Rel" not in options.quantity:
            val = sig_json["systematic"][options.quantity][int(options.tagbin)-1]

        valUp = sig_json["systematics"][options.quantity+"Up"][int(options.tagbin)-1]
        valDn = sig_json["systematics"][options.quantity+"Dn"][int(options.tagbin)-1]

        print "\n BUILDING RESULT", json_path, "xlabel", x_label, "ylabel", y_label
        print val, valUp, valDn

        thisResult = result(x_label, y_label, val, valUp, valDn)                                   
        self.results.append(thisResult)

    # fix the x label
    def build_x_slices(self):
        thisSlice = slice_container("x")
        for result in self.results:
            thisSlice.add_result2slice(result, result.x_label, result.y_label)

        return thisSlice

    # fix the y label
    def build_y_slices(self):
        thisSlice = slice_container("y")
        for result in self.results:
            thisSlice.add_result2slice(result, result.y_label, result.x_label)

        return thisSlice
    
    def graphs_from_slice_container(self, sliceCont):
        
        #sort the slices before building graphs in terms of the floating variable
        slice_name = sliceCont.name
        sliceCont.sort()
        
        graph_dict = {}
        for key in sliceCont.val.keys():
            #build empty lists for building tgraph arrays
            float_var  = [] 
            val_list   = [] 
            valUp_list = []
            valDn_list = []  

            #the first entry in the pair is the floating label (mass, ctau0)  
            for pair in sliceCont.val[key]: float_var.append(float(pair[0]))            
            #the second entry in the pair is the quantity were drawing

            for pair in sliceCont.val[key]: 
                val_change = 0 if pair[1] == 0 else math.log10((pair[1]*100))
                val_change = 0 if val_change < 0 else val_change
                val_list.append(val_change)
            for pair in sliceCont.valUp[key]: 
                val_change = 0 if pair[1] == 0 else math.log10((pair[1]*100))
                val_change = 0 if val_change < 0 else val_change
                valUp_list.append(val_change)
            for pair in sliceCont.valDn[key]:
                val_change = 0 if pair[1] == 0 else (-1 * math.log10((pair[1]*100)))
                val_change = 0 if val_change > 0 else val_change
                valDn_list.append(val_change)

            print "valList", val_list
            print "valListUp", valUp_list
            print "valListDn", valDn_list

            # build a graph for each key and add it to the graph dictionary
            val_graph   = rt.TGraph( len(float_var), array.array("d",float_var),  array.array("d", val_list))
            val_graph.SetName("%s_%s_val_graph" % (slice_name, key))
            valUp_graph = rt.TGraph( len(float_var), array.array("d",float_var),  array.array("d", valUp_list))
            valUp_graph.SetName("%s_%s_valUp_graph" % (slice_name, key))
            valDn_graph = rt.TGraph( len(float_var), array.array("d",float_var),  array.array("d", valDn_list))
            valDn_graph.SetName("%s_%s_valDn_graph" % (slice_name, key))

            # print "\n LOADING DICTIONARY", key
            # print "len float var", len(float_var)
            # print "valList", val_list
            # print "valListUp array", array.array("d", valUp_list)
            # print "valListDn array", array.array("d", valDn_list)
            
            # print "float_var", float_var

            graph_dict[key, "val"]   = val_graph
            graph_dict[key, "valUp"] = valUp_graph
            graph_dict[key, "valDn"] = valDn_graph
            
        return graph_dict

if __name__ == '__main__':
    (options, args) = parser.parse_args()
    parser.print_help()


    myInterpretation = interpretation()
    result_list = []
    option_list = options.result_list 

    print "Combine Result Root Files provided...Parsing lines..."
    strip_lines = map( lambda x: x.rstrip("\n") , open(option_list, "r").readlines())

    for line in strip_lines: 
        x_label = line.split("/")[-1].split("_")[2]
        y_label = line.split("/")[-1].split("_")[-1].split("mm")[0]

        print "x_label", x_label, "y_label", y_label
        #x_label = line.split("XX")[1].split("YY")[0]
        #y_label = line.split("YY")[0].split(".")[0]
        result_list.append((line, x_label, y_label))

    for result_tuple in result_list:
        (path, x_label, y_label) = result_tuple
        myInterpretation.add_result(path, x_label, y_label)

    output_file = rt.TFile(options.quantity + "_" + str(options.tagbin) + "tags_graphs.root", "RECREATE")
    xslices     = myInterpretation.build_x_slices()
    xgraph_dict = myInterpretation.graphs_from_slice_container(xslices)
    output_file.cd()

    for key in xgraph_dict.keys(): 
        xgraph_dict[key].Write()
        #xgraph_dict[key].Print()

    print "Taking y slices of limits...."
    #  y slices
    yslices = myInterpretation.build_y_slices()
    ygraph_dict = myInterpretation.graphs_from_slice_container(yslices)
    output_file.cd()
    for key in ygraph_dict.keys(): 
        ygraph_dict[key].Write()
        #ygraph_dict[key].Print()
