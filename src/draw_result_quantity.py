from  optparse  import OptionParser
import ROOT as rt
import rootlogon, json, CMS_lumi
parser = OptionParser()

parser.add_option("-s", "--siglist", dest="sig_list",
                  help="list of signal result jsons to be analyzed",default="xx4j_result.list",
                  action="store",type="string")

parser.add_option( "--ctaulabel", dest="ctau_label",
                  help="label for the ctau scan",default="mm",
                  action="store",type="string")


parser.add_option( "--draw_options", dest="draw_options",
                  help="graph draw options",default="LP",
                  action="store",type="string")


parser.add_option( "--ylabel", dest="ylabel",
                  help="label for the y axis of the scan",default="log_{10}(HLT % Error)",
                  action="store",type="string")

parser.add_option( "-f", dest="path",
                  help="path to the combine root file",default="nTagHLTSYSRel_graphs.root",
                  action="store",type="string")

parser.add_option( "--ymin", dest="ymin",
                  help="y axis min ",default=-100,
                  action="store",type="float")

parser.add_option( "--ymax", dest="ymax",
                  help="y axis max ",default=100,
                  action="store",type="float")


parser.add_option( "--dsusy", dest="dsusy",
                  help="y axis max ",default=False,
                  action="store_true",)


rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)
(options, args) = parser.parse_args()
color_scheme = {}
if not options.dsusy:
    # XX4J COLOR SCHEME
    # mass 
    color_scheme["x_50"]   = rt.kBlack
    color_scheme["x_100"]  = rt.kBlue+2
    color_scheme["x_200"]  = rt.kBlue
    color_scheme["x_300"]  = rt.kBlue-3
    color_scheme["x_500"]  = rt.kBlue-5
    color_scheme["x_700"]  = rt.kBlue-9
    color_scheme["x_1000"] = rt.kViolet-9
    color_scheme["x_1500"] = rt.kViolet-5
    color_scheme["x_2000"] = rt.kViolet
    color_scheme["x_3000"] = rt.kViolet+2
# ctau
    color_scheme["y_1"]    = rt.kBlue+2
    color_scheme["y_3"]    = rt.kBlue
    color_scheme["y_10"]   = rt.kBlue-5
    color_scheme["y_30"]   = rt.kBlue-9
    color_scheme["y_100"]  = rt.kRed-9
    color_scheme["y_300"]  = rt.kRed-5
    color_scheme["y_1000"] = rt.kRed
    color_scheme["y_2000"] = rt.kRed+2
    color_scheme["label"] = "X^{0} Mass [GeV]" 
    color_scheme["mass_label"]  = "m_{X^{0}}" 
else: 
    # DSUSY COLOR SCHEME
    color_scheme["x_200"]   = rt.kBlack
    color_scheme["x_300"]  = rt.kBlue+2
    color_scheme["x_400"]  = rt.kBlue
    color_scheme["x_500"]  = rt.kBlue-3
    color_scheme["x_600"]  = rt.kBlue-5
    color_scheme["x_700"]  = rt.kBlue-9
    color_scheme["x_800"] = rt.kViolet-9
    color_scheme["x_900"] = rt.kViolet-5
    color_scheme["x_1000"] = rt.kViolet
    color_scheme["x_1200"] = rt.kViolet+2
    color_scheme["label"]  = "#tilde{t} Mass [GeV]" 
    color_scheme["mass_label"]  = "m_{#tilde{t}}" 

# #NEUTRALINO COLOR SCHEME
# color_scheme["x_400"]  = rt.kBlue
# color_scheme["x_800"]  = rt.kBlue-9
# color_scheme["x_1200"]  = rt.kViolet-9
# color_scheme["x_1600"]  = rt.kViolet-3
# # ctau
# color_scheme["y_100"]    = rt.kBlue
# color_scheme["y_300"]   = rt.kBlue-9
# color_scheme["y_1000"]  = rt.kRed-9
# color_scheme["y_10000"] = rt.kRed
# color_scheme["label"]  = "#tilde{#chi}^{0} Mass [GeV]" 
# color_scheme["mass_label"]  = "m_{#tilde{#chi}^{0}}" 

def get_graph_with_full_range(graphs_to_draw):
    nMax = 0
    keyMax = None

    for key in graphs_to_draw.keys():
        graph = graphs_to_draw[key]
        
        if graph.GetN() > nMax:
            graphMax = graph
            nMax     = graph.GetN()
            keyMax   = key

    return keyMax
    

if __name__ == '__main__':
    rootlogon.style()
    rt.gStyle.SetPadRightMargin(0.2)

    (options, args) = parser.parse_args()

    input_file = rt.TFile(options.path,"READ")

    x_canvas  = rt.TCanvas("x_slice_canvas", "x", 1024,768)
    y_canvas  = rt.TCanvas("y_slice_canvas", "y", 1024, 768)

    signal_json_list_file = open(options.sig_list,"r")
    sig_json_lines        = map(lambda x: x.rstrip('\n') , signal_json_list_file.readlines())
    sig_jsons             = []

    x_graph_labels              = [] 
    y_graph_labels              = [] 
    
    # read in each json
    for line in sig_json_lines:
        print "checking json", line
        sig_json_data = open(line).read()
        sig_json      = json.loads(sig_json_data)

        #from the json determine the label for the limit root result
        x_label   = sig_json["x_limit_label"]
        y_label   = sig_json["y_limit_label"]

        print "x label", x_label, "y label", y_label

        #build the graph lists 
        if "x_%s" % x_label not in x_graph_labels:
            x_graph_labels.append("x_%s" % x_label)

        if "y_%s" % y_label not in y_graph_labels:
            y_graph_labels.append("y_%s" % y_label)

    x_graphs_to_draw = {}
    y_graphs_to_draw = {}
    x_graphs_to_drawUp = {}
    y_graphs_to_drawUp = {}
    x_graphs_to_drawDn = {}
    y_graphs_to_drawDn = {}

    #build the graphs
    print "looping over graph labels"
    for graph_label in x_graph_labels+y_graph_labels:
        type_labels = ["val","valUp","valDn"]
        print "looping over value types"
        for type_label in type_labels:
            graph_name = "%s_%s_graph" % (graph_label, type_label)
            print "searching for graph name", graph_name

            graph = input_file.Get(graph_name)
            if graph == None:
                print "GRAPH NOT FOUND....continuing"
                continue 
            graph.Print()
            #set the color based on the scheme 
            if graph_label in color_scheme.keys():
                graph.SetLineColor(color_scheme[graph_label])
                graph.SetMarkerColor(color_scheme[graph_label])
            else:
                print "graph does not have color sceheme defined: ", graph_label
                
            graph.SetMarkerStyle(25)
            graph.SetMarkerSize(1.1)
            graph.SetLineWidth(3)

            # put the graphs into the appropriate dictionary for draphing
            print "placing graph label in dictionary", graph_label, type_label 
            if graph_label in x_graph_labels:
                if type_label == "valUp": 
                    x_graphs_to_drawUp[(graph_label, type_label)] = graph
                elif type_label == "valDn":
                    x_graphs_to_drawDn[(graph_label, type_label)] = graph
                else:
                    x_graphs_to_draw[(graph_label, type_label)] = graph
            else: #otherwise its a y slice
                if type_label == "valUp": 
                    y_graphs_to_drawUp[(graph_label, type_label)] = graph
                elif type_label == "valDn":
                    y_graphs_to_drawDn[(graph_label, type_label)] = graph
                else:
                    y_graphs_to_draw[(graph_label, type_label)] = graph

    print "keys for y graph val", y_graphs_to_draw.keys()
    print "keys for y graph valUp", y_graphs_to_drawUp.keys()
    print "keys for y graph valDn", y_graphs_to_drawDn.keys()

    print "keys for x graph val", x_graphs_to_draw.keys()
    print "keys for x graph valUp", x_graphs_to_drawUp.keys()
    print "keys for x graph valDn", x_graphs_to_drawDn.keys()

    ##sort the labels
    sort_lambda = lambda x: int(x.split("_")[1])
    x_graph_labels = sorted(x_graph_labels, key=sort_lambda)
    y_graph_labels = sorted(y_graph_labels, key=sort_lambda)
    print x_graph_labels
    print y_graph_labels
    
    #draw the x slices
    x_canvas.cd()
    first_key    = get_graph_with_full_range(x_graphs_to_draw)
    draw_first_x = None

    if first_key[1] == "val":
        draw_first_x = x_graphs_to_draw[first_key]
    if first_key[1] == "valUp":
        draw_first_x = x_graphs_to_drawUp[first_key]
    if first_key[1] == "valDn":
        draw_first_x = x_graphs_to_drawDn[first_key]

    # draw_firstDn = x_graphs_to_drawDn[first_key]
    draw_first_x.SetMinimum(-50)
    draw_first_x.SetMaximum(50)
    draw_first_x.GetXaxis().SetRangeUser(1,2000)
    draw_first_x.GetYaxis().SetRangeUser(-50,50)


    #draw_first_x.Draw("A"+options.draw_options)
    draw_first_x.Draw("AL")                      
    #draw_first_x.Draw("A")

    # draw_first_x.SetMinimum(-2)
    # draw_first_x.SetMaximum(2)
    # draw_first_x.GetXaxis().SetRangeUser(1,2000)
    # draw_first_x.GetYaxis().SetRangeUser(-2,2)


    draw_first_x.GetXaxis().SetTitle("c#tau_{0} [%s]" % options.ctau_label)
    draw_first_x.GetYaxis().SetTitle(options.ylabel)
    #for key in x_graphs_to_draw.keys(): x_graphs_to_draw[key].Draw("same"+options.draw_options)        
    for key in x_graphs_to_drawUp.keys(): x_graphs_to_drawUp[key].Draw("same"+options.draw_options)        
    for key in x_graphs_to_drawDn.keys(): x_graphs_to_drawDn[key].Draw("same"+options.draw_options)        

    #x_canvas.SetLogy(1)
    x_canvas.SetLogx(1)
    x_canvas.SetGridy(1)
    x_canvas.SetGridx(1)
    x_canvas.Modified()    
    #deal with the legend
    x_leg = rt.TLegend(.81, .45, .99, .91)    
    #x_leg_th = rt.TLegend(.53, .41, .72, .63)

    print "x graph labels", x_graph_labels
    for label in x_graph_labels:
        if (label,"val") not in x_graphs_to_draw.keys(): continue
        graph = x_graphs_to_draw[(label,"val")]
        x_leg.AddEntry(graph, "%s=%s" % (color_scheme["mass_label"], label.split("_")[1]), "pl")

    x_leg.SetFillColor(0)
    x_leg.SetLineColor(0)
    x_leg.Draw('same')
    CMS_lumi.CMS_lumi(x_canvas, 4, 0)


    # draw_first_x.SetMinimum(-2)
    # draw_first_x.SetMaximum(2)
    # draw_first_x.GetXaxis().SetRangeUser(2,2E3+5)
    # draw_first_x.GetYaxis().SetRangeUser(-2,2)
    # x_canvas.Modified()
    
    #draw the y slices
    y_canvas.cd()

    first_key = get_graph_with_full_range(y_graphs_to_draw)
    draw_first_y = None 

    if first_key[1] == "val":
        draw_first_y = y_graphs_to_draw[first_key]
    if first_key[1] == "valUp":
        draw_first_y = y_graphs_to_drawUp[first_key]
    if first_key[1] == "valDn":
        draw_first_y = y_graphs_to_drawDn[first_key]


    draw_first_options = "AP" if first_key != "val"else "A"
    draw_first_y.Draw(draw_first_options+options.draw_options)

    draw_first_y.SetMinimum(-2)
    draw_first_y.SetMaximum(2)
    draw_first_y.GetXaxis().SetRangeUser(1,1E6)
    draw_first_y.GetYaxis().SetRangeUser(-50,50)

    draw_first_y.GetXaxis().SetTitle(color_scheme["label"])
    draw_first_y.GetYaxis().SetTitle(options.ylabel)

    for key in y_graphs_to_draw.keys(): y_graphs_to_draw[key].Draw("same"+options.draw_options)
    for key in y_graphs_to_drawUp.keys(): y_graphs_to_drawUp[key].Draw("same"+options.draw_options)
    for key in y_graphs_to_drawDn.keys(): 
        y_graphs_to_drawDn[key].Print()
        y_graphs_to_drawDn[key].Draw("same" +options.draw_options)

    y_canvas.SetLogx(1)
    y_canvas.SetGridy(1)
    y_canvas.SetGridx(1)

    print "y graph labels", y_graph_labels
    #deal with the legend
    y_leg = rt.TLegend(.81, .45, .99, .91)    
    #y_leg_th = rt.TLegend(.53, .41, .72, .63)
    included = []
    for label in y_graph_labels:
        if (label,"val") not in y_graphs_to_draw.keys(): continue
        graph = y_graphs_to_draw[(label,"val")]
        y_leg.AddEntry(graph, "c#tau=%s %s" % (label.split("_")[1],options.ctau_label), "pl")

    y_leg.SetFillColor(0)
    y_leg.SetLineColor(0)
    y_leg.Draw('same')    

    # zero_line_y = rt.TLine(draw_first_y.GetXaxis().GetXmin(), 0, draw_first_y.GetXaxis().GetXmax(), 0)
    # zero_line_y.SetLineWidth(3)
    # #zero_line.SetLineColor(rt.Black)
    # zero_line_y.Draw("same")

    CMS_lumi.CMS_lumi(y_canvas, 4, 0)

    x_canvas.cd()
    draw_first_x.SetMinimum(-50)
    draw_first_x.SetMaximum(50)
    draw_first_x.GetXaxis().SetRangeUser(1,2000)
    draw_first_x.GetYaxis().SetRangeUser(options.ymin, options.ymax)
    x_canvas.Update()
    x_canvas.Modified()

    zero_line_x = rt.TLine(1, 0, 2000, 0)
    zero_line_x.SetLineWidth(7)
    zero_line_x.SetLineColor(rt.kBlack)
    zero_line_x.Draw("same")


    y_canvas.cd()
    draw_first_y.SetMinimum(-50)
    draw_first_y.SetMaximum(50)
    draw_first_y.GetXaxis().SetRangeUser(50,3000)
#    draw_first_y.GetYaxis().SetRangeUser(options.ymin, options.ymax)
#    draw_first_y.GetYaxis().SetRangeUser(-2.5, 2.5)
    draw_first_y.GetYaxis().SetRangeUser(-50, 50)

    zero_line_y = rt.TLine(50, 0, 3000, 0)
    zero_line_y.SetLineWidth(6)
    #zero_line.SetLineColor(rt.kBlack)
    zero_line_y.Draw("same")


    x_canvas.Update()
    x_canvas.Modified()

    y_canvas.Update()
    y_canvas.Modified()

    raw_input("RAW INPUT")

    sig_name = options.sig_list.split(".list")[0] + "_" + options.path.split(".root")[0]
    quantity_type = options.path.split(".root")[0]
    # x_canvas.SaveAs("%s_%sxslice.pdf" % (sig_name, quantity_type))
    # y_canvas.SaveAs("%s_%s_yslice.pdf" % (sig_name, quantity_type))
