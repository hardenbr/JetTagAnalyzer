import json, sys, math

x_unit = "m_{X}~[GeV]"
y_unit = "c\\tau_0 [mm]"

json_names = sys.argv[1:]

class column:
    def __init__(self, json):
        self.label                = json["label"]
        self.isData               = "data" in self.label.lower()
        self.isSignal             = json["isSignal"]
        # nominal values
        self.nEventsAnalyzed      = float(json["nEventsAnalyzed"])
        self.nEventsInFile        = float(json["nEventsInFile"])
        self.nEventsPassTrigger   = float(json["nEventsPassTriggerSelection"])
        self.nEventsPassSelection = float(json["nEventsPassSelection"])
        self.nTagPred             = json["nTagPred"]
        self.nTagTru              = json["nTagTrue"]
        self.norm = 0;
        
        #stats errors on values
        self.nEventsPassTriggerE   = (math.sqrt(json["nEventsPassTriggerSelection"]))
        self.nEventsPassSelectionE = (math.sqrt(json["nEventsPassSelection"]))
        self.nTagPredE             = []
        for ii in self.nTagPred: self.nTagPredE.append(math.sqrt(ii))
        self.nTagTruE             = []
        for ii in self.nTagTru: self.nTagTruE.append(math.sqrt(ii))
        
        self.x_limit_label        = json["x_limit_label"]
        self.y_limit_label        = json["y_limit_label"]
        self.xsec                 = json["xsec"]
        self.nEventsNoSelection   = 100 #self.xsec * 2502 

    def get_observed_array(self):
        str_array = []
        
        if self.isData:
            str_array.append("N/A")
            str_array.append(("%2.0f" % (self.nEventsPassTrigger)) + "")
            str_array.append(("%2.0f" % (self.nEventsPassSelection)) + "")
            str_array.append(("%2.0f" % self.nTagTru[0]) + "")
            #str_array.append("%5i" % (self.nEventsInFile))
        else:
            str_array.append("%5i" % (self.nEventsAnalyzed))
            str_array.append(("%2.1f \pm %2.2f" % (self.nEventsPassTrigger, self.nEventsPassTriggerE)) + "\\%")
            str_array.append(("%2.1f \pm %2.2f" % (self.nEventsPassSelection, self.nEventsPassSelectionE)) + "\\%")
            str_array.append(("%2.1f \pm %2.2f" % (self.nTagTru[0]*self.norm, math.sqrt(self.nTagTru[0]) * self.norm )) + "\\%")

        # blind the data
        if self.isData:
            str_array.append("Blinded")
            str_array.append("Blinded")
            str_array.append("Blinded")
        else:
            str_array.append(("%2.1f \pm %2.2f" % (sum(self.nTagTru[1:])*self.norm, math.sqrt(sum(self.nTagTru[1:])) * self.norm)) + "\\%")
            str_array.append(("%2.1f \pm %2.2f" % (sum(self.nTagTru[2:])*self.norm, math.sqrt(sum(self.nTagTru[2:])) * self.norm)) + "\\%")
            str_array.append(("%2.1f \pm %2.2f" % (sum(self.nTagTru[3:])*self.norm, math.sqrt(sum(self.nTagTru[3:])) * self.norm)) + "\\%")
            
        # str_array.append("%s", self.nTagTrue[1])
        # str_array.append("%s", self.nTagTrue[1])
        print str_array
        return str_array

    def get_pred_array(self):
        str_array = []
        
        if self.isData:
            str_array.append("N/A")
        else:
            str_array.append("%5i" % (self.nEventsAnalyzed))
        str_array.append("%2.1f" % self.nEventsPassTrigger)
        str_array.append("%2.1f" % self.nEventsPassSelection)
        str_array.append("%2.1f" % (self.nTagPred[0]))

        print "NTAGPRED", self.nTagPred
        print "NTAGPRED EXACT1", (self.nTagPred[0])
        print "NTAGPRED SUM1", sum(self.nTagPred[1:])
        print "NTAGPRED SUM2", sum(self.nTagPred[2:])
        
        str_array.append("%2.1f" % sum(self.nTagPred[1:]))
        str_array.append("%2.1f" % sum(self.nTagPred[2:]))
        str_array.append("%2.1f" % sum(self.nTagPred[3:]))
        print str_array,"\n"
        return str_array


    def get_error_array(self):
        str_array = []
        
        if self.isData:
            str_array.append("N/A")
        else:
            str_array.append(0)
        str_array.append("%2.1f" % self.nEventsPassTriggerE)
        str_array.append("%2.1f" % self.nEventsPassSelectionE)
        str_array.append("%2.1f" % (self.nTagPredE[0]))

        print "NTAGPRED", self.nTagPred
        print "NTAGPRED EXACT1", (self.nTagPred[0])
        print "NTAGPRED SUM1", sum(self.nTagPred[1:])
        print "NTAGPRED SUM2", sum(self.nTagPred[2:])
        
        str_array.append("%2.1f" % sum(self.nTagPred[1:]))
        str_array.append("%2.1f" % sum(self.nTagPred[2:]))
        str_array.append("%2.1f" % sum(self.nTagPred[3:]))
        print str_array,"\n"
        return str_array


        
    #normalize the signal and background MC by the appropriate value
    def normalize_mc(self):
        self.norm = 1. / (self.nEventsInFile if self.isSignal else self.nEventsAnalyzed) * self.nEventsNoSelection
        #errors
        self.nEventsInFile =  self.nEventsNoSelection
        self.nEventsPassSelection = float(self.nEventsPassSelection) * self.norm
        self.nEventsPassTrigger = float(self.nEventsPassTrigger) * self.norm

        # dont normalize the tags as they need to be normalized during the sums later
        # ----------
        # for ntag in range(len(self.nTagPred)):
        #     self.nTagPred[ntag] = float(self.nTagPred[ntag]) * self.norm
        # for ntag in range(len(self.nTagTru)):
        #     self.nTagTru[ntag] = float(self.nTagTru[ntag]) * self.norm

        # nominal values 
        self.nEventsInFile =  self.nEventsNoSelection
        self.nEventsPassSelectionE = float(self.nEventsPassSelectionE) * self.norm
        self.nEventsPassTriggerE = float(self.nEventsPassTriggerE) * self.norm
        # for ntag in range(len(self.nTagPredE)):
        #     self.nTagPredE[ntag] = float(self.nTagPredE[ntag]) * self.norm
        # for ntag in range(len(self.nTagTruE)):
        #     self.nTagTruE[ntag] = float(self.nTagTruE[ntag]) * self.norm
            
    def add_column(self, column2):
        self.nEventsInFile        += column2.nEventsInFile
        self.nEventsPassSelection += column2.nEventsPassSelection
        for ntag in self.nTagPred:
            self.nTagPred[ntag] += column2.nTagPred[ntag]
        for ntag in self.nTagTrue:
            self.nTagTrue[ntag] += column2.nTagTrue[ntag]

class table:
    def __init__(self, name):
        self.name      = name
        self.columns   = []
        # names of the cuts  considered 
        self.row_names = [ "Total Events", "Trigger", "Event Sel.", "= 1 Tag","2+ Tags", "3+ Tags", "4+ Tags"]
        #storage by cut
        self.rows      = [] 
    
    # turn the columns into rows for latex
    def parse_to_rows(self):        
        xlabelrow = [x_unit] 
        # get the names of the samples
        for col in self.columns:
            xlabelrow.append(col.x_limit_label)
        self.rows.append(xlabelrow)

        ylabelrow = [y_unit] 
        for col in self.columns:
            ylabelrow.append(col.y_limit_label)
        self.rows.append(ylabelrow)

        #parse out the rows from the columns in pairs of obs, pred
        for rowNum in range(len(self.row_names)):
            row = [self.row_names[rowNum]]                                    
            for column in self.columns:
                obs_array = column.get_observed_array()
                pred_array = column.get_pred_array()
                row.append((obs_array[rowNum]))
                #row.append((obs_array[rowNum],pred_array[rowNum]))
            print row
            self.rows.append(row)

        return self.rows

    # add a column into the table
    def add_column(self, col): 
        self.columns.append(col)
                           
columns = []
qcd_columns = []
signal_columns = []
#parse the json and build the columns
for file_name in json_names:
    json_data  = open(file_name, "r").read()
    thisJson   = json.loads(json_data)
    thisColumn = column(thisJson)

    if not thisColumn.isData: thisColumn.normalize_mc() 
        
    # check if the result is a qcd column
    if "qcd" in thisColumn.label.lower():
        qcd_columns.append(thisColumn)
    elif "data" in thisColumn.label.lower():
        columns.append(thisColumn)
    else:
        signal_columns.append(thisColumn)

# qcd_total_column = qcd_columns[0]
# for column in qcd_columns[1:]: qcd_total_column.add_column(column)
# columns.append(qcd_total_column)

output_table = table("Baseline")

#add the columns into the table
for column in columns: output_table.add_column(column)
#output_table.add_column(qcd_total_column)
for column in signal_columns: output_table.add_column(column)

output_rows = output_table.parse_to_rows()

for ii in output_rows: print ii
    
#now build the latex output
table_string = "\\begin{tabular}{|"
#write the top column row |c|c|c|c|
for ii in range(len(output_table.columns)+1): 
    table_string += "c|"

table_string += "} \n \hline \n"
table_string += "\\multicolumn{%i}{|c|}{\\textbf{Displaced SUSY}} \n \hline \n" % (len(output_table.columns)+1)

# use the first of the output rows to label the table
# labelrow = output_rows[0]
# for ii in labelrow{-1: table_string += ii + " & "
# table_string += "\\\\ \n \hline \n"

# use the first of the output rows to label the table
for row in output_rows:
    for item in row[:-1]: table_string += "\t %s &" %  item 
    index = output_rows.index(row)
    table_string += "\t %s \\\\ \n \hline \n" % row[-1]
    if index == 2: table_string += "\hline \n"

table_string += "\\end{tabular}"

print "\n\n\n\n", table_string
