import os,re 
import pandas as pd
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='compine depth file')

parser.add_argument('-i', '--dir', help='depth dirpath')
parser.add_argument('-o','--output_file', help='output file.')


args = parser.parse_args()
dir = args.dir
out_file = args.output_file


 
#dir = "./"
#out_file = "test"
file = os.listdir(dir)
file = [i for i in file if ".depth" in i]
file_1 = pd.read_csv(file[0],sep = "\t")
file_1 = file_1[["contigName","totalAvgDepth"]]
file_1["totalAvgDepth"] = file_1["totalAvgDepth"]/file_1["totalAvgDepth"].sum()
#sum_scaler = lambda x:(x)/sum(x)
file_1 = file_1.rename(columns={'totalAvgDepth':re.sub(".depth","",file[0])})
file.pop(0)
for i in range(len(file)):
        file_2 = pd.read_csv(file[i],sep = "\t")
        file_2 = file_2[["contigName","totalAvgDepth"]]
        file_2["totalAvgDepth"] = file_2["totalAvgDepth"]/file_2["totalAvgDepth"].sum()
        file_2 = file_2.rename(columns={'totalAvgDepth':re.sub(".depth","",file[i])})
        file_1 = pd.merge(file_1,file_2,on="contigName")


file_1.to_csv(out_file+"_contig.profile",sep = "\t",index=False)


def nth_repl(s, sub, repl, nth):
    find = s.find(sub)
    # if find is not p1 we have found at least one match for the substring
    i = find != -1
    # loop util we find the nth or we find no match
    while find != -1 and i != nth:
        # find + 1 means we start at the last match start index + 1
        find = s.find(sub, find + 1)
        i += 1
    # if i  is equal to nth we found nth matches so replace
    if i == nth:
        return s[:find]+repl+s[find + len(sub):]
    return s
def re_sub(x):
        list_n = []
        for i in x:
                i = re.sub(r'^([^_]*\_[^_]*)_.*',r'\1',i)
                list_n.append(i)
        return list_n

file_1["contigName"] = re_sub(file_1["contigName"])
file_2 = file_1.groupby(file_1['contigName']).median()

file_2.to_csv(out_file+"_SGB.profile",sep = "\t")
