#python replaceGeneIDwithGeneName.py ./mergedtranscripts/merged.gtf ./mergedtranscripts/mergedID.gtf

import sys

def main():
    outfile = open(sys.argv[2],'w')
    for line in open(sys.argv[1],'r'):
        newlinelist = line.strip().split("\t")
        oldattributes = newlinelist[8].split("; ")
        attributeDict = {}
        for attribute in oldattributes:
            attributeDict[attribute.split(" ")[0]]=attribute.split(" ")[1]
        newattributes = oldattributes
        newattributes[0] = "gene_id "+attributeDict["gene_name"]
        newlinelist[8] = "; ".join(newattributes)
        newline = "\t".join(newlinelist)
        #print newline
        outfile.write(newline+"\n")

if __name__=="__main__":
    main()
