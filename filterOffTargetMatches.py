import sys
def main(): 
	blatFile=sys.argv[1]
	blatFileExtensionless=blatFile[0:blatFile.index(".")]
	outputFile=blatFileExtensionless+"Targeted.psl"

	geneStart, geneStop= 43776250, 43795340
	with open(blatFile, "r") as inFile, open(outputFile, "w") as outFile: 
		count=0
		for line in inFile:
			lineList=line.split()
			chrm, strand=lineList[13], lineList[8]
			chrmStart, chrmStop, readLen= map(int,[lineList[15], lineList[16], lineList[10]])
			if (chrm=="chr15"): 
				if (strand=="+" and chrmStart>geneStart-readLen and chrmStop<geneStop+readLen): 
					count+=1
					outFile.write(line)
				elif(strand=="-" and chrmStart<geneStop+readLen and chrmStop>geneStart-readLen): 
					count+=1
					outFile.write(line)
		print "count", count

main()