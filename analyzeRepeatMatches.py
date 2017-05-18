import sys 

def main():
	totalReads, duplicateReads = {}, {}
	blatFile=sys.argv[1]
	blatFileExtenionless=blatFile[0:blatFile.index(".")]
	outputFile= blatFileExtenionless+ "RepeatMatches.txt"
	#print outputFile
	duplicateReads=analyzeMatches(blatFile, duplicateReads, totalReads)
	#print duplicateReads
	with open (outputFile, "w") as file: 
		for element in duplicateReads:
			line=""
			line=element + "\t"+ str(len(duplicateReads[element])) 
			for index in range(0, len(duplicateReads[element])): 
				line+="\t"+str(duplicateReads[element][index]) 
			line.rstrip()
			file.write(line+"\n")





def analyzeMatches(blatFile, duplicateReads, totalReads): 
	with open (blatFile, "r") as file:
		duplicateReadNumber=0
		for line in file: 
			readName, numMatches, strand =line.split()[9], line.split()[0], line.split()[8]
			if(readName not in totalReads): 
				totalReads[readName]=[numMatches]
			else:
				duplicateReadNumber+=1
				totalReads[readName].append(numMatches)
	for element in totalReads: 
		if (len(totalReads[element])!=1): 
			duplicateReads[element]=totalReads[element]

	#print "duplicateReadNumber", duplicateReadNumber
	return duplicateReads	


main()
