import sys

def main(): 
	largestMatch={}
	repeatMatches=sys.argv[2]
	blatFile=sys.argv[1]
	blatFileExtensionless=blatFile[0:blatFile.index(".")]
	outputFile= blatFileExtensionless+ "NoDup.psl"
	with open(repeatMatches, "r") as file: 
		for line in file: 
			lineList=line.split()
			readName, matches= lineList[0], lineList[1:]
			matches = [int(element) for element in matches]
			#print "matches", matches
			largestMatch[readName]=max(matches[1:])
			#print "largestMatch", largestMatch[readName] 

	with open(blatFile, "r") as file1, open(outputFile, "w") as file2:
		count=0
		for line in file1: 
			lineList1=line.split()
			readName=lineList1[9]
			numMatches=int(lineList1[0])
			#print len(largestMatch)
			if (readName not in largestMatch):
				file2.write(line)
			elif(numMatches == largestMatch[readName]):
				count+=1
				file2.write(line)
				
		print count



def getLargest(matches): 
	largest=0 
	numberMatches=len(matches)
	for i in xrange(1,numberMatches): 
		if (matches[i]>largest): 
			largest=matches[i]
	return largest

main()