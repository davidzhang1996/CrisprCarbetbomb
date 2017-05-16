import sys

def main ():
    filename= sys.argv[1]
    print filename
    misalignedQNames, unexpectedlyModifiedQNames=[],[]
    sizeIndel= {-15:0, -14:0, -13:0, -12:0, -11:0, -10:0, -9:0, -8:0, -7:0, -6:0, -5:0, -4:0, -3:0, -2:0, -1:0, 0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0}
    caseOccurences = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0, 18:0, 19:0, 20:0, 21:0, 22:0, 23:0, 24:0, 25:0, 26:0, 27:0, 28:0, 29:0, 30:0, 31:0, 32:0, 33:0, 34:0}
    gRNAChrPos= [[43776408, 43776428], [43776468, 43776488], [43777048, 43777068], [43776907, 43776887]]
    with open(filename, "r") as file:
        for line in file:
            misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences, = analyzeLine(line, gRNAChrPos, misalignedQNames, unexpectedlyModifiedQNames,  sizeIndel, caseOccurences)

    print "sizeIndel"
    print sizeIndel

    print 'caseOccurences'
    print caseOccurences

    print 'misalignedQNames'
    print misalignedQNames

def analyzeLine(line, gRNAChrPos, misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences):
    lineList=line.split()
    #print "lineList 18"
    #print lineList[18]
    strand, qName, blockSize, qStart, tStart, qEnd, tEnd=lineList[8], lineList[9], map(int, lineList[18].split(",")[:-1]), map(int, lineList[19].split(",")[:-1]), map(int, lineList[20].split(",")[:-1]), [], []
    numBlocks, numQStarts, numTStarts= len(blockSize), len(qStart), len(tStart)
    spliceRange, indelRange= 15,15

    #print "qStart", "tStart", "blockSize", "qEnd", "tEnd"
    #print qStart, tStart, blockSize, qEnd, tEnd

    if (numBlocks!= numQStarts and numQStarts!= numTStarts):
        return misalignedQNames.append(qName), unexpectedlyModifiedQNames, sizeIndel, caseOccurences
    for i in range (0, numBlocks):
        tEnd.append(tStart[i]+blockSize[i])
        qEnd.append(qStart[i]+blockSize[i])
        
        
    i=0
    while(i<numQStarts): #Helps remove artifacts before the start of first guide RNA
        if (tStart[i]< gRNAChrPos[0][0]-15):
            qStart, tStart, blockSize= qStart[1:], tStart[1:], blockSize[1:]
            numQStarts, numTStarts, numBlocks= numQStarts-1, numTStarts-1, numBlocks-1
        else:
            i+=1

    if (numBlocks==1 or strand == "-"): #If there is only one block then there was no editing of genome. Also exception if read starts in the middle of the exome
        return misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences
            
    else: 
        for i in range(0, numBlocks-1):
            if ((tEnd)>(gRNAChrPos[i][0]-spliceRange) and (tEnd)<(gRNAChrPos[i][1]+spliceRange) and (tStart[i+1]) > (gRNAChrPos[i][0]-spliceRange) and (tStart[i+1]) < (gRNAChrPos[i][1]+spliceRange)): #gRNA1-gRNA2 del
                caseOccurences[11]+=1
            elif (qEnd == qStart[i+1] and tEnd==tStart[i+1]):
                return misalignedQNames.append(qName), unexpectedlyModifiedQNames, sizeIndel, caseOccurences
            elif (qEnd[i] == qStart[i+1] and tEnd[i] < tStart[i+1] and (tStart[i+1]-tEnd[i])<indelRange): #deletion
                print "deletion"
                caseOccurences[1]+=1
                sizeIndel[-(tStart[i+1]-tEnd[i])]+=1
            elif (qEnd[i]<qStart[i+1] and tEnd[i]==tStart[i+1] and (qStart[i+1]-qEnd[i])<indelRange): #insertion
                print "insertion"
                caseOccurences[1]+=1
                sizeIndel[-(qStart[i+1]-qEnd[i])]+=1 
            elif (qEnd[i]<qStart[i+1] and tEnd[i]<tStart[i+1] and (qStart[i+1]-qEnd[i])<indelRange and (tStart[i+1]-tEnd[i])<indelRange): 
                print "indel"
                caseOccurences[1]+=1
                sizeIndel[0]+=1

    return misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences



main()  


