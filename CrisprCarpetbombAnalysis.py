import sys

def main ():
    lineCount=1 
    blatfilename= sys.argv[1]
    gRNAfilename= sys.argv[2]

    misalignedQNames, unexpectedlyModifiedQNames=[],[]
    sizeIndel= {-15:0, -14:0, -13:0, -12:0, -11:0, -10:0, -9:0, -8:0, -7:0, -6:0, -5:0, -4:0, -3:0, -2:0, -1:0, 0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0}
    caseOccurences = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0, 17:0, 18:0, 19:0, 20:0, 21:0, 22:0, 23:0, 24:0, 25:0, 26:0, 27:0, 28:0, 29:0, 30:0, 31:0, 32:0, 33:0, 34:0}
    gRNAChrPos= getgRNAS(gRNAfilename)
    #print gRNAChrPos
    with open(blatfilename, "r") as file:
        for line in file:

            #print "line", lineCount
            misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences, = analyzeLine(line, gRNAChrPos, misalignedQNames, unexpectedlyModifiedQNames,  sizeIndel, caseOccurences)
            lineCount+=1

    print "sizeIndel"
    print sizeIndel

    print 'caseOccurences'
    print caseOccurences

    #print 'misalignedQNames'
    #print misalignedQNames

def analyzeLine(line, gRNAChrPos, misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences):
    lineList=line.split()
    strand, qName, blockSize, qStart, tStart, qEnd, tEnd=lineList[8], lineList[9], map(int, lineList[18].split(",")[:-1]), map(int, lineList[19].split(",")[:-1]), map(int, lineList[20].split(",")[:-1]), [], []
    numBlocks, numQStarts, numTStarts, numgRNAs= len(blockSize), len(qStart), len(tStart), len(gRNAChrPos)
    spliceRange, indelRange= 20,15

    #print "qStart", "tStart", "blockSize", "qEnd", "tEnd"
    #print qStart, tStart, blockSize, qEnd, tEnd

    if (numBlocks!= numQStarts and numQStarts!= numTStarts):
        return misalignedQNames.append(qName), unexpectedlyModifiedQNames, sizeIndel, caseOccurences
    for i in range (0, numBlocks):
        tEnd.append(tStart[i]+blockSize[i])
        qEnd.append(qStart[i]+blockSize[i])
        
    #print "blockSize"
    #print blockSize

    #print "tStart"
    #print tStart



    #print "gRNAChrPos"
    #print gRNAChrPos

    i=0
    while(i<numQStarts): #Helps remove artifacts before the start of first guide RNA
        if (tEnd[i]< gRNAChrPos[0][0]-15):
            #print "tEnd1", tEnd
            #print "numQStarts1", numQStarts


            qStart, tStart, blockSize = qStart[1:], tStart[1:], blockSize[1:]
            qEnd, tEnd= qEnd[1:], tEnd[1:]
            numQStarts, numTStarts, numBlocks= numQStarts-1, numTStarts-1, numBlocks-1

            #print "tEnd2", tEnd
            #print "numQStarts2", numQStarts
        else:
            i+=1

    #print "blockSize"
    #print blockSize

    if (numBlocks==1 or strand == "-"): #If there is only one block then there was no editing of genome. Also exception if read starts in the middle of the exome
        return misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences
            
    else: 
        for i in range(0, numBlocks-1):
            nearbygRNAs= determineNearbygRNAs(tEnd[i], tStart[i+1], gRNAChrPos, numgRNAs, spliceRange)
            if (nearbygRNAs[0] != None and nearbygRNAs[1] != None):
                if (nearbygRNAs[0] != nearbygRNAs[1]): #gRNA1-gRNA2 del
                    if (nearbygRNAs[1]-nearbygRNAs[0]==1): 
                        caseOccurences[11]+=1 
                    elif(nearbygRNAs[1]-nearbygRNAs[0]==2): 
                        caseOccurences[12]+=1 
                    elif(nearbygRNAs[1]-nearbygRNAs[0]==3): 
                        caseOccurences[13]+=1

                elif (qEnd[i] == qStart[i+1] and tEnd[i]==tStart[i+1]):
                    return misalignedQNames.append(qName), unexpectedlyModifiedQNames, sizeIndel, caseOccurences
                elif (qEnd[i] == qStart[i+1] and tEnd[i] < tStart[i+1] and (tStart[i+1]-tEnd[i])<indelRange): #deletion
                    #print "deletion"
                    caseOccurences= findCase(nearbygRNAs[])
                    if (nearbygRNAs[0]==0): 
                        caseOccurences[1]+=1
                    elif(nearbygRNAs[0]==1): 
                        caseOccurences[2]+=1
                    elif(nearbygRNAs[0]==2): 
                        caseOccurences[3]+=1 
                    elif(nearbygRNAs[0]==3): 
                        caseOccurences[4]+=1
                    sizeIndel[-(tStart[i+1]-tEnd[i])]+=1
                elif (qEnd[i]<qStart[i+1] and tEnd[i]==tStart[i+1] and (qStart[i+1]-qEnd[i])<indelRange): #insertion
                    #print "insertion"
                    if (nearbygRNAs[0]==0): 
                        caseOccurences[1]+=1
                    elif(nearbygRNAs[0]==1): 
                        caseOccurences[2]+=1
                    sizeIndel[(qStart[i+1]-qEnd[i])]+=1 
                elif (qEnd[i]<qStart[i+1] and tEnd[i]<tStart[i+1] and (qStart[i+1]-qEnd[i])<indelRange and (tStart[i+1]-tEnd[i])<indelRange): 
                    if (nearbygRNAs[0]==0): 
                        caseOccurences[1]+=1
                    elif(nearbygRNAs[0]==1): 
                        caseOccurences[2]+=1
                    caseOccurences[1]+=1
                    sizeIndel[0]+=1

    return misalignedQNames, unexpectedlyModifiedQNames, sizeIndel, caseOccurences

def determineNearbygRNAs(tEndSite, tStartSite, gRNAChrPos, numgRNAs, spliceRange): 
    startgRNA, stopgRNA= None,None 

    for gInd in range (0, numgRNAs): 
        if (tEndSite > gRNAChrPos[gInd][0]-spliceRange and tEndSite < gRNAChrPos[gInd][1]+spliceRange):
            startgRNA=gInd
        if(tStartSite > gRNAChrPos[gInd][0]-spliceRange and tStartSite < gRNAChrPos[gInd][1]+spliceRange): 
            stopgRNA=gInd 
    return [startgRNA, stopgRNA]

def getgRNAS(gRNAfilename): 
    gRNAChrPos=[]
    with open (gRNAfilename, "r") as file: 
        for line in file: 
            line=map(int, line.split()[1:3])
            gRNAChrPos.append([line[0], line[1]])
    return gRNAChrPos


main()  


