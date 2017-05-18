import sys

def main ():
    lineCount=1 
    blatfile1= sys.argv[1]
    blatfile2=sys.argv[2]
    gRNAfilename= sys.argv[3]

    matchInfo = {}
    sizeIndel= {-15:0, -14:0, -13:0, -12:0, -11:0, -10:0, -9:0, -8:0, -7:0, -6:0, -5:0, -4:0, -3:0, -2:0, -1:0, 0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0}
    caseDict = {}
    pairedEndStrandMM=0 
    pairedEndCaseMM=0 

    gRNAChrPos= getgRNAS(gRNAfilename)
    #print gRNAChrPos
    with open(blatfile1, "r") as file1, open(blatfile2, "r") as file2:
        for line in file1:
            lineList=line.split()
            readName=lineList[9]
            matchInfo[readName]=[lineList]
        for line in file2:
            lineList=line.split()
            readName=lineList[9] 
            if(readName in matchInfo): 
                matchInfo[readName].append(lineList)
            #else: 
            #    matchInfo[readName]=[lineList] Will add this feature later!

    for line in matchInfo:
        #print "matchInfo[line]", matchInfo[line]
        if (len(matchInfo[line])==2):
            readOne= matchInfo[line][0]
            readTwo=matchInfo[line][1]
            case, sizeIndel, pairedEndStrandMM, pairedEndCaseMM= analyzeLine(gRNAChrPos, readOne, readTwo, sizeIndel, pairedEndStrandMM, pairedEndCaseMM)
            caseDict= updateCaseList(case, caseDict) 
    #print "pairedEndStrandMM", pairedEndStrandMM, "pairedEndCaseMM", pairedEndCaseMM
    for element in caseDict:
        print element, caseDict[element]
    print "pairedEndCaseMM", pairedEndCaseMM

    #print "sizeIndel"
    #print sizeIndel

    #print 'caseOccurences'
    #print caseOccurences

    #print 'misalignedQNames'
    #print misalignedQNames

def analyzeLine(gRNAChrPos, readOne, readTwo, sizeIndel, pairedEndStrandMM, pairedEndCaseMM):
    r1List=readOne
    r2List=readTwo
    r1Strand, r1QName, r1BlockSize, r1QStart, r1TStart, r1QEnd, r1TEnd=r1List[8], r1List[9], map(int, r1List[18].split(",")[:-1]), map(int, r1List[19].split(",")[:-1]), map(int, r1List[20].split(",")[:-1]), [], []
    r2Strand, r2QName, r2BlockSize, r2QStart, r2TStart, r2QEnd, r2TEnd=r2List[8], r2List[9], map(int, r2List[18].split(",")[:-1]), map(int, r2List[19].split(",")[:-1]), map(int, r2List[20].split(",")[:-1]), [], []
    r1NumBlocks, r1NumQStarts, r1NumTStarts = len(r1BlockSize), len(r1QStart), len(r1TStart)
    r2NumBlocks, r2NumQStarts, r2NumTStarts = len(r2BlockSize), len(r2QStart), len(r2TStart) 
    numgRNAs = len(gRNAChrPos)

    caseTotal=[]
    spliceRange, indelRange= 20,15

    if (r1Strand == r2Strand):
        return caseTotal, sizeIndel, pairedEndStrandMM+1, pairedEndCaseMM
    elif (r1Strand == "-" and r2Strand == "+"): 
        for i in range (0, r1NumBlocks):
            r1TEnd.append(r1TStart[i]-r1BlockSize[i])
        for i in range(0, r2NumBlocks): 
            r2TEnd.append(r2TStart[i]+r2BlockSize[i])  
    elif(r1Strand == "+" and r2Strand == "-"): 
        for i in range (0, r1NumBlocks):
            r1TEnd.append(r1TStart[i]+r1BlockSize[i])
            r1QEnd.append(r1QStart[i]+r1BlockSize[i])

        for i in range(0, r2NumBlocks): 
            r2TEnd.append(r2TStart[i]-r2BlockSize[i])
            r2QEnd.append(r2QStart[i]+r2BlockSize[i])

            

    if (r1NumBlocks == r1NumQStarts and r1NumQStarts == r1NumTStarts and r2NumBlocks==r2NumQStarts and r2NumQStarts==r2NumTStarts):    
        case1=[]
        case2=[]

        #TODO figure out way of updating cases such that it factors in reads that start near guide RNAs

        case1= determineNearbygRNAs(r1TEnd, r1TStart, gRNAChrPos, numgRNAs, spliceRange, r1Strand)
        case2= determineNearbygRNAs(r2TEnd, r2TStart, gRNAChrPos, numgRNAs, spliceRange, r2Strand)    
        
        if (r1Strand=="+" and r2Strand=="-"): 
            caseTotal, pairedEndCaseMM=getCaseTotal(case1, case2, numgRNAs, pairedEndCaseMM)
        elif(r1Strand=="-" and r2Strand=="+"): 
            caseTotal, pairedEndCaseMM=getCaseTotal(case2, case1, numgRNAs, pairedEndCaseMM)

        """if (nearbygRNAs[0] != None and nearbygRNAs[1] != None):

                if (qEnd[i] == qStart[i+1] and tEnd[i]==tStart[i+1]):
                    return misalignedQNames.append(qName), unexpectedlyModifiedQNames, sizeIndel, caseOccurences
                elif (qEnd[i] == qStart[i+1] and tEnd[i] < tStart[i+1] and (tStart[i+1]-tEnd[i])<indelRange): #deletion
                    #print "deletion"
                    sizeIndel[-(tStart[i+1]-tEnd[i])]+=1
                elif (qEnd[i]<qStart[i+1] and tEnd[i]==tStart[i+1] and (qStart[i+1]-qEnd[i])<indelRange): #insertion
                    #print "insertion"
                    sizeIndel[(qStart[i+1]-qEnd[i])]+=1 
                elif (qEnd[i]<qStart[i+1] and tEnd[i]<tStart[i+1] and (qStart[i+1]-qEnd[i])<indelRange and (tStart[i+1]-tEnd[i])<indelRange): 
                    caseOccurences[1]+=1
                    sizeIndel[0]+=1
        """
    return caseTotal, 0, pairedEndStrandMM, pairedEndCaseMM

def getCaseTotal(case1, case2, numgRNAs, pairedEndCaseMM): 
    #print case1, case2
    case1count4=case1.count(4)
    case2count4=case2.count(4)
    case=[0]*numgRNAs
    """for i in range(0,numgRNAs):
        if (case1[i]==4 and case2[i]!=4):
            case[i]=case2[i]
        elif(case1[i]!=4 and case2[i]==4):
            case[i]=case1[i]
        elif(case1[i]==4 and case2[i]==4): 
            case[i]=4
        elif(case1[i]!=4 and case2[i]!=4): 
            pairedEndCaseMM+=1
            if(case1count4>case2count4): 
                case[i]=case2[i]
            else: 
                case[i]=case1[i]
    """
    for i in range(0,numgRNAs): 
        case[i]=case1[i]
    return case, pairedEndCaseMM




def updateCaseList(case, caseDict): 
    caseString=''.join(str(e) for e in case)
    if(caseString not in caseDict): 
        caseDict[caseString]=1
    else: 
        caseDict[caseString]+=1
    return caseDict 

def determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange, strand): 
    gRNAsEdits=[1]*numgRNAs
    startgRNA, stopgRNA= 0,0 
    numTEnd, numTStart=len(tEnd), len(tStart)
    flankedPositions=[]

    for gInd in range (0, numgRNAs): 
        for i in range(0, numTStart-1): 
            tEndSite, tStartSite=tEnd[i], tStart[i+1]
            flanked=False

            if (strand=="+"):
                if (tEndSite > gRNAChrPos[gInd][0]-spliceRange and tEndSite < gRNAChrPos[gInd][1]+spliceRange):
                    gRNAsEdits[gInd]=3
                    flanked=True
                if(tStartSite > gRNAChrPos[gInd][0]-spliceRange and tStartSite < gRNAChrPos[gInd][1]+spliceRange): 
                    if (flanked): gRNAsEdits[gInd]=2
                    else: gRNAsEdits[gInd]=3
                if (gRNAChrPos[gInd][0]-spliceRange > tEndSite  and gRNAChrPos[gInd][1]+spliceRange < tStartSite): #Deleted gRNA 
                    gRNAsEdits[gInd]=0
            elif(strand=="-"): 
                if(tEndSite < gRNAChrPos[numgRNAs-gInd-1][1]+spliceRange and tEndSite > gRNAChrPos[numgRNAs-gInd-1][0]-spliceRange): 
                    gRNAsEdits[numgRNAs-gInd-1]=3
                    flanked=True

                if(tStartSite < gRNAChrPos[numgRNAs-gInd-1][1]+spliceRange and tStartSite > gRNAChrPos[numgRNAs-gInd-1][0]-spliceRange): 
                    if (flanked): gRNAsEdits[numgRNAs-gInd-1]=2
                    else: gRNAsEdits[numgRNAs-gInd-1]=3
                if (gRNAChrPos[numgRNAs-gInd-1][1]+spliceRange < tEndSite  and gRNAChrPos[numgRNAs-gInd-1][0]-spliceRange > tStartSite): #Deleted gRNA 
                    gRNAsEdits[numgRNAs-gInd-1]=0

                 

 
    for i in range(0, numgRNAs): #Determine which reads weren't covered
        if (strand=="+"):
            if (tEnd[numTEnd-1] > gRNAChrPos[i][0]-spliceRange and tEnd[numTEnd-1] < gRNAChrPos[i][1]+spliceRange):
                gRNAsEdits[i]=3
            if (tEnd[numTEnd-1] < gRNAChrPos[i][0]-spliceRange): 
                gRNAsEdits[i]=4 
        elif(strand=="-"): 
            if (tEnd[numTEnd-1] < gRNAChrPos[i][0]-spliceRange and tEnd[numTEnd-1] > gRNAChrPos[i][1]+spliceRange):
                gRNAsEdits[i]=3
            if (tEnd[numTEnd-1] > gRNAChrPos[i][0]-spliceRange): 
                gRNAsEdits[i]=4 


    return gRNAsEdits

def getgRNAS(gRNAfilename): 
    gRNAChrPos=[]
    with open (gRNAfilename, "r") as file: 
        for line in file: 
            line=map(int, line.split()[1:3])
            gRNAChrPos.append([line[0], line[1]])
    return gRNAChrPos

if __name__ == "__main__":
   main()



