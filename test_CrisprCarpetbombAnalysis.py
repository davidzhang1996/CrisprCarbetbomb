import unittest
import CrisprCarpetbombAnalysis as crispr 
from CrisprCarpetbombAnalysis import getCase

class TestGetCasesMethod(unittest.TestCase):

    def test_spliceNormal(self):
        edits, strand, numgRNAs=[[0,3],[4,4]], "+", 5
        self.assertEqual(getCase(edits, strand, numgRNAs), [3,0,0,3,2])

class TestDetermineNearbygRNAS(unittest.TestCase): 
    def test_oneIndel(self): 
        tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange=[50,90,170], [0,70,100], [[40,60],[120,140]], 2, 15
        output, strand= [2,1], "+"
        self.assertEqual(crispr.determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange, strand), output)

    def test_oneSplice(self): 
        tEnd, tStart=[50,300], [0,110]
        gRNAChrPos, numgRNAs, spliceRange = [[40,60],[120,140],[200,220]], 3, 15
        output, strand=[3,3,1], "+"
        self.assertEqual(crispr.determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange, strand), output)



    def test_multipleSpliceMultipleIndel(self):
        tEnd, tStart,  =[50,155,200,250], [0,105,160,210] 
        gRNAChrPos, numgRNAs, spliceRange = [[40,60],[100,120], [160,180], [200, 220]], 4, 15
        output,strand=[3,3,2,2], "+"
        self.assertEqual(crispr.determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange, strand), output)
    
    def test_nonconsequetiveSplices(self):
        tEnd, tStart,  =[100,150, 350], [0,110, 290] 
        gRNAChrPos, numgRNAs, spliceRange = [[40,60],[100,120], [160,180], [200, 220] , [300, 320]], 5, 15
        output,strand=[1,2,3,0,3],"+"
        self.assertEqual(crispr.determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange, strand), output)
    
    def test_noncomprehensive(self): 
        tEnd, tStart,  =[100,150], [0,110] 
        gRNAChrPos, numgRNAs, spliceRange = [[40,60],[100,120], [160,180], [200, 220] , [300, 320]], 5, 15
        output,strand=[1,2,3,4,4],"+"
        self.assertEqual(crispr.determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange, strand), output)

    def test_bigDeletion(self): 
        tEnd, tStart,  =[100,140], [0,110] 
        gRNAChrPos, numgRNAs, spliceRange = [[40,60],[100,120], [160,180], [200, 220] , [300, 320]], 5, 15
        output,strand=[1,2,4,4,4],"+"
        self.assertEqual(crispr.determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange, strand), output)        

    def test_reverseSplicing(self): 
        tEnd, tStart,  =[250,140], [500,170] 
        gRNAChrPos, numgRNAs, spliceRange = [[40,60],[100,120], [160,180], [240, 260] , [300, 320]], 5, 15
        output,strand=[4,4,3,3,1],"-"
        self.assertEqual(crispr.determineNearbygRNAs(tEnd, tStart, gRNAChrPos, numgRNAs, spliceRange,strand), output)        

if __name__ == '__main__':
    unittest.main()