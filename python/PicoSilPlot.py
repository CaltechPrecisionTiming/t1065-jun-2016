from optparse import OptionParser
import ROOT as rt
from array import array
import numpy as np
from math import sqrt
import sys

if __name__ == '__main__':
    parser = OptionParser()

    parser.add_option('-c','--cut',dest="cut",type="string",default='',
                  help="cut to apply before averaging, e.g. to look only at event 101, do event==101; ring0int = 1 central hexagon, ring1int = 6 central hexagons, ring2int = 18 central hexagons")
    
    (options,args) = parser.parse_args()

    f = rt.TFile.Open(args[0])

    tree = f.Get("pulse")

    empty = rt.TH2D("empty","empty",100, -6, 6, 100, -6, 6)
    
    h2p = rt.TH2Poly()
    h2p.SetName("picosil")
    h2p.SetTitle("picosil")

    xVertex = array('d',[1, 0.5, -0.5, -1, -0.5, 0.5])
    yVertex = array('d',[0, -sqrt(3.)/2.,-sqrt(3.)/2., 0, sqrt(3.)/2., sqrt(3.)/2.])

    centerPositions = [[-3, -sqrt(3)], [-3, 0], [-3, sqrt(3)],
                       [-1.5, -5*sqrt(3)/2], [-1.5, -3*sqrt(3)/2], [-1.5, -sqrt(3)/2], [-1.5, sqrt(3)/2], [-1.5, 3*sqrt(3)/2], [-1.5, 5*sqrt(3)/2],
                       [0, -3*sqrt(3)], [0, -2*sqrt(3)], [0, -1*sqrt(3)],
                       [0, 0],
                       [0, sqrt(3)], [0, 2*sqrt(3)], [0, 3*sqrt(3)],
                       [1.5, -5*sqrt(3)/2], [1.5, -3*sqrt(3)/2], [1.5, -sqrt(3)/2], [1.5, sqrt(3)/2], [1.5, 3*sqrt(3)/2], [1.5, 5*sqrt(3)/2],
                       [3, -sqrt(3)], [3, 0], [3, sqrt(3)]
                       ]


    mapArrayToCenterPos = {
        29: [-3, 0], 
        16: [-3, -sqrt(3)], 
        24: [-1.5, 5*sqrt(3)/2],   
        6: [-1.5, sqrt(3)/2],
        5: [-1.5, -sqrt(3)/2],  
        15: [-1.5, -3*sqrt(3)/2],       
        21: [-1.5, -5*sqrt(3)/2],  
        25: [0, 3*sqrt(3)],        
        7: [0, 2*sqrt(3)],   
        31: [0, sqrt(3)],      
        34: [0,0],        
        30: [0,-sqrt(3)],      
        14: [0,-2*sqrt(3)],   
        22: [0,-3*sqrt(3)],        
        28: [1.5, 5*sqrt(3)/2],   
        3: [1.5, 3*sqrt(3)/2], 
        11: [1.5, sqrt(3)/2],       
        4: [1.5, -sqrt(3)/2],     
        13: [1.5, -3*sqrt(3)/2], 
        23: [1.5, -5*sqrt(3)/2],       
        20: [3, sqrt(3)],   
        19: [3, 0],
        12: [3, -sqrt(3)]
        }

    ring0Index = []
    ring1Index = []
    ring2Index = []
    for key, val in mapArrayToCenterPos.iteritems():
        if np.linalg.norm(np.array(val)) <= 0: ring0Index.append(key)
        if np.linalg.norm(np.array(val)) <= sqrt(3): ring1Index.append(key)
        if np.linalg.norm(np.array(val)) <= 2*sqrt(3): ring2Index.append(key)


    for centerPos in centerPositions:
        x2 = array('d',[xi+centerPos[0] for xi in xVertex])
        y2 = array('d',[yi+centerPos[1] for yi in yVertex])
        h2p.AddBin(6,x2,y2)


    options.cut = options.cut.replace('ring0int','+'.join(['int[%i]'%index for index in ring0Index]))
    options.cut = options.cut.replace('ring1int','+'.join(['int[%i]'%index for index in ring1Index]))
    options.cut = options.cut.replace('ring2int','+'.join(['int[%i]'%index for index in ring2Index]))
                                    
    nevents = tree.Draw('>>elist',options.cut,'entrylist')
        
    elist = rt.gDirectory.Get('elist')
    
    entry = -1
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)
        for key, val in mapArrayToCenterPos.iteritems():
            h2p.Fill(val[0],val[1], tree.amp[key])
        #print tree.event
        
    h2p.Scale(1.0/nevents)
    
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetPaintTextFormat("3.3f")

    c = rt.TCanvas('c','c',500,500)
    c.SetRightMargin(0.2)
        
    empty.Draw("")
    h2p.GetZaxis().SetTitle("Amplitude [V]")
    h2p.GetZaxis().SetTitleOffset(2)
    h2p.SetMaximum(0.1)
    h2p.Draw("colztextsame")
    
    l = rt.TLatex()
    l.SetTextAlign(22)
    l.SetTextSize(0.03)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.5,0.94,"cut: %s #rightarrow %i events"%(options.cut,nevents))
    
    c.Print('picosil.pdf')
    c.Print('picosil.C')
    
    f.Close()
