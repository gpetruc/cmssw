import ROOT
ROOT.gROOT.SetBatch(True)

from collections import defaultdict

from sys import argv
for filename in argv[1:]:
    tfile = ROOT.TFile.Open(filename)
    tree = tfile.Get("pixelOccupancyByLumi/tree")
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        print "%d %d %d %d %d" % (tree.run, tree.lumi, tree.detid, tree.tpm, tree.roc);
    tfile.Close()

