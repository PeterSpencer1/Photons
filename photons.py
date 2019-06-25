import  ROOT
import time
import argparse
import numpy as np
from math import pi
from array import array

max_num = 1000

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

nPhoton = array('i', [0])
photonPt = array('f', np.zeros(max_num, dtype=float))
photonPhi = array('f', np.zeros(max_num, dtype=float))
photonEta = array('f', np.zeros(max_num, dtype=float))


# load FWlite python libraries
from DataFormats.FWLite import Handle, Events
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--test', help = 'Only go over the first 100 events for testing', action = 'store_true')
args = parser.parse_args()

def writeTree(inputFile):

        print('######## Creating branches ########')

        #Create a new ROOT fill
        output = ROOT.TFile('Photons.root', 'RECREATE')

        #Create a new ROOT TTree
        eventTree = ROOT.TTree('eventTree', 'List of different variables in events from VBF_HToInvisible dataset')

        eventTree.Branch('nPhoton', nPhoton, 'nPhoton/I')
        eventTree.Branch('photonPhi', photonPhi, 'photonPhi/F')
        eventTree.Branch('photonPt', photonPt, 'photonPt/F')
        eventTree.Branch('photonEta', photonEta, 'photonEta/F')

        photons, photonLabel = Handle('std::vector<pat::Photon>'), 'slimmedPhotons'

        events = Events(inputFile)

        print('Took the input file successfully')

        for i, event in enumerate(events):

                event.getByLabel(photonLabel, photons)

                photons_ = photons.product()

                nPhoton[0] = len(photons_)

                for i, ph in enumerate(photons_):

                        photonPt[i] = ph.pt()
                        photonPhi[i] = ph.phi()
                        photonEta[i] = ph.eta()

                eventTree.Fill()

        #Save the output root file
        output.Write()

if __name__ == '__main__':

        inputFile = '/cms/data/store/user/dsperka/lowmassdiphoton/ggh_m10/MiniAOD_1393595.64.root'

        writeTree(inputFile)
