import ROOT
import time
import argparse
import numpy as np
from math import pi
from array import array
import glob
import os

max_num = 1000
                
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

#Create a new ROOT fill
output = ROOT.TFile('Photons.root', 'RECREATE')
#Create a new ROOT TTree
eventTree = ROOT.TTree('eventTree', 'List of different variables in events from VBF_HToInvisible dataset')
eventTree.Branch('nPhoton', nPhoton, 'nPhoton/I')
eventTree.Branch('photonPhi', photonPhi, 'photonPhi[nPhoton]/F')
eventTree.Branch('photonPt', photonPt, 'photonPt[nPhoton]/F')
eventTree.Branch('photonEta', photonEta, 'photonEta[nPhoton]/F')
photons, photonLabel = Handle('std::vector<pat::Photon>'), 'slimmedPhotons'


def writeTree(inputFile):

        print('######## Creating branches ########')

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

dirname = '/cms/data/store/user/dsperka/lowmassdiphoton/ggh_m10'

if __name__ == '__main__':

        for filename in os.listdir(dirname):

                if 'MiniAOD' in filename:

                        writeTree(os.path.join(dirname,filename))

print('Part 1 Finished')

def readTree(inputFile):

        f = ROOT.TFile.Open(inputFile,"update")

        num_photons=ROOT.TH1F('numphotons', 'Number of photons', 4, 0, 4)
        num_photons.GetXaxis().SetTitle('Number of photons')
        num_photons.GetYaxis().SetTitle('Number of Events')


        photon_phi=ROOT.TH1F('photonphi', 'Phi', 48, -3., 3.)
        photon_phi.GetXaxis().SetTitle('Phi')
        photon_phi.GetYaxis().SetTitle('Number of Events')

        photon_pt=ROOT.TH1F('photonpt', 'Pt', 75, 0., 150.)
        photon_pt.GetXaxis().SetTitle('Pt')
        photon_pt.GetYaxis().SetTitle('Number of Events')

        photon_eta=ROOT.TH1F('photoneta', 'Eta', 48, -3., 3.)
        photon_eta.GetXaxis().SetTitle('Eta')
        photon_eta.GetYaxis().SetTitle('Number of Events')

        event_count_before = 0
        event_count_after = 0

        for event in f.eventTree:

                # Reading the branches of eventTree             
                nPhoton = event.nPhoton
                photonPhi = event.photonPhi
                photonPt = event.photonPt
                photonEta = event.photonEta

                num_photons.Fill(nPhoton)
                for Phi in photonPhi:

                        photon_phi.Fill(Phi)

                for Pt in photonPt:

                        photon_pt.Fill(Pt)

                for Eta in photonEta:

                        photon_eta.Fill(Eta)

        num_photons.Draw()
        photon_phi.Draw()
        photon_pt.Draw()
        photon_eta.Draw()
        f.Write()
if __name__ == '__main__':

        inputFile = 'Photons.root'
        readTree(inputFile)
