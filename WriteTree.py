import ROOT
import time
import argparse
import numpy as np
from math import pi
from array import array
import glob
import os

max_num = 1000

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

nPhoton = array('i', [0])
photonPt = array('f', np.zeros(max_num, dtype=float))
photonPhi = array('f', np.zeros(max_num, dtype=float))
photonEta = array('f', np.zeros(max_num, dtype=float))
pdgId = array('i', np.zeros(max_num, dtype=int))
nParticles = array('i', [0])
particleStatus = array('i', np.zeros(max_num, dtype=int))
nMothers = array('i', np.zeros(max_num, dtype=int))
mothers = array('i', np.zeros(max_num, dtype=int))
photonEnergy = array('f', np.zeros(max_num, dtype=float))
photonPx = array('f', np.zeros(max_num, dtype=float))
photonPy = array('f', np.zeros(max_num, dtype=float))
photonPz = array('f', np.zeros(max_num, dtype=float))

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--test', help = 'Only go over the first 100 events for testing', action = 'store_true')
args = parser.parse_args()

#Create a new ROOT fill
output = ROOT.TFile('Photongraphs.root', 'RECREATE')
#Create a new ROOT TTree
eventTree = ROOT.TTree('eventTree', 'List of different variables in events from VBF_HToInvisible dataset')
eventTree.Branch('nPhoton', nPhoton, 'nPhoton/I')
eventTree.Branch('photonPhi', photonPhi, 'photonPhi[nPhoton]/F')
eventTree.Branch('photonPt', photonPt, 'photonPt[nPhoton]/F')
eventTree.Branch('photonEta', photonEta, 'photonEta[nPhoton]/F')
eventTree.Branch('nParticles', nParticles, 'nParticles/I')
eventTree.Branch('pdgId', pdgId, 'pdgId[nParticles]/I')
eventTree.Branch('particleStatus', particleStatus, 'particleStatus[nParticles]/I')
eventTree.Branch('nMothers', nMothers, 'nMothers[nParticles]/I')
eventTree.Branch('mothers', mothers, 'mothers[nParticles]/I')
eventTree.Branch('photonEnergy', photonEnergy, 'photonEnergy[nPhoton]/F')
eventTree.Branch('photonPx', photonPx, 'photonPx[nPhoton]/F')
eventTree.Branch('photonPy', photonPy, 'photonPy[nPhoton]/F')
eventTree.Branch('photonPz', photonPz, 'photonPz[nPhoton]/F')

photons, photonLabel = Handle('std::vector<pat::Photon>'), 'slimmedPhotons'
genParticles, genParticlesLabel = Handle('std::vector<reco::GenParticle>'), 'prunedGenParticles'
triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")


def writeTree(inputFile):

        print('######## Creating branches ########')

        events = Events(inputFile)

        print('Took the input file successfully')

        for i, event in enumerate(events):

                event.getByLabel(photonLabel, photons)
                event.getByLabel(genParticlesLabel, genParticles)
                event.getByLabel(triggerBitLabel, triggerBits)

                photons_ = photons.product()
                genParticles_ = genParticles.product()

                nPhoton[0] = len(photons_)
                nParticles[0] = len(genParticles_)

                for i, prt in enumerate(genParticles_):

                        pdgId[i] = prt.pdgId()
                        particleStatus[i] = prt.status()
                        nMothers[i] = prt.numberOfMothers()

                        if i > 1:

                                mothers[i] = prt.mother(0).pdgId()

                for i, ph in enumerate(photons_):

                        photonPt[i] = ph.pt()
                        photonPhi[i] = ph.phi()
                        photonEta[i] = ph.eta()
                        photonEnergy[i] = ph.energy()
                        photonPx[i] = ph.px()
                        photonPy[i] = ph.py()
                        photonPz[i] = ph.pz()

                eventTree.Fill()

        print('Number of trigger paths: %d' % triggerBits.product().size())

        names = event.object().triggerNames(triggerBits.product())

        for i in range(triggerBits.product().size()):

                print('Trigger ', names.triggerNames()[i], ('PASS' if triggerBits.product().accept(i) else 'FAIL'))

dirname = '/cms/data/store/user/dsperka/lowmassdiphoton/ggh_m10'

if __name__ == '__main__':

        for filename in os.listdir(dirname):

                if 'MiniAOD' in filename:

                        writeTree(os.path.join(dirname,filename))

        output.Write()
        output.Close()
