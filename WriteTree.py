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
photonEnergyg = array('f', np.zeros(max_num, dtype=float))
photonPxg = array('f', np.zeros(max_num, dtype=float))
photonPyg = array('f', np.zeros(max_num, dtype=float))
photonPzg = array('f', np.zeros(max_num, dtype=float))
numPhotons_gen = array('i', [0])
photonEtag = array('f', np.zeros(max_num, dtype=float))
photonPhig = array('f', np.zeros(max_num, dtype=float))
photonPtg = array('f', np.zeros(max_num, dtype=float))
nPhoton_L1 = array('i', [0])
l1photonPt = array('f', np.zeros(max_num, dtype=float))
l1photonPhi = array('f', np.zeros(max_num, dtype=float))
l1photonEta = array('f', np.zeros(max_num, dtype=float))
l1photonEnergy = array('f', np.zeros(max_num, dtype=float))
l1photonPx = array('f', np.zeros(max_num, dtype=float))
l1photonPy = array('f', np.zeros(max_num, dtype=float))
l1photonPz = array('f', np.zeros(max_num, dtype=float))
hltPhoton20 = array('i', [0])
hltPhoton33 = array('i', [0])
hltP20H = array('i', [0])
hltP30H = array('i', [0])
hltDP30 = array('i', [0])

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
eventTree.Branch('photonEnergyg', photonEnergyg, 'photonEnergyg[nParticles]/F')
eventTree.Branch('photonPxg', photonPxg, 'photonPxg[nParticles]/F')
eventTree.Branch('photonPyg', photonPyg, 'photonPyg[nParticles]/F')
eventTree.Branch('photonPzg', photonPzg, 'photonPzg[nParticles]/F')
eventTree.Branch('numPhotons_gen', numPhotons_gen, 'numPhotons_gen/I')
eventTree.Branch('numPhotons_gen', numPhotons_gen, 'numPhotons_gen/I')
eventTree.Branch('photonEnergyg', photonEnergyg, 'photonEnergyg[numPhotons_gen]/F')
eventTree.Branch('photonPxg', photonPxg, 'photonPxg[numPhotons_gen]/F')
eventTree.Branch('photonPyg', photonPyg, 'photonPyg[numPhotons_gen]/F')
eventTree.Branch('photonPzg', photonPzg, 'photonPzg[numPhotons_gen]/F')
eventTree.Branch('photonEtag', photonEtag, 'photonEtag[numPhotons_gen]/F')
eventTree.Branch('photonPhig', photonPhig, 'photonPhig[numPhotons_gen]/F')
eventTree.Branch('photonPtg', photonPtg, 'photonPtg[numPhotons_gen]/F')
eventTree.Branch('nPhoton_L1', nPhoton_L1, 'nPhoton_L1/I')
eventTree.Branch('l1photonPhi', l1photonPhi, 'l1photonPhi[nPhoton_L1]/F')
eventTree.Branch('l1photonPt', l1photonPt, 'l1photonPt[nPhoton_L1]/F')
eventTree.Branch('l1photonEta', l1photonEta, 'l1photonEta[nPhoton_L1]/F')
eventTree.Branch('l1photonEnergy', l1photonEnergy, 'l1photonEnergy[nPhoton_L1]/F')
eventTree.Branch('l1photonPx', l1photonPx, 'l1photonPx[nPhoton_L1]/F')
eventTree.Branch('l1photonPy', l1photonPy, 'l1photonPy[nPhoton_L1]/F')
eventTree.Branch('l1photonPz', l1photonPz, 'l1photonPz[nPhoton_L1]/F')
eventTree.Branch('hltPhoton20', hltPhoton20, 'hltPhoton20/I')
eventTree.Branch('hltPhoton33', hltPhoton33, 'hltPhoton33/I')
eventTree.Branch('hltP20H', hltP20H, 'hltP20H/I')
eventTree.Branch('hltP30H', hltP30H, 'hltP30H/I')
eventTree.Branch('hltDP30', hltDP30, 'hltDP30/I')

photons, photonLabel = Handle('std::vector<pat::Photon>'), 'slimmedPhotons'
genParticles, genParticlesLabel = Handle('std::vector<reco::GenParticle>'), 'prunedGenParticles'
triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults","","HLT")
photonL1, photonL1Label = Handle("BXVector<l1t::EGamma>"), "caloStage2Digis:EGamma"


def writeTree(inputFile):

        print('######## Creating branches ########')

        events = Events(inputFile)

        print('Took the input file successfully')

        for i, event in enumerate(events):

                event.getByLabel(photonLabel, photons)
                event.getByLabel(genParticlesLabel, genParticles)
                event.getByLabel(triggerBitLabel, triggerBits)
                event.getByLabel(photonL1Label, photonL1)
                
                numPhotons_gen[0] = 0
                
                triggerBits_ = triggerBits.product()
                photons_ = photons.product()
                genParticles_ = genParticles.product()

                nPhoton[0] = len(photons_)
                nParticles[0] = len(genParticles_)

                names = event.object().triggerNames(triggerBits_)

                for i in range(triggerBits_.size()):

                        if names.triggerNames()[i]=='HLT_Photon20_v2':
                                if triggerBits_.accept(i):
                                        hltPhoton20[0] = 1
                                else:
                                        hltPhoton20[0] = 0

                        if names.triggerNames()[i]=='HLT_Photon33_v5':
                                if triggerBits_.accept(i):
                                        hltPhoton33[0] = 1
                                else:
                                        hltPhoton33[0] = 0

                        if names.triggerNames()[i]=='HLT_Photon20_HoverELoose_v10':
                                if triggerBits_.accept(i):
                                        hltP20H[0] = 1
                                else:
                                        hltP20H[0] = 0

                        if names.triggerNames()[i]=='HLT_Photon30_HoverELoose_v10':
                                if triggerBits_.accept(i):
                                        hltP30H[0] = 1
                                else:
                                        hltP30H[0] = 0

                        if names.triggerNames()[i]=='HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_v2':
                                if triggerBits_.accept(i):
                                        hltDP30[0] = 1
                                else:
                                        hltDP30[0] = 0
                                        
                for i, prt in enumerate(genParticles_):

                        pdgId[i] = prt.pdgId()
                        particleStatus[i] = prt.status()
                        nMothers[i] = prt.numberOfMothers()

                        if prt.pdgId()==22:

                                photonEnergyg[i] = prt.energy()
                                photonPxg[i] = prt.px()
                                photonPyg[i] = prt.py()
                                photonPzg[i] = prt.pz()
                                numPhotons_gen[0] += 1
                                #print('in if loop')
                                
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

                bxVector_photon = photonL1.product()

                bx=0

                nPhoton_L1[0] = bxVector_photon.size(bx)

                for i in range(bxVector_photon.size(bx)):

                        photon = bxVector_photon.at(bx, i)

                        l1photonPt[i] = photon.pt()
                        l1photonEta[i] = photon.eta()
                        l1photonPhi[i] = photon.phi()
                        l1photonEnergy[i] = photon.energy()
                        l1photonPx[i] = photon.px()
                        l1photonPy[i] = photon.py()
                        l1photonPz[i] = photon.pz()
                        
                eventTree.Fill()

        print('Number of trigger paths: %d' % triggerBits.product().size())

        names = event.object().triggerNames(triggerBits.product())

        for i in range(triggerBits.product().size()):

                print('Trigger ', names.triggerNames()[i], ('PASS' if triggerBits.product().accept(i) else 'FAIL'))

#/cms/data/store/user/dsperka/lowmassdiphoton/ggh_m10

dirname = '/home/peterspencer/CMSWork/ggh_m30'

if __name__ == '__main__':

        for filename in os.listdir(dirname):

                if 'MiniAOD' in filename:

                        writeTree(os.path.join(dirname,filename))

        output.Write()
        output.Close()
