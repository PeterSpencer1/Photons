import ROOT
import math
from array import array

def readTree(inputFile):

        f = ROOT.TFile.Open(inputFile,"update")

        num_photons=ROOT.TH1F('numphotons', 'Number of photons', 4, 0, 4)
        num_photons.GetXaxis().SetTitle('Number of photons')
        num_photons.GetYaxis().SetTitle('Number of Events')

        photon_phi=ROOT.TH1F('photonphi', 'Phi', 150, -3., 3.)
        photon_phi.GetXaxis().SetTitle('Phi')
        photon_phi.GetYaxis().SetTitle('Number of Events')

        photon_pt=ROOT.TH1F('photonpt', 'Pt', 150, 0., 150.)
        photon_pt.GetXaxis().SetTitle('Pt')
        photon_pt.GetYaxis().SetTitle('Number of Events')

        photon_eta=ROOT.TH1F('photoneta', 'Eta', 150, -3., 3.)
        photon_eta.GetXaxis().SetTitle('Eta')
        photon_eta.GetYaxis().SetTitle('Number of Events')

        invariant_mass=ROOT.TH1F('invariantmass', 'Invariant Mass Of 2 Photons', 50, 50, 150)
        invariant_mass.GetXaxis().SetTitle('Mass (GeV)')
        invariant_mass.GetYaxis().SetTitle('Number of Events')

        invariant_massg=ROOT.TH1F('invariantmassg', 'Invariant Mass Of 2 Photons (Gen)', 50, 7, 13)
        invariant_massg.GetXaxis().SetTitle('Mass (GeV)')
        invariant_massg.GetYaxis().SetTitle('Number of Events')
        
        photon_phig=ROOT.TH1F('photonphig', 'Photon Phi Difference (Gen)', 75, 0, 4.)
        photon_phig.GetXaxis().SetTitle('Phi')
        photon_phig.GetYaxis().SetTitle('Number of Events')

        photon_etag=ROOT.TH1F('photonetag', 'Photon Eta Difference (Gen)', 75, 0, 4.)
        photon_etag.GetXaxis().SetTitle('Eta')
        photon_etag.GetYaxis().SetTitle('Number of Events')

        photon_etadiff=ROOT.TH1F('photonetadiff', 'Photon Eta Difference', 75, 0, 4.)
        photon_etadiff.GetXaxis().SetTitle('Eta')
        photon_etadiff.GetYaxis().SetTitle('Number of Events')

        photon_ptg=ROOT.TH1F('photonptg', 'Photon Pt (Gen)', 50, 0, 50.)
        photon_ptg.GetXaxis().SetTitle('Pt')
        photon_ptg.GetYaxis().SetTitle('Number of Events')

        l1photon_phi=ROOT.TH1F('l1photonphi', 'Phi', 50, -3., 3.)
        l1photon_phi.GetXaxis().SetTitle('Phi (L1)')
        l1photon_phi.GetYaxis().SetTitle('Number of Events')

        l1photon_pt=ROOT.TH1F('l1photonpt', 'Pt', 50, 0., 50.)
        l1photon_pt.GetXaxis().SetTitle('Pt (L1)')
        l1photon_pt.GetYaxis().SetTitle('Number of Events')

        l1photon_eta=ROOT.TH1F('l1photoneta', 'Eta', 50, -3., 3.)
        l1photon_eta.GetXaxis().SetTitle('Eta (L1)')
        l1photon_eta.GetYaxis().SetTitle('Number of Events')

        l1invariant_mass=ROOT.TH1F('l1invariantmass', 'Invariant Mass Of 2 Photons (L1)', 50, 0, 60)
        l1invariant_mass.GetXaxis().SetTitle('Mass (GeV)')
        l1invariant_mass.GetYaxis().SetTitle('Number of Events')

        l1photon_etadiff=ROOT.TH1F('l1photonetadiff', 'Photon Eta Difference (L1)', 75, 0, 4.)
        l1photon_etadiff.GetXaxis().SetTitle('Eta')
        l1photon_etadiff.GetYaxis().SetTitle('Number of Events')

        numPhoton_L1=ROOT.TH1F('l1numphotons', 'Number of Photons (L1)', 10, 0, 10)
        numPhoton_L1.GetXaxis().SetTitle('Number Of Photons')
        numPhoton_L1.GetYaxis().SetTitle('Number of Events')

        matchingPhotons_Pt=ROOT.TH1F('matchingPhotons_Pt','Pt of Matching Photons (L1 vs Gen)', 50, 0, 50.)
        matchingPhotons_Pt.GetXaxis().SetTitle('Pt')
        matchingPhotons_Pt.GetYaxis().SetTitle('Number of Events')
        matchingPhotons_Energy=ROOT.TH1F('matchingPhotons_Energy','Energy of Matching Photons (L1 vs Gen)', 50, 0, 40)
        matchingPhotons_Energy.GetXaxis().SetTitle('Energy')
        matchingPhotons_Energy.GetYaxis().SetTitle('Number of Events')

        num_matchingPhotons=ROOT.TH1F('num_matchingPhotons', 'Number of Matching Photons (L1 vs Gen)', 10, 0, 10)
        num_matchingPhotons.GetXaxis().SetTitle('Number of Matching Photons')
        num_matchingPhotons.GetYaxis().SetTitle('Number Of Events')

        event_count_before = 0
        event_count_after = 0

        counterDictBef = {}
        counterDictAft = {}

        for pt in range(5, 65, 3):

                counterDictBef[pt] = 0
                counterDictAft[pt] = 0
        
        for event in f.eventTree:

                #reading the branches of eventTree             
                nPhoton = event.nPhoton
                photonPhi = event.photonPhi
                photonPt = event.photonPt
                photonEta = event.photonEta
                photonPx = event.photonPx
                photonPy = event.photonPy
                photonPz = event.photonPz
                photonEnergy = event.photonEnergy
                photonPxg = event.photonPxg
                photonPyg = event.photonPyg
                photonPzg = event.photonPzg
                photonEnergyg = event.photonEnergyg
                pdgId = event.pdgId
                numPhotons_gen = event.numPhotons_gen
                photonPhig = event.photonPhig
                photonEtag = event.photonEtag
                photonPtg = event.photonPtg
                nPhoton_L1 = event.nPhoton_L1
                l1photonPhi = event.l1photonPhi
                l1photonPt = event.l1photonPt
                l1photonEta = event.l1photonEta
                l1photonPx = event.l1photonPx
                l1photonPy = event.l1photonPy
                l1photonPz = event.l1photonPz
                l1photonEnergy = event.l1photonEnergy
                
                matchingPhotonsPt = []
                matchingPhotonsEnergy = []
                matchingPhotons = 0
                
                for pt in range(5, 65, 3):
                        if len(photonPt) != 0:
                                if pt < photonPt[0] < pt+3:
                                        counterDictBef[pt] += 1
                                        if hltPhoton20 == 1:
                                                counterDictAft[pt] += 1
                
                for i in range(numPhotons_gen):

                        for j in range(nPhoton_L1):

                                eta_diff = l1photonEta[j] - photonEtag[i]
                                phi_diff = l1photonPhi[j] - photonPhig[i]
                                deltaR = math.sqrt(eta_diff**2 + phi_diff**2)

                                if deltaR < .5:

                                        matchingPhotons += 1

                                if matchingPhotons==2:

                                        matchingPhotonsPt.append(l1photonPt[j])
                                        matchingPhotonsEnergy.append(l1photonEnergy[j])

                                for mEnergy in matchingPhotonsEnergy:

                                        matchingPhotons_Energy.Fill(mEnergy)

                                for mPt in matchingPhotonsPt:

                                        matchingPhotons_Pt.Fill(mPt)

                num_matchingPhotons.Fill(matchingPhotons)
                
                if nPhoton==2:

                        px = photonPx[0] + photonPx[1]
                        py = photonPy[0] + photonPy[1]
                        pz = photonPz[0] + photonPz[1]
                        energy = photonEnergy[0] + photonEnergy[1]
                        invariantMass = energy**2 - px**2 - py**2 - pz**2
                        invariant_mass.Fill(invariantMass)

                photonEnergygen = [value for value in photonEnergyg if value != 0]
                photonPxgen = [value for value in photonPxg if value != 0]
                photonPygen = [value for value in photonPyg if value != 0]
                photonPzgen = [value for value in photonPzg if value != 0]
                
                if nPhoton_L1==2:

                        px = l1photonPx[0] + l1photonPx[1]
                        py = l1photonPy[0] + l1photonPy[1]
                        pz = l1photonPz[0] + l1photonPz[1]
                        energy = l1photonEnergy[0] + l1photonEnergy[1]
                        l1invariantMass = math.sqrt(energy**2 - px**2 - py**2 - pz**2)
                        l1invariant_mass.Fill(l1invariantMass)
                        l1etadiff = abs(l1photonEta[0] - l1photonEta[1])
                        l1photon_etadiff.Fill(l1etadiff)

                numPhoton_L1.Fill(nPhoton_L1)

                pt_lower_10 = False

                for pt in photonPtg:

                        photon_ptg.Fill(pt)

                        if pt < 10:

                                pt_lower_10 = True

                if pt_lower_10: continue

                if numPhotons_gen==2:

                        pxg = photonPxg[0] + photonPxg[1]
                        pyg = photonPyg[0] + photonPyg[1]
                        pzg = photonPzg[0] + photonPzg[1]
                        energyg = photonEnergyg[0] + photonEnergyg[1]
                        invariantMassg = math.sqrt(energyg**2 - pxg**2 - pyg**2 - pzg**2)
                        invariant_massg.Fill(invariantMassg)
                        etadiff = abs(photonEtag[0] - photonEtag[1])
                        phidiff = abs(photonPhig[0] - photonPhig[1])
                        photon_etag.Fill(etadiff)
                        photon_phig.Fill(phidiff)

                num_photons.Fill(nPhoton)

                for Phi in photonPhi:

                        photon_phi.Fill(Phi)

                for Pt in photonPt:

                        photon_pt.Fill(Pt)

                for Eta in photonEta:

                        photon_eta.Fill(Eta)
                        
                for Phi in l1photonPhi:

                        l1photon_phi.Fill(Phi)

                for Pt in l1photonPt:

                        l1photon_pt.Fill(Pt)

                for Eta in l1photonEta:

                        l1photon_eta.Fill(Eta)
         
        ratio = []

        for i in counterDictBef.keys():

                ratio.append(float(counterDictAft[i])/counterDictBef[i])

        x_axis = []

        for num in range(5, 65, 3):

                x_axis.append(float(num))

        ratio = array('f', ratio)
        x_axis = array('f', x_axis)

        n = len(x_axis)
        ptAccept=ROOT.TGraph(n, x_axis, ratio)
        ptAccept.Draw('AC')
        ptAccept.Write('ptAccept')
        
        num_photons.Draw()
        photon_phi.Draw()
        photon_pt.Draw()
        photon_eta.Draw()
        invariant_mass.Draw()
        invariant_massg.Draw()
        photon_phig.Draw()
        photon_etag.Draw()
        photon_etadiff.Draw()
        l1photon_phi.Draw()
        l1photon_pt.Draw()
        l1photon_eta.Draw()
        l1invariant_mass.Draw()
        l1photon_etadiff.Draw()
        numPhoton_L1.Draw()
        matchingPhotons_Pt.Draw()
        matchingPhotons_Energy.Draw()
        num_matchingPhotons.Draw()
        f.Write()

if __name__ == '__main__':

        inputFile = 'Photongraphs.root'
        readTree(inputFile)
