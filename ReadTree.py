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

        recoMatchingPhotons_Pt=ROOT.TH1F('recoMatchingPhotons_Pt','Pt of Matching Photons (Reco vs Gen)', 50, 0, 50.)
        recoMatchingPhotons_Pt.GetXaxis().SetTitle('Pt')
        recoMatchingPhotons_Pt.GetYaxis().SetTitle('Number of Events')

        recoMatchingPhotons_Energy=ROOT.TH1F('recoMatchingPhotonsEnergy','Energy of Matching Photons (Reco vs Gen)', 50, 0, 40)
        recoMatchingPhotons_Energy.GetXaxis().SetTitle('Energy')
        recoMatchingPhotons_Energy.GetYaxis().SetTitle('Number of Events')

        recoNum_matchingPhotons=ROOT.TH1F('recoNummatchingPhotons', 'Number of Matching Photons (Reco vs Gen)', 10, 0, 10)
        recoNum_matchingPhotons.GetXaxis().SetTitle('Number of Matching Photons')
        recoNum_matchingPhotons.GetYaxis().SetTitle('Number Of Events')

        etaSpecial=ROOT.TH1F('etaSpecial', 'Difference in Eta Divided By Gen', 50, -2, 2)
        etaSpecial.GetXaxis().SetTitle('Eta Difference')
        etaSpecial.GetYaxis().SetTitle('Number Of Events')

        phiSpecial=ROOT.TH1F('phiSpecial', 'Difference in Phi Divided By Gen', 50, -2, 2)
        phiSpecial.GetXaxis().SetTitle('Phi Difference')
        phiSpecial.GetYaxis().SetTitle('Number Of Events')

        energySpecial=ROOT.TH1F('energySpecial', 'Difference in Energy Divided By Gen', 50, -2, 3)
        energySpecial.GetXaxis().SetTitle('Energy Difference')
        energySpecial.GetYaxis().SetTitle('Number Of Events')

        delta_r=ROOT.TH1F('deltar', 'Delta R', 20, 0, 10)
        delta_r.GetXaxis().SetTitle('Delta R')
        delta_r.GetYaxis().SetTitle('Number Of Events')

        min_r=ROOT.TH1F('minr', 'Lowest Value of Delta R', 40, 0, 10)
        min_r.GetXaxis().SetTitle('Delta R')
        min_r.GetYaxis().SetTitle('Number Of Events')

        delta_reco=ROOT.TH1F('deltareco', 'Delta R (Reco)', 20, 0, 10)
        delta_reco.GetXaxis().SetTitle('Delta R')
        delta_reco.GetYaxis().SetTitle('Number Of Events')

        min_reco=ROOT.TH1F('minreco', 'Lowest Value of Delta R (Reco)', 40, 0, 10)
        min_reco.GetXaxis().SetTitle('Delta R')
        min_reco.GetYaxis().SetTitle('Number Of Events')
        
        event_count_before = 0
        event_count_after = 0

        counterDictBef = {}
        counterDictAft = {}

        counterDictBef20 = {}
        counterDictAft20 = {}

        for pt in range(5, 65):

                counterDictBef20[pt] = 0
                counterDictAft20[pt] = 0

        counterDictBef33 = {}
        counterDictAft33 = {}

        for pt in range(5, 65):

                counterDictBef33[pt] = 0
                counterDictAft33[pt] = 0

        counterDictBefP20 = {}
        counterDictAftP20 = {}

        for pt in range(5, 65):

                counterDictBefP20[pt] = 0
                counterDictAftP20[pt] = 0

        counterDictBefP30 = {}
        counterDictAftP30 = {}

        for pt in range(5, 65):

                counterDictBefP30[pt] = 0
                counterDictAftP30[pt] = 0

        counterDictBefDP30 = {}
        counterDictAftDP30 = {}

        for pt in range(5, 65):

                counterDictBefDP30[pt] = 0
                counterDictAftDP30[pt] = 0
        
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
                hltPhoton20 = event.hltPhoton20
                hltPhoton33 = event.hltPhoton33
                hltP20H = event.hltP20H
                hltP30H = event.hltP30H
                hltDP30 = event.hltDP30
                
                matchingPhotonsPt = []
                matchingPhotonsEnergy = []
                matchingPhotons = 0
                index = []
                
                for pt in range(5, 65):

                        if len(photonPt) != 0:

                                if pt < photonPt[0] < pt+1:
                                        counterDictBef20[pt] += 1
                                        if hltPhoton20 == 1:
                                                counterDictAft20[pt] += 1

                for pt in range(5, 65):

                        if len(photonPt) != 0:

                                if pt < photonPt[0] < pt+1:
                                        counterDictBef33[pt] += 1
                                        if hltPhoton33 == 1:
                                                counterDictAft33[pt] += 1

                for pt in range(5, 65):

                        if len(photonPt) != 0:

                                if pt < photonPt[0] < pt+1:
                                        counterDictBefP20[pt] += 1
                                        if hltP20H == 1:
                                                counterDictAftP20[pt] += 1

                for pt in range(5, 65):

                        if len(photonPt) != 0:

                                if pt < photonPt[0] < pt+1:
                                        counterDictBefP30[pt] += 1
                                        if hltP30H == 1:
                                                counterDictAftP30[pt] += 1

                for pt in range(5, 65):

                        if len(photonPt) != 0:

                                if pt < photonPt[0] < pt+1:
                                        counterDictBefDP30[pt] += 1
                                        if hltDP30 == 1:
                                                counterDictAftDP30[pt] += 1
                
                for i in range(numPhotons_gen):

                        for j in range(nPhoton_L1):

                                eta_diff = l1photonEta[j] - photonEtag[i]
                                phi_diff = l1photonPhi[j] - photonPhig[i]
                                energy_diff = l1photonEnergy[j] - photonEnergyg[i]
                                deltaR = math.sqrt(eta_diff**2 + phi_diff**2)
                                delta_r.Fill(deltaR)

                                if deltaR < .5:

                                        matchingPhotons += 1
                                        index.append(j)
                                        eta_difference = eta_diff/photonEtag[i]
                                        etaSpecial.Fill(eta_difference)
                                        phi_difference = phi_diff/photonPhig[i]
                                        phiSpecial.Fill(phi_difference)
                                        energy_difference = energy_diff/photonEnergyg[i]
                                        energySpecial.Fill(energy_difference)
                                        #delta_r.Fill(deltaR)
                                        
                if matchingPhotons==0:

                        for h in range(nPhoton_L1):

                                delta_R = []

                                for u in range(numPhotons_gen):

                                        if photonPtg[u] < 2: continue

                                        eta_diff = l1photonEta[h] - photonEtag[u]
                                        phi_diff = l1photonPhi[h] - photonPhig[u]
                                        deltaR = math.sqrt(eta_diff**2 + phi_diff**2)
                                        delta_R.append(deltaR)

                                if len(delta_R) != 0:

                                        r = min(delta_R)

                                        min_r.Fill(r)
                                       
                if matchingPhotons==2:

                        for j in index:

                                matchingPhotonsPt.append(l1photonPt[j])
                                matchingPhotonsEnergy.append(l1photonEnergy[j])


                for mEnergy in matchingPhotonsEnergy:

                        matchingPhotons_Energy.Fill(mEnergy)

                for mPt in matchingPhotonsPt:

                        matchingPhotons_Pt.Fill(mPt)

                num_matchingPhotons.Fill(matchingPhotons)
                
                recoMatchingPhotonsPt = []
                recoMatchingPhotonsEnergy = []
                recoMatchingPhotons = 0
                recoIndex = []
                
                for o in range(numPhotons_gen):

                        for k in range(nPhoton):

                                eta_diff = photonEta[k] - photonEtag[o]
                                phi_diff = photonPhi[k] - photonPhig[o]
                                deltaR = math.sqrt(eta_diff**2 + phi_diff**2)
                                energy_diff = photonEnergy[k] - photonEnergyg[o]
                                delta_reco.Fill(deltaR)

                                if deltaR < .5:

                                        recoMatchingPhotons += 1
                                        recoIndex.append(k)
                                        eta_difference = eta_diff/photonEtag[o]
                                        etaSpecialReco.Fill(eta_difference)
                                        phi_difference = phi_diff/photonPhig[o]
                                        phiSpecialReco.Fill(phi_difference)
                                        energy_difference = energy_diff/photonEnergyg[o]
                                        energySpecialReco.Fill(energy_difference)

                if matchingPhotons==0:

                        for h in range(nPhoton):

                                delta_Reco = []

                                for u in range(numPhotons_gen):

                                        if photonPtg[u] < 2: continue

                                        eta_diff = photonEta[h] - photonEtag[u]
                                        phi_diff = photonPhi[h] - photonPhig[u]
                                        deltaReco = math.sqrt(eta_diff**2 + phi_diff**2)
                                        delta_Reco.append(deltaR)

                                if len(delta_Reco) != 0:

                                        r = min(delta_Reco)

                                        min_reco.Fill(r)
                                        
                if recoMatchingPhotons==2:

                        for k in recoIndex:

                                recoMatchingPhotonsPt.append(photonPt[k])
                                recoMatchingPhotonsEnergy.append(photonEnergy[k])


                for rmEnergy in recoMatchingPhotonsEnergy:

                        recoMatchingPhotons_Energy.Fill(rmEnergy)

                for rmPt in recoMatchingPhotonsPt:

                        recoMatchingPhotons_Pt.Fill(rmPt)

                recoNum_matchingPhotons.Fill(recoMatchingPhotons)
                
                if nPhoton==2:

                        px = photonPx[0] + photonPx[1]
                        py = photonPy[0] + photonPy[1]
                        pz = photonPz[0] + photonPz[1]
                        energy = photonEnergy[0] + photonEnergy[1]
                        invariantMass = energy**2 - px**2 - py**2 - pz**2
                        invariant_mass.Fill(invariantMass)
                        etadifference = abs(photonEta[0] - photonEta[1])
                        photon_etadiff.Fill(etadifference)

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
                        
                for PhiL1 in l1photonPhi:

                        l1photon_phi.Fill(PhiL1)

                for PtL1 in l1photonPt:

                        l1photon_pt.Fill(PtL1)

                for EtaL1 in l1photonEta:

                        l1photon_eta.Fill(EtaL1)
         
        ratio20 = []

        for i in sorted(counterDictAft20.keys()):

                ratio20.append(float(counterDictAft20[i])/counterDictBef20[i])

        x_axis20 = []

        for num in range(5, 65):

                x_axis20.append(float(num))

        ratio20 = array('f', ratio20)
        x_axis20 = array('f', x_axis20)

        n20 = len(x_axis20)
        ptAccept20=ROOT.TGraph(n20, x_axis20, ratio20)
        ptAccept20.GetXaxis().SetTitle('Photon Pt')
        ptAccept20.GetYaxis().SetTitle('Acceptance')
        ptAccept20.Draw('AC')
        ptAccept20.Write('ptAccept20')

        ratio33 = []

        for i in sorted(counterDictAft33.keys()):

                ratio33.append(float(counterDictAft33[i])/counterDictBef33[i])

        x_axis33 = []

        for num in range(5, 65):

                x_axis33.append(float(num))

        ratio33 = array('f', ratio33)
        x_axis33 = array('f', x_axis33)

        n33 = len(x_axis33)
        ptAccept33=ROOT.TGraph(n20, x_axis33, ratio33)
        ptAccept33.GetXaxis().SetTitle('Photon Pt')
        ptAccept33.GetYaxis().SetTitle('Acceptance')
        ptAccept33.Draw('AC')
        ptAccept33.Write('ptAccept33')
        
        ratioP20 = []

        for i in sorted(counterDictAftP20.keys()):

                ratioP20.append(float(counterDictAftP20[i])/counterDictBefP20[i])

        x_axisP20 = []

        for num in range(5, 65):

                x_axisP20.append(float(num))

        ratioP20 = array('f', ratioP20)
        x_axisP20 = array('f', x_axisP20)

        nP20 = len(x_axisP20)
        ptAcceptP20=ROOT.TGraph(nP20, x_axisP20, ratioP20)
        ptAcceptP20.GetXaxis().SetTitle('Photon Pt')
        ptAcceptP20.GetYaxis().SetTitle('Acceptance')
        ptAcceptP20.Draw('AC')
        ptAcceptP20.Write('ptAcceptP20')

        ratioP30 = []

        for i in sorted(counterDictAftP30.keys()):

                ratioP30.append(float(counterDictAftP30[i])/counterDictBefP30[i])

        x_axisP30 = []

        for num in range(5, 65):

                x_axisP30.append(float(num))

        ratioP30 = array('f', ratioP30)
        x_axisP30 = array('f', x_axisP30)

        nP30 = len(x_axisP30)
        ptAcceptP30=ROOT.TGraph(nP30, x_axisP30, ratioP30)
        ptAcceptP30.GetXaxis().SetTitle('Photon Pt')
        ptAcceptP30.GetYaxis().SetTitle('Acceptance')
        ptAcceptP30.Draw('AC')
        ptAcceptP30.Write('ptAcceptP30')
        
        ratioDP30 = []

        for i in sorted(counterDictAftDP30.keys()):

                ratioDP30.append(float(counterDictAftDP30[i])/counterDictBefDP30[i])

        x_axisDP30 = []

        for num in range(5, 65):

                x_axisDP30.append(float(num))

        ratioDP30 = array('f', ratioDP30)
        x_axisDP30 = array('f', x_axisDP30)

        nDP30 = len(x_axisDP30)
        ptAcceptDP30=ROOT.TGraph(nDP30, x_axisDP30, ratioDP30)
        ptAcceptDP30.GetXaxis().SetTitle('Photon Pt')
        ptAcceptDP30.GetYaxis().SetTitle('Acceptance')
        ptAcceptDP30.Draw('AC')
        ptAcceptDP30.Write('ptAcceptDP30')
        
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
        recoMatchingPhotons_Pt.Draw()
        recoMatchingPhotons_Energy.Draw()
        recoNum_matchingPhotons.Draw()
        f.Write()

if __name__ == '__main__':

        inputFile = 'Photongraphs.root'
        readTree(inputFile)
