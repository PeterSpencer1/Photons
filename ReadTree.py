import ROOT
import math
from array import array

def drawTriggerEff(inputFile, trigger):

        f = ROOT.TFile.Open(inputFile, 'UPDATE')
        if trigger == 'hltPhoton20':

                photon_array = array('f', [0., 10., 15., 17., 18., 19., 20., 21., 22., 23., 25., 27.,  30., 40., 50., 70.])
                triggerTitle = 'HLT_Photon20_v2'

        elif trigger == 'hltPhoton33':

                photon_array = array('f', [0., 10., 20., 25., 28., 30., 31., 32., 33., 34., 35., 36., 38., 41., 45., 50., 60.])
                triggerTitle = 'HLT_Photon33_v5'

        elif trigger == 'hltP20H':

                photon_array = array('f', [0., 10., 15., 17., 18., 19., 20., 21., 22., 23., 25., 27.,  30., 40., 50., 70.])
                triggerTitle = 'HLT_Photon20_HoverELoose_v10'

        elif trigger == 'hltP30H':

                photon_array = array('f', [0., 10., 20., 25., 27., 28., 29., 30., 31., 32., 33., 35., 37., 40., 50., 60., 80.])
                triggerTitle = 'HLT_Photon30_HoverELoose_v10'

        elif trigger == 'hltDP30':

                photon_array = array('f', [0., 10., 20., 25., 27., 28., 29., 30., 31., 32., 33., 35., 37., 40., 50., 60., 80.])
                triggerTitle = 'HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_v2'

        photonHist = ROOT.TH1F('photonHist', 'photonHist', len(photon_array)-1, photon_array)

        photonHistTrigger = ROOT.TH1F('photonHistTrigger', 'photonHistTrigger', len(photon_array)-1, photon_array)

        triggerCuts = trigger + ' == 1'

        f.eventTree.Draw('photonPt[0]>>photonHist')
        f.eventTree.Draw('photonPt[0]>>photonHistTrigger', triggerCuts, '')

        if ROOT.TEfficiency.CheckConsistency(photonHistTrigger, photonHist):

                photonEff = ROOT.TEfficiency(photonHistTrigger, photonHist)
                photonEff.SetTitle(triggerTitle + " Efficiency;Leading Photon Pt (GeV);Trigger Efficiency")
                photonEff.Write('photonEff' + trigger + '_photon')

                print('Efficiency graph constructed')

        f.Write()
        f.Close()

def readTree(inputFile):

        f = ROOT.TFile.Open(inputFile,"update")

        num_photons=ROOT.TH1F('numphotons', 'Number of photons', 10, 0, 10)
        num_photons.GetXaxis().SetTitle('Number of photons')
        num_photons.GetYaxis().SetTitle('Number of Events')

        photon_phi=ROOT.TH1F('photonphi', 'Phi', 150, -3., 3.)
        photon_phi.GetXaxis().SetTitle('Phi')
        photon_phi.GetYaxis().SetTitle('Number of Events')

        photon_PtTotal=ROOT.TH1F('totalphotonpt', 'Total Photon Pt', 80, 0., 80.)
        photon_PtTotal.GetXaxis().SetTitle('Pt (GeV)')
        photon_PtTotal.GetYaxis().SetTitle('Number of Events')

        leadingPhoton_pt=ROOT.TH1F('leadingphotonpt', 'Leading Photon Pt', 80, 0., 80.)
        leadingPhoton_pt.GetXaxis().SetTitle('Pt (GeV)')
        leadingPhoton_pt.GetYaxis().SetTitle('Number of Events')

        trailingPhoton_pt=ROOT.TH1F('trailingphotonpt', 'Trailing Photon Pt', 80, 0., 80.)
        trailingPhoton_pt.GetXaxis().SetTitle('Pt (GeV)')
        trailingPhoton_pt.GetYaxis().SetTitle('Number of Events')

        photon_eta=ROOT.TH1F('photoneta', 'Eta', 150, -3., 3.)
        photon_eta.GetXaxis().SetTitle('Eta')
        photon_eta.GetYaxis().SetTitle('Number of Events')

        invariant_mass=ROOT.TH1F('invariantmass', 'Invariant Mass Of 2 Photons', 50, 20, 40)
        invariant_mass.GetXaxis().SetTitle('Mass (GeV)')
        invariant_mass.GetYaxis().SetTitle('Number of Events')

        invariant_massg=ROOT.TH1F('invariantmassg', 'Invariant Mass Of 2 Photons (Gen)', 50, 29.99, 30.01)
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

        photon_ptg=ROOT.TH1F('photonptg', 'Photon Pt (Gen)', 50, 0, 75.)
        photon_ptg.GetXaxis().SetTitle('Pt')
        photon_ptg.GetYaxis().SetTitle('Number of Events')

        l1photon_phi=ROOT.TH1F('l1photonphi', 'Phi', 50, -3., 3.)
        l1photon_phi.GetXaxis().SetTitle('Phi (L1)')
        l1photon_phi.GetYaxis().SetTitle('Number of Events')

        l1photon_pt=ROOT.TH1F('l1photonpt', 'Pt', 50, 0., 75.)
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

        numPhoton_L1=ROOT.TH1F('l1numphotons', 'Number of Photons (L1)', 15, 0, 15)
        numPhoton_L1.GetXaxis().SetTitle('Number Of Photons')
        numPhoton_L1.GetYaxis().SetTitle('Number of Events')

        matchingPhotons_Pt=ROOT.TH1F('matchingPhotonsPt','Pt of Matching Photons (L1 vs Gen)', 50, 0, 75.)
        matchingPhotons_Pt.GetXaxis().SetTitle('Pt')
        matchingPhotons_Pt.GetYaxis().SetTitle('Number of Events')

        matchingPhotons_Energy=ROOT.TH1F('matchingPhotonsEnergy','Energy of Matching Photons (L1 vs Gen)', 50, 0, 50)
        matchingPhotons_Energy.GetXaxis().SetTitle('Energy')
        matchingPhotons_Energy.GetYaxis().SetTitle('Number of Events')

        num_matchingPhotons=ROOT.TH1F('nummatchingPhotons', 'Number of Matching Photons (L1 vs Gen)', 10, 0, 10)
        num_matchingPhotons.GetXaxis().SetTitle('Number of Matching Photons')
        num_matchingPhotons.GetYaxis().SetTitle('Number Of Events')

        recoMatchingPhotons_Pt=ROOT.TH1F('recoMatchingPhotons_Pt','Pt of Matching Photons (Reco vs Gen)', 50, 0, 75.)
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

        etaSpecialReco=ROOT.TH1F('etaSpecialReco', 'Difference in Eta Divided By Gen', 50, -2, 2)
        etaSpecialReco.GetXaxis().SetTitle('Eta Resolution')
        etaSpecialReco.GetYaxis().SetTitle('Number Of Events')

        phiSpecialReco=ROOT.TH1F('phiSpecialReco', 'Difference in Phi Divided By Gen', 50, -2, 2)
        phiSpecialReco.GetXaxis().SetTitle('Phi Resolution')
        phiSpecialReco.GetYaxis().SetTitle('Number Of Events')

        energySpecialReco=ROOT.TH1F('energySpecialReco', 'Difference in Energy Divided By Gen', 50, -2, 3)
        energySpecialReco.GetXaxis().SetTitle('Energy Resolution')
        energySpecialReco.GetYaxis().SetTitle('Number Of Events')
        
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
        
        massSquared=ROOT.TH1F('masssquared', 'Single Photon Mass Squared', 50, -.001, 0.001)
        massSquared.GetXaxis().SetTitle('Mass (GeV)')
        massSquared.GetYaxis().SetTitle('Number of Events')

        failedPhoton_Pt=ROOT.TH1F('failedphotonpt', 'Photons that failed ID cuts', 50, 0, 50)
        failedPhoton_Pt.GetXaxis().SetTitle('Pt (GeV)')
        failedPhoton_Pt.GetYaxis().SetTitle('Number of Events')

        failedPhoton_Energy=ROOT.TH1F('failedphotonenergy', 'Photons that failed ID cuts', 50, 0, 50)
        failedPhoton_Energy.GetXaxis().SetTitle('Energy (GeV)')
        failedPhoton_Energy.GetYaxis().SetTitle('Number of Events')

        delta_reco_l1=ROOT.TH1F('deltarecol1', 'Delta R (Reco vs L1)', 20, 0, 10)
        delta_reco_l1.GetXaxis().SetTitle('Delta R')
        delta_reco_l1.GetYaxis().SetTitle('Number of Events')

        energySpecialRecoL1=ROOT.TH1F('energyspecialrecol1', 'Difference in Reco and L1 Energy Divided by L1 Energy', 50, -2, 3)
        energySpecialRecoL1.GetXaxis().SetTitle('Energy Resolution')
        energySpecialRecoL1.GetYaxis().SetTitle('Number of Events')

        min_L1_reco=ROOT.TH1F('minrecol1', 'Minimum Value of R (Reco vs L1)', 40, 0, 10)
        min_L1_reco.GetXaxis().SetTitle('Minimum Value of R')
        min_L1_reco.GetYaxis().SetTitle('Number of Events')

        l1recoMatchingPhotons_Energy=ROOT.TH1F('l1recomatchingphotonsenergy', 'Matching Photons Energy (L1 vs Reco)', 50, 0, 80)
        l1recoMatchingPhotons_Energy.GetXaxis().SetTitle('Energy (GeV)')
        l1recoMatchingPhotons_Energy.GetYaxis().SetTitle('Number of Events')

        l1recoMatchingPhotons_Pt=ROOT.TH1F('l1recomatchingphotonspt', 'Matching Photons Pt (L1 vs Reco)', 50, 0, 75)
        l1recoMatchingPhotons_Pt.GetXaxis().SetTitle('Pt (GeV)')
        l1recoMatchingPhotons_Pt.GetYaxis().SetTitle('Number of Events')

        l1recoNum_matchingPhotons=ROOT.TH1F('l1reconummatchingphotons', 'Number of Matching Photons (L1 vs Reco)', 6, 0, 6)
        l1recoNum_matchingPhotons.GetXaxis().SetTitle('Number of Matching Photons')
        l1recoNum_matchingPhotons.GetYaxis().SetTitle('Number of Events')

        invariant_mass_L1_Reco=ROOT.TH1F('invariantmassl1reco', 'Invariant Mass of two Photons (Matched L1 vs Reco)', 100, 20, 40)
        invariant_mass_L1_Reco.GetXaxis().SetTitle('Invariant Mass (GeV)')
        invariant_mass_L1_Reco.GetYaxis().SetTitle('Number of Events')

        recoGenMassResolution=ROOT.TH1F('recogenmassresolution', 'Mass Resolution (Reco vs Gen)', 30, -3, 3)
        recoGenMassResolution.GetXaxis().SetTitle('Mass Resolution')
        recoGenMassResolution.GetYaxis().SetTitle('Number of Events')

        photon20l1 = 0
        photon33l1 = 0
        totalphotons = 0
        hoverphotonl1 = 0
        
        for event in f.eventTree:

                #reading the branches of eventTree             
                nPhoton = event.nPhoton
                nParticles = event.nParticles
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
                failedPhotonPt = event.failedPhotonPt
                failedPhotonEnergy = event.failedPhotonEnergy
                failedPhotons = event.failedPhotons
                
                totalphotons += 1

                if nPhoton_L1 >= 1:

                        if l1photonEnergy[0] > 15: #and abs(l1photonEta[0]) < 2.5:

                                photon20l1 += 1

                        if l1photonEnergy[0] > 26: #and abs(l1photonEta[0]) < 2.5:

                                photon33l1 += 1

                        if l1photonEnergy[0] > 10: #and abs(l1photonEta[0]) < 2.5:

                                hoverphotonl1 += 1
                
                l1genMatchingPhotonsPt = []
                l1genMatchingPhotonsEnergy = []
                l1genMatchingPhotonsPx = []
                l1genMatchingPhotonsPy = []
                l1genMatchingPhotonsPz = []
                l1genMatchingPhotons = 0
                l1genIndex = []
                
                for i in range(numPhotons_gen):

                        for j in range(nPhoton_L1):

                                eta_diff = l1photonEta[j] - photonEtag[i]
                                phi_diff = abs(l1photonPhi[j] - photonPhig[i])
                                energy_diff = l1photonEnergy[j] - photonEnergyg[i]

                                if phi_diff > math.pi:
                                        phi_diff2 = 2*math.pi - phi_diff
                                        deltaR = math.sqrt(eta_diff**2 + phi_diff2**2)
                                        delta_r.Fill(deltaR)
                                else:
                                        deltaR = math.sqrt(eta_diff**2 + phi_diff**2)
                                        delta_r.Fill(deltaR)

                                if deltaR < .5:

                                        l1genMatchingPhotons += 1
                                        l1genIndex.append(j)
                                        eta_difference = eta_diff/photonEtag[i]
                                        etaSpecial.Fill(eta_difference)
                                        phi_difference = phi_diff/photonPhig[i]
                                        phiSpecial.Fill(phi_difference)
                                        energy_difference = energy_diff/photonEnergyg[i]
                                        energySpecial.Fill(energy_difference)
                                        
                if l1genMatchingPhotons==0:

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
                                       
                if l1genMatchingPhotons==2:

                        for j in l1genIndex:

                                l1genMatchingPhotonsPt.append(l1photonPt[j])
                                l1genMatchingPhotonsEnergy.append(l1photonEnergy[j])
                                l1genMatchingPhotonsPx.append(l1photonPx[j])
                                l1genMatchingPhotonsPy.append(l1photonPy[j])
                                l1genMatchingPhotonsPz.append(l1photonPz[j])

                for mEnergy in l1genMatchingPhotonsEnergy:

                        matchingPhotons_Energy.Fill(mEnergy)

                for mPt in l1genMatchingPhotonsPt:

                        matchingPhotons_Pt.Fill(mPt)
                        
                num_matchingPhotons.Fill(l1genMatchingPhotons)

                if l1genMatchingPhotons==2:

                        px = l1genMatchingPhotonsPx[0] + l1genMatchingPhotonsPx[1]
                        py = l1genMatchingPhotonsPy[0] + l1genMatchingPhotonsPy[1]
                        pz = l1genMatchingPhotonsPz[0] + l1genMatchingPhotonsPz[1]
                        energy = l1genMatchingPhotonsEnergy[0] + l1genMatchingPhotonsEnergy[1]
                        l1invariantMassSquared = energy**2 - px**2 - py**2 - pz**2

                        if l1invariantMassSquared > 0:
                                l1invariantMass = math.sqrt(l1invariantMassSquared)
                                l1invariant_mass.Fill(l1invariantMass)

                        #if l1invariantMass < 0:

                                #print(px, py, pz, energy)

                        #l1etadiff = abs(l1photonEta[0] - l1photonEta[1])
                        #l1photon_etadiff.Fill(l1etadiff)
                
                recoMatchingPhotonsPt = []
                recoMatchingPhotonsEnergy = []
                recoMatchingPhotons = 0
                recoIndex = []
                recoMatchingPhotonsPx = []
                recoMatchingPhotonsPy = []
                recoMatchingPhotonsPz = []
                
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

                if recoMatchingPhotons==0:

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
                                recoMatchingPhotonsPx.append(photonPx[k])
                                recoMatchingPhotonsPy.append(photonPy[k])
                                recoMatchingPhotonsPz.append(photonPz[k])

                for rmEnergy in recoMatchingPhotonsEnergy:

                        recoMatchingPhotons_Energy.Fill(rmEnergy)

                for rmPt in recoMatchingPhotonsPt:

                        recoMatchingPhotons_Pt.Fill(rmPt)

                recoNum_matchingPhotons.Fill(recoMatchingPhotons)
                
                if recoMatchingPhotons==2:

                        px = recoMatchingPhotonsPx[0] + recoMatchingPhotonsPx[1]
                        py = recoMatchingPhotonsPy[0] + recoMatchingPhotonsPy[1]
                        pz = recoMatchingPhotonsPz[0] + recoMatchingPhotonsPz[1]
                        energy = recoMatchingPhotonsEnergy[0] + recoMatchingPhotonsEnergy[1]
                        invariantMassSquared = energy**2 - px**2 - py**2 - pz**2

                        #tried to get rid of negative values, didn't work
                        #px = recoMatchingPhotonsPx[0]
                        #py = recoMatchingPhotonsPy[0]
                        #pz = recoMatchingPhotonsPz[0]
                        #pMag = math.sqrt(px**2 + py**2 +pz**2)
                        #px2 = recoMatchingPhotonsPx[1]
                        #py2 = recoMatchingPhotonsPy[1]
                        #pz2 = recoMatchingPhotonsPz[1]
                        #pMag2 = math.sqrt(px2**2 + py2**2 + pz2**2)
                        #dotproduct = px*px2 + py*py2 + pz*pz2
                        #energy = recoMatchingPhotonsEnergy[0]
                        #energy2 = recoMatchingPhotonsEnergy[1]
                        #invariantmass = math.sqrt(2*((energy*energy2) - dotproduct))
                        #invariant_mass.Fill(invariantmass)

                        if invariantMassSquared > 0:

                                invariantMass = math.sqrt(invariantMassSquared)
                                invariant_mass.Fill(invariantMass)
                        #etadifference = abs(photonEta[0] - photonEta[1])
                        #photon_etadiff.Fill(etadifference)

                for i in range(len(recoMatchingPhotonsPx)):

                        px = recoMatchingPhotonsPx[i]
                        py = recoMatchingPhotonsPy[i]
                        pz = recoMatchingPhotonsPz[i]
                        energy = recoMatchingPhotonsEnergy[i]
                        singleMass = energy**2 - px**2 - py**2 - pz**2
                        massSquared.Fill(singleMass)
                                
                photonEnergygen = [value for value in photonEnergyg if value != 0]
                photonPxgen = [value for value in photonPxg if value != 0]
                photonPygen = [value for value in photonPyg if value != 0]
                photonPzgen = [value for value in photonPzg if value != 0]

                pt_lower_10 = False

                for pt in photonPtg:

                        photon_ptg.Fill(pt)

                        if pt < 10:

                                pt_lower_10 = True

                if pt_lower_10: continue

                genMass = []
                recoMass = []
                l1Mass = []        
                        
                if numPhotons_gen==2:

                        pxg = photonPxg[0] + photonPxg[1]
                        pyg = photonPyg[0] + photonPyg[1]
                        pzg = photonPzg[0] + photonPzg[1]
                        energyg = photonEnergyg[0] + photonEnergyg[1]
                        invariantMassg = math.sqrt(energyg**2 - pxg**2 - pyg**2 - pzg**2)
                        invariant_massg.Fill(invariantMassg)
                        genMass.append(invariantMassg)
                        etadiff = abs(photonEtag[0] - photonEtag[1])
                        phidiff = abs(photonPhig[0] - photonPhig[1])
                        photon_etag.Fill(etadiff)
                        photon_phig.Fill(phidiff)

                if nPhoton==2:

                        pxr = photonPx[0] + photonPx[1]
                        pyr = photonPy[0] + photonPy[1]
                        pzr = photonPz[0] + photonPz[1]
                        energyr = photonEnergy[0] + photonEnergy[1]
                        invariantMassr = math.sqrt(energyr**2 - pxr**2 - pyr**2 - pzr**2)
                        recoMass.append(invariantMassr)

                if nPhoton_L1==2:

                        pxl1 = l1photonPx[0] + l1photonPx[1]
                        pyl1 = l1photonPy[0] + l1photonPy[1]
                        pzl1 = l1photonPz[0] + l1photonPz[1]
                        energyl1 = l1photonEnergy[0] + l1photonEnergy[1]
                        invariantMassl1 = math.sqrt(energyl1**2 - pxl1**2 - pyl1**2 - pzl1**2)
                        l1Mass.append(invariantMassl1)

                for i in range(len(recoMass)):

                        for j in range(len(genMass)):

                                resolution1 = recoMass[i] - genMass[j]
                                resolution2 = resolution1/genMass[j]
                                recoGenMassResolution.Fill(resolution2)

                for o in range(len(l1Mass)):

                        for k in range(len(genMass)):

                                resolution3 = l1Mass[o] - genMass[k]
                                resolution4 = resolution3/genMass[k]
                                l1GenMassResolution.Fill(resolution4)
                                
                num_photons.Fill(nPhoton)
                numPhoton_L1.Fill(nPhoton_L1)

                l1recoMatchingPhotonsPt = []
                l1recoMatchingPhotonsEnergy = []
                l1recoMatchingPhotons = 0
                l1recoIndex = []
                l1recoMatchingPhotonsPx = []
                l1recoMatchingPhotonsPy = []
                l1recoMatchingPhotonsPz = []

                for o in range(nPhoton_L1):

                        for k in range(nPhoton):

                                eta_diff = photonEta[k] - l1photonEta[o]
                                phi_diff = photonPhi[k] - l1photonPhi[o]
                                deltaR = math.sqrt(eta_diff**2 + phi_diff**2)
                                energy_diff = photonEnergy[k] - l1photonEnergy[o]
                                delta_reco_l1.Fill(deltaR)

                                if deltaR < .5:

                                        l1recoMatchingPhotons += 1
                                        l1recoIndex.append(k)
                                        energy_difference = energy_diff/l1photonEnergy[o]
                                        energySpecialRecoL1.Fill(energy_difference)

                if l1recoMatchingPhotons==0:

                        for h in range(nPhoton):

                                delta_Reco = []

                                for u in range(nPhoton_L1):

                                        if l1photonPt[u] < 2: continue

                                        eta_diff = photonEta[h] - l1photonEta[u]
                                        phi_diff = photonPhi[h] - l1photonPhi[u]
                                        deltaReco = math.sqrt(eta_diff**2 + phi_diff**2)
                                        delta_Reco.append(deltaR)

                                if len(delta_Reco) != 0:

                                        r = min(delta_Reco)

                                        min_L1_reco.Fill(r)

                if l1recoMatchingPhotons==2:

                        for k in l1recoIndex:

                                l1recoMatchingPhotonsPt.append(photonPt[k])
                                l1recoMatchingPhotonsEnergy.append(photonEnergy[k])
                                l1recoMatchingPhotonsPx.append(photonPx[k])
                                l1recoMatchingPhotonsPy.append(photonPy[k])
                                l1recoMatchingPhotonsPz.append(photonPz[k])

                for rmEnergy in l1recoMatchingPhotonsEnergy:

                        l1recoMatchingPhotons_Energy.Fill(rmEnergy)

                for rmPt in l1recoMatchingPhotonsPt:

                        l1recoMatchingPhotons_Pt.Fill(rmPt)

                l1recoNum_matchingPhotons.Fill(l1recoMatchingPhotons)
                
                if l1recoMatchingPhotons==2:

                        px = l1recoMatchingPhotonsPx[0] + l1recoMatchingPhotonsPx[1]
                        py = l1recoMatchingPhotonsPy[0] + l1recoMatchingPhotonsPy[1]
                        pz = l1recoMatchingPhotonsPz[0] + l1recoMatchingPhotonsPz[1]
                        energy = l1recoMatchingPhotonsEnergy[0] + l1recoMatchingPhotonsEnergy[1]
                        invariantMass_L1_RecoS = energy**2 - px**2 - py**2 - pz**2

                        if invariantMass_L1_RecoS > 0:

                                invariantMass_L1_Reco = math.sqrt(energy**2 - px**2 - py**2 - pz**2)
                                invariant_mass_L1_Reco.Fill(invariantMass_L1_Reco)
                                
                if nPhoton != 0:

                        if hltPhoton20 == 1:

                                leadingPhoton_Pt = photonPt[0]
                                leadingPhoton_pt.Fill(leadingPhoton_Pt)

                        if nPhoton > 1:

                                trailingPhoton_Pt = photonPt[1]
                                trailingPhoton_pt.Fill(trailingPhoton_Pt)

                for Pt in photonPt:

                        photon_PtTotal.Fill(Pt)
                        
                for Phi in photonPhi:

                        photon_phi.Fill(Phi)

                for Eta in photonEta:

                        photon_eta.Fill(Eta)
                        
                for PhiL1 in l1photonPhi:

                        l1photon_phi.Fill(PhiL1)

                for PtL1 in l1photonPt:

                        l1photon_pt.Fill(PtL1)

                for EtaL1 in l1photonEta:

                        l1photon_eta.Fill(EtaL1)
                        
                if failedPhotons != 0:

                        for i in range(failedPhotons):

                                leadingfailedPt = failedPhotonPt[i]
                                leadingfailedEnergy = failedPhotonEnergy[i]

                                failedPhoton_Pt.Fill(leadingfailedPt)
                                failedPhoton_Energy.Fill(leadingfailedEnergy)
                                
        num_photons.Draw()
        photon_phi.Draw()
        photon_PtTotal.Draw()
        leadingPhoton_pt.Draw()
        trailingPhoton_pt.Draw()
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
        failedPhoton_Pt.Draw()
        failedPhoton_Energy.Draw()
        f.Write()
        photon20percent = float(photon20l1)/totalphotons
        photon33percent = float(photon33l1)/totalphotons
        photonhoverpercent = float(hoverphotonl1)/totalphotons
        print('Total amount of events: %d' % totalphotons)
        print('Events that passed 20 L1 trigger: %d' % photon20l1)
        print('Ratio for L1 20: %f' % photon20percent)
        print('Events that passes 33 L1 trigger: %d' % photon33l1)
        print('Ratio for L1 33: %f' % photon33percent)
        print('Events that passed hover triggers: %d' % hoverphotonl1)
        print('Ratio for L1 hover: %f' % photonhoverpercent)

if __name__ == '__main__':

        inputFile = 'Photongraphs.root'
        readTree(inputFile)
        #Plots trigger graphs
        #trigger = 'hltPhoton20'
        #drawTriggerEff(inputFile, trigger)
        #trigger2 = 'hltPhoton33'
        #drawTriggerEff(inputFile, trigger2)
        #trigger3 = 'hltP20H'
        #drawTriggerEff(inputFile, trigger3)
        #trigger4 = 'hltP30H'
        #drawTriggerEff(inputFile, trigger4)
        #trigger5 = 'hltDP30'
        #drawTriggerEff(inputFile, trigger5)
