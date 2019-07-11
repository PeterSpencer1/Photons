import ROOT
import math
from array import array

def deltaPhi(phi1, phi2):

        PHI = math.fabs(phi1-phi2)
        if (PHI<=3.14159265):
                return PHI;
        else:
                return 2*3.14159265-PHI;

def deltaR(eta1, eta2, phi1, phi2):

        return math.sqrt((eta2-eta1)*(eta2-eta1)+deltaPhi(phi1,phi2)*deltaPhi(phi1,phi2))


def vectorSumMass(px1, py1, pz1, px2, py2, pz2):

        E1 = math.sqrt(px1*px1 + py1*py1 + pz1*pz1);
        E2 = math.sqrt(px2*px2 + py2*py2 + pz2*pz2);
        cosTheta = (px1*px2 + py1*py2 + pz1*pz2)/ (E1*E2);
        return math.sqrt(2*E1*E2*(1-cosTheta));

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

        photon_phig=ROOT.TH1F('photonphig', 'Photon Phi Difference (Gen)', 75, 0, 4.)
        photon_phig.GetXaxis().SetTitle('Phi')
        photon_phig.GetYaxis().SetTitle('Number of Events')

        photon_etag=ROOT.TH1F('photonetag', 'Photon Eta Difference (Gen)', 75, 0, 4.)
        photon_etag.GetXaxis().SetTitle('Eta')
        photon_etag.GetYaxis().SetTitle('Number of Events')

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

        recoMatchingPhotons_Energy=ROOT.TH1F('recoMatchingPhotonsEnergy','Energy of Matching Photons (Reco vs Gen)', 50, 0, 50)
        recoMatchingPhotons_Energy.GetXaxis().SetTitle('Energy')
        recoMatchingPhotons_Energy.GetYaxis().SetTitle('Number of Events')

        recoNum_matchingPhotons=ROOT.TH1F('recoNummatchingPhotons', 'Number of Matching Photons (Reco vs Gen)', 10, 0, 10)
        recoNum_matchingPhotons.GetXaxis().SetTitle('Number of Matching Photons')
        recoNum_matchingPhotons.GetYaxis().SetTitle('Number Of Events')

        min_r=ROOT.TH1F('minr', 'Lowest Value of Delta R (L1 vs Gen)', 40, 0, 10)
        min_r.GetXaxis().SetTitle('Delta R')
        min_r.GetYaxis().SetTitle('Number Of Events')

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

        recoGenMassResolution=ROOT.TH1F('recogenmassresolution', 'Mass Resolution (Reco vs Gen)', 100, -.5, .5)
        recoGenMassResolution.GetXaxis().SetTitle('Mass Resolution (Reco vs Gen)')
        recoGenMassResolution.GetYaxis().SetTitle('Number of Events')

        l1GenMassResolution=ROOT.TH1F('l1genmassresolution', 'Mass Resolution (L1 vs Gen)', 100, -.8, .8)
        l1GenMassResolution.GetXaxis().SetTitle('Mass Resolution (L1 vs Gen)')
        l1GenMassResolution.GetYaxis().SetTitle('Number of Events')

        l1recoMassResolution=ROOT.TH1F('l1recomassresolution', 'Mass Resolution (L1 vs Reco)', 50, -1, 1)
        l1recoMassResolution.GetXaxis().SetTitle('Mass Resolution (L1 vs Reco)')
        l1recoMassResolution.GetYaxis().SetTitle('Number of Events')

        gen_invariant_mass=ROOT.TH1F('gen_invariantmass', 'Invariant Mass Of 2 Gen Photons', 100, 9.99, 10.01)
        gen_invariant_mass.GetXaxis().SetTitle('Invariant Mass at Gen (GeV)')
        gen_invariant_mass.GetYaxis().SetTitle('Number of Events')

        l1_invariant_mass=ROOT.TH1F('l1_invariantmass', 'Invariant Mass Of 2 L1 Photons', 150, 0, 20)
        l1_invariant_mass.GetXaxis().SetTitle('Invariant Mass at L1 (GeV)')
        l1_invariant_mass.GetYaxis().SetTitle('Number of Events')

        reco_invariant_mass=ROOT.TH1F('reco_invariantmass', 'Invariant Mass of 2 Reco Photons', 150, 6, 14)
        reco_invariant_mass.GetXaxis().SetTitle('Invariant Mass at Reco (GeV)')
        reco_invariant_mass.GetYaxis().SetTitle('Number of Events')

        l1gendeltar=ROOT.TH1F('l1gendeltar', 'Delta R (L1 vs Gen)', 50, 0, 7)
        l1gendeltar.GetXaxis().SetTitle('Delta R')
        l1gendeltar.GetYaxis().SetTitle('Number of Events')

        recogendeltar=ROOT.TH1F('recogendeltar', 'Delta R (Reco vs Gen)', 50, 0, 7)
        recogendeltar.GetXaxis().SetTitle('Delta R')
        recogendeltar.GetYaxis().SetTitle('Number of Events')

        l1recodeltar=ROOT.TH1F('l1recodeltar', 'Delta R (L1 vs Reco)', 50, 0, 7)
        l1recodeltar.GetXaxis().SetTitle('Delta R')
        l1recodeltar.GetYaxis().SetTitle('Number of Events')

        l1genEnergyResolution=ROOT.TH1F('l1genenergyresolution', 'Energy Resolution (L1 vs Gen)', 150, -1, 1)
        l1genEnergyResolution.GetXaxis().SetTitle('Energy Resolution')
        l1genEnergyResolution.GetYaxis().SetTitle('Number of Events')

        recogenEnergyResolution=ROOT.TH1F('recogenenergyresolution', 'Energy Resolution (Reco vs Gen)', 150, -.4, .4)
        recogenEnergyResolution.GetXaxis().SetTitle('Energy Resolution')
        recogenEnergyResolution.GetYaxis().SetTitle('Number of Events')

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
                photonEnergy = event.photonEnergy
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
                l1photonEnergy = event.l1photonEnergy
                hltPhoton20 = event.hltPhoton20
                hltPhoton33 = event.hltPhoton33
                hltP20H = event.hltP20H
                hltP30H = event.hltP30H
                hltDP30 = event.hltDP30
                failedPhotonPt = event.failedPhotonPt
                failedPhotonEnergy = event.failedPhotonEnergy
                failedPhotons = event.failedPhotons
                hltInvariantMass = event.hltInvariantMass
                genInvariantMass = event.genInvariantMass
                recoInvariantMass = event.recoInvariantMass
                l1InvariantMass = event.l1InvariantMass

                #reco_invariant_mass.Fill(recoInvariantMass)
                #l1_invariant_mass.Fill(l1InvariantMass)
                #gen_invariant_mass.Fill(genInvariantMass)

                #if genInvariantMass != 0.:

                        #if recoInvariantMass != 0.:
                        
                                #recomingen = recoInvariantMass - genInvariantMass
                                #genrecoresolution = recomingen/genInvariantMass
                                #recoGenMassResolution.Fill(genrecoresolution)

                #if genInvariantMass != 0:

                        #if l1InvariantMass != 0:

                                #l1mingen = l1InvariantMass - genInvariantMass
                                #genl1resolution = l1mingen/genInvariantMass
                                #l1GenMassResolution.Fill(genl1resolution)

                #if l1InvariantMass != 0:

                        #if recoInvariantMass != 0:

                                #l1minreco = l1InvariantMass - recoInvariantMass
                                #recol1resolution = l1minreco/recoInvariantMass
                                #l1recoMassResolution.Fill(recol1resolution)

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
                l1genMatchingPhotons = 0
                l1Index = []
                genIndex = []

                for i in range(numPhotons_gen):

                        if event.particleStatus[i] == 21 and event.mothers[i] == 2212:
                                continue

                        for j in range(nPhoton_L1):

                                l1_deltaR = deltaR(event.photonEtag[i], event.l1photonEta[j], event.photonPhig[i], event.l1photonPhi[j])

                                l1gendeltar.Fill(l1_deltaR)

                                energydiff = l1photonEnergy[j] - photonEnergyg[i]

                                if l1_deltaR < .3:

                                        if i not in genIndex and j not in l1Index:

                                                l1genMatchingPhotons += 1
                                                genIndex.append(i)
                                                l1Index.append(j)
                                                energyresolution = energydiff/photonEnergyg[i]
                                                l1genEnergyResolution.Fill(energyresolution)

                if l1genMatchingPhotons==0:

                        for h in range(nPhoton_L1):

                                delta_R = []

                                for u in range(numPhotons_gen):

                                        if photonPtg[u] < 2: continue

                                        eta_difference = l1photonEta[h] - photonEtag[u]
                                        phi_diff = l1photonPhi[h] - photonPhig[u]
                                        deltaRa = math.sqrt(eta_difference**2 + phi_diff**2)
                                        delta_R.append(deltaRa)

                                if len(delta_R) != 0:

                                        r = min(delta_R)

                                        min_r.Fill(r)

                if l1genMatchingPhotons==2:
                        #make a diphoton mass from gen
                        gen_gg_mass = vectorSumMass(event.photonPxg[genIndex[0]], event.photonPyg[genIndex[0]], event.photonPzg[genIndex[0]], event.photonPxg[genIndex[1]], event.photonPyg[genIndex[1]], event.photonPzg[genIndex[1]])
                        #make a diphoton mass from l1
                        l1_gg_mass = vectorSumMass(event.l1photonPx[l1Index[0]], event.l1photonPy[l1Index[0]], event.l1photonPz[l1Index[0]], event.l1photonPx[l1Index[1]], event.l1photonPy[l1Index[1]], event.l1photonPz[l1Index[1]])

                        gen_invariant_mass.Fill(gen_gg_mass)
                        l1_invariant_mass.Fill(l1_gg_mass)
                        l1GenMassResolution.Fill((l1_gg_mass-gen_gg_mass)/gen_gg_mass)

                for mEnergy in l1genMatchingPhotonsEnergy:

                        matchingPhotons_Energy.Fill(mEnergy)

                for mPt in l1genMatchingPhotonsPt:

                        matchingPhotons_Pt.Fill(mPt)

                num_matchingPhotons.Fill(l1genMatchingPhotons)

                recoMatchingPhotonsPt = []
                recoMatchingPhotonsEnergy = []
                recogenMatchingPhotons = 0
                recoIndex = []
                recogenIndex = []

                for o in range(numPhotons_gen):

                        if event.particleStatus[o] == 21 and event.mothers[o] == 2212:

                                continue

                        for k in range(nPhoton):

                                reco_deltaR = deltaR(event.photonEtag[o], event.photonEta[k], event.photonPhig[o], event.photonPhi[k])

                                recogendeltar.Fill(reco_deltaR)

                                energyDiff = photonEnergy[k] - photonEnergyg[o]

                                if reco_deltaR < .3:

                                        if o not in recogenIndex and k not in recoIndex:

                                                recogenMatchingPhotons += 1
                                                recogenIndex.append(o)
                                                recoIndex.append(k)
                                                energyresolution = energyDiff/photonEnergyg[o]
                                                recogenEnergyResolution.Fill(energyresolution)

                if recogenMatchingPhotons==0:

                        for h in range(nPhoton):

                                delta_Reco = []

                                for u in range(numPhotons_gen):

                                        if photonPtg[u] < 2: continue

                                        eta_diff = photonEta[h] - photonEtag[u]
                                        phi_diff = photonPhi[h] - photonPhig[u]
                                        deltaReco = math.sqrt(eta_diff**2 + phi_diff**2)
                                        delta_Reco.append(deltaReco)

                                if len(delta_Reco) != 0:

                                        r = min(delta_Reco)

                                        min_reco.Fill(r)

                if recogenMatchingPhotons==2:

                        #make a diphoton mass from gen
                        gen_gg_mass = vectorSumMass(event.photonPxg[recogenIndex[0]], event.photonPyg[recogenIndex[0]], event.photonPzg[recogenIndex[0]], event.photonPxg[recogenIndex[1]], event.photonPyg[recogenIndex[1]], event.photonPzg[recogenIndex[1]])
                        #make a diphoton mass from l1
                        reco_gg_mass = vectorSumMass(event.photonPx[recoIndex[0]], event.photonPy[recoIndex[0]], event.photonPz[recoIndex[0]], event.photonPx[recoIndex[1]], event.photonPy[recoIndex[1]], event.photonPz[recoIndex[1]])
                        reco_invariant_mass.Fill(reco_gg_mass)
                        recoGenMassResolution.Fill((reco_gg_mass-gen_gg_mass)/gen_gg_mass)

                for rmEnergy in recoMatchingPhotonsEnergy:

                        recoMatchingPhotons_Energy.Fill(rmEnergy)

                for rmPt in recoMatchingPhotonsPt:

                        recoMatchingPhotons_Pt.Fill(rmPt)

                recoNum_matchingPhotons.Fill(recogenMatchingPhotons)
                
                photonEnergygen = [value for value in photonEnergyg if value != 0]

                #pt_lower_10 = False

                #for pt in photonPtg:

                        #photon_ptg.Fill(pt)

                        #if pt < 10:

                                #pt_lower_10 = True

                #if pt_lower_10: continue

                num_photons.Fill(nPhoton)
                numPhoton_L1.Fill(nPhoton_L1)

                l1recoMatchingPhotonsPt = []
                l1recoMatchingPhotonsEnergy = []
                l1recoMatchingPhotons = 0
                l1recoIndex = []
                recol1Index = []

                for o in range(nPhoton_L1):

                        for k in range(nPhoton):

                                l1_deltaR = deltaR(event.photonEta[k], event.l1photonEta[o], event.photonPhi[k], event.l1photonPhi[o])

                                l1recodeltar.Fill(l1_deltaR)

                                if l1_deltaR < .3:

                                        if k not in l1recoIndex and o not in recol1Index:
                                                l1recoMatchingPhotons += 1
                                                l1recoIndex.append(k)
                                                recol1Index.append(o)

                if l1recoMatchingPhotons==0:

                        for h in range(nPhoton):

                                delta_Reco = []

                                for u in range(nPhoton_L1):
                                
                                        if l1photonPt[u] < 2: continue

                                        eta_diff = photonEta[h] - l1photonEta[u]
                                        phi_diff = photonPhi[h] - l1photonPhi[u]
                                        deltaReco = math.sqrt(eta_diff**2 + phi_diff**2)
                                        delta_Reco.append(deltaReco)

                                if len(delta_Reco) != 0:

                                        r = min(delta_Reco)

                                        min_L1_reco.Fill(r)

                if l1recoMatchingPhotons==2:

                        reco_gg_mass = vectorSumMass(event.photonPx[l1recoIndex[0]], event.photonPy[l1recoIndex[0]], event.photonPz[l1recoIndex[0]], event.photonPx[l1recoIndex[1]], event.photonPy[l1recoIndex[1]], event.photonPz[l1recoIndex[1]])
                        l1_gg_mass = vectorSumMass(event.l1photonPx[recol1Index[0]], event.l1photonPy[recol1Index[0]], event.l1photonPz[recol1Index[0]], event.l1photonPx[recol1Index[1]], event.l1photonPy[recol1Index[1]], event.l1photonPz[recol1Index[1]])
                        l1recoMassResolution.Fill((l1_gg_mass-reco_gg_mass)/reco_gg_mass)

                for rmEnergy in l1recoMatchingPhotonsEnergy:

                        l1recoMatchingPhotons_Energy.Fill(rmEnergy)

                for rmPt in l1recoMatchingPhotonsPt:

                        l1recoMatchingPhotons_Pt.Fill(rmPt)

                l1recoNum_matchingPhotons.Fill(l1recoMatchingPhotons)

                if nPhoton != 0:

                        if nPhoton > 0:

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

                for Phil1 in l1photonPhi:

                        l1photon_phi.Fill(Phil1)

                for Ptl1 in l1photonPt:

                        l1photon_pt.Fill(Ptl1)

                for Etal1 in l1photonEta:

                        l1photon_eta.Fill(Etal1)

                if failedPhotons != 0:

                        for i in range(failedPhotons):

                                leadingfailedPt = failedPhotonPt[0]
                                leadingfailedEnergy = failedPhotonEnergy[0]

                                failedPhoton_Pt.Fill(leadingfailedPt)
                                failedPhoton_Energy.Fill(leadingfailedEnergy)

        num_photons.Draw()
        photon_phi.Draw()
        photon_PtTotal.Draw()
        leadingPhoton_pt.Draw()
        trailingPhoton_pt.Draw()
        photon_eta.Draw()
        photon_phig.Draw()
        photon_etag.Draw()
        l1photon_phi.Draw()
        l1photon_pt.Draw()
        l1photon_eta.Draw()
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
        #photon20percent = float(photon20l1)/totalphotons
        #photon33percent = float(photon33l1)/totalphotons
        #photonhoverpercent = float(hoverphotonl1)/totalphotons
        #print('Total amount of events: %d' % totalphotons)
        #print('Events that passed 20 L1 trigger: %d' % photon20l1)
        #print('Ratio for L1 20: %f' % photon20percent)
        #print('Events that passes 33 L1 trigger: %d' % photon33l1)
        #print('Ratio for L1 33: %f' % photon33percent)
        #print('Events that passed hover triggers: %d' % hoverphotonl1)
        #print('Ratio for L1 hover: %f' % photonhoverpercent)

if __name__ == '__main__':

        inputFile = 'Photongraphs.root'
        readTree(inputFile)
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


