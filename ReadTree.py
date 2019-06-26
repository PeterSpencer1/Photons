import ROOT

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

        event_count_before = 0
        event_count_after = 0

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

                if nPhoton==2:

                        px = photonPx[0] + photonPx[1]
                        py = photonPy[0] + photonPy[1]
                        pz = photonPz[0] + photonPz[1]
                        energy = photonEnergy[0] + photonEnergy[1]
                        invariantMass = energy**2 - px**2 - py**2 - pz**2
                        invariant_mass.Fill(invariantMass)

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
        invariant_mass.Draw()
        f.Write()

if __name__ == '__main__':

        inputFile = 'Photongraphs.root'
        readTree(inputFile)
