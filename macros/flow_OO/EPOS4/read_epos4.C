void read_epos4() {
    TString dirName = "simOO";
    TSystemDirectory dir("dir", dirName);
    TList *fInList = dir.GetListOfFiles();

    if (!fInList) return;

    TIter next(fInList);
    TSystemFile *file;
    int counterJpsi = 0;

    TH1F *hMass = new TH1F("hMass", ";#it{M} (GeV/#it{c}^{2});Counts", 200, 2, 4);
    TH1F *hPt = new TH1F("hPt", ";#it{p}_{T} (GeV/#it{c});Counts", 20, 0, 10);
    TH1F *hEta = new TH1F("hEta", ";#it{#eta};Counts", 20, -5, 5);

    while ((file = (TSystemFile*)next())) {
        TString fInName = file->GetName();
        if (file->IsDirectory()) continue;
        if (!fInName.EndsWith(".root")) continue;

        TString fullName = dirName + "/" + fInName;
        std::cout << fullName << std::endl;
        TFile *fIn = TFile::Open(fullName);
        if (!fIn || fIn->IsZombie()) continue;

        TTree *treeEposHead = (TTree*) fIn->Get("teposhead");
        TTree *treeEposEvent = (TTree*) fIn->Get("teposevent");
        if (!treeEposHead || !treeEposEvent) continue;

        Long64_t nEntries = treeEposEvent->GetEntries();
        int nEvents, nParticles;
        int particleId[99999];
        float fMass[99999], fPx[99999], fPy[99999], fPz[99999];

        treeEposEvent->SetBranchAddress("nev", &nEvents);
        treeEposEvent->SetBranchAddress("np", &nParticles);
        treeEposEvent->SetBranchAddress("id", particleId);
        treeEposEvent->SetBranchAddress("e", fMass);
        treeEposEvent->SetBranchAddress("px", fPx);
        treeEposEvent->SetBranchAddress("py", fPy);
        treeEposEvent->SetBranchAddress("pz", fPz);

        for (int iEntry = 0;iEntry < nEntries;iEntry++) {
            treeEposEvent->GetEntry(iEntry);
            for (int iParticle = 0;iParticle < nParticles;iParticle++) {
                if (particleId[iParticle] == 441) {
                    float mass = fMass[iParticle];
                    float px = fPx[iParticle];
                    float py = fPy[iParticle];
                    float pz = fPz[iParticle];
                    float pt = TMath::Sqrt(px*px + py*py);
                    float ptot = TMath::Sqrt(px*px + py*py + pz*pz);
                    float eta = 0.5*TMath::Log((ptot + pz) / (ptot - pz));
                    hMass->Fill(mass);
                    hPt->Fill(pt);
                    hEta->Fill(eta);
                    counterJpsi++;
                }
            }
        }
        fIn->Close();
        delete fIn;
    }

    std::cout << "N. Jpsi = " << counterJpsi << std::endl;

    TCanvas *canvas = new TCanvas("canvas", "", 1800, 600);
    canvas->Divide(3, 1);
    canvas->cd(1); hMass->Draw("H");
    canvas->cd(2); hPt->Draw("H");
    canvas->cd(3); hEta->Draw("H");
}
/////////////////////////////
void read_epos4_from_AO2D() {

    int nFiles = 6;

    int counterJpsi = 0;
    for (int iFile = 0;iFile < nFiles;iFile++) {
        TFile *fIn = TFile::Open(Form("AO2D_%d.root", iFile));
        fIn->ls();
        TTree *treeMcCollision = (TTree*) fIn -> Get("DF_1/O2mccollision_001");
        std::cout << "n. Events: " << treeMcCollision->GetEntries() << std::endl;
        //treeMcCollision->Print();

        float posX, posY, posZ, impactParameter, eventPlaneAngle;

        treeMcCollision->SetBranchAddress("fPosX", &posX);
        treeMcCollision->SetBranchAddress("fPosY", &posY);
        treeMcCollision->SetBranchAddress("fPosZ", &posZ);
        treeMcCollision->SetBranchAddress("fImpactParameter", &impactParameter);
        treeMcCollision->SetBranchAddress("fEventPlaneAngle", &eventPlaneAngle);

        TTree *treeMcParticle = (TTree*) fIn -> Get("DF_1/O2mcparticle_001");
        std::cout << "n. Particles: " << treeMcParticle->GetEntries() << std::endl;
        //treeMcParticle->Print();
        
        int indexMcCollisions, pdgCode, statusCode;

        treeMcParticle->SetBranchAddress("fIndexMcCollisions", &indexMcCollisions);
        treeMcParticle->SetBranchAddress("fPdgCode", &pdgCode);
        treeMcParticle->SetBranchAddress("fStatusCode", &statusCode);

        Long64_t nColls = treeMcCollision->GetEntries();
        for (int iColl = 0;iColl < nColls;iColl++) {
            treeMcCollision->GetEntry(iColl);
            //Printf("[%d],[%f, %f, %f], %f, %f", iColl, posX, posY, posZ, impactParameter, eventPlaneAngle);
        }

        Long64_t nParts = treeMcParticle->GetEntries();
        for (int iPart = 0;iPart < nParts;iPart++) {
            treeMcParticle->GetEntry(iPart);
            if (pdgCode == 421) { // 421 D0
                counterJpsi++;
                std::cout << "Jpsi in the tree!" << std::endl;
                std::cout << "Status code: " << statusCode << " -> " << indexMcCollisions << std::endl;
                treeMcCollision->GetEntry(indexMcCollisions);
                //Printf("[%f, %f, %f], %f, %f", posX, posY, posZ, impactParameter, eventPlaneAngle);
            }
        }
    }
    std::cout << "N. J/psi: " << counterJpsi << std::endl;
    
}



