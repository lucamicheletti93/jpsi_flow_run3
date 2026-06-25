void read_epos4() {

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