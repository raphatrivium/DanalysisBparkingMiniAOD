//***********************************************************************************
// MAIN FUNCTION
// ------------ method called for each event  ------------
void DstarD0TTree::analyze{	
	//Calling Secondary Functions
	RecDstar(iEvent,iSetup,RecVtx); //Reconstruction of D*
	GenDstarInfo(iEvent,iSetup); //Stores information from D0 and its products.
	GenDstarMesonBInfo(iEvent,iSetup);	
	//Fill the tree
	data->Fill();
}//End analyze
//***********************************************************************************
void DstarD0TTree::RecDstar{
   			//Cuts
				//Fill vectors
				Dsmass.push_back(dS_p4.M());				
}//End RecDstar
//***********************************************************************************
void DstarD0TTree::GenDstarInfo{
			//Fill vector
			MCDsmass.push_back(p.mass());			
}//End GenDstarInfo
//***********************************************************************************
void DstarD0TTree::GenDstarMesonBInfo{
			// D* fromm B0 or B+-		
			FlagDstarfromB.push_back(1);			
}//End GenDstarMesonBInfo
//***********************************************************************************
void DstarD0TTree::initialize(){
	//Clear Vectors
}
//***********************************************************************************
void DstarD0TTree::beginJob(){
		
	data->Branch("Dsmass",&Dsmass);	
	data->Branch("MCDsmass",&MCDsmass);
}
