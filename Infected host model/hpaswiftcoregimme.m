%% load the transcriptomic data of the infected tissues
run trancriptomics_data_integration.m;
%% core reactions important for lungs and gastric cells
run organs_core_reactions.m
%% pneumonia GEM reconstruction
%% map the expression data only to the genes present in the lung cell line models
mpi1=SPlung;
model1=modelA549final;
hparemove=setdiff(splunggenes,model1.genes);
hpaonly=matches(splunggenes,hparemove);
hpainclude=find(hpaonly==1);
hpaincludeids=setdiff(1:length(splunggenes),hpainclude')';
%% repeat the same process for other infection points
lungdata=mpi1(:,1);
lowerThs = prctile(lungdata(hpaincludeids,:),30,'all');
%upperThs=prctile(hpaconsensusdata(hpaincludeids,:),70,'all');
%%
expressionData1.gene=splunggenes(hpaincludeids,:);
expressionData1.value=lungdata(hpaincludeids,:);
expressionRxns1=mapExpressionToReactions(modelA549final,expressionData1,'minSum');
%% important core reactions
rxns_add={'DATPtn'
'DGTPtn'
'DTTPtn'
'NDPK5n'
'NDPK8n'
'DM_atp_c_'
'biomass_reaction'};
corem1=findRxnIDs(modelA549final,rxns_add);
%threshold_lb1=prctile(mpi1(:,3),30);
expressionRxns1(corem1)=5;
tissueModel_sp_0 = GIMME(modelA549final, expressionRxns1,lowerThs);
model3=tissueModel_sp_0;
minimalmedianutrient.m
tissueModel_sp_0=model3;
save('tissueModel_sp_0','tissueModel_sp_0');
%% gastric cancer GEM reconstruction
%% map the expression data only to the genes present in the cell line models
mpi2=HPcolon;
model2=modelHK74final;
hparemove=setdiff(hpcolongenes,model2.genes);
hpaonly=matches(hpcolongenes,hparemove);
hpainclude=find(hpaonly==1);
hpaincludeids=setdiff(1:length(hpcolongenes),hpainclude')';
%% repeat the same process for other infection points
colondata=mpi2(:,1);
lowerThs2 = prctile(colondata(hpaincludeids,:),15,'all');
%upperThs=prctile(hpaconsensusdata(hpaincludeids,:),70,'all');
expressionData2.gene=hpcolongenes(hpaincludeids,:);
expressionData2.value=colondata(hpaincludeids,:);
expressionRxns2=mapExpressionToReactions(modelHK74final,expressionData2,'minSum');
corem1=findRxnIDs(modelHK74final,rxns_add);
%threshold_lb=prctile(mpi2(:,3),20);
expressionRxns2(corem1)=5;
tissueModel_hp_0 = GIMME(modelHK74final, expressionRxns2, lowerThs2);
model3=tissueModel_hp_0;
minimalmedianutrient.m
tissueModel_hp_0=model3;
save('tissueModel_hp_0','tissueModel_hp_0');