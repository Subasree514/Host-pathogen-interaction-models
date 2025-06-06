load('model_final.mat');
load('tissueModel_hp_0.mat');
coreRxnInd = findRxnIDs(model_final,tissueModel_hp_0.rxns);
find(coreRxnInd==0)
tissueModel_hp_0.rxns(ans)
coreRxnInd = findRxnIDs(model_final,setdiff(tissueModel_hp_0.rxns,ans));
tissueModel = fastcore(model_final, coreRxnInd)
optimizeCbModel(tissueModel)
model3=tissueModel;
memhuman
tissueModel_hp_0=model3;