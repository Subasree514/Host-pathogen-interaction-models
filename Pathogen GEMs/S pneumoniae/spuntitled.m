%%
load('SPD39_host.mat')
load('RSmodel_Recon3D1_2023_ecoli.mat')
modelrs=RSmodel_Recon3D1_2023_ecoli;
modelload=SPD39;
%%
%ids=findMetIDs(modelrs,{'glc_D[c]';'gln_L[c]';'trp_L[c]';'tyr_L[c]';'rib_D[c]';'cys_L[c]';'val_L[c]';'tyr_L[c]';'trp_L[c]';'ser_L[c]';'rib_D[c]';'glu_L-1[c]';'hcys_L[c]'});
%modelrs.mets(ids)={'glc-D[c]';'gln-L[c]';'trp-L[c]';'tyr-L[c]';'rib-D[c]';'cys-L[c]';'val-L[c]';'tyr-L[c]';'trp-L[c]';'ser-L[c]';'rib-D[c]';'glu-L[c]';'hcys-L[c]'};
%%
modelload = addReaction(modelload,'ura_tr','reactionFormula','ura[e] <=> ura[c]');
modelload = addReaction(modelload,'sucr_tr','reactionFormula','sucr[e] <=> sucr[c]');
%modelload = addReaction(modelload,'no2_tr','reactionFormula','no2[e] <=> no2[c]');
modelload = addReaction(modelload,'lac-D_tr','reactionFormula','lac_D[e] <=> lac_D[c]');
%modelload = addReaction(modelload,'o2s_tr','reactionFormula','o2s[e] <=> o2s[c]');
%modelload=addExchangeRxn(modelload,{'o2s[e]'});
%%
modelrs = removeRxns(modelrs,{'RS_43';'RS_115'});
modelload=addExchangeRxn(modelload,{'o2s[e]'});
rsexchanges={'EXh2s';'EXh2o2';'EX_o2s[e]'};
rsids=findRxnIDs(modelload,rsexchanges);
modelload.lb(rsids)=-1000;
%%
intersect(modelload.rxns,modelrs.rxns);
modelrs = removeRxns(modelrs,ans);
[modelload, removedRxnInd1, keptRxnInd1] = checkDuplicateRxn(modelload,'FR');
[model_1_rs] = mergeTwoModels(modelload,modelrs);
optimizeCbModel(model_1_rs)
model_final_rs=model_1_rs;
optimizeCbModel(SPD39)
[model_1_rs1, removedRxnInd1, keptRxnInd1] = checkDuplicateRxn(model_1_rs,'S');
[model_1_rs2, removedRxnInd2, keptRxnInd2] = checkDuplicateRxn(model_1_rs,'FR');
SPD39_rs_host=model_1_rs2;
SPD39_rs_host.lb(806)=-1000;
optimizeCbModel(SPD39)
optimizeCbModel(SPD39_rs_host)
save('SPD39_rs_host','SPD39_rs_host');
