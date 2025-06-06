%load('bacteria_RS.mat')
%find(modelCW2.c==1);
%modelCW2.ub(ans)=0.033;
%save('modelCW2','modelCW2');
%modelCOV318 = addSinkReactions(modelCOV318,biomassmets([11;15;16;18;41;32;33;35;36;37;34;38]));
%find(modelcervixhpa.c==1);
%modelfallopianhpa.ub(ans)=0.02;
%save('modelcervixhpa.mat','modelcervixhpa')
%load('model_RSModule.mat')
%load('RSmodel.mat')
%model_final=model;
%RSmodelnew = removeRxns(RSmodel,{%'DM_2_mthio_et-adduct[c]'
%'DM_co3-2[c]'
%'DM_s2-[c]'
%'EX_ura-adduct[e]'
%'RS_11'
%'RS_158'});
%'RS_159'
%'RS_160'
%'RS_164'});
%'RS_185'
%'RS_193'});
%'sink_6-Benzylaminopurine-copper(II)[c]'
%'sink_Alanylhistidine-cu[c]'
%'sink_Cimetidine-copper(II)[c]'
%'sink_co3-[c]'
%'sink_o-2[c]'});
%RSmodel = removeRxns(RSmodel,RSmodel.rxns(151:200));
%load('RSmodelchecknew.mat')
%load('RSmodel_Recon3D1_2023.mat')
modelload=iIT341;
modelrs=RSmodel_Recon3D1_2023_ecoli;
%%
ids=findMetIDs(modelrs,'ser_L[c]');
modelrs.mets(ids)={'ser__L[c]'}
%% glutamate
ids=findMetIDs(modelrs,'glu_L-1[c]');
modelrs.mets(ids)={'glu__L[c]'};
%% glutamine
ids=findMetIDs(modelrs,'gln_L[c]');
modelrs.mets(ids)={'gln__L[c]'};
%% glucose
ids=findMetIDs(modelrs,'glc_D[c]');
modelrs.mets(ids)={'glc__D[c]'};
%%
ids=findMetIDs(modelrs,'rib_D[c]');
modelrs.mets(ids)={'rib__D[c]'};
%%
ids=findMetIDs(modelrs,'trp_L[c]');
modelrs.mets(ids)={'trp__L[c]'};
%%
ids=findMetIDs(modelrs,'tyr_L[c]');
modelrs.mets(ids)={'tyr__L[c]'};
%%
ids=findMetIDs(modelrs,'val_L[c]');
modelrs.mets(ids)={'val__L[c]'};
%%
modelrs = addReaction(modelrs,'ura_tr','reactionFormula','ura[e] <=> ura[c]');
%modelrs = addReaction(modelrs,'h2o2_tr','reactionFormula','h2o2[e] <=> h2o2[c]');
%modelrs = addReaction(modelrs,'h2s_tr','reactionFormula','h2s[e] <=> h2s[c]');
%modelrs = addReaction(modelrs,'o2s_tr','reactionFormula','o2s[e] <=> o2s[c]');

%%
intersect(modelrs.rxns,modelload.rxns);
modelload = removeRxns(modelload,ans);
[modelload, removedRxnInd1, keptRxnInd1] = checkDuplicateRxn(modelload,'FR');
[model_1_rs] = mergeTwoModels(modelload,modelrs);
%run dmemfinal.m
%[model_1_rs, removedRxnInd, keptRxnInd] = checkDuplicateRxn(model_1_rs,'FR')
optimizeCbModel(model_1_rs)
model_final_rs=model_1_rs;
%rxnids_exes3=findRxnIDs(model_final_rs,{'EX_o2[e]';'EX_h2o[e]';'EX_h[e]';'EX_co3[e]'});
%model_final_rs.lb(rxnids_exes3)=[-1000;-1.207674955;-1.143436078;0];
%rxnids_exes3=findRxnIDs(model_final_rs,{'EX_h2o_e';'EX_h2_e'});
%model_final_rs.lb(rxnids_exes3)=[-1.207674955;-1.143436078];
%save('modelsinglecellovaryharvetta_rsnew','modelsinglecellovaryharvetta_rs');
%%
%rxnids_exes3=findRxnIDs(model_final_rs,{'EX_o2(e)';'EX_h2o(e)';'EX_h(e)';'EX_o2s(e)';'EX_h2o2(e)'});
%model_final_rs.lb(rxnids_exes3)=[-1;-1;-1;-1;-1];
%%
[model_1_rs1, removedRxnInd1, keptRxnInd1] = checkDuplicateRxn(model_1_rs,'S');
[model_1_rs2, removedRxnInd2, keptRxnInd2] = checkDuplicateRxn(model_1_rs,'FR');
iIT341_rs=model_1_rs2;
optimizeCbModel(iIT341_rs)
%%
%%
iIT341_rs = removeRxns(iIT341_rs,{'RS_60';'RS_43';'RS_115'});
%modelrs=addExchangeRxn(modelrs,{'o2s[e]';'h2s[e]'});
%modelrs.lb(655)=-1000;
rsexchanges={'EX_o2s[e]';'EX_h2s[e]';'EX_h2o2[e]'};
rsids=findRxnIDs(iIT341_rs,rsexchanges);
iIT341_rs.lb(rsids)=-1000;
iIT341_rs.ub(rsids)=1000;
%removedRxnInd3rxns=setdiff(RSmodel.rxns,model_1_rs.rxns);
%RSmodel_Recon3D = removeRxns(RSmodel,[removedRxnInd3rxns;model_1_rs.rxns(removedRxnInd1);model_1_rs.rxns(removedRxnInd2)]);
%%
%load('ihuman_modified_RSmodel.mat')
%[removedRxnInd3rxns;model_1_rs.rxns(removedRxnInd1);model_1_rs.rxns(removedRxnInd2)]
%RSmodel_ihuman = removeRxns(RSmodel,ans);
%save('RSmodel_ihuman','RSmodel_ihuman');
%%

