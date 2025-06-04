% RS integrated host model
load('model_hp6_rs.mat')
model1=model_hp6_rs;
model1 = addExchangeRxn(model1,{'no[e]'
'lac_L[e]'
'adn[e]'
'k[e]'
'Lkynr[e]'
'for[e]'
'succ[e]'
'gthrd[e]'
'urate[e]'
'etoh[e]'
'ac[e]'
'o2s[e]'
'caro[e]'
'ascb_L[e]'
'avite1[e]'});
mem_tr={'no[e] <=> no[c]'
'lac_L[e] <=> lac_L[c]'
'adn[e] <=> adn[c]'
'k[e] <=> k[c]'
'Lkynr[e] <=> Lkynr[c]'
'for[e] <=> for[c]'
'succ[e] <=> succ[c]'
'gthrd[e] <=> gthrd[c]'
'urate[e] <=> urate[c]'
'etoh[e] <=> etoh[c]'
'ac[e] <=> ac[c]'
'o2s[e] <=> o2s[c]'
'caro[e] <=> caro[c]'
'ascb_L[e] <=> ascb_L[c]'
'avite1[e] <=> avite1[c]'};
rxnNamesnew={'U1';'U2';'U3';'U4';'U5';'U6';'U7';'U8';'U9';'U10';'U11';'U12';'U13';'U14';'U15';'U16';'U17';'U18';'U19';'U20';'U21';'U22';'U23';'U24';'U25';'U26';'U27';'U28';'U29';'U30';'U31';'U32';'U33';'U34';'U35';'U36';'U37';'U38';'U39';'U40';'U41';'U42';'U43';'U44';'U45';'U46';'U47';'U48';'U49';'U50'};
memsubsytem={'Transport, extracellular'};
for i = 1:15
model1 = addReaction(model1,rxnNamesnew{i,1},'reactionFormula',mem_tr{i,1},'subSystem',memsubsytem{1,1});
end
model2=model1;
[model1, removedRxnInd1, keptRxnInd1] = checkDuplicateRxn(model1,'S');
[model1, removedRxnInd2, keptRxnInd2] = checkDuplicateRxn(model1,'FR');
blockexchanges={'EX_no[e]'
'EX_lac_L[e]'
'EX_adn[e]'
'EX_k[e]'
'EX_Lkynr[e]'
'EX_for[e]'
'EX_succ[e]'
'EX_gthrd[e]'
'EX_urate[e]'
'EX_etoh[e]'
'EX_ac[e]'
'EX_o2s[e]'
'EX_caro[e]'
'EX_ascb_L[e]'
'EX_avite1[e]'};
findRxnIDs(model1,blockexchanges);
model1.lb(ans)=0;
load('iIT341_rsmin.mat')
%%
rxnlist=contains(model1.rxns,'EX_');
find(rxnlist==1);
model1.lb(ans)=-1000;
modelHost=model1;
%%
%%
models={};
nameTagsModels={};
models{1,1}=iIT341_rsmin;
bioID={};
bioID{1,1}=iIT341_rsmin.rxns(find(strncmp(iIT341_rsmin.rxns,'BIOMASS_HP_published',7)));
nameTagsModels{1,1}=strcat('HP_');
nameTagHost='Gastric_';
%reducedDietConstraints={};
[modelJoint] = createMultipleSpeciesModel(models,'nameTagsModels',nameTagsModels,'modelHost',modelHost,'nameTagHost',nameTagHost);
rxnlist=contains(modelJoint.rxns,'HP_IEX_');
find(rxnlist==1);
modelJoint.lb(ans)=-1000;
%%
%modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('HP_IEX_',modelJoint.rxns)),0,'l');
%modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('HP_IEX_',modelJoint.rxns)),1000,'u');
%%
modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('Gastric_EX_',modelJoint.rxns)),0,'l');
modelJoint=changeRxnBounds(modelJoint,'Gastric_EX_o2[e]b',-2,'l');
modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('Gastric_IEX',modelJoint.rxns)),1000,'u');
%%
sp_IEX_exchanges={'HP_IEX_ala_D[u]tr'
'HP_IEX_ala_L[u]tr'
'HP_IEX_arg_L[u]tr'
'HP_IEX_his_L[u]tr'
'HP_IEX_ile_L[u]tr'
'HP_IEX_leu_L[u]tr'
'HP_IEX_met_L[u]tr'
'HP_IEX_val_L[u]tr'
'HP_IEX_pi[u]tr'
'HP_IEX_so4[u]tr'
'HP_IEX_fe2[u]tr'
'HP_IEX_fe3[u]tr'
'HP_IEX_thm[u]tr'
'HP_IEX_h2o[u]tr'
'HP_IEX_h[u]tr'
'HP_IEX_pime[u]tr'
'HP_IEX_glc_D[u]tr'
'HP_IEX_o2[u]tr'};
SP_IEX_exchanges_lb=[-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-1
-0.001];
spids=findRxnIDs(modelJoint,sp_IEX_exchanges);
modelJoint.lb(spids)=1*SP_IEX_exchanges_lb;
%%
%ids=findRxnIDs(modelJoint,'SP_RS_139');
%modelJoint.lb(ids)=0;
%%
lunglikeexchanges={'Gastric_EX_arg_L[e]b'
'Gastric_EX_Lcystin[e]b'
'Gastric_EX_his_L[e]b'
'Gastric_EX_ile_L[e]b'
'Gastric_EX_leu_L[e]b'
'Gastric_EX_lys_L[e]b'
'Gastric_EX_met_L[e]b'
'Gastric_EX_phe_L[e]b'
'Gastric_EX_thr_L[e]b'
'Gastric_EX_trp_L[e]b'
'Gastric_EX_tyr_L[e]b'
'Gastric_EX_val_L[e]b'
'Gastric_EX_btn[e]b'
'Gastric_EX_chol[e]b'
'Gastric_EX_pnto_R[e]b'
'Gastric_EX_fol[e]b'
'Gastric_EX_ncam[e]b'
'Gastric_EX_pydxn[e]b'
'Gastric_EX_ribflv[e]b'
'Gastric_EX_thm[e]b'
'Gastric_EX_inost[e]b'
'Gastric_EX_ca2[e]b'
'Gastric_EX_mg2[e]b'
'Gastric_EX_k[e]b'
'Gastric_EX_hco3[e]b'
'Gastric_EX_na1[e]b'
'Gastric_EX_pi[e]b'
'Gastric_EX_glc_D[e]b'
'Gastric_EX_o2[e]b'
'Gastric_EX_cl[e]b'
'Gastric_EX_so4[e]b'};
bloodbounds=[0.007854014
0.00403395
0.004072988
0.015662349
0.015662349
0.015726769
0.003972188
0.007891414
0.015915457
0.001547336
0.007861179
0.015850276
3.23419E-05
5.63672E-05
1.65438E-05
1.78944E-05
6.46837E-05
3.86834E-05
2.09878E-06
2.34167E-05
8.76824E-05
0.014218764
0.006424926
0.042087542
0.206679894
0.008005782
0.04384119
2
0.981506585
0.006424926
1.139885955];
bloodids=findRxnIDs(modelJoint,lunglikeexchanges);
modelJoint.lb(bloodids)=-1*bloodbounds;
optimizeCbModel(modelJoint,'max')
%%
modelHP6new=modelJoint;
save('modelHP6new','modelHP6new');