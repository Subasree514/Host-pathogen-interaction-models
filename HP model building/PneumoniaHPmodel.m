% RS integrated host model. The same procedure applies to the other infection points host tissue models
load('model_sp0_rs.mat')
model1=model_sp0_rs;
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
%% remove duplicate reactions
[model1, removedRxnInd1, keptRxnInd1] = checkDuplicateRxn(model1,'S');
[model1, removedRxnInd2, keptRxnInd2] = checkDuplicateRxn(model1,'FR');

%% update exchanges of original host model
rxnlist=contains(model1.rxns,'EX_');
find(rxnlist==1);
model1.lb(ans)=-1000;
modelHost=model1;

%% load the RS integrated pathogen
load('SPD39_rsmin_host.mat');
findRxnIDs(SPD39_rsmin_host,{'Exbtn'});
SPD39_rsmin_host.rxns(ans)={'EXbtn'};

%% default methodology to integrate th host and pathogen models
models={};
nameTagsModels={};
models{1,1}=SPD39_rsmin_host;
bioID={};
bioID{1,1}=SPD39_rsmin_host.rxns(find(strncmp(SPD39_rsmin_host.rxns,'bioP',7)));
nameTagsModels{1,1}=strcat('SP_');
nameTagHost='Lung_';
%reducedDietConstraints={};
[modelJoint] = createMultipleSpeciesModel(models,'nameTagsModels',nameTagsModels,'modelHost',modelHost,'nameTagHost',nameTagHost);

%% Nutreint uptake for the host and pathogen
%% constrain the uptake reactions as lunglike media for the pathogen
modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('SP_IEX_',modelJoint.rxns)),0,'l');
modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('SP_IEX_',modelJoint.rxns)),1000,'u');
sp_IEX_exchanges={'SP_IEX_acgam[u]tr'
'SP_IEX_ala_L[u]tr'
'SP_IEX_arg_L[u]tr'
'SP_IEX_asn_L[u]tr'
'SP_IEX_asp_L[u]tr'
'SP_IEX_cys_L[u]tr'
'SP_IEX_gln_L[u]tr'
'SP_IEX_glu_L[u]tr'
'SP_IEX_gly[u]tr'
'SP_IEX_his_L[u]tr'
'SP_IEX_ile_L[u]tr'
'SP_IEX_leu_L[u]tr'
'SP_IEX_lys_L[u]tr'
'SP_IEX_met_L[u]tr'
'SP_IEX_phe_L[u]tr'
'SP_IEX_pro_L[u]tr'
'SP_IEX_ser_L[u]tr'
'SP_IEX_thr_L[u]tr'
'SP_IEX_trp_L[u]tr'
'SP_IEX_val_L[u]tr'
'SP_IEX_mg2[u]tr'
'SP_IEX_Cl[u]tr'
'SP_IEX_ca2[u]tr'
'SP_IEX_pi[u]tr'
'SP_IEX_k[u]tr'
'SP_IEX_nh4[u]tr'
'SP_IEX_cit[u]tr'
'SP_IEX_na1[u]tr'
'SP_IEX_ac[u]tr'
'SP_IEX_pyr[u]tr'
'SP_IEX_so4[u]tr'
'SP_IEX_zn2[u]tr'
'SP_IEX_ade[u]tr'
'SP_IEX_gua[u]tr'
'SP_IEX_ura[u]tr'
'SP_IEX_xan[u]tr'
'SP_IEX_thm[u]tr'
'SP_IEX_ribflv[u]tr'
'SP_IEX_nac[u]tr'
'SP_IEX_pydam[u]tr'
'SP_IEX_pydxn[u]tr'
'SP_IEX_fol[u]tr'
'SP_IEX_chol[u]tr'
'SP_IEX_cyncblm[u]tr'
'SP_IEX_4abz[u]tr'
'SP_IEX_lipoate[u]tr'
'SP_IEX_pnto_R[u]tr'
'SP_IEX_btn[u]tr'
'SP_IEX_o2[u]tr'};
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
-1];
spids=findRxnIDs(modelJoint,sp_IEX_exchanges);
modelJoint.lb(spids)=1*SP_IEX_exchanges_lb;

%% set the uptake of nutreints from blood as minimal media metabolites
modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('Lung_EX_',modelJoint.rxns)),0,'l');
modelJoint = changeRxnBounds(modelJoint,modelJoint.rxns(strmatch('Lung_EX_',modelJoint.rxns)),1000,'u');
bloodexchanges={'Lung_EX_arg_L[e]b'
'Lung_EX_Lcystin[e]b'
'Lung_EX_his_L[e]b'
'Lung_EX_ile_L[e]b'
'Lung_EX_leu_L[e]b'
'Lung_EX_lys_L[e]b'
'Lung_EX_met_L[e]b'
'Lung_EX_phe_L[e]b'
'Lung_EX_thr_L[e]b'
'Lung_EX_trp_L[e]b'
'Lung_EX_tyr_L[e]b'
'Lung_EX_val_L[e]b'
'Lung_EX_btn[e]b'
'Lung_EX_chol[e]b'
'Lung_EX_pnto_R[e]b'
'Lung_EX_fol[e]b'
'Lung_EX_ncam[e]b'
'Lung_EX_pydxn[e]b'
'Lung_EX_ribflv[e]b'
'Lung_EX_thm[e]b'
'Lung_EX_inost[e]b'
'Lung_EX_ca2[e]b'
'Lung_EX_mg2[e]b'
'Lung_EX_k[e]b'
'Lung_EX_hco3[e]b'
'Gastric_EX_pi[e]b'
'Gastric_EX_glc_D[e]b'
'Gastric_EX_o2[e]b'
'Gastric_EX_cl[e]b'
'Gastric_EX_so4[e]b'
'Gastric_EX_na1[e]b'};
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
bloodids=findRxnIDs(modelJoint,bloodexchanges);
modelJoint.lb(bloodids)=-1*bloodbounds;
optimizeCbModel(modelJoint,'max')
modelSP0new=modelJoint;
save('modelSP0new','modelSP0new');