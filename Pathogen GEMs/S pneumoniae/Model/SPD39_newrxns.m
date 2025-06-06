clear
load('SPD39_modified.mat')
model2=SPD39;
%% add additional reactions to the model
rxns_new={'achms[c] + trdrd[c] + tsul[c] <=> ac[c] + h[c] + hcys-L[c] + so3[c] + trdox[c] '
'cys-L[c] + h2o[c] -> h2s[c] + nh4[c] + pyr[c]'
'acser[c] + trdrd[c] + tsul[c] <=> ac[c] + cys-L[c] + h[c] + so3[c] + trdox[c] '
'2 fe2[c] + h[c] + nad[c] <=> 2 fe3[c] + nadh[c] '
'4 fe2[c] + 4 h[c] + o2[c] -> 4 fe3[c] + 2 h2o[c] '
'cgly[c] + glu-L[c] -> gthrd[c] + h2o[c]'
'h2o[c] + lgt-S[c] -> gthrd[c] + h[c] + lac-D[c] '
'h2o2[c] + h[c] + nadh[c] <=> 2 h2o[c] + nad[c] '
'no2[c] + no3[e] -> no2[e] + no3[c] '
'h[c] + nadph[c] + no3[c] -> h2o[c] + nadp[c] + no2[c] '
'atp[c] + trdrd[c] -> datp[c] + h2o[c] + trdox[c] '
'gtp[c] + trdrd[c] -> dgtp[c] + h2o[c] + trdox[c] '
'2 h[c] + 2 o2s[c] -> h2o2[c] + o2[c] '
'h2o2[c] + trdrd[c] -> 2 h2o[c] + trdox[c] '
'2 gthrd[c] + tsul[c] <=> gthox[c] + h2s[c] + so3[c] '};
rxnNamesnew={'AHSERL3'
'CYSDS'
'CYSS3r'
'FE2DH'
'FERO'
'GGTA'
'GLYOX'
'H202D'
'NO3t7'
'NTRARy'
'RNTR1'
'RNTR2'
'SPODM'
'THIORDXi'
'TSULST'};
rxns_newsubsytem={'Methionine and cysteine metabolism'
'Methionine and cysteine metabolism'
'Methionine and cysteine metabolism'
'Transport, extracellular'
'Vitamin B2 metabolism'
'Glutathione metabolism'
'Pyruvate metabolism'
'Miscellaneous'
'Transport, extracellular'
'Nitrogen metabolism'
'Nucleotide interconversion'
'Nucleotide interconversion'
'ROS detoxification'
'ROS detoxification'
'Sulfur metabolism'};
for i = 1:length(rxnNamesnew)
model2 = addReaction(model2,rxnNamesnew{i,1},'reactionFormula',rxns_new{i,1},'subSystem',rxns_newsubsytem{i,1});
end

%% add transport reactions to the following metabolites 
model2 = addReaction(model2,'no2_tr','reactionFormula','no2[e] <=> no2[c]');
model2 = addReaction(model2,'lac-D_tr','reactionFormula','lac-D[e] <=> lac-D[c]');
model2 = addReaction(model2,'o2s_tr','reactionFormula','o2s[e] <=> o2s[c]');

%% add exchange to the superoxide reaction
model2=addExchangeRxn(model2,{'o2s[e]'});
rsexchanges={'EX_o2s[e]'};
rsids=findRxnIDs(model2,rsexchanges);
model2.lb(rsids)=-1000;
SPD39=model2;