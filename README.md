# Host-pathogen-interaction-models
Codes used to build host-pathogen interaction models based on Recon 3D and Dual RNA-seq of pathogen and host. Two cases were considered. First case was pneumonia and the second case was gastric cancer. The RNA seq data of the infected host tissues belonging to the two cases were taken from the following sources 
1. Aprianto, Rieza, et al. "Time-resolved dual RNA-seq reveals extensive rewiring of lung epithelial and pneumococcal transcriptomes during early infection." _Genome biology_ 17 (2016): 1-16.
2. Sharafutdinov, Irshad, et al. "Early and late genome-wide gastric epithelial transcriptome response during infection with the human carcinogen Helicobacter pylori." _Cell Insight_ 1.3 (2022): 100032.

## Description of the codes in the following folders
### 1. Pathogen GEMs
   Two pathogen genome-scale metabolic models were used in this project. They are _Streptococcus pneumoniae_ D39 and _Helicobacter pylori_. The models were taken from the following sources.
   1. Pedram, Narges, Hamid Rashedi, and Ehsan Motamedian. "A systematic strategy using a reconstructed genome-scale             metabolic network for pathogen Streptococcus pneumoniae D39 to find novel potential drug targets." Pathogens and           Disease 78.6 (2020): ftaa051.
   2. Thiele, Ines, et al. "Expanded metabolic reconstruction of helicobacter pylori (i it341 gsm/gpr): an in silico             genome-scale characterization of single-and double-deletion mutants." Journal of bacteriology 187.16 (2005): 5818-         5830.

   These models were further modified by adding a few additional reactions and removing thermodynamically infeasible cycles. Further the nomenclature of the metabolites were changed according to the nomenclature of the Rexon 3D derived host GEMs. These models were constrained using the lunglike media for _Streptococcus pneumoniae_ D39 and a minimal media for _Helicobacter pylori_ to support their growth respectively. The nutrient media were taken from the sources of their respective GEMs.
   
For the _H.pylori_ GEM: [...Host-pathogen-interaction-models/Pathogen GEMs/H pylori/Model](https://github.com/Subasree514/Host-pathogen-interaction-models/tree/main/Pathogen%20GEMs/H%20pylori/Model) <br>
hp_model_modifications.m - Updated the reaction bounds and nomenclature <br>
iit341_newrxns.m - Added a few additional reactions

[...Host-pathogen-interaction-models/Pathogen GEMs/H pylori/TICs removal/](https://github.com/Subasree514/Host-pathogen-interaction-models/tree/main/Pathogen%20GEMs/H%20pylori/TICs%20removal) <br>
iit341_changes.m - Removed TICs and update the nomenclature

For the _S.pneumoniae_ model: [...Host-pathogen-interaction-models/Pathogen GEMs/S pneumoniae/Model](https://github.com/Subasree514/Host-pathogen-interaction-models/tree/main/Pathogen%20GEMs/S%20pneumoniae/Model) <br>
SPD39_newrxns - Added a few additional reactions <br>
sp_lunglikemedia.m - Updated the constraints of the lunglike media <br>
SPmetsmod.m - Updated the nomenclature

[......Host-pathogen-interaction-models/Pathogen GEMs/S pneumoniae/TICs removal/](https://github.com/Subasree514/Host-pathogen-interaction-models/tree/main/Pathogen%20GEMs/S%20pneumoniae/TICs%20removal) <br>
sp_removetics.m - Removed TICs

### 2. Infected host model
     Following steps were taken to build the infected host models from the CCLE cancer cell line GEMs as both the dual RNA seq experiments were done on the cell lines. The detailed methodology common to both the scenarious is shown in the figure below.
   ![image-url](https://github.com/Subasree514/Host-pathogen-interaction-models/blob/main/hp.png)
- The first step to reconstruct the two cancer cell line GEMs were same as the one followed in S. Sridhar et al.,_Metabolic Engineering_,2023 (doi.org/10.1016/j.ymben.2023.08.006). The code is available in https://github.com/Subasree514/Building-cancer-specific-genome-scale-models.
- The second step was to map the transcriptomic data of the infected tissues to the cancer cell line GEMs and reconstruct models using GIMME algorithm. The code is available in https://github.com/Subasree514/Host-pathogen-interaction-models/blob/main/Infected%20host%20model/Step1_gimme.m. The GIMME derived models could not be used for flux sampling analysis. So an additional step was also considered. 
- The final step was to get the reactions from the GIMME models and use them to reconstruct the infected host tissue GEMs with the modified Recon 3D model through FASTCORE algorithm, given in https://github.com/Subasree514/Host-pathogen-interaction-models/blob/main/Infected%20host%20model/Step2_fastcore.m

### 3. HP model building
      Codes used to reconstruct the host-pathogen interaction models for pneumonia and gastric cancer.
   1. For pneumonia, codes to integrate the RS integrated host and pathogen models along with the inclusion of                   constraints for lunglike media to the pathogen and minimal media to the host are available in https://github.com/Subasree514/Host-pathogen-interaction-models/blob/main/HP%20model%20building/PneumoniaHPmodel.m
   2. For gastric cancer, codes to integrate the RS integrated host and pathogen models along with the inclusion of              minimal media constraints for each of the model are available in https://github.com/Subasree514/Host-pathogen-interaction-models/blob/main/HP%20model%20building/GastriccancerHPmodel.m

### Additional details
RSmodel_Recon3D1_2023_ecoli.mat - Reconstructed the reactive species reactions model (RS model) from the human RS model (S. Sridhar et al.,_Metabolic Engineering_,2023) by retaining only the reactive species and reactions that are relevant to be present in a prokaryote _Escherichia coli_. This RS model is integrated to the two pathogen GEMs discussed in this project. We created two scenarios, one is the basal RS levels at the zeroth time point and the normal levels at the later time points. For the zeroth time point, the bounds of the RS reactions are adjusted such that the biomass of the RS integrated pathogen model is same to the biomass of the pathogen model. The bounds of the RS reactions are adjusted such that the biomass of the RS integrated pathogen model is only slighly greater (~1.5 times) than the biomass of the pathogen model for the later time points. <br>

RS model taken from S. Sridhar et al.,_Metabolic Engineering_,2023 is intgegrated to the host tissue GEMs. For the host GEMs at the zeroth time point in both the cases, the bounds of the RS reactions integrated to the host tissue models are kept at basal levels (may be -1 to 1) and [-500, 500] for the second time point and full levels [-1000, 1000] to the final time point of infection.

##### Important note
Duplicate reactions and TICs generated by the RS reactions are removed in the RS integrated GEMs using the default cobratoolbox functions like checkDuplicateRxn and findMinNull respectively.
