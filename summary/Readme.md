### Quality Control
- [x]  Generate Raw mSet File
- [x]  Extract median intensities, check for outliers
- [x]  Assess bisulfite conversion efficiency
- [x]  Sex Check
- [x]  Relatedness check using SNP probes
- [x]  p-filter mSet
- [x]  Cell type deconvolution
  - [x]   CETs reference
  - [x]   Internal three-cell reference
- [x]  Epigenetic clock test
- [x]  Dasen normalisation
  - [x]  Normalise all together
  - [x]  Normalise CNG and PFC regions seperately

### Data Exploration
- [x] Principal Component generation
- [x] PC correlation assessment with known potential sources of variation
   If significant confounding effects apparent 
- [x] ComBat batch correction


### Epigenome Wide Association Analysis
***Neuropathology Hypothesis 1: There is an associated methylation signature common across all lewy body positive cases***
Outcome : Binary disease status, controls vs. all LB positive cases
- [ ] Linear mixed model OR Logistic model
    - [ ] Random effects cross regional analysis
    - [ ] Linear individual region analysis
    
***Neuropathology Hypothesis 2: There is an associated methylation signature with linear Braak LB pathology staging***    
Outcome: Lewy Body Braak Stage    
- [x] Linear mixed model OR Logistic model
    - [x] Random effects cross regional analysis
    - [x] Linear individual region analysis
    
    
***Clinical Hypothesis 1: There is are associated methylation signatures across the differing LBD subtypes***     
Outcome: Three group comparison of control : DLB : PDD
- [x] Regress out confounding variable effects
  - [x] ANOVA with Tukey's HSD post hoc individual regional analysis
  - [ ] Cross regional analysis??

**Differentially Methylated Region Analysis**


**Meta Analysis and Validation**


**Network Analysis**


**Machine Learning Based Differentiation**


**Machine Learning Based Differentiation**
