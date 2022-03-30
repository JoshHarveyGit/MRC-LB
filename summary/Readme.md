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
- [ ] Flow diagram of sample retention

### Data Exploration
- [x] Principal Component generation
- [x] PC correlation assessment with known potential sources of variation \
   If significant confounding effects apparent 
- [ ] ComBat batch correction: Techincal effects, plate & BB

### Epigenome Wide Association Analysis

***Pre-analysis***
- [ ] Median absolute deviation based filtering of non-variable probes
- [ ] Sex chromosome filtering!

***Neuropathology Hypothesis 1: There is an associated methylation signature common across all lewy body positive cases***
Outcome : Binary disease status, controls vs. all LB positive cases
- [ ] Linear mixed model OR Logistic model
    - [ ] Random effects cross regional analysis
    - [ ] Linear individual region analysis
    
***Neuropathology Hypothesis 2: There is an associated methylation signature with linear Braak LB pathology staging***    
Outcome: Lewy Body Braak Stage    
- [ ] Linear mixed model OR Logistic model
    - [ ] Random effects cross regional analysis
    - [ ] Linear individual region analysis
    
    
***Clinical Hypothesis 1: There is are associated methylation signatures across the differing LBD subtypes***     
Outcome: Three group comparison of control : DLB : PDD
- [ ] Regress out confounding variable effects, extract residuals
  - [ ] ANOVA with Tukey's HSD post hoc individual regional analysis
  - [ ] Cross regional analysis??

***Clinical Hypothesis 2: There is are associated methylation signatures across the differing LBD subtypes, regardless of coincident AD pathology*** \
Outcome: Three group comparison of control/AD : DLB : PDD
- [ ] Regress out confounding variable effects, extract residuals
  - [ ] ANOVA with Tukey's HSD post hoc individual regional analysis
  - [ ] Cross regional analysis??

### Differentially Methylated Region Analysis
- [ ] Comb-p based DMR identification 
  - [ ] Neuropath Hypothesis 1
  - [ ] Neuropath Hypothesis 2
  - [ ] Clinical Hypothesis 1
  - [ ] Clinical Hypothesis 2

### Meta Analysis and Validation
- [ ] Correlate effects size with other potential datasets:

For independent datasets with similar analysis strategy: meta analysis
- [ ] 




### ML based differentiation


### mQTL




