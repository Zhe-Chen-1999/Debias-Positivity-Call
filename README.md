# Debias-Positivity-Call
Determining vaccine responders using single-cell assays and paired control samples

In vaccine studies, determining whether participants mount an immune response is a critical objective. Cellular immune responses are commonly assessed using single-cell assays such as intracellular cytokine staining (ICS) and B-cell phenotyping (BCP). These assays help profile immune cell phenotypes and cytokine production before and after vaccination. However, assay variations (e.g., batch effects) between time points can result in misclassification of vaccine responders.

This repository provides an implementation of a statistical framework that addresses this challenge by:

1. Incorporating paired control data to account for batch effects.

2. Calculating two types of adjusted p-values:
  - Maximally Adjusted P-value: Conservative adjustment ensuring validity over all plausible batch effects.
  - Minimally Adjusted P-value: Minimal adjustment ensuring the adjusted p-value is not falsified by the control data.

These p-values balance type-I error rates and statistical power in the presence of assay variation.
