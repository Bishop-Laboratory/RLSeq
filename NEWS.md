# RLSeq 1.3.0

Bioc 3.16 release with the addition of a noise analysis function.

# RLSeq 1.2.0

Same as 1.1.0 but released on Bioc 3.15

# RLSeq 1.1.0

Bioc 3.14 release

# RLSeq 1.0.1

Updated with RLHub update. Added 118 new samples to the database. 

# RLSeq 1.0.0

First official release on bioconductor -- [link](https://bioconductor.org/packages/release/bioc/html/RLSeq.html).

# RLSeq 0.99.11

Pre-release. This version is prior to the first
official release (v1.0.0), anticipated in Oct. 2021.

## Package

* Base class
  - RLRanges stores R-loop data and RLSeq results.
* Workflow
  -The wrapper function RLSeq() performs all analysis
  steps. However, individual functions are also accessible.
* Model
  - The model for predicting sample quality label is
  updated to the latest version, incorporating 231 samples
  in training. This model showed high accuracy (.9304) on a 
  test set of 115 samples. 
