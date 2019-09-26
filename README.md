**The R package for ‘mirPLS: a partial linear structure identifier method for miRNA data.’** 

**Example data:**

- x.rda: A 300 subjects * 300 miRNAs matrix. The first 6 columns are 6 miRNAs have non-linear associations, the next 14 columns are 14 miRNAs have non-linear associations and the others are miRNAs have no associations.

- y.rda: A vector indicating whether a subject is a case or a control (there are 200 cases and 100 controls). 1s represents cases and 0s represents controls.

- true_label.rda: A vector indicating the true clustering lables of 200 cases.

- The cancer datasets used in the paper are available in the Cancer Genome Atlas (TCGA): https://portal.gdc.cancer.gov/

**Code:**

- cv_mirPLS: determine regularization parameters in mirPLS.

- mirPLS: select the miRNAs linearly or non-linearly associated to outcomes.

- predict_mirPLS: predicting the regression mean of new data for mirPLS.

- sp_clustering: perform the famous spectral clustering algorithms.


**Tutorial:**

The code has been tested on the following systems: macOS 10.14.6 and Windows 10, with R 3.5.2 installed.

Load sample data x.rda and y.rda or download miRNA data from TCGA.

Run cv_mirPLS to determine the regularization parameters: lambda1 and lambda2. We suggest to search lambda1 and lambda2 in the ranges of (0, 1) and (0, 10), respectively. With lambda1 and lambda2, run mirPLS to select miRNAs linearly or non-linearly associated with disease outcomes. Finally, run sp_clustering to conduct spectral clustering on patients using the selected non-linearly associated miRNAs. In the whole process, cv_mirPLS are time consuming. It may take hours to run the whole process for a typical miRNA data from TCGA.
