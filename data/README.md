The complete data (no missing values) for the [mouse challenge](http://www.crm.umontreal.ca/2016/Genetics16/competition_e.php) are in mouse.csv.
There are 615 rows and the header is:
```
"id","cen","strain.id","sex","geno","zyg","litter","IMPC_HEM_001_001","IMPC_HEM_002_001","IMPC_HEM_003_001","IMPC_HEM_004_001","IMPC_HEM_005_001","IMPC_HEM_006_001","IMPC_HEM_007_001","IMPC_HEM_008_001","IMPC_HEM_019_001","IMPC_HEM_027_001","IMPC_HEM_029_001","IMPC_HEM_030_001","IMPC_HEM_031_001","IMPC_HEM_032_001","IMPC_HEM_033_001","IMPC_HEM_034_001","IMPC_HEM_035_001","IMPC_HEM_036_001","IMPC_HEM_037_001","IMPC_HEM_038_001","IMPC_HEM_039_001","IMPC_HEM_040_001"
```

Specifically: 
- "id" mice have each a unique id 
- "cen" in the competition they were specifically selected to come from the same centre (10)
- "strain.id" in the competition data they endep up being all of the same strain (35)
- "sex" male/female
- **"geno"** describes the gene knockout, and specifically 0 = wild-type mice, while the other values describe the gene knockout
- "zyg" is always 1 in the selected data
- "litter" is the id of the litter each mouse is part of 
- the rest of the columns are phenotypes, described as:
```
290 IMPC_HEM_001_001                        White blood cell count        numerical Hematology
291 IMPC_HEM_002_001                          Red blood cell count        numerical Hematology
292 IMPC_HEM_003_001                                    Hemoglobin        numerical Hematology
293 IMPC_HEM_004_001                                    Hematocrit        numerical Hematology
294 IMPC_HEM_005_001                              Mean cell volume        numerical Hematology
295 IMPC_HEM_006_001                   Mean corpuscular hemoglobin        numerical Hematology
296 IMPC_HEM_007_001            Mean cell hemoglobin concentration        numerical Hematology
297 IMPC_HEM_008_001                                Platelet count        numerical Hematology
298 IMPC_HEM_019_001                          Mean platelet volume        numerical Hematology
299 IMPC_HEM_027_001             Red blood cell distribution width        numerical Hematology
503 IMPC_HEM_029_001                 Neutrophil differential count        numerical Hematology
504 IMPC_HEM_030_001                         Neutrophil cell count        numerical Hematology
505 IMPC_HEM_031_001                 Lymphocyte differential count        numerical Hematology
506 IMPC_HEM_032_001                         Lymphocyte cell count        numerical Hematology
507 IMPC_HEM_033_001                   Monocyte differential count        numerical Hematology
508 IMPC_HEM_034_001                           Monocyte cell count        numerical Hematology
509 IMPC_HEM_035_001                 Eosinophil differential count        numerical Hematology
510 IMPC_HEM_036_001                         Eosinophil cell count        numerical Hematology
511 IMPC_HEM_037_001                           Basophil cell count        numerical Hematology
512 IMPC_HEM_038_001                   Basophil differential count        numerical Hematology
513 IMPC_HEM_039_001              Large Unstained Cell (LUC) count        numerical Hematology
514 IMPC_HEM_040_001 Large Unstained Cell (LUC) differential count        numerical Hematology
```

Some of the variables are deterministically related to each other, e.g. the White blood cell count is the sum of the counts of neutrophils, lymphocytes, monocytes, eosinophils, basophils and large unstained cells.

Just for completeness, these two files are from the original challenge:
- impc_quantitative_phen_submatrix.RData (that contains more than the Hematology phenotypes from centre 10 used in the challenge), notably:
-- pmet contains the description of the different phenotypes
-- final contains the data used in the competition
- extract_data_mice.R

The batch-preprocessing/ folder contains some scripts for creating datasets by randomly selecting 3 of each phenotypes and 2 of the gene knockouts (+ wild type mice), and creating two types of files:
-  real-{DATASET_COUNTER}.csv - just a subset of the original data for these 3 phenotypes (renamed v1, v2, v3) and the 2 knockouts (+ wild type) which have an index in the newly created "experiment" column (1=wild-type, 2 and 3 genetically modified)
-  nan-{DATASET_COUNTER}.csv in which a randomly selected phenotype in a randomly selected genetically modified mice is set to NAN

There are a few examples of these files in the folder so that it is maybe clearer.
