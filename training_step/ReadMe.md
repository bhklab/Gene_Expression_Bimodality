## Training step to run RLOBICO

To run the training step:
- The user needs to navigate to the "training_step" directory.
- Inside that directory, there is an R script file named "run_1drug_model_withSpaces.R". This file is used to run the modeling step for a drug. 
- Arguments needed to run the script:
  * Drug name
  * Dataset (e.g. CTRPv2)
  * K_CV: the number of folds for cross validation
  * numOfSolutions: the number od solutions produced by mRMRe for feature selection.
  * outputPath: output directory
  * dictionary: a named array of drugs. Used to give coded names to drugs that have special chars or spaces
  * data: data needed for different annotations
