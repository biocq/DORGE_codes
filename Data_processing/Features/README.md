# Codes used for training feature profile generation

## Codes to process feature profile can be found at:

https://figshare.com/s/0e7dc381deb8893ee8f6

## Dependency

Java v1.8.0

Java SDK: Sun Java SE Development Kit 7 or higher (Links of JDK 8 https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) is needed. Setting JAVA_HOME may also be needed to launch Java at any folders(see http://www.sajeconsultants.com/how-to-set-java_home-on-mac-os-x/?utm_source=rss&utm_medium=rss&utm_campaign=how-to-set-java_home-on-mac-os-x).

As a convenience way, Bioconda users may also try 'conda install -c bioconda java-jdk' to install Java easily.

## Steps

1. Extract training feature set from the complete feature profile
```
cd Features
Rscript feature_subset_generation.R
cd ..
```

  Input: All_features.csv (Complete feature profile including those not included in the final feature profile)

  Output: All_features_training.csv (Feature profile for DORGE model training)
