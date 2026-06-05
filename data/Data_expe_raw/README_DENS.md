RDG README File Template --- General --- Version: 0.1 (2026-06-03)

This README file was generated on 2026-06-03 by Raphaël Royauté.

Last updated: 2026-06-03.

# GENERAL INFORMATION:

## Dataset title: "Predicting Energy Assimilation from Resource Quantity and Quality: A Mechanistic DEB Approach for Aporrectodea caliginosa"

## DOI:

# Contact email: raphael.royaute\@inrae.fr

# METHODOLOGICAL INFORMATION

## Environmental/experimental conditions:

Monitoring of earthworm growth and reproduction from hatching under different population densities and feeding regimes in laboratory microcosms.

## Methods for processing the data:

All processing steps are available at the [following website](https://lisagllt.github.io/Acali-DEB/) which explains all data analyses steps.

# DATA & FILE OVERVIEW

All data are stored under the `Data_DENS_growth.csv` & `Data_DENS_repro.csv` files.

## Variables names :

The `Data_DENS_growth.csv` file contains 13 columns and 775 rows. The columns defintions and units are as follows:

-   `date` : Date of the measurement (DD-MM-YYY format)
-   `Heure` : Hour of the measurement (HH-MM-SS format)
-   `t` : Number of days since the start of the experiment (range: 0-154)
-   `Lot` : Two batches of earthworms were use and started at different times (A & B).
-   `Condition` : Condition names which refers to the number of individual in the cosm and the amount of food given.
-   `Nb_rep` : Replicate number (range: 1-4).
-   `Nb_vdt` : Number of earthworms in the corresponding cosm (range: 1-7).
-   `Food` : Amount of horse dung (dry equivalent) given in each cosm per 14 days (range: 0.4-3).
-   `w` : Fresh weight of the earthworm (mg, range 4-1499).
-   `Diapause` : Number of earthworm in diapause (range 0-0).
-   `Status` :
    -   “J” : Juvenile
    -   “A” : Adulte
    -   “D” : Dead
-   `Keep` : Selection criteria (Yes-No). Selects data from cosm where no earthworm died.
-   `Comment` : Any comment about the measurement.

The `Data_DENS_repro.csv` file contains 11 columns and 126 rows. The columns defintions and units are as follows:

-   `date` : Date of the measurement (DD-MM-YYY format)
-   `t` : Number of days since the start of the experiment (range: 0-154)
-   `Lot` : Two batches of earthworms were use and started at different times (A & B).
-   `Condition` : Condition names which refers to the number of individual in the cosm and the amount of food given.
-   `Nb_rep` : Replicate number (range: 1-4).
-   `Nb_vdt` : Number of earthworms in the corresponding cosm (range: 1-7).
-   `Food` : Amount of horse dung (dry equivalent) given in each cosm per 14 days (range: 0.4-3).
-   `Nb_cocoons` : Number of cocoon collected in the cosm which is the number of cocoon produced by the earthworm in one cosm during the last 28 days (range: 0-49).
-   `Nb_cocoons_ind` : `Nb_cocoons` divided by `Nb_vdt` (range: 0-18).
-   `Reproduction` : Mean cumulated reproduction per individual (range: 0-36).
-   `Keep` : Selection criteria (Yes-No). Selects data from cosm where no earthworm died.

## Missing data codes

Missing data are represented by empty cells which are transformed to `NA` by default during data processing with R.
