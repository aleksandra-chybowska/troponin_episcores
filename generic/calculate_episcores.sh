#!/bin/bash

trainEpiScores() {
    echo "Calculate EpiScores for a set of traits - procedure begins"
    SETTINGS=$1 #"/home/shirin/Projects/Sequencing/Data_albicanis/VCF"

    for TRAIT in $SETTINGS/*
    do
        echo "\nTrait: $TRAIT"
        
        Rscript --vanilla /Users/shirin/Projects/R/troponin_episcores/generic/data_prep.R --settings $TRAIT
        # Rscript --vanilla ~/Cluster_Filespace/Marioni_Group/Ola/Code/troponin_episcores/generic/data_prep.R --settings $TRAIT
    done

    echo "\nProcedure finished."
}

trainEpiScores /Users/shirin/Projects/R/troponin_episcores/generic/settings/test_settings