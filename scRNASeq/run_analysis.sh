#!/bin/bash

####################################################################
# Modify this variables by adding the path to your singularity image
container="ENTER CONTAINER PATH"
####################################################################

#mkdir data
#mkdir procdata
#mkdir figures

echo ${container}


singularity exec ${container} jupytext --sync 01_processing.py
singularity exec ${container} jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 01_processing.ipynb

singularity exec ${container} jupytext --sync 02_DE.R
singularity exec ${container} jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 02_DE.R

singularity exec ${container} jupytext --sync 03_projections.py
singularity exec ${container} jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 03_projections.ipynb

singularity exec ${container} jupytext --sync 04_dotscore.py
singularity exec ${container} jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 04_dotscore.ipynb

singularity exec ${container} jupytext --sync 05_varfigures.py
singularity exec ${container} jupyter nbconvert --to notebook --inplace --execute --ExecutePreprocessor.timeout=0 04_varfigures.ipynb
