# DeepMirCut

# Setup
All versions of python packages needed to run DeepMirCut can be found in the deepMirCut_env.yml file.  The environment can be easily installed using conda.
```sh
$ conda env create -f deepMirCut_env.yml file
```

DeepMirCut works best when a dot-bracket fold and bpRNA structure array are provided.  Programs which fold RNA's and annotate thier loop structures can be found at the following links:

* [The ViennaRNA Package](https://www.tbi.univie.ac.at/RNA/) - The ViennaRNA Package consists of several programs which are used to fold and compare RNA structures.   A program called RNAfold may be used to obtain a dot-bracket structure array.
* [bpRNA](https://github.com/hendrixlab/bpRNA) -  The bpRNA script is a tool which is able to annotate the features of RNA secondary structures.

# Data Preparation

A script to prepare sequences for analysis will be uploaded shortly...

# Predicting miRNA Cleavage Sites

The deepMirCut_predict.py script may be used for Dicer and Drosha cleavage site prediction.  deepMirCut_predict can base it's predictions on a single model using the -m option or a list of models with the -L option.   There are three options for the type of information to predict on which may be set using the  --input_setting flag: 1) sequence Only, 2) Sequence and dot-bracket structure array, 3) Sequence, dot-bracket, and bpRNA context array.

### Examples

Predicting with sequence, dot-bracket structure array, and bpRNA context:
```sh
$ python3 deepMirCut_predict.py examples.txt -m seqBPRNA.model -o output --input_setting 2
```

Predicting with sequence and dot-bracket structure array:
```sh
$ python3 deepMirCut_predict.py examples.txt -m seqFold.model -o output --input_setting 1
```

Predicting with sequence only:
```sh
$ python3 deepMirCut_predict.py examples.txt -m seqOnly.model -o output --input_setting 0
```

Predicting with sequence, dot-bracket structure array, and bpRNA context using an ensemble:
```sh
$ python3 deepMirCut_predict.py examples.txt -L seqBPRNA_ensemble_list.txt -o output --input_setting 2
```
Note: you may have to edit seqBPRNA_ensemble_list.txt to change the directory for each of the models.

# Visualizing Results

If deepMircut_predict.py is used with the -d option, it will output the decision values output by the algorithm for each nucleotide along the sequence.  Decision values can be graphed using either plot_predictions.py for a 2D visualization or plot_predictions_3d.py for a 3D visualization.

### Examples

Producing a 2D Graph of output:
```sh
$ python3 deepMirCut_predict.py examples.txt -m seqBPRNA.model -o output --input_setting 2 -d
$ python ../Scripts/plot_predictions.py  output_classification_DVs.txt examples.txt ex4857
```

Producing a 3D Graph of output:
```sh
$ python3 deepMirCut_predict.py examples.txt -m seqBPRNA.model -o output --input_setting 2 -d
$ python ../Scripts/plot_predictions.py  output_classification_DVs.txt examples.txt ex4857
```

# Training Models



# Testing Performance

