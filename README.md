# DeepMirCut

# Setup
All versions of python packages needed to run DeepMirCut can be found in the deepMirCut_env.yml file.  The environment can be easily installed using conda.
```sh
$ conda env create -f deepMirCut_env.yml file
```

DeepMirCut works best when a dot-bracket fold and bpRNA structure array are provided.  Please install the following:

* [The ViennaRNA Package](https://www.tbi.univie.ac.at/RNA/) - The ViennaRNA Package includes RNAfold, a program that folds sequences and returns their dot-bracket structure array.
* [bpRNA](https://github.com/hendrixlab/bpRNA) -  The bpRNA script annotates the features of RNA secondary structures.
* [Graph.pm](https://metacpan.org/pod/distribution/Graph/lib/Graph.pod) - A dependency needed to run bpRNA, available on CPAN

After installing these dependencies run the following script to set up deep-mir-cut.  

```sh
$ bash setup.sh
```

# Data Preparation

The prepareData.pl script takes in a fasta file and uses RNAfold and bpRNA to find the dot-bracket and bpRNA structure array for each sequence. 

```sh
$ perl prepareData.pl examples.fa prepared_examples.txt
```

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
$ python plot_predictions.py  output_classification_DVs.txt examples.txt ex445
```

Producing a 3D Graph of output:
```sh
$ python3 deepMirCut_predict.py examples.txt -m seqBPRNA.model -o output --input_setting 2 -d
$ python plot_predictions_3d.py  output_classification_DVs.txt examples.txt ex445
```

# Training Models

The trainModel.py script was used to train the models contained in this github repository. 

### Examples

```sh
python3 trainModel.py trainSet.txt -s validationSet.txt -m trained_model.model -o output --embedding_layer_output 32 --embedding_dropout 0.417 --bi_lstm1_units 64 --bi_lstm2_units 160 --learning_rate 0.00357 --epsilon 1.34896288259165e-07 --input_setting 2
```

# Testing Performance

The testModel.py script may be used to test the performance of the trained model on a test or validation set.

```sh
python3 testModel.py testSet.txt -m trained_model.model -o output --input_setting 2
```

# Datasets

Sample datasets and scripts used to produce them can be found in the [DeepMirCut_Data repository](https://github.com/JimBell/deepMirCut_data)
