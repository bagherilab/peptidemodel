# Modulating and evaluating receptor promiscuity through directed evolution and modeling

<sub>LAST UPDATED: 2017-03-25</sub>

## Content overview

- __`scripts`__ : `.pbs` scripts for submitting jobs on [QUEST](http://www.it.northwestern.edu/research/user-services/quest/index.html)
- __`lib`__ : library functions
- __`pipeline`__ : main code for simulation, regression, and analysis

PBS scripts for running grouping exhaustive regression of receptors and single
amino acid variants. Make sure Quest is setup as follows (starred files for shuffling only):

## Running on QUEST

Pipeline was written specifically for running on Northwestern's high performance computing core [QUEST](http://www.it.northwestern.edu/research/user-services/quest/index.html).

#### File structure

The working directory on QUEST contains the following:

```bash
`~/Matlab/`
    logs/ -- standard out and error logs
    results/ -- results stored here
        Vectorizations.mat
        Data_R_*.mat
        Data_P_*.mat
        Data_Shuffle.mat
    *.m (all Matlab files)
    *.pbs (all submission scripts)
```

## Pipeline

`run_regress_P.pbs`

By Property regression.
There are 255 vectorizations, completed in sets of 25, for 9 different responses.
Run using:

```bash
msub -t runP[1-11] run_regress_P.pbs
```

`run_regress_R.pbs`

By Residue regression.
There are 8191 vectorizations, completed in sets of 25, for 9 different responses.
Run using:

```bash
msub -t runR[1-328] run_regress_R.pbs
```

`run_merge_R.pbs` / `run_merge_P.pbs`

Merges results into a single Matlab file.
Run only after previous jobs are completed.

```bash
msub -t mergeP[1-7] run_merge_P.pbs
msub -t mergeR[1-7] run_merge_R.pbs
```

`run_shuffle_P.pbs`

By Property shuffled regression.
There are 255 vectorizations, completed in sets of 25, for 51 reps.
Run using:

```bash
msub -t shuffleP[1-11] run_shuffle_P.pbs
```

`run_shuffle_R.pbs`

By Residue shuffled regression.
There are 8191 vectorizations, completed in sets of 25, for 51 reps.
Run using:

```bash
msub -t shuffleR[1-328] run_shuffle_R.pbs
```
