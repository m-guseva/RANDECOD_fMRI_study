# RANDECOD fMRI Study 

This repository contains the analysis code and behavioral dataset from our fMRI study investigating the neural correlates of human randomization. 

## Dataset

The mat file `behavioralData.mat` contains each participant's age, sex and responses during the random sequence generation task across all six runs.

Description of columns in `behavioralData.response` of each participant:
- *runNumber*, *blockNumber*, *trialNumber*
- *timingData*: Indicates whether the button was pressed during the allotted time period ('on time') or not ('too late')
- *keys* : 1/2=left button press, 3/4=right button press
- *fixationStart*, *fixationEnd*: Timestamps of start and end of fixation phase in each trial
- *choiceStart*, *choiceTimestamp*: Time stamps of beginning of choice phase and actual button press respectively
- *feedbackStart*, *feedbackEnd*: Timestamps of start and end of feedback phase in each trial
- *dur_choice*: Time elapsed between stimulus onset (*choiceStart*) to button press (*choiceTimestamp*)
- *condition*: "ER": Explicit Randomness, "FC": Free Choice, "MC": Mental Coin Toss
- *headsPosition*: Position of heads, left or right (didn't change throughout the experiment within subject)
- *focusResponse* and *instructionsResponse*: Responses to the 5 point likert-scaled questions pertaining to the focus and instruction adherence (respectively) in each run

This mat file also includes a matrix of pseudorandomly generated binary sequences (*randseq*), which was created with the `createPseudoRandSeq.m` function. It creates 28 sequences of length 540, replicating the sequence length in the experiment, 28 representing the average samples size per condition. 

The excel sheet `FreeFormAnswers` lists all the answers to the post-experiment questionnaire. The internalID corresponds to the internalID in `behavioralData.mat`. The second sheet in this file shows the text used to prompt the free-form responses.

## Behavioral Analysis
Scripts include:
- `LikertAnalysis.m`: Descriptive Analysis of responses to focus and instruction adherence questions
- `ProportionValueScript.m`: Calculation of proportion values and comparison between conditions
- `ReactionTimes.m`: Calculation of reaction times and comparison between conditions
- `RunLengthScript.m`: Calculation of run lengths and comparison between conditions
- `RandomnessMeasures.m`;`markovOrderBIC.m`: Calculation of conditional entropy and optimal Markov order values and comparison between conditions

## Imaging Analysis

This folder contains the preprocessing script of the imaging data as well as the code for the univariate analysis (first&second level) and multivariate analysis (first level, searchlight decoding, second level).
