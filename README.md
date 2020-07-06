# ssylla

Under Construction.

Ssylla is an ensemble protein secondary structure prediction method that combines a template-based tool with a non-template-based tool.

To use Ssylla, you need to following tools:
Nnessy: available at https://github.com/gkrieg/nnessy/
Porter: available at https://github.com/mircare/Porter5/

The simplest way to use Ssylla is to take the output files from Nnessy and estimate their accuracy with nnessyaccest.py, with usage python3 nnessyaccest.py \< Nnessy probabilities output file> < Nnessy percent identity output file ends in i  > < "3state" or "8state" >

The output is the estimated accuracy of the Nnessy prediction. If this estimated accuracy exceeds 0.8, use the Nnessy prediction. Otherwise, use the prediction from Porter.

A wrapper function will soon be added that will run Nnessy for you and then run Porter if needed and return the hybrid prediction.

