# ssylla

Ssylla is an ensemble protein secondary structure prediction method that combines a template-based tool with a non-template-based tool.

To use Ssylla, you need to following tools:
Nnessy: available at https://github.com/gkrieg/nnessy/
Porter: available at https://github.com/mircare/Porter5/

A wrapper function is provided in runhybrid that will take the Nnessy and Porter executables and run a hybrid prediction with runhybrid <input file name> <desired output location>
  
The wrapper function runhybridfromfile takes in the secondary structure prediction of a tool and hybridizes it with Nnessy's prediction. It is run with runhybridfromfile <input fasta for Nnessy> <output location for Nnessy> <secondary structure prediction from second tool>

The datasets for the paper based on this work is available at https://github.com/gkrieg/nnessy/ under the datasets folder
