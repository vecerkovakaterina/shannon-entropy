# Shannon Entropy Calculator for Multiple Sequence Alignment

## Table of Contents

-   [Introduction](#introduction)
-   [Features](#features)
-   [Requirements](#requirements)
-   [Usage](#usage)
-   [Output](#output)

## Introduction {#introduction}

This Python script calculates Shannon entropy for a given multiple sequence alignment (MSA). Shannon entropy is a measure of the diversity or information content at each position in the alignment. This tool can be useful for analyzing the conservation or variability of positions in biological sequences, such as DNA or protein alignments.

### Understanding Shannon Entropy Values

Here's how to interpret the values:

-   **Higher Entropy Values**: Positions with higher entropy values indicate greater variability among the sequences at that position. In biological terms, this suggests lower conservation, as there is more diversity in the amino acids or nucleotides found at that position.

-   **Lower Entropy Values**: Conversely, positions with lower entropy values indicate higher conservation. These positions have a more limited set of amino acids or nucleotides, suggesting that they play a critical role in the structure or function of the biological sequence. Such positions are often associated with functionally or structurally important regions in proteins or genetic elements in DNA.

The entropy values typically range from 0 to a maximum value, with higher values indicating more diversity and lower values indicating more conservation. It's essential to analyze the entropy values in the context of your specific biological problem and domain knowledge.

## Example Interpretation

For instance, if you calculate the Shannon entropy for a protein alignment and find that a certain position has an entropy value of 0.2, this suggests that this position is relatively conserved, as there is limited amino acid diversity at that site. On the other hand, if another position has an entropy value of 1.0, it indicates high variability, implying that this position is not conserved and may not have a critical functional role.

Remember that the interpretation of entropy values should always be done with consideration of the specific biological context and the goals of your analysis.

## Features {#features}

-   Calculate Shannon entropy for each position in an MSA.
-   Generate an IGV track for graphical representation of entropy profiles.
-   Easily customizable for different input formats and options.

## Requirements {#requirements}

-   Python 3.x
-   BioPython
-   click
-   pandas

You can install the required libraries using `conda`:

``` bash
conda create --name <env> --file requirements.txt
```

## Usage {#usage}

To calculate Shannon entropy, use this command:

``` bash
python shannon_entropy.py INPUT_FILENAME
                          [[clustal|phylip|stockholm|fasta]] [OUTPUT_FILENAME]
                          [STUDIED_GENOME_ACCESSION] [MAX_FRACTION_GAPS]
                          [WINDOW_SIZE] [CHROMOSOME] [WIG_START_INDEX]
```

Arguments:
- **[INPUT_FILENAME]**: Input file with multiple sequence alignment.
- **[MSA_FORMAT]**: File format of the MSA file, can be either clustal, phylip, stockholm or fasta.
- **[OUTPUT_FILENAME]**: Output file name.
- **[STUDIED_GENOME_ACCESSION]**: If calculating Shannon entropy for a specific genome, provide the accession number.
- **[MAX_FRACTION_GAPS]**: Allowed fraction of gaps in each MSA position.
- **[WINDOW_SIZE]**: Size of the moving average window.
- **[CHROMOSOME]**: Chromosome ID for the resulting .wig file.
- **[WIG_START_INDEX]**: Starting index for the resulting .wig file.



