import click
from pathlib import Path
from Bio import AlignIO
import pandas as pd
import math
import statistics


def msa_to_list(msa_alignio):
    msa = []
    for aligned_sequence in msa_alignio:
        msa.append([aligned_sequence.description] + list(str(aligned_sequence.seq)))
    return msa


def msa_list_to_df(msa_list):
    df = pd.DataFrame(msa_list).set_index(0, drop=True)
    return df


def msa_to_df(msa):
    return msa_list_to_df(msa_to_list(msa))


def drop_gaps(df, max_frac_gaps):
    max_gaps = round((max_frac_gaps) * (len(df.index)))
    df = df.loc[:, ((df == "-").sum(axis=0) <= max_gaps)]
    return df


def calculate_shannon_entropy(df):
    shannon_entropy_msa = []

    for col in df.columns:
        bases_list = df.loc[:, col].values.tolist()
        shannon_entropy_msa.append([col, shannon_entropy_position(bases_list)])

    return shannon_entropy_msa


def shannon_entropy_position(bases):
    unique_bases = list(set(bases))
    entropies = []

    for base in unique_bases:
        n_i = bases.count(base)
        P_i = n_i / float(len(bases))
        entropy_i = P_i * (math.log(P_i, 2))
        entropies.append(entropy_i)

    shannon_entropy = -(sum(entropies))

    return shannon_entropy


def rolling_window_average(entropy_list, window_size):
    entropy_df = pd.DataFrame(entropy_list)
    temp_list = entropy_df[1].values.tolist()

    for i in range(len(temp_list) - window_size + 1):
        window = temp_list[i : i + window_size]
        middle_index = int((len(window) - 1) / 2)
        window_mean = statistics.fmean(window)
        temp_list[middle_index + i] = window_mean

    entropy_df[1] = temp_list
    entropy_list = entropy_df.values.tolist()

    return entropy_list


def get_output_filename(input_filename, output_filename):
    if output_filename == "":
        output_filename = Path(input_filename).stem + ".wig"
    return output_filename


def write_shannon_entropy_to_wig_file(
    entropy_list, file_in, file_out, chromosome, start_index
):
    output_filename = get_output_filename(file_in, file_out)
    with open(output_filename, "w") as wig:
        wig.write("variableStep chrom=" + chromosome + "\n")
        for index, entropy in entropy_list:
            wig.write(f"{int(index + start_index)} {entropy}\n")


@click.command()
@click.argument("input_filename", type=click.Path(exists=True), required=True)
@click.argument(
    "filetype",
    type=click.Choice(
        ["clustal", "phylip", "stockholm", "fasta"], case_sensitive=False
    ),
    default="clustal",
)
@click.argument("output_filename", type=click.Path(), default="")
@click.argument("studied_genome_accession", type=str, required=False)
@click.argument("max_fraction_gaps", type=float, default=0.0)
@click.argument("window_size", type=int, default=1)
@click.argument("chromosome", type=str, required=False)
@click.argument("wig_start_index", type=int, required=False)
def output_entropy_wig(
    input_filename,
    filetype,
    output_filename,
    studied_genome_accession,
    max_fraction_gaps,
    window_size,
    chromosome,
    wig_start_index,
):
    # parse alignment
    alignment = AlignIO.read(input_filename, filetype)

    # create alignment dataframe
    msa_df = msa_to_df(alignment)

    # drop positions with a gap in a studied genome
    if studied_genome_accession:
        studied_genome_row = pd.DataFrame(msa_df.loc[studied_genome_accession, :]).T
        studied_genome_row = drop_gaps(studied_genome_row, 0.0)
        gapless_positions = studied_genome_row.columns
        msa_df = msa_df.loc[:, gapless_positions]

        # reindex columns
        msa_df.columns = range(1, len(msa_df.columns) + 1)

    # drop positions with more than max_fraction_gaps
    msa_df = drop_gaps(msa_df, max_fraction_gaps)

    # calculate Shannon entropy for each position of the MSA
    shannon_entropy_list = calculate_shannon_entropy(msa_df)

    if window_size > 1:
        assert window_size % 2 == 1, "Window size must be an odd number"

        # average Shannon entropy in rolling window
        shannon_entropy_list = rolling_window_average(shannon_entropy_list, window_size)

    write_shannon_entropy_to_wig_file(
        shannon_entropy_list,
        input_filename,
        output_filename,
        chromosome,
        wig_start_index,
    )


output_entropy_wig()
