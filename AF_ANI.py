import os

import pandas as pd

query_CDS_file = "1.Arcanobacterium_hemolyticum_DSM_20595_GCF_000092365.1_ASM9236v1_cds_from_genomic.fna"
query_Protein_file = "1.Arcanobacterium_hemolyticum_DSM_20595_GCF_000092365.1_ASM9236v1_protein.faa"

# Get all files inside directory
CDS_sequences_dir = os.listdir('Fasta_files/CDS/')
CDS_sequences_dir.remove(query_CDS_file)

# Get all files inside directory
Protein_sequences_dir = os.listdir('Fasta_files/Protein/')
Protein_sequences_dir.remove(query_Protein_file)

Blastn_results_path = 'Blastn_results/'
Blastn_results_dir = os.listdir(Blastn_results_path)
Blastp_results_path = 'Blastp_results/'
Blastp_results_dir = os.listdir(Blastp_results_path)

blast_headers_file = 'Fasta_files/blast_headers.txt'

# Prepare dataframe
columns = ["AAI", "AF", "ANI"]
rows = ["1vs2", "1vs3", "1vs4", "1vs5"]
blast_df = pd.DataFrame([[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                        columns=columns, index=rows)
print(blast_df)


def run_blastn():
    for CDS_sequence in CDS_sequences_dir:
        print("Starting Blastn")
        command = f"blastn -query Fasta_files/CDS/{query_CDS_file} " \
                  f"-subject Fasta_files/CDS/{CDS_sequence} " \
                  f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart" \
                  f" qend sstart send evalue bitscore qlen qcovs slen'" \
                  f" -out {Blastn_results_path}1vs{CDS_sequence[0]}.txt"
        os.system(command)
        print("End Blastn")


def run_blastp():
    for Prot_sequence in Protein_sequences_dir:
        print("Starting Blastp")
        command = f"blastp -query Fasta_files/Protein/{query_Protein_file} " \
                  f"-subject Fasta_files/Protein/{Prot_sequence} " \
                  f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart" \
                  f" qend sstart send evalue bitscore qlen qcovs slen'" \
                  f" -out {Blastp_results_path}1vs{Prot_sequence[0]}.txt"
        os.system(command)
        print("End Blastp")


def add_headers():
    for blastn_result in Blastn_results_dir:
        os.system(f"cat {blast_headers_file} Blastn_results/{blastn_result} "
                  f"> Blastn_results/headers_{blastn_result}")
        os.system(f"mv Blastn_results/headers_{blastn_result} Blastn_results/{blastn_result}")
    for blastp_result in Blastp_results_dir:
        os.system(f"cat {blast_headers_file} Blastp_results/{blastp_result} "
                  f"> Blastp_results/headers_{blastp_result}")
        os.system(f"mv Blastp_results/headers_{blastp_result} Blastp_results/{blastp_result}")


def calculate_statistics():
    for blastn_result in Blastn_results_dir:
        file_name, _ = os.path.splitext(blastn_result)
        blastn_result_data, blastn_result_BBH = get_bbh(Blastn_results_path,
                                                        blastn_result, 70, 70)
        # Calculate AF
        AF = blastn_result_BBH["query_length"].sum() / blastn_result_data["query_length"].sum()
        blast_df.loc[f"{file_name}", "AF"] = round(AF, 2)
        # Calculate ANI
        blastn_result_BBH["BBH_product"] = (blastn_result_BBH["aln_length"] / 100) \
                                           * blastn_result_BBH["pct_identity"]
        ANI = blastn_result_BBH["BBH_product"].sum() / blastn_result_BBH["query_length"].sum()
        blast_df.loc[f"{file_name}", "ANI"] = round(ANI, 2)
    for blastp_result in Blastp_results_dir:
        file_name, _ = os.path.splitext(blastp_result)
        blastn_result_data, blastp_result_BBH = get_bbh(Blastp_results_path,
                                                        blastp_result, 30, 80)
        # Calculate AAI
        blastp_result_BBH["BBH_product"] = (blastp_result_BBH["aln_length"] / 100) \
                                           * blastp_result_BBH["pct_identity"]
        AAI = blastp_result_BBH["BBH_product"].sum() / blastp_result_BBH["query_length"].sum()
        blast_df.loc[f"{file_name}", "AAI"] = round(AAI, 2)
    print(blast_df)
    blast_df.to_csv("AF_ANI_AAI.tsv", header=True, index=True, sep="\t")


def get_bbh(directory, blast_result, identity, coverage):
    # Read file
    blast_result_data = pd.read_csv(f'{directory}{blast_result}', delimiter="\t")
    # Calculate subject coverage
    print("Calculating subject coverage")
    blast_result_data["subject_coverage"] \
        = round((blast_result_data["aln_length"] / blast_result_data["subject_length"]) * 100)
    # Filter by identity
    print(f"Filtering by identity >= {identity}%")
    id_sup = blast_result_data["pct_identity"] >= identity
    blastn_result_BBH = blast_result_data[id_sup]
    # Filter by coverage
    print(f"Dropping rows not covered by at least {coverage}% for smallest gene")
    for index, row in blastn_result_BBH.iterrows():
        if query_not_covered(row, coverage) \
                or subject_not_covered(row, coverage) \
                or both_not_covered(row, coverage):
            print(row["subject_length"], row["subject_coverage"],
                  row["query_length"], row["query_coverage"])
            print("dropping row")
            blastn_result_BBH.drop(index, inplace=True)
    return blast_result_data, blastn_result_BBH


def query_not_covered(row, coverage):
    return row["query_length"] < row["subject_length"] and row["query_coverage"] < coverage


def subject_not_covered(row, coverage):
    return row["query_length"] > row["subject_length"] and row["subject_coverage"] < coverage


def both_not_covered(row, coverage):
    return row["query_length"] == row["subject_length"] \
           and row["subject_coverage"] < coverage \
           and row["query_coverage"] < coverage


if __name__ == "__main__":
    # run_blastn()
    # run_blastp()
    # add_headers()
    calculate_statistics()
