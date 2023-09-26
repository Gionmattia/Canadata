import os
import argparse

# Arguments: path to folder containing files, extension of the files to modify.
# 
# NEED TO UPDATE THE SCRIPT TO THE NEW NAMES OF FILES
index_dict = {

    # Reps 1, one sample is missing

    "GTAGTCAG-CCTAACAG": "Rep1_VCctrl_input",
    "GTACCTTG-TACGGTCT": "Rep1_VCctrl_poly",

    "TGTAGCCA-CTCAGAAG": "Rep1_BPctrl_poly",

    "CAGTGAAG-GAGTGTGT": "Rep1_VCanoi_input",
    "GCGATAGT-ACGCAGTA": "Rep1_VCanoi_poly",

    "CACCTTAC-TAGTGCCA": "Rep1_BPanoi_input",
    "TTCTCGAC-ACCAAGCA": "Rep1_BPanoi_poly",

    # Reps 2

    "CCAATAGG-CTGGTCAT": "Rep2_VCctrl_input",
    "AGCCAAGT-ACAGTGAC": "Rep2_VCctrl_poly",

    "AAGCGCAT-CATGAGCA": "Rep2_BPctrl_input",
    "CCACTTCT-ACACGAGA": "Rep2_BPctrl_poly",

    "CTCGTCTT-ACGAACGA": "Rep2_VCanoi_input",
    "CAGAGTGT-GTCTGAGT": "Rep2_VCanoi_poly",

    "GATTCAGC-GCACACAA": "Rep2_BPanoi_input",
    "GAGTCTCT-CCTCATCT": "Rep2_BPanoi_poly",

    # Reps 3

    "ATCTTCGG-GTCTCATC": "Rep3_VCctrl_input",
    "TCGCGATA-ACCGCTAT": "Rep3_VCctrl_poly",

    "GCTACGTT-CGGTTGTT": "Rep3_BPctrl_input",
    "TCTTCTGC-AGACGCTA": "Rep3_BPctrl_poly",

    "AAGGCGTT-CGGAGTAT": "Rep3_VCanoi_input",
    "AAGTCCGT-CTGTGGTA": "Rep3_VCanoi_poly",

    "AAGCCACA-TGGACCAT": "Rep3_BPanoi_input",
    "CCTATGGT-CGGCATTA": "Rep3_BPanoi_poly",
    
    # Reps 4

    "ATGTAGCG-AAGGCTCT": "Rep4_VCctrl_input",
    "TCTAGCTG-GCCTTCTT": "Rep4_VCctrl_poly",

    "GTTAGACG-CCGCTTAA": "Rep4_BPctrl_input",
    "TTGCTGGA-ACCATGTC": "Rep4_BPctrl_poly",

    "GTCTTGCA-CACCATGA": "Rep4_VCanoi_input",
    "TCCTTAGC-TAGCAGGA": "Rep4_VCanoi_poly",

    "ACATCCTG-CTTGCTAG": "Rep4_BPanoi_input",
    "GTTCGGTT-TGCCTCAA": "Rep4_BPanoi_poly",
}

def rename_files (folder_path, extension, output_path):

    # Iterates through the files in the given path

    for filename in os.listdir(folder_path):

        # Checks if the file matches the specified extension

        if extension in filename:

            # Retrieves the nucleotide index (the index applied by the sequencing lab)

            index_start = 14
            index_end = 31
            X = filename[index_start:index_end]

            # Checks if the index is among the keys of the conversion dictionary (index_dict)

            try:
                index_dict[X]
            except KeyError:
                raise CustomException('NUCLEOTIDE INDEX NOT PRESENT')
            
            # Gets the "human readable name" for the file from the conversion dictionary and stores it in Y
        
            Y = index_dict.get(X)

            # Checks whether it's the technical replicate 1 or 2. Specific to these files.

            BOL = "4_1" in filename
            rep_number = "_1" if BOL else "_2"

            # Assembles the new name

            new_name = Y + "_" + X + "_" + rep_number + extension

            # Reconstructs the old and the new path
            
            old_path = os.path.join(folder_path, filename)
            new_path = os.path.join(output_path, new_name)
            command = "cp " + old_path + " " + new_path
            os.system(command)

            # THIS NEEDS TO BE CHANGED: Now it renames, I want it to save to a different location

            # os.rename(old_path, new_path)


# TO DO:
# Then it should be tested against some dummy files (txt files) -> TXT file should have a sequence inside to double-check name changed correctly.
# The above test should be run from the terminal
# If everything goes smoothly, then it can be run on the actual files on svr03.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="renames the Canadata files accordingly")
    parser.add_argument("folder_path", type = str, help = "Path to the folder containing the data")
    parser.add_argument("extension", type = str, help = "extension of the files to rename (ex. '.fastq.gz', '.html'...)")
    parser.add_argument("output_path", type = str, help = "outhput folder path")
    args = parser.parse_args()

    print("Folder currently searched ", args.folder_path)
    print("Extension currently selected ", args.extension)

    # AND NOW YOU CALL THE RENAME FUNCTION FROM HERE
    rename_files(args.folder_path, args.extension, args.output_path)

    #df = check_and_assign_header(df) 
    #save_with_new_name(df, args.Path)