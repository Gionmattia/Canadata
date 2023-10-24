import argparse
import os
import pandas as pd


def input_files_list(path, keyword):
    '''
    Reads all the files in a input directory and creates a list of their names
    '''

    folder_files = os.listdir(path)
    file_list = []
    for each in folder_files:
        # Files with "missing_data" in their names are created by the "tximport_working.R" script,
        # but they are of no use here
        if "missing_data" in each:
            continue
        # Only include files that contain the keyword
        elif keyword.lower() not in each.lower():
            continue
        else:
            bits = each.split("/")
            name = bits[-1]
            file_list.append(name)
    file_list.sort()

    return file_list

def tsv_merger(file_paths):
    '''
    Extract gene counts from .tsv files and merges them into a single dataframe based on the "Name" column.
    '''
    try:
        # Initialize the base DataFrame with the first file in the list
        print(f"Processing {file_paths[0]}.")
        base_df = pd.read_csv(file_paths[0], sep="\t")
        base_df.insert(0, 'Name', base_df.index) 
        #base_df["Name"] = base_df.index

        # Loop through the remaining files and merge them into the base DataFrame
        for file_path in file_paths[1:]:
            temp_df = pd.read_csv(file_path, sep="\t")
            if "Name" not in temp_df.columns:
                print(f"Processing {file_path}.")
                temp_df["Name"] = temp_df.index
                base_df = pd.merge(base_df, temp_df, on="Name", how="outer")
                continue


        return base_df

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def tsv_saver(df, filename):
    '''
    Saves the merged df in the current directory
    '''

    full_path = f"{filename}.tsv"
    df.to_csv(full_path, sep="\t", index=False)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="merges the gene_count.tsv files into a single tsv for Anota2Seq analysis")
    parser.add_argument("input_dir", type = str, help = "Path to the directory containing the input files")
    parser.add_argument("keyword", type=str, default='', help="Keyword to filter files to be merged")
    parser.add_argument("output_name", type = str, help = "Output file name. The extension will be .tsv by default")

    args = parser.parse_args()

    files = input_files_list(args.input_dir, args.keyword)
    full_file_paths = [os.path.join(args.input_dir, file_name) for file_name in files]
    merged_tsv = tsv_merger(full_file_paths)
    tsv_saver(merged_tsv, args.output_name)

# NEED TO TEST THE BEHAVIOUR OF THIS NEW FUNCTION
