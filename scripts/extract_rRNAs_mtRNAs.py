import argparse
import Bio
from Bio import SeqIO


def insert_newline(string):
    """
    Formats a string so that every 80 characters a new_line characters is inserted.
    """
    
    return "\n".join(string[i:i+80] for i in range(0, len(string), 80))



def rRNA_fasta_filter2(path, output_file):
    """
    Parses through a fasta file, extracting all the reads listed as rRNA or ribosomal RNA to a separate fasta (fna) file
    
    """
    
    with open(output_file, "w") as outfile:
        
        for seq_record in SeqIO.parse(path, "fasta"):
            line = seq_record.description
            words = line.split(",")
            
            # You can modify what to extract by changing the words listed here in these blocks!
    
            if "rRNA" in words[-1]:
        
                descriptor = "\n" + ">" + line
        
                sequence = str(seq_record.seq)
                sequence = insert_newline(sequence)
            
                read = descriptor + "\n" + sequence
                
                outfile.write(read)
        
            elif "ribosomal RNA" in words[-1]:

                descriptor = "\n" + ">" + line
                
                sequence = str(seq_record.seq)
                sequence = insert_newline(sequence)
                
                read = descriptor + "\n" + sequence
                
                outfile.write(read)

            elif "Mt_tRNA" in words[-1]:
                
                descriptor = "\n" + ">" + line
                
                sequence = str(seq_record.seq)
                sequence = insert_newline(sequence)
                
                read = descriptor + "\n" + sequence
                
                outfile.write(read)

                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="extracts all the reads listed as 'rRNA' or 'ribosomal RNA'")
    parser.add_argument("input", type = str, help = "Path to the fasta file")
    parser.add_argument("output", type = str, help = "Path to the desired outputfile. Can specicify a new file to be created, but always remember the extension")
    args = parser.parse_args()

    # AND NOW YOU CALL THE RENAME FUNCTION FROM HERE
    rRNA_fasta_filter2(args.input, args.output) 