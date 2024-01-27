# Author: Daniell McCann
# Date: 12/23-??
# Gene annotation pipeline using Blastx and Uniprot
# code written in python that can be run through the command line w/ seqauence as argument

## Author notes:
# consider adding user feedback and prompts to inform user of progress on pipeline when reading files, running blast, or quering uniprot (try-except blocks?)
# conisder adding check for imported biopython before code executes
# consider adding check for gene sequence file existence before code executes

'''
#Install biopython using pip to access biopython libraries for Bio.Blast module
# if pip is not installed on computer, need to install pip manually through command line as follows:
```bash
   pip install biopython 
   ```
'''

### Python code: yayy ###

#Added the `argparse` module to parse command-line arguments. will nee
import argparse
import requests
# using biopython library module bio.blast to import necessary libraries
from Bio.Blast import NCBIWWW, NCBIXML

## Function to perform BLAST called run_blast
# create function called run_blast
# arguements are sequence and output_file
def run_blast(sequence, output_file):
    result_handle = NCBIWWW.qblast("blastx", "swissprot", sequence)
    # open file as write to write in the result_handle that ran blastx
    with open(output_file, "w") as out_handle:
        out_handle.write(result_handle.read())

## Function to parse BLAST results
# create function called parse_blast_results
# arguement is blast_file
def parse_blast_results(blast_file):
    # open blast_file as result_handle
    result_handle = open(blast_file)
    # parse the result_handle using NCVIXML as blast_records
    blast_records = NCBIXML.parse(result_handle)

    # Empty list to store gene IDs from BLAST to be used later to compare to inprot query
    gene_ids_from_blast = []

    # Loop through paresed blast_records and print relevant information and store gene_id's in emply list []
    for record in blast_records:
        for alignment in record.alignments:
            gene_id = alignment.hit_id
            # add gene_id's to empty list 
            gene_ids_from_blast.append(gene_id)
            # print relevant information 
            print(f"Gene ID: {gene_id}, E-value: {alignment.hsps[0].expect}")

    return gene_ids_from_blast

## Function to query UniProt for annotations
# create function called query_uniprot with gene_id as arguement
def query_uniprot(gene_id):
    # Construct the URL for the UniProt API using the provided gene_id from ^^
    url = f"https://www.uniprot.org/uniprot/{gene_id}.xml"
    # Send an HTTP GET request to the UniProt API
    response = requests.get(url)
    # Return the text content of the response (XML data)
    return response.text


if __name__ == "__main__":
    # Command-line argument parser
    # Create arg parse object w/ description
    parser = argparse.ArgumentParser(description="Gene Annotation Pipeline")
    # positional argument for gene sequence file for path file to be FASTA formatting
    parser.add_argument("gene_sequence", help="Path to the gene sequence in FASTA format")
    args = parser.parse_args()

    # Read gene sequence from the file
    # reads gene file conent as gene_sequence
    with open(args.gene_sequence, "r") as gene_file:
        gene_sequence = gene_file.read()

    # Run BLAST and save results to a file
    # calls fun_blast function (woohoo) with gene_sequence as input and outputting to file as blast_results.xml
    run_blast(gene_sequence, "blast_results.xml")

    # Parse BLAST results and get gene IDs
    # calls parse_blast_results function (woohoo) to extract gene_ids from blast_results.xml file
    gene_ids_from_blast = parse_blast_results("blast_results.xml")

    # Query UniProt for annotations using gene IDs from BLAST
    # for loop to iterate gene_ids from blast results and query uniprot calling query_uniprot function (woohoo)
    for gene_id in gene_ids_from_blast:
        # call query uni_prot function
        uniprot_xml = query_uniprot(gene_id)
        # print results!
        print(f"Gene Annotation for {gene_id}:\n{uniprot_xml}")
```


# the above script ^^ can be run from the command line by providing the path to the gene sequence file as an argument. For example:
  ```bash
  python script_name.py input_gene_sequence.fasta
  ```

# the above script ^^ can be run in RStudio using system/system2 providing the path to the gene sequence file as an argument. For example:
'''
system("python script_name.py input_gene_sequence.fasta")
'''
