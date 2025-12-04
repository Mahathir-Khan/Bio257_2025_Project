import os # import os library

input_directory = "/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data/trimmed_fastq_for_RNAseq/trimming_reports" # assign directory

txt_files = [f for f in os.listdir(input_directory) if f.endswith(".txt")] # create list of files
# assign headers for reads to collect in table
headers = [
    "Total initial reads", # total initial reads
    "Surviving reads", # surviving reads
    "Discarded reads", # discarded reads
    "Reads with adapters:", # reads with adapters
    "Reads written:", # reads written
    "Total basepairs processed:", # total basepairs processed
    "Quality-trimmed:", # quality trimmed
    "Total written(filtered):" # total written
]
# create output file
output_file = "/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data/trimmed_fastq_for_RNAseq/trimming_reports/trimming_su$
with open(output_file, "w") as out:
    if not txt_files: # check to make sure fastq files exist
        print("no fastq files in current directory")
    else: # if exists, continue
        header_row = ["File Name"] + headers # assign row for headers
        print("\t".join(header_row)) # print each header tab separated
        out.write("\t".join(header_row) + "\n") # write to output file 

        for file in txt_files: 
            print("Generating report on " + str(file)) # print statement to keep track of sample 
            data = {} # create dictionary for data
            try:
                with open(os.path.join(input_directory, file), 'r') as f: # try to open each file
                    for line in f: # iterate through file
                        for key in headers:  # select each key in headers
                            if line.strip().startswith(key): # strip, select line starting with kep
                                parts = line.strip().split(key, 1) # strip and split based on key
                                if len(parts) > 1:
                                    value = parts[1].strip() 
                                    data[key] = value # assign values to data
                                    break
            except FileNotFoundError:
                print(str(file) + " not found") # check for if file not found
                continue

            row = [file] + [data.get(key, "") for key in headers] # add each key to row
            print("\t".join(row)) # print rows
            out.write("\t".join(row) + "\n") # write to output file

      
