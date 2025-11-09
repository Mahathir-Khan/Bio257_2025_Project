import os

input_directory = "/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data/trimmed_fastq_for_RNAseq/trimming_reports"

txt_files = [f for f in os.listdir(input_directory) if f.endswith(".txt")]

headers = [
    "Total initial reads",
    "Surviving reads",
    "Discarded reads",
    "Reads with adapters:",
    "Reads written:",
    "Total basepairs processed:",
    "Quality-trimmed:",
    "Total written(filtered):"
]
output_file = "/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data/trimmed_fastq_for_RNAseq/trimming_reports/trimming_su$
with open(output_file, "w") as out:
    if not txt_files:
        print("no fastq files in current directory")
    else:
        header_row = ["File Name"] + headers
        print("\t".join(header_row))
        out.write("\t".join(header_row) + "\n")

        for file in txt_files:
            print("Generating report on " + str(file))
            data = {}
            try:
                with open(os.path.join(input_directory, file), 'r') as f:
                    for line in f:
                        for key in headers:
                            if line.strip().startswith(key):
                                parts = line.strip().split(key, 1)
                                if len(parts) > 1:
                                    value = parts[1].strip()
                                    data[key] = value
                                    break
            except FileNotFoundError:
                print(str(file) + " not found")
                continue

            row = [file] + [data.get(key, "") for key in headers]
            print("\t".join(row))
            out.write("\t".join(row) + "\n")

      
