import sys

matrix_file_tsv = sys.argv[1]
tissues_file_tsv = sys.argv[2]
output_matrix_tsv = sys.argv[3]
output_tissues_csv = sys.argv[4]

tissues = {}
with open(tissues_file_tsv) as f, open(output_tissues_csv, "w") as out:
    # Skip the header
    next(f)
    for line in f:
        taxa_name, tissues_csv = line.strip().split("\t")
        tissues[taxa_name] = tissues_csv.split(",")
        if len(tissues[taxa_name]) > 1:
            # Write out the expanded tissues to the output file
            i=0
            for tissue in tissues[taxa_name]:
                out.write(f"{taxa_name}_{i},{tissue}\n")
                i += 1
        else:
            out.write(f"{taxa_name},{tissues[taxa_name][0]}\n")

with open(matrix_file_tsv) as f, open(output_matrix_tsv, "w") as out:
    header = f.readline().strip().split("\t")
    out.write("\t".join(header) + "\n")
    for line in f:
        taxa_name, *values = line.strip().split("\t")
        if taxa_name in tissues:
            if len(tissues[taxa_name]) > 1:
                for i, tissue in enumerate(tissues[taxa_name]):
                    out.write(f"{taxa_name}_{i}\t" + "\t".join(values) + "\n")
            else:
                out.write(f"{taxa_name}\t" + "\t".join(values) + "\n")
        else:
            # If the taxa name is not found in the tissues file, give an error
            sys.stderr.write(f"Error: {taxa_name} not found in tissues file\n")

