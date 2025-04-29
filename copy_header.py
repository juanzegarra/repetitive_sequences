import os 

def seqheader(file, newfile):
    fixedlines = []
    current_header = None
    copy_number = 1
    first_sequence = False

    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                current_header = line
                copy_number = 1
                first_sequence = True
                fixedlines.append(f"{current_header}_copy_{copy_number}")

            else:
                if not first_sequence:
                    copy_number += 1
                    fixedlines.append(f"{current_header}_copy_{copy_number}")

                fixedlines.append(line)
                first_sequence = False

    with open(newfile, 'w') as f_out:
        for item in fixedlines:
            f_out.write(f"{item}\n")
def copy_header(arquivo):
    arquivo = input(str("Nomeie o arquivo fasta"))
    abspath = os.get.abspath(f"{arquivo}")
    working_path = os.getcwd()
    novo_arquivo = os.path.join(f"{working_path}", f"correct_{arquivo}\n")
    seqheader(f"{arquivo}", f"{novo_arquivo}")
copy_header()