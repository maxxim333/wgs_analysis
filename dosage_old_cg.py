#### PART 1: Libraries and packages ####
import subprocess
from matplotlib.lines import Line2D
import math
import random
import datetime

# Use pip to install the package if they are not installed already
def install_program(package_name):
    try:
        subprocess.check_call(["/Library/Frameworks/Python.framework/Versions/3.10/bin/pip3.10", "install", package_name])
        print(f"Successfully installed {package_name}")
    except subprocess.CalledProcessError:
        print(f"Failed to install {package_name}")


package_name="pysam"
install_program(package_name)

package_name="natsort"
install_program(package_name)

package_name="mplcursors"
install_program(package_name)

import pysam
from natsort import natsorted
import mplcursors
import matplotlib.pyplot as plt

current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

#### PART 2: Extracting data ####

### Format explanation: This is how my VCF file is organized (header and first line as example):
### #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NBX2Y382
### 1	125271	rs377100675	C	T	.	PASS	.	GT:RC:AC:GP:DS	1/1:0:1:6.34842e-06,0.0858881,0.914106:1.9141
### The 1/1 means that the variant is homozygous for the alternative allele and 0:1 means that there was zero reads of the reference allele and 1 of the alternative
### I filtered the data to have only the 1/1 or 1/1 variants. This is my input file

vcf_file = "/Users/maxxim333/Desktop/projects/haplogroupsY/filtered_only_alternate.vcf"  # Replace "your_vcf_file.vcf" with the path to your VCF file
vcf = pysam.VariantFile(vcf_file)

# Initialize lists to store chromosome numbers, positions, and p-values
chromosomes = []
positions = []
dosages = []
rs = []

# Iterate over variants in the VCF. Extract the info on chromosomes, positions and ID of the variants
for variant in vcf:
    chromosomes.append(variant.chrom)
    positions.append(variant.pos)
    rs.append(variant.id)
vcf.close()

# Open the VCF file again to extract the information of dosage (hetero or homozygous)
with open(vcf_file, 'r') as f:
    # Iterate over each line in the file
    for line in f:
        # Skip header lines
        if line.startswith('#'):
            continue

        # Split the line by tab to extract fields
        fields = line.strip().split('\t')

        # Extract the "NBX2Y382" column value
        nbx2y382_column = fields[-1]  # Assuming "NBX2Y382" column is the last column


        # Extract the last number (dosage) from the "NBX2Y382" column
        dosage = float(nbx2y382_column.split(':')[-1])

        # Append the dosage to the list
        dosages.append(dosage)


### Here, I load the Human Genome Reference assembly 37 file to extract clinical significance and number of scientific papers mentioning the variant (evidence)
vcf_file2 = "/Users/maxxim333/Desktop/projects/haplogroupsY/clin_var_1by1/ass37_all.vcf"  # Replace "your_vcf_file.vcf" with the path to your VCF file
clinical_sign = ['not_found'] * len(rs)  # Initialize clinical_sign array with empty strings
evidence = [0] * len(rs)  # Initialize clinical_sign array with empty strings


total_lines = sum(1 for line in open(vcf_file2))  # Count total lines in the file

processed_lines = 0  # Counter for processed lines

counter = 0

rs_test = ['not_found'] * len(rs)  # Initialize clinical_sign array with empty strings
position_test = ['not_found'] * len(rs)  # Initialize clinical_sign array with empty strings

### The following block of code is basically to perform a join operation with some specificities on all the information and create a document containing the info of variant ID, clinical significance, and the position on genome
### This block of code is very heavy but it stores the info on a separate file and can be commented out once its done once
'''
with open(vcf_file2, 'r') as f:

        # Iterate over each line in the file
    for line in f:
        counter = counter+1
        if counter < 200000000:
            processed_lines += 1
            # Calculate progress percentage
            progress_percent = (processed_lines / total_lines) * 100


            # Split the line by tab to extract fields
            fields = line.strip().split('\t')

            # Check if the fields list has at least 3 columns
            if len(fields) < 3:
                print("skiped")  # Skip the loop if less than 3 columns are found

            else:
                rs_column = fields[1]  # Column containing the string
                clinical_column = fields[3]  # Column containing clinical_sign value
                evidence_column = len(fields[-1].split(","))
                if rs_column in rs:  # Check if string is in rs and clinical_sign is not empty

                    index = rs.index(rs_column)  # Get index of string in rs array
                    clinical_sign[index] = clinical_column  # Assign clinical_sign value to corresponding index
                    evidence[index] = evidence_column

                    rs_test[index] = rs[index]
                    position_test[index] = positions[index]


                # Draw the progress bar
                bar_length = 50
                filled_length = int(bar_length * processed_lines // total_lines)
                bar = '=' * filled_length + '-' * (bar_length - filled_length)
                print(f'\rProgress: [{bar}] {progress_percent:.2f}%', end='', flush=True)




# Generate filename with current time and date
output_file = f"clinical_sign_{current_time}.txt"

# Write clinical_sign array to the file
with open(output_file, 'w') as outfile:
    for item, item2, item3, item4 in zip(clinical_sign, rs_test, position_test, evidence):
        outfile.write("%s," % item)
        outfile.write("%s," % str(item2))
        outfile.write("%s," % str(item3))
        outfile.write("%s\n" % str(item4))

print(f"Clinical sign array stored in file: {output_file}")
'''

### Open the file generated by the code above. The clinical_sign_2024-03-23_21-15-57.txt string can be changed to the most current file generated by the previous block of code like this: #input_file = "{}".format(output_file)  # Provide the filename here
input_file = "/Users/maxxim333/Desktop/projects/haplogroupsY/clinical_sign_2024-03-23_21-15-57.txt"


clinical_sign = []
evidence = []
rs_vcf2 = []

# Read the content of the file into the clinical_sign array
with open(input_file, 'r') as infile:
    for line in infile:
        clinical_sign.append(str(line.strip().split(",")[0]))
        rs_vcf2.append(line.strip().split(",")[1])
        evidence.append(line.strip().split(",")[3])
print("Clinical sign array loaded from file.")


# Assuming `chromosomes`, `positions`, `dosages`, `clinical_sign`, and `evidence` are defined arrays. Create genomic positions array with all the info needed
genomic_positions = []
for chrom, pos, sig, ev, rs, dose in zip(chromosomes, positions, clinical_sign, evidence, rs, dosages):
    genomic_positions.append((chrom, pos, sig, ev, rs, dose))

# Sort the combined genomic positions list
sorted_positions = natsorted(genomic_positions, key=lambda x: (x[0], x[1], x[2], x[3], x[4]))

# Extract sorted chromosomes and positions
sorted_chromosomes = [chrom for chrom, _, _, _, _, _ in sorted_positions]
sorted_positions_pos = [pos for _, pos, _, _, _, _ in sorted_positions]
sorted_significances = [sig for _, _, sig, _, _, _ in sorted_positions]
sorted_evidences = [ev for _, _, _, ev, _, _ in sorted_positions]
sorted_rs = [rs for _, _, _, _, rs, _ in sorted_positions]
sorted_dosage = [dose for _, _, _, _, _, dose in sorted_positions]


#### Part 3: Plotting ###
plt.figure(figsize=(10, 6))

# Define color based on pathogenicity
colors = ['red' if sign == 'pathogenic' else 'yellow' if sign== 'risk-factor' else 'green' if sign == 'protective' else 'grey' if sign == 'drug-response' else 'orange' if sign== 'likely-pathogenic' or sign == 'pathogenic-likely-pathogenic' else 'blue' for sign in sorted_significances]

# Define size of dots based on evidence
sizes = [(math.log((float(e))+2))*50 for e in sorted_evidences]
dosages_alpha = [0.15 if float(e) < 0.5 else 0.3 if float(e) < 1  else 0.85 if float(e) < 1.5 else 1 for e in sorted_dosage]

### Create Scatter plot
plt.scatter(sorted_positions_pos, sorted_chromosomes, marker='.', c=colors, s=sizes, alpha=dosages_alpha)

#Add labels
plt.xlabel('Genomic Position')
plt.ylabel('Chromosome')
plt.title('Genome-Wide Distribution of Variants')

# Create a dictionary to store legend labels and corresponding colors
legend_labels_colors = {
    'Pathogenic': 'red',
    'Risk Factor': 'yellow',
    'Protective': 'green',
    'Affects Drug Response': 'grey',
    'Likely Pathogenic': 'orange',
    'Probably Benign': 'blue'
    # Add more categories and colors as needed
}

# Create legend handles and labels based on legend_labels_colors dictionary
legend_handles = [Line2D([0], [0], marker='o', color='w', markersize=10, markerfacecolor=color) for color in legend_labels_colors.values()]
legend_labels = legend_labels_colors.keys()

# Add legend
plt.legend(legend_handles, legend_labels, title='Pathogenicity', loc='best')

# Define a function to display labels when hovering over points
def on_hover(sel):
    index = sel.target.index
    if index is not None:
        rs_label = sorted_rs[index]
        sel.annotation.set_text(f'RS ID: {rs_label}')

# Connect the hovering function to the scatter plot
mplcursors.cursor(hover=True).connect("add", on_hover)

plt.show()
plt.savefig("combined_positions{}".format(current_time))


#### Part 4: Graph of individual Chromosomes ####

# Get unique chromosomes
unique_chromosomes = sorted(set(chromosomes))
print(unique_chromosomes)

# Plotting individual graphs for each chromosome
for chrom in unique_chromosomes:
    # Filter positions and dosages for the current chromosome
    ## A little trick is employed to duplicate data if the value of dosage is above 1.5, meaning that the variant is in homozygosis
    if chrom != "X":
        sorted_positions_pos_ind = [pos if dose < 1.5 else [pos, pos] for ch, pos, _, _, _, dose in sorted_positions  if ch == chrom]
        sorted_significances_ind = [sig if dose < 1.5 else [sig, sig] for ch, _, sig, _, _, dose in sorted_positions if ch == chrom]
        sorted_evidences_ind = [ev if dose < 1.5 else [ev, ev] for ch, _, _, ev, _, dose in sorted_positions if ch == chrom]
        sorted_rs_ind = [rs if dose < 1.5 else [rs, rs] for ch, _, _, _, rs, dose in sorted_positions if ch == chrom]
        sorted_dosage_ind = [1 if dose < 1.5 else [1, 2] for ch, _, _, _, _, dose in sorted_positions if ch == chrom]

    else:
        sorted_positions_pos_ind = [pos for ch, pos, _, _, _, dose in sorted_positions  if ch == chrom]
        sorted_significances_ind = [sig  for ch, _, sig, _, _, dose in sorted_positions if ch == chrom]
        sorted_evidences_ind = [ev  for ch, _, _, ev, _, dose in sorted_positions if ch == chrom]
        sorted_rs_ind = [rs  for ch, _, _, _, rs, dose in sorted_positions if ch == chrom]
        sorted_dosage_ind = [1 for ch, _, _, _, _, dose in sorted_positions if ch == chrom]


    if chrom != "X":  #This basically flattens the array. Chr X functions a little bit differently. There is never homozygosis
        sorted_positions_pos_ind = [val if isinstance(val, int) else val for sublist in sorted_positions_pos_ind for val in
                             ([sublist] if isinstance(sublist, int) else sublist)]
        sorted_significances_ind = [val if isinstance(val, str) else val for sublist in sorted_significances_ind for val in
                             ([sublist] if isinstance(sublist, str) else sublist)]
        sorted_evidences_ind = [val if isinstance(val, (int, str)) else val for sublist in sorted_evidences_ind for val in
                         ([sublist] if isinstance(sublist, (int, str)) else sublist)]
        sorted_rs_ind = [val if isinstance(val, str) else val for sublist in sorted_rs_ind for val in
                         ([sublist] if isinstance(sublist, str) else sublist)]
        sorted_dosage_ind = [val if isinstance(val, int) else val for sublist in sorted_dosage_ind for val in
                             ([sublist] if isinstance(sublist, int) else sublist)]


    # Plotting
    plt.figure(figsize=(18, 4))

    # Define color based on pathogenicity
    colors = ['red' if sign == 'pathogenic' else 'yellow' if sign == 'risk-factor' else 'green' if sign == 'protective' else 'grey' if sign == 'drug-response' else 'orange' if sign == 'likely-pathogenic' or sign == 'pathogenic-likely-pathogenic' else 'blue' for sign in sorted_significances_ind]

    sizes = [(math.log((float(e)) + 2)) * 50 for e in sorted_evidences_ind]  # Multiply by 10 for better visualization
    sign_alpha = [1 if sign == 'pathogenic' else 0.5 if sign == 'risk-factor' else 0.004 if sign == 'benign' else 0.75 if sign == 'likely-pathogenic' or sign == 'pathogenic-likely-pathogenic' else 0.002 for sign in sorted_significances_ind]

    fake_ys = [0] * len(sorted_positions_pos_ind)

    plt.xlabel('Genomic Position')
    plt.ylabel('Dosage')

    plt.title(f'Dosage Distribution of Variants on Chromosome {chrom}')

    # Create a dictionary to store legend labels and corresponding colors
    legend_labels_colors = {
        'Pathogenic': 'red',
        'Risk Factor': 'yellow',
        'Protective': 'green',
        'Affects Drug Response': 'grey',
        'Likely Pathogenic': 'orange',
        'Probably Benign': 'blue'
        # Add more categories and colors as needed
    }

    # Create legend handles and labels based on legend_labels_colors dictionary
    legend_handles = [Line2D([0], [0], marker='o', color='w', markersize=6, markerfacecolor=color) for color in
                      legend_labels_colors.values()]
    legend_labels = legend_labels_colors.keys()

    # Add legend
    plt.legend(legend_handles, legend_labels, title='Pathogenicity', loc='right')


    # Define a function to display labels when hovering over points
    def on_hover(sel):
        index = sel.target.index
        if index is not None:
            rs_label = sorted_rs_ind[index]
            sel.annotation.set_text(f'RS ID: {rs_label}')
    mplcursors.cursor(hover=True).connect("add", on_hover)

    #Define graph scales
    plt.ylim(0.5, 2.5)
    plt.xlim(0, 350000000)
    plt.yticks([1, 2], ["Heterozygote", "Homozygote"])


    plt.scatter(sorted_positions_pos_ind, sorted_dosage_ind, marker='.', c=colors, s=sizes, alpha=sign_alpha)

    plt.show()
    plt.savefig("individual_chr{}_constantscale".format(chrom))