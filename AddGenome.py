import sys, os, gzip, time, shutil
import subprocess, pickle
import multiprocessing
import urllib.request


#####################
#####   USAGE   #####
#####################

def command():
        printUsage = False
        for arg in sys.argv:
                if (arg == '-h') or (arg == '-help') or (arg == '--help'): printUsage = True
        if (len(sys.argv) < 2) or (printUsage):
                sys.stderr.write("\nAddGenome version 3.0.1\n")
                sys.stdout.write("adds a new genome to TargetRNA3 that can be searched for sRNA:target interactions\n\n")
                sys.stdout.write("Usage: python AddGenome.py -a assembly_accession [options]\n")
                sys.stdout.write("Usage: python AddGenome.py -g genome_dir [options]\n\n")
                sys.stdout.write("*****   Required arguments   *****\n\n")
                sys.stdout.write("\tEither the -a flag or the -g flag is required\n")
                sys.stdout.write("\t\t\t\tbut not both and not neither\n\n")
                sys.stdout.write("\t-a STRING\tAsssembly accession identifier for\n")
                sys.stdout.write("\t\t\t\ta genome assembly, e.g.,\n")
                sys.stdout.write("\t\t\t\tGCF_000005845.2 (for E. coli)\n")
                sys.stdout.write("\t\t\t\tGCF_000006765.1 (for P. aeruginosa)\n")
                sys.stdout.write("\t-g STRING\tPath to directory containing genome information\n")
                sys.stdout.write("\t\t\t\tincluding files genomic.fna, protein.faa,\n")
                sys.stdout.write("\t\t\t\trna_from_genomic.fna, and feature_table.txt\n")
                sys.stdout.write("\t\t\t\t(files may be gzipped or not)\n\n")
                sys.stdout.write("*****   Optional arguments   *****\n\n")
                sys.stdout.write("\t-n_threads INT\tnumber of threads\n")
                sys.stdout.write("\t\t\t\t(default is based on self-identification\n")
                sys.stdout.write("\t\t\t\tof number of processors)\n")
                sys.stdout.write("\t-h\t\tprint USAGE and DESCRIPTION, ignore all other flags\n")
                sys.stdout.write("\t-help\t\tprint USAGE and DESCRIPTION, ignore all other flags\n")
                sys.stdout.write("\n")
                sys.exit(1)


#########################
#####   VARIABLES   #####
#########################

# Download genome files
GENOME_LOCATION = 'Genomes/'
ASSEMBLY_ACCESSION = ''
ASSEMBLY_BACTERIA_FILE_LOCAL = 'DataFiles/assembly_summary.txt.bacteria.gz'
ASSEMBLY_ARCHAEA_FILE_LOCAL = 'DataFiles/assembly_summary.txt.archaea.gz'
ASSEMBLY_BACTERIA_FILE_CLOUD = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
ASSEMBLY_ARCHAEA_FILE_CLOUD = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'

# Closest relatives
DB_16S = 'DataFiles/16S.fna'
GENOME_TO_ID_FILENAME = 'DataFiles/genome_IDs.pickle'
RELATIVES_FILENAME = 'closest_relatives.txt'
BLAST_RESULTS = 250

# Target homologs
DB_FAA = 'DataFiles/combined.faa'
NUM_RELATIVES = 10
TARGET_GENOME_MAP_FILE = 'DataFiles/gene_genome.pickle'
FILTER_FILE = 'blast_filter.txt'
TARGET_HOMOLOGS_FILE = 'mRNA.homologs'

# Target accessibility
RNAPLFOLD_DIR = 'RNAplfold_results/'

# COMMAND LINE ARGUMENTS
GENOME_DIR = ''
NUM_THREADS = int(max(multiprocessing.cpu_count()/2, 1))


#########################
#####   FUNCTIONS   #####
#########################

def arguments():
        global GENOME_DIR, ASSEMBLY_ACCESSION, NUM_THREADS
        for i in range(1, len(sys.argv)):
                if (sys.argv[i] == '-a'): ASSEMBLY_ACCESSION = sys.argv[i+1]
                elif (sys.argv[i] == '-g'): GENOME_DIR = sys.argv[i+1]
                elif (sys.argv[i] == '-n_threads'): NUM_THREADS = int(sys.argv[i+1])
        if (not GENOME_DIR.endswith('/')): GENOME_DIR += '/'
        if (len(GENOME_DIR) <= 1) and (len(ASSEMBLY_ACCESSION) <= 1):
                error('As a command line argument, either the -a flag or the -g flag must be provided appropriately.')
        if (len(GENOME_DIR) > 1) and (not os.path.isdir(GENOME_DIR)): error('As a command line argument, the flag -g is expected to be followed by the path to a directory containing genome files')


def setParameters(gen_dir, a_accession, n_threads):
        global GENOME_DIR, ASSEMBLY_ACCESION, NUM_THREADS
        GENOME_DIR = gen_dir
        if (not GENOME_DIR.endswith('/')): GENOME_DIR += '/'
        ASSEMBLY_ACCESSION = a_accession
        NUM_THREADS = n_threads


def error(s):
        sys.stderr.write('\nError - ' + s + '\n\n')
        sys.exit(1)


# Helper function for downloading genome files
def readIn_FTP_path(IN_FILENAME):
        if (not os.path.exists(IN_FILENAME)): return ''
        if (IN_FILENAME.endswith('.gz')): in_file = gzip.open(IN_FILENAME, 'rt')
        else: in_file = open(IN_FILENAME, 'r')
        line = in_file.readline()  # Ignore header line
        line = in_file.readline()  # Ignore header line
        line = in_file.readline()
        while (line != ''):
                parse_line = line.split('\t')
                accession = parse_line[0]
                ftp_path = parse_line[19]
                if (accession == ASSEMBLY_ACCESSION): return ftp_path
                line = in_file.readline()
        in_file.close()
        return ''


# Helper function for downloading genome files
def get_FTP():
        ftp = readIn_FTP_path(ASSEMBLY_BACTERIA_FILE_LOCAL)  # Bacteria locally
        if (len(ftp) == 0): ftp = readIn_FTP_path(ASSEMBLY_ARCHAEA_FILE_LOCAL)  # Archaea locally
        if (len(ftp) == 0):  # Bacteria cloud
                TEMP_FILENAME = str(time.time()) + '.bacteria'
                try: urllib.request.urlretrieve(ASSEMBLY_BACTERIA_FILE_CLOUD, TEMP_FILENAME)
                except: pass
                ftp = readIn_FTP_path(TEMP_FILENAME)
                if (os.path.exists(TEMP_FILENAME)): os.remove(TEMP_FILENAME)
        if (len(ftp) == 0):  # Archaea cloud
                TEMP_FILENAME = str(time.time()) + '.archaea'
                try: urllib.request.urlretrieve(ASSEMBLY_ARCHAEA_FILE_CLOUD, TEMP_FILENAME)
                except: pass
                ftp = readIn_FTP_path(TEMP_FILENAME)
                if (os.path.exists(TEMP_FILENAME)): os.remove(TEMP_FILENAME)
        if (len(ftp) == 0):
                error('Unable to determine the location online at RefSeq with genome files for ' + ASSEMBLY_ACCESSION + '\n' + 'As an alternative, you can download the genome files manually and re-run this program with the flag -g to point to the directory containing the downloaded genome files.')
        return ftp


# Main function for downloading genome files
def read_in_genome_files(ftp):
        if (not os.path.exists(GENOME_LOCATION)): os.makedirs(GENOME_LOCATION)
        if (not os.path.exists(GENOME_LOCATION + ASSEMBLY_ACCESSION)): os.makedirs(GENOME_LOCATION + ASSEMBLY_ACCESSION)
        base_filename = ftp[(ftp.rfind('/')+1):]
        try:
                urllib.request.urlretrieve(ftp + '/' + base_filename + '_protein.faa.gz', GENOME_LOCATION + ASSEMBLY_ACCESSION + '/' + base_filename + '_protein.faa.gz')
                urllib.request.urlretrieve(ftp + '/' + base_filename + '_genomic.fna.gz', GENOME_LOCATION + ASSEMBLY_ACCESSION + '/' + base_filename + '_genomic.fna.gz')
                urllib.request.urlretrieve(ftp + '/' + base_filename + '_feature_table.txt.gz', GENOME_LOCATION + ASSEMBLY_ACCESSION + '/' + base_filename + '_feature_table.txt.gz')
                urllib.request.urlretrieve(ftp + '/' + base_filename + '_genomic.gff.gz', GENOME_LOCATION + ASSEMBLY_ACCESSION + '/' + base_filename + '_genomic.gff.gz')
                urllib.request.urlretrieve(ftp + '/' + base_filename + '_rna_from_genomic.fna.gz', GENOME_LOCATION + ASSEMBLY_ACCESSION + '/' + base_filename + '_rna_from_genomic.fna.gz')
        except: error('Unable to download the genome files for ' + ASSEMBLY_ACCESSION + '\n' + 'As an alternative, you can download the genome files manually and re-run this program with the flag -g to point to the directory containing the downloaded genome files.')


# Helper function for computing closest_relatives.txt
def map_genome_to_id():
        if (not os.path.exists(GENOME_TO_ID_FILENAME)):
                error('Could not read in file ' + GENOME_TO_ID_FILENAME)
        with gzip.open(GENOME_TO_ID_FILENAME, 'rb') as in_file:
                genome_IDs = pickle.load(in_file)

        genome_to_id, id_to_genome, accession_to_id, id_to_accession = {}, {}, {}, {}
        for genome, ID, accession in genome_IDs:
                genome_to_id[genome] = ID
                id_to_genome[ID] = genome
                accession_to_id[accession] = ID
                id_to_accession[ID] = accession
        return genome_to_id, id_to_genome, accession_to_id, id_to_accession


# Helper function for computing closest_relatives.txt
def get_16S_sequence(GENOME_DIR):
        filelist = os.listdir(GENOME_DIR)
        for f in filelist:
                if (f.endswith('_rna_from_genomic.fna.gz') or f.endswith('_rna_from_genomic.fna')):
                        if (f.endswith('.gz')): in_file = gzip.open(GENOME_DIR + f, 'rt')
                        else: in_file = open(GENOME_DIR + f, 'r')
                        line = in_file.readline()
                        while (line != '') and ('16S ribosomal RNA' not in line) and ('16S Ribosomal RNA' not in line) and ('ribosomal RNA-16S' not in line): line = in_file.readline()
                        seq = ''
                        line = in_file.readline().strip()
                        while (line != '') and (not line.startswith('>')):
                                seq += line
                                line = in_file.readline().strip()
                        in_file.close()
                        return seq
        error('Could not extract 16s ribosomal RNA sequence from file ' + GENOME_DIR + '*_rna_from_genomic.fna')
        return ''


# Main function for computing closest_relatives.txt
def determine_closest_relatives(GENOME_DIR):
        FILENAME = str(time.time())
        TEMP_FILENAME1 = FILENAME + '.tmp1'
        TEMP_FILENAME2 = FILENAME + '.tmp2'
        genome_to_id, id_to_genome, accession_to_id, id_to_accession = map_genome_to_id()

        # Get accession
        accession = '???'
        filelist = os.listdir(GENOME_DIR)
        for f in filelist:
                if ((f.endswith('_genomic.fna.gz')) or (f.endswith('_genomic.fna'))) and ('_rna_from_' not in f):
                        if (f.endswith('.gz')): in_file = gzip.open(GENOME_DIR + f, 'rt')
                        else: in_file = open(GENOME_DIR + f, 'r')
                        accession = in_file.readline().split()[0][1:]
                        in_file.close()

        # BLAST 16S sequence
        seq_16S = get_16S_sequence(GENOME_DIR)
        with open(TEMP_FILENAME1, 'w') as out_file:
                out_file.write('>16S' + '\n' + seq_16S + '\n')
        p = subprocess.run(['blastn', '-db', DB_16S, '-query', TEMP_FILENAME1, '-evalue', '0.01', '-max_target_seqs', str(BLAST_RESULTS+1), '-num_threads', str(NUM_THREADS), '-outfmt', '6 qseqid sseqid evalue bitscore', '-out', TEMP_FILENAME2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Read in BLAST results
        replicates = {}
        with open(TEMP_FILENAME2, 'r') as in_file:
                with open(GENOME_DIR + RELATIVES_FILENAME, 'w') as out_file:
                        line = in_file.readline().strip()
                        while (line != ''):
                                parse_line = line.split('\t')
                                ID, bitscore = parse_line[1], parse_line[3]
                                if (ID not in replicates) and ((accession not in accession_to_id) or (ID != accession_to_id[accession])):
                                        out_file.write(ID + '\t' + id_to_genome[ID] + '\t' + bitscore + '\n')
                                        replicates[ID] = True
                                line = in_file.readline().strip()

        if os.path.exists(TEMP_FILENAME1): os.remove(TEMP_FILENAME1)
        if os.path.exists(TEMP_FILENAME2): os.remove(TEMP_FILENAME2)


# Helper function for computing target homologs
def map_target_ID_to_genome():
        if (os.path.exists(TARGET_GENOME_MAP_FILE)):
                with gzip.open(TARGET_GENOME_MAP_FILE, 'rb') as in_file:
                        return pickle.load(in_file)
        else: error('Could not read in file ' + TARGET_GENOME_MAP_FILE)


# Helper function for computing target homologs
def get_gene_mappings(GENOME_DIR):
        filelist = os.listdir(GENOME_DIR)
        target_ID_to_name = {}  # Map mRNA ID to its name
        for f in filelist:
                if (f.endswith('_feature_table.txt.gz') or f.endswith('_feature_table.txt')):
                        if (f.endswith('.gz')): in_file = gzip.open(GENOME_DIR + f, 'rt')
                        else: in_file = open(GENOME_DIR + f, 'r')
                        line = in_file.readline()  # Ignore header line
                        line = in_file.readline()
                        while (line != ''):
                                parse_line = line.split('\t')
                                if (parse_line[0] == 'CDS'):
                                        accession, mRNA_ID, mRNA_name, mRNA_synonym = parse_line[6], parse_line[10], parse_line[14], parse_line[16]
                                        if (mRNA_name.strip() == ''): mRNA_name = mRNA_synonym
                                        mRNA_name = mRNA_name.replace('/', '_')
                                        target_ID_to_name[mRNA_ID] = mRNA_name
                                line = in_file.readline()
                        in_file.close()
        return target_ID_to_name


# Helper function for computing target homologs
def read_in_closest_relatives(GENOME_DIR):
        if (not os.path.exists(GENOME_DIR + RELATIVES_FILENAME)):
                error('Could not read in file ' + GENOME_DIR + RELATIVES_FILENAME)
        relatives = []
        with open(GENOME_DIR + RELATIVES_FILENAME, 'r') as in_file:
                line = in_file.readline()
                while (line != ''):
                        ID, d, bitscore = line.strip().split()
                        relatives.append((ID, d))
                        line = in_file.readline()
        relatives = relatives[:NUM_RELATIVES]  # Only keep closest relatives
        return relatives


# Helper function for computing target homologs
def unzip_protein_file(GENOME_DIR):
        protein_file = ''
        filelist = os.listdir(GENOME_DIR)
        for f in filelist:
                if (f.endswith('_protein.faa.gz')):
                        protein_file = f[:-3]
                        with gzip.open(GENOME_DIR + f, 'rb') as f_in:
                                with open(protein_file, 'wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                if (f.endswith('_protein.faa')):
                        protein_file = f
                        shutil.copy2(GENOME_DIR + f, protein_file)
        return protein_file


# Helper function for computing target homologs
# GET THE IDs OF ALL THE PROTEINS IN ALL THE CLOSE RELATIVES AND STORE THESE
# IN A FILE THAT CAN BE USED TO FILTER THE LARGE BLAST DATABASE.
def create_blast_filter_for_relatives(GENOME_DIR, target_ID_to_genome):
        genome_to_id, id_to_genome, accession_to_id, id_to_accession = map_genome_to_id()
        genome_ID_to_target_ID = {}
        for target_id in target_ID_to_genome:
                genome_ID = target_ID_to_genome[target_id]
                if (genome_ID not in genome_ID_to_target_ID): genome_ID_to_target_ID[genome_ID] = {}
                genome_ID_to_target_ID[genome_ID][target_id] = True

        # Output IDs of all proteins of all relatives to Blast filter file
        relatives = read_in_closest_relatives(GENOME_DIR)
        with open(FILTER_FILE, 'w') as out_file:
                for ID, d in relatives:
                        accession = id_to_accession[ID]
                        targets = genome_ID_to_target_ID[accession]
                        for t in targets: out_file.write(t + '\n')

        # Format Blast filter file
        p = subprocess.run(['blastdb_aliastool', '-seqid_file_in', FILTER_FILE], stdout=subprocess.PIPE, stderr=subprocess.PIPE)


# Helper function for computing target homologs
# PERFORM BLAST WITHOUT USING MULTIPLE PROCESSORS
def blast_targets_against_database_SINGLE_PROCESS(protein_file):
        BLAST_OUTPUT_FILE = str(time.time()) + '.blast'
        p = subprocess.run(['blastp', '-db', DB_FAA, '-query', protein_file, '-outfmt', '6 qseqid sseqid evalue bitscore', '-out', BLAST_OUTPUT_FILE, '-evalue', '0.01', '-max_target_seqs', '100', '-num_threads', str(NUM_THREADS), '-seqidlist', FILTER_FILE + '.bsl'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if (os.path.exists(protein_file)): os.remove(protein_file)
        if (os.path.exists(FILTER_FILE)): os.remove(FILTER_FILE)
        if (os.path.exists(FILTER_FILE + '.bsl')): os.remove(FILTER_FILE + '.bsl')
        return BLAST_OUTPUT_FILE


# Helper function for computing target homologs
def run_blast(target_subset):
        accession_file = target_subset[0][1] + '.faa'
        with open(accession_file, 'w') as out_file:
                for target in target_subset:
                        GENOME_DIR, accession, protein_seq = target
                        out_file.write('>' + accession + '\n' + protein_seq + '\n')
        p = subprocess.run(['blastp', '-db', DB_FAA, '-query', accession_file, '-outfmt', '6 qseqid sseqid evalue bitscore', '-out', GENOME_DIR + accession + '.blast.xyz', '-evalue', '0.01', '-max_target_seqs', '100', '-num_threads', str(NUM_THREADS), '-seqidlist', FILTER_FILE + '.bsl'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if (os.path.exists(accession_file)): os.remove(accession_file)


# Helper function for computing target homologs
def blast_MULTI_PROCESS(targets):
        with multiprocessing.Pool() as pool: pool.map(run_blast, targets)


# Helper function for computing target homologs
# PERFORM BLAST POSSIBLY USING MULTIPLE PROCESSORS
def blast_targets_against_database(GENOME_DIR, protein_file):
        if (NUM_THREADS == 1): return blast_targets_against_database_SINGLE_PROCESS(protein_file)
        BLAST_OUTPUT_FILE = str(time.time()) + '.blast'

        # Get list of targets
        NUM_BATCHES = 5
        targets = []
        for index in range(NUM_BATCHES): targets.append([])
        index = 0
        with open(protein_file, 'r') as in_file:
                line = in_file.readline()
                while (line != ''):
                        if (line.startswith('>')):
                                accession = line.split()[0][1:]
                                seq = ''
                                line = in_file.readline()
                                while (line != '') and (not line.startswith('>')):
                                        seq += line.strip()
                                        line = in_file.readline()
                                targets[index].append((GENOME_DIR, accession, seq))
                                index = (index + 1) % NUM_BATCHES
        blast_MULTI_PROCESS(targets)
        if (os.path.exists(protein_file)): os.remove(protein_file)
        if (os.path.exists(FILTER_FILE)): os.remove(FILTER_FILE)
        if (os.path.exists(FILTER_FILE + '.bsl')): os.remove(FILTER_FILE + '.bsl')

        # Combine blast output files
        with open(BLAST_OUTPUT_FILE, 'wb') as out_file:
                filelist = os.listdir(GENOME_DIR)
                for f in filelist:
                        if (f.endswith('.blast.xyz')):
                                with open(GENOME_DIR + f, 'rb') as in_file: shutil.copyfileobj(in_file, out_file)
                                os.remove(GENOME_DIR + f)
        return BLAST_OUTPUT_FILE


# Main function for computing target homologs
def create_target_homolog_file(GENOME_DIR, BLAST_OUTPUT_FILE, target_ID_to_genome):
        # Read in BLAST results
        mRNA_homologs = {}
        target_ID_to_name = get_gene_mappings(GENOME_DIR)
        with open(BLAST_OUTPUT_FILE, 'r') as in_file:
                line = in_file.readline().strip()
                while (line != ''):
                        mRNA_ID, homolog, evalue, bitscore = line.split()
                        if (homolog.startswith('ref|')): homolog = homolog[4:]
                        if (homolog.endswith('|')): homolog = homolog[:-1]
                        mRNA_name = target_ID_to_name[mRNA_ID]
                        homolog_ID = target_ID_to_genome[homolog]
                        if (mRNA_name not in mRNA_homologs): mRNA_homologs[mRNA_name] = {}
                        if (len(mRNA_homologs[mRNA_name]) < NUM_RELATIVES):
                                mRNA_homologs[mRNA_name][homolog_ID] = True
                        line = in_file.readline().strip()
        with gzip.open(GENOME_DIR + TARGET_HOMOLOGS_FILE, 'wb') as out_file:
                pickle.dump(mRNA_homologs, out_file)
        if (os.path.exists(BLAST_OUTPUT_FILE)): os.remove(BLAST_OUTPUT_FILE)


# Helper function for computing target accessibility
def readInGenes(GENOME_DIR):
        UPSTREAM, DOWNSTREAM = 200, 100
        genome, genes = {}, []
        filelist = os.listdir(GENOME_DIR)
        for f in filelist:
                # Read in replicons
                if (f.endswith('_genomic.fna.gz') or (f.endswith('_genomic.fna'))) and ('_rna_from_' not in f):
                        if (f.endswith('.gz')): in_file = gzip.open(GENOME_DIR + f, 'rt')
                        else: in_file = open(GENOME_DIR + f, 'r')
                        line = in_file.readline().strip()
                        while (line != ''):
                                if (line.startswith('>')):
                                        ID = line.split()[0][1:]
                                        replicon = ''
                                line = in_file.readline().strip()
                                while (line != '') and (not line.startswith('>')):
                                        replicon += line
                                        line = in_file.readline().strip()
                                genome[ID] = replicon.upper()
                        in_file.close()

                # Read in genes
                if (f.endswith('_feature_table.txt.gz') or f.endswith('_feature_table.txt')):
                        gene_types = {'CDS', 'rRNA', 'tRNA', 'ncRNA'}
                        if (f.endswith('.gz')): in_file = gzip.open(GENOME_DIR + f, 'rt')
                        else: in_file = open(GENOME_DIR + f, 'r')
                        line = in_file.readline()  # Ignore header line
                        line = in_file.readline()
                        while (line != ''):
                                parse_line = line.split('\t')
                                gene_type = parse_line[0]
                                if (gene_type in gene_types):
                                        accession, start, stop, strand = parse_line[6:10]
                                        name, synonym = parse_line[14], parse_line[16]
                                        genes.append((gene_type, accession, name, synonym, start, stop, strand))
                                line = in_file.readline()
                        in_file.close()

        # Get target info needed by RNAplfold
        targets = []
        count = 1
        for g in genes:
                gene_type, accession, mRNA_name, mRNA_synonym, start, stop, strand = g[0], g[1], g[2], g[3], int(g[4]), int(g[5]), g[6]
                if (mRNA_name.strip() == ''): mRNA_name = mRNA_synonym
                mRNA_name = mRNA_name.replace('/', '_')
                if (gene_type == 'CDS'):
                        if (strand == '+'): mRNA_sequence = genome[accession][max(start-UPSTREAM-1,0):min(start+DOWNSTREAM-1,len(genome[accession]))]
                        else: mRNA_sequence = reverse(complement(genome[accession][max(stop-DOWNSTREAM,0):min(stop+UPSTREAM,len(genome[accession]))]))
                        targets.append((accession, mRNA_name, mRNA_sequence, str(count)))
                        count += 1
        return targets


# Helper function for computing target accessibility
def reverse(s):
        return s[::-1]


# Helper function for computing target accessibility
def complement(s):
        s = s.replace("C", "?")
        s = s.replace("G", "C")
        s = s.replace("?", "G")
        s = s.replace("A", "?")
        s = s.replace("T", "A")
        s = s.replace("?", "T")
        return s


# Helper function for computing target accessibility
def renameAndRemoveOutputFiles(prefix, name):
        os.rename(prefix + '_openen', GENOME_DIR + RNAPLFOLD_DIR + name + '_openen')
        if (os.path.exists(prefix + '_dp.ps')): os.remove(prefix + '_dp.ps')
        if (os.path.exists(prefix + '_uplex')): os.remove(prefix + '_uplex')


# Helper function for computing target accessibility
def run_RNAplfold(target):
        accession, mRNA_name, mRNA_sequence, count = target
        with open(GENOME_DIR + RNAPLFOLD_DIR + accession + '____' + mRNA_name + '.fa', 'w') as out_file:
                out_file.write('>' + accession + '____' + mRNA_name + '\n' + mRNA_sequence + '\n')
        p = subprocess.run(['RNAplfold', '-u', '40', '-O', '--plex_output', '--auto-id', '--id-prefix', 'TR' + count], input=mRNA_sequence.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if (p.returncode != 0) or (len(p.stderr.decode()) > 0):
                error('Problem executing RNAplfold:\t' + str(p.stderr.decode()) + '\n')
        renameAndRemoveOutputFiles('TR' + count + '_0001', accession + '____' + mRNA_name)


# Helper function for computing target accessibility
def RNAplfold_SINGLE_PROCESS(targets):
        for t in targets: run_RNAplfold(t)


# Helper function for computing target accessibility
def RNAplfold_MULTI_PROCESS(targets):
        with multiprocessing.Pool() as pool: pool.map(run_RNAplfold, targets)


# Helper function for computing target accessibility
# CREATE BINARY VERSIONS OF ACCESSIBILITY FILES USING RNAplex
def create_binary_files():
        p = subprocess.run(['RNAplex', '-a', GENOME_DIR + RNAPLFOLD_DIR, '-k'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        filelist = os.listdir(GENOME_DIR + RNAPLFOLD_DIR)
        for f in filelist:
                if (f.endswith('_openen')): os.remove(GENOME_DIR + RNAPLFOLD_DIR + f)


# Main function for computing target accessibility
def determine_target_accessibility(targets):
        if (NUM_THREADS == 1): RNAplfold_SINGLE_PROCESS(targets)
        else: RNAplfold_MULTI_PROCESS(targets)
        create_binary_files()


##############################
##########   MAIN   ##########
##############################

if __name__ == "__main__":
        command()
        arguments()

        # Download genome files
        if (len(ASSEMBLY_ACCESSION) > 1):
                sys.stdout.write('Downloading genome files...\n')
                ftp = get_FTP()
                read_in_genome_files(ftp)
                GENOME_DIR = GENOME_LOCATION + ASSEMBLY_ACCESSION
                if (not GENOME_DIR.endswith('/')): GENOME_DIR += '/'

        # Determine closest relatives via 16S
        sys.stdout.write('Computing closest relatives...\n')
        determine_closest_relatives(GENOME_DIR)

        # Determine target homologs
        sys.stdout.write('Identifying homologs...\n')
        target_ID_to_genome = map_target_ID_to_genome()
        create_blast_filter_for_relatives(GENOME_DIR, target_ID_to_genome)
        protein_file = unzip_protein_file(GENOME_DIR)
        BLAST_OUTPUT_FILE = blast_targets_against_database(GENOME_DIR, protein_file)
        create_target_homolog_file(GENOME_DIR, BLAST_OUTPUT_FILE, target_ID_to_genome)

        # Determine target structural accessibility
        sys.stdout.write('Calculating accessibility of RNA structures...\n')
        if (not os.path.exists(GENOME_DIR + RNAPLFOLD_DIR)): os.makedirs(GENOME_DIR + RNAPLFOLD_DIR)
        targets = readInGenes(GENOME_DIR)
        determine_target_accessibility(targets)
        sys.stdout.write('DONE.' + '\t\t' + 'New files have been generated in '  + GENOME_DIR + '\n\n')

