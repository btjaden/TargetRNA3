import sys, os, gzip, time, shutil
import subprocess, pickle
import multiprocessing
from functools import partial
import numpy as np
import scipy.stats as stats
import pandas as pd
import warnings; warnings.filterwarnings("ignore")


#####################
#####   USAGE   #####
#####################

def command():
        printUsage = False
        for arg in sys.argv:
                if (arg == '-h') or (arg == '-help') or (arg == '--help'): printUsage = True
        if (len(sys.argv) < 5) or (printUsage):
                sys.stderr.write("\nTargetRNA3 version 3.0.1\n")
                sys.stdout.write("TargetRNA3 predicts targets of a small regulatory RNA\n\n")
                sys.stdout.write("Usage: python TargetRNA3.py -s *.fa -g genome_dir [options]\n\n")
                sys.stdout.write("*****   Required arguments   *****\n\n")
                sys.stdout.write("\t-s STRING\tFile in FASTA format containing sRNA sequence\n")
                sys.stdout.write("\t-g STRING\tPath to directory containing genome information\n")
                sys.stdout.write("\t\t\t\tincluding files genomic.fna, portein.faa,\n")
                sys.stdout.write("\t\t\t\trna_from_genomic.fna, and feature_table.txt\n")
                sys.stdout.write("\t\t\t\t(files may be gzipped or not)\n\n")
                sys.stdout.write("*****   Optional arguments   *****\n\n")
                sys.stdout.write("\t-o STRING\tFile to which results should be output\n")
                sys.stdout.write("\t\t\t\t(default is standard out)\n")
                sys.stdout.write("\t-prob FLOAT\tProbability above which a target is predicted\n")
                sys.stdout.write("\t\t\t\tas a target\n")
                sys.stdout.write("\t\t\t\t(default is 0.5)\n")
                sys.stdout.write("\t-pval FLOAT\tP-value below which a target is predicted\n")
                sys.stdout.write("\t\t\t\tas a target\n")
                sys.stdout.write("\t\t\t\t(default is 0.05)\n")
                sys.stdout.write("\t-n_threads INT\tnumber of threads\n")
                sys.stdout.write("\t\t\t\t(default is based on self-identification\n")
                sys.stdout.write("\t\t\t\tof number of processors)\n")
                sys.stdout.write("\t-db STRING\tPath to BLAST database: DataFiles/combined.fna\n")
                sys.stdout.write("\t-model STRING\tPath to ML model file: DataFiles/model.pickle\n")
                sys.stdout.write("\t-h\t\tprint USAGE and DESCRIPTION, ignore all other flags\n")
                sys.stdout.write("\t-help\t\tprint USAGE and DESCRIPTION, ignore all other flags\n")
                sys.stdout.write("\n")
                sys.exit(1)


#########################
#####   VARIABLES   #####
#########################

MRNA_HOMOLOGS_FILE = 'mRNA.homologs'
RNAPLFOLD_DIR = 'RNAplfold_results/'

# COMMAND LINE ARGUMENTS
SRNA_FILENAME = ''
GENOME_DIR = ''
OUTPUT_FILE = ''
PROBABILITY = 0.5
PVALUE = 0.05
NUM_THREADS = int(max(multiprocessing.cpu_count()/2, 1))
DB = 'DataFiles/combined.fna'
MODEL_FILE = 'DataFiles/model.pickle'


#########################
#####   FUNCTIONS   #####
#########################

def arguments():
        global SRNA_FILENAME, GENOME_DIR, OUTPUT_FILE, PROBABILITY, PVALUE, NUM_THREADS, DB, MODEL_FILE
        for i in range(1, len(sys.argv)):
                if (sys.argv[i] == '-s'): SRNA_FILENAME = sys.argv[i+1]
                elif (sys.argv[i] == '-g'): GENOME_DIR = sys.argv[i+1]
                elif (sys.argv[i] == '-o'): OUTPUT_FILE = sys.argv[i+1]
                elif (sys.argv[i] == '-prob'): PROBABILITY = float(sys.argv[i+1])
                elif (sys.argv[i] == '-pval'): PVALUE = float(sys.argv[i+1])
                elif (sys.argv[i] == '-n_threads'): NUM_THREADS = int(sys.argv[i+1])
                elif (sys.argv[i] == '-DB'): DB = sys.argv[i+1]
                elif (sys.argv[i] == '-model'): MODEL_FILE = sys.argv[i+1]
        if (not GENOME_DIR.endswith('/')): GENOME_DIR += '/'
        if (len(GENOME_DIR) <= 1) or (not os.path.isdir(GENOME_DIR)): error('As a command line argument, the flag -g is expected to be followed by the path to a directory containing genome files')
        if (len(SRNA_FILENAME) < 1): error('As a command line argument, the flag -s is expected to be followed by a path to a file containing a sRNA sequence in FASTA format')


def setParameters(sRNA_file, gen_dir, out_filename, prob, pval, n_threads, db, model):
        global SRNA_FILENAME, GENOME_DIR, OUTPUT_FILE, PROBABILITY, PVALUE, NUM_THREADS, DB, MODEL_FILE
        SRNA_FILENAME = sRNA_file
        GENOME_DIR = gen_dir
        if (not GENOME_DIR.endswith('/')): GENOME_DIR += '/'
        OUTPUT_FILE = out_filename
        PROBABILITY = prob
        PVALUE = pval
        NUM_THREADS = n_threads
        DB = db
        MODEL_FILE = model


def error(s):
        sys.stderr.write('\nError - ' + s + '\n\n')
        sys.exit(1)


def read_in_sRNA():
        with open(SRNA_FILENAME, 'r') as in_file:
                line = in_file.readline()
                if (not line.startswith('>')): error(SRNA_FILENAME + ' not in FASTA format. Expecting first line of file to start with ">".')
                sRNA_name = line.split()[0][1:]
                sRNA_sequence = ''
                line = in_file.readline().strip()
                while (line != ''):
                        sRNA_sequence += line
                        line = in_file.readline().strip()
        sRNA_sequence = sRNA_sequence.upper().replace('U', 'T')
        return sRNA_name, sRNA_sequence


def read_in_genes(GENOME_DIR):
        if (not os.path.isdir(GENOME_DIR)): error(GENOME_DIR + ' is not a valid directory.\n')
        filelist = os.listdir(GENOME_DIR)
        UPSTREAM, DOWNSTREAM = 200, 100

        # Read in replicons
        genome = {}
        for f in filelist:
                if ((f.endswith('_genomic.fna.gz') or (f.endswith('_genomic.fna'))) and ('_rna_from_' not in f)):
                        if (f.endswith('.gz')): in_file = gzip.open(GENOME_DIR + f, 'rt')
                        else: in_file = open(GENOME_DIR + f, 'r')
                        line = in_file.readline().strip()
                        while (line != ''):
                                if (line.startswith('>')):
                                        accession = line.split()[0][1:]
                                        replicon = ''
                                line = in_file.readline().strip()
                                while (line != '') and (not line.startswith('>')):
                                        replicon += line
                                        line = in_file.readline().strip()
                                genome[accession] = replicon.upper()
                        in_file.close()
        if (len(genome) == 0): error('Could not read in genome file - *_genomic.fna.gz or *_genomic.fna')

        # Read in genes
        genes = []
        for f in filelist:
                if (f.endswith('_feature_table.txt.gz')) or (f.endswith('_feature_table.txt')):
                        gene_types = {'CDS', 'rRNA', 'tRNA', 'ncRNA'}
                        if (f.endswith('.gz')): in_file = gzip.open(GENOME_DIR + f, 'rt')
                        else: in_file = open(GENOME_DIR + f, 'r')
                        count = 0
                        line = in_file.readline()  # Ignore header line
                        line = in_file.readline()
                        while (line != ''):
                                parse_line = line.split('\t')
                                gene_type = parse_line[0]
                                if (gene_type in gene_types):
                                        assembly = parse_line[2]
                                        accession, start, stop, strand, ID = parse_line[6:11]
                                        product, name, synonym = parse_line[13], parse_line[14], parse_line[16]
                                        if (name.strip() == ''): name = synonym
                                        name = name.replace('/', '_')
                                        mRNA_sequence = ''
                                        if (gene_type == 'CDS'):
                                                if (strand == '+'): mRNA_sequence = genome[accession][max(int(start)-UPSTREAM-1,0):min(int(start)+DOWNSTREAM-1,len(genome[accession]))]
                                                else: mRNA_sequence = reverse(complement(genome[accession][max(int(stop)-DOWNSTREAM,0):min(int(stop)+UPSTREAM,len(genome[accession]))]))
                                                count += 1
                                        genes.append((gene_type, accession, name, synonym, start, stop, strand, product, mRNA_sequence, str(count), ID))
                                line = in_file.readline()
                        in_file.close()
        if (len(genes) == 0): error('Could not read in gene file - *_feature_table.txt.gz or *_feature_table.txt')
        return genome, genes


def reverse(s):
        return s[::-1]


def complement(s):
        s = s.replace("C", "?")
        s = s.replace("G", "C")
        s = s.replace("?", "G")
        s = s.replace("A", "?")
        s = s.replace("T", "A")
        s = s.replace("?", "T")
        return s


def getSeeds(sRNA_sequence, genes):
        SEED_LENGTH = 7
        seeds = {}
        sRNA_seq = reverse(complement(sRNA_sequence))

        # Create dictionary of sRNA seed sequences of given length
        sRNA_seeds = {}
        for i in range(len(sRNA_seq)-SEED_LENGTH+1):
                sRNA_seeds[sRNA_seq[i:(i+SEED_LENGTH)]] = True

        # Determine if mRNA has matching seed to sRNA
        for g in genes:
                if (g[0] != 'CDS'): continue
                mRNA_name, mRNA_seq = g[2], g[8]
                seeds[mRNA_name] = 0  # Default is no matching seed
                for i in range(len(mRNA_seq)-SEED_LENGTH+1):
                        if (mRNA_seq[i:(i+SEED_LENGTH)] in sRNA_seeds): seeds[mRNA_name] = 1
        return seeds


# READ IN INFO ABOUT THE DISTANCE FROM THE START OF EACH GENE TO ITS UPSTREAM GENE
def getGeneDistances(genome, genes):
        # For each gene, get the distance to its upstream gene
        # on the same strand or on any strand
        distances_same_strand, distances_any_strand = {}, {}
        is_upstream_same_strand, is_upstream_any_strand = {}, {}
        for accession in genome:
                genes_accession = []
                for g in genes:
                        if (g[1] == accession): genes_accession.append(g)
                for i in range(len(genes_accession)):
                        g = genes_accession[i]
                        REPLICON_LENGTH = len(genome[accession])
                        start, stop, strand = int(g[4]), int(g[5]), g[6]
                        if (strand == '+'):
                                j = i - 1
                                while (j >= 0) and (genes_accession[j][6] != strand): j -= 1
                                if (j < 0): upstream_same_strand = start  # First plus strand gene
                                else: upstream_same_strand = start - int(genes_accession[j][5]) - 1
                                if (i-1 < 0): upstream_any_strand = start  # First plus strand gene
                                else: upstream_any_strand = start - int(genes_accession[i-1][5]) - 1
                        elif (strand == '-'):
                                j = i + 1
                                while (j < len(genes_accession)) and (genes_accession[j][6] != strand): j += 1
                                if (j >= len(genes_accession)): upstream_same_strand = REPLICON_LENGTH - stop - 1  # Last minus strand gene
                                else: upstream_same_strand = int(genes_accession[j][4]) - stop - 1
                                if (i+1 >= len(genes_accession)): upstream_any_strand = REPLICON_LENGTH - stop - 1  # Last minus strand gene
                                else: upstream_any_strand = int(genes_accession[i+1][4]) - stop - 1
                        distances_same_strand[g[2]] = abs(upstream_same_strand)
                        distances_any_strand[g[2]] = abs(upstream_any_strand)
                        is_upstream_same_strand[g[2]] = 1
                        is_upstream_any_strand[g[2]] = 1
                        if (upstream_same_strand < 0): is_upstream_same_strand[g[2]] = 0
                        if (upstream_any_strand < 0): is_upstream_any_strand[g[2]] = 0
        return is_upstream_same_strand, distances_same_strand, is_upstream_any_strand, distances_any_strand


# HELPER FUNCTION FOR COMPUTING HOMOLOG INFO. COMPUTES HOMOLOGS OF SRNA.
def get_sRNA_homologs(SRNA_FILENAME, genome):
        FILENAME = str(time.time())
        SRNA_HOMOLOGS_FILENAME = FILENAME + '.sRNA_blast'

        # Restrict (self) replicons
        RESTRICT_FILENAME = FILENAME + '.restrict'
        with open(RESTRICT_FILENAME, 'w') as out_file:
                for accession in genome: out_file.write(accession + '\n')
        p = subprocess.run(['./blastdb_aliastool', '-seqid_file_in', RESTRICT_FILENAME], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # BLAST sRNA sequence
        p = subprocess.run(['./blastn', '-db', DB, '-query', SRNA_FILENAME, '-outfmt', '6 qseqid sseqid evalue bitscore qstart qend', '-out', SRNA_HOMOLOGS_FILENAME, '-evalue', '0.01', '-max_target_seqs', '100', '-num_threads', str(NUM_THREADS), '-negative_seqidlist', RESTRICT_FILENAME + '.bsl'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Determine homologs of sRNA
        sRNA_homologs = {}
        sRNA_conservation = {}
        with open(SRNA_HOMOLOGS_FILENAME, 'r') as in_file:
                line = in_file.readline()
                while (line != ''):
                        parse_line = line.split('\t')
                        accession = parse_line[1]
                        if (accession.startswith('ref|')): accession = accession[4:]
                        if (accession.endswith('|')): accession = accession[:-1]
                        sRNA_homologs[accession] = True
                        if (accession not in sRNA_conservation):
                                qstart, qend = int(parse_line[4]), int(parse_line[5])
                                sRNA_conservation[accession] = (qstart, qend)
                        line = in_file.readline()

        # Clean up temp files
        if (os.path.exists(RESTRICT_FILENAME)): os.remove(RESTRICT_FILENAME)
        if (os.path.exists(RESTRICT_FILENAME + '.bsl')): os.remove(RESTRICT_FILENAME + '.bsl')
        if (os.path.exists(SRNA_HOMOLOGS_FILENAME)): os.remove(SRNA_HOMOLOGS_FILENAME)
        return sRNA_homologs, sRNA_conservation


def get_homologs(GENOME_DIR, SRNA_FILENAME, genome, genes):

        sRNA_homologs, sRNA_conservation = get_sRNA_homologs(SRNA_FILENAME, genome)  # sRNA homologs
        if (not os.path.exists(GENOME_DIR + MRNA_HOMOLOGS_FILE)): error('Could not locate file ' + GENOME_DIR + MRNA_HOMOLOGS_FILE)
        with gzip.open(GENOME_DIR + MRNA_HOMOLOGS_FILE, 'rb') as in_file:  # mRNA homologs
                mRNA_homologs = pickle.load(in_file)

        # Determine number of homologs for each sRNA:mRNA pair
        pairs_homologs = {}
        for m in mRNA_homologs:
                if (m not in pairs_homologs): pairs_homologs[m] = {}
                for homolog in mRNA_homologs[m]:
                        if (homolog in sRNA_homologs): pairs_homologs[m][homolog] = True
        return sRNA_homologs, mRNA_homologs, pairs_homologs, sRNA_conservation


# HELPER FUNCTION FOR COMPUTING INTERACTION ENERGIES. COMPUTES SRNA ACCESSIBILITY.
def determine_sRNA_accessibility(GENOME_DIR, sRNA_name, sRNA_sequence):
        TIME_STR = str(time.time())
        WINDOW_SIZE = min(70, len(sRNA_sequence))
        p = subprocess.run(['nice', './RNAplfold', '-u', '40', '-O', '--plex_output', '-W', str(WINDOW_SIZE), '--auto-id', '--id-prefix', TIME_STR], input=sRNA_sequence.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if (p.returncode != 0) or (len(p.stderr.decode()) > 0):
                sys.stderr.write('ERROR executing RNAplfold:\t' + str(p.stderr.decode()) + '\n')
        p = subprocess.run(['./RNAplex', '-a', '.', '-k'], input=sRNA_sequence.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Clean up
        shutil.move(TIME_STR + '_0001_openen_bin', GENOME_DIR + RNAPLFOLD_DIR + sRNA_name + '_openen_bin')
        os.remove(TIME_STR + '_0001_dp.ps')
        os.remove(TIME_STR + '_0001_uplex')
        os.remove(TIME_STR + '_0001_openen')


# HELPER FUNCTION FOR COMPUTING INTERACTION ENERGIES. COMPUTES ENERGIES USING RNAPLEX.
def run_RNAplex(GENOME_DIR, SRNA_FILENAME, f):
        p = subprocess.run(['nice', './RNAplex', '-f', '0', '-q', SRNA_FILENAME, '-t', GENOME_DIR + RNAPLFOLD_DIR + f, '-a', GENOME_DIR + RNAPLFOLD_DIR, '-b'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        structure_info = p.stdout.decode().strip().split('\n')[2]
        return (f, structure_info)


def get_interaction_energies(GENOME_DIR, SRNA_FILENAME, sRNA_name, sRNA_sequence):
        if (not os.path.exists(GENOME_DIR + RNAPLFOLD_DIR)):
                error('Could not find directory ' + GENOME_DIR + RNAPLFOLD_DIR + ' containing RNAplfold results.')
        determine_sRNA_accessibility(GENOME_DIR, sRNA_name, sRNA_sequence)

        # Get list of mRNA files
        mRNA_files = []
        filelist = os.listdir(GENOME_DIR + RNAPLFOLD_DIR)
        for f in filelist:
                if (f.endswith('.fa')): mRNA_files.append(f)

        # Execute RNAplex to determine interaction energy
        energies, structure_data = {}, {}
        with multiprocessing.Pool() as pool: energy_list = pool.map(partial(run_RNAplex, GENOME_DIR, SRNA_FILENAME), mRNA_files)
        for f, structure_info in energy_list:
                energy = float(structure_info.split()[-1][1:-1])
                parse_filename = f.split('____')
                accession = parse_filename[0]
                mRNA_name = parse_filename[1][:-3]
                energies[mRNA_name] = energy
                structure_data[accession + '____' + mRNA_name] = structure_info
        os.remove(GENOME_DIR + RNAPLFOLD_DIR + sRNA_name + '_openen_bin')
        return energies, structure_data


def generate_features(genes, sRNA_length, is_upstream_same_strand, distances_same_strand, is_upstream_any_strand, distances_any_strand, num_sRNA_homologs, mRNA_homologs, pairs_homologs, seeds, energies):
        NUM_FEATURES = 9
        X = np.zeros((len(seeds), NUM_FEATURES))
        previous = {}
        i = 0
        for g in genes:
                if (g[0] != 'CDS'): continue
                name = g[2]
                if (name in previous):
                        X[previous[name], 0] = int(g[5]) - int(g[4]) + 1  # mRNA length
                        continue
                X[i, 0] = int(g[5]) - int(g[4]) + 1  # mRNA length
                if (name in is_upstream_same_strand): X[i, 1] = is_upstream_same_strand[name]
                if (name in distances_same_strand): X[i, 2] = distances_same_strand[name]
                if (name in is_upstream_any_strand): X[i, 3] = is_upstream_any_strand[name]
                if (name in distances_any_strand): X[i, 4] = distances_any_strand[name]
                if (name in pairs_homologs):
                        if (len(pairs_homologs[name]) > 0): X[i, 5] = 1
                        X[i, 6] = len(pairs_homologs[name])
                if (name in seeds): X[i, 7] = seeds[name]  # Seed
                if (name in energies): X[i, 8] = energies[name]
                previous[name] = i
                i += 1
        return X


def predictions(X):
        with open(MODEL_FILE, 'rb') as in_file:
                model_info = pickle.load(in_file)
        model, scaler = model_info['model'], model_info['scaler']
        X_scaled = scaler.transform(X)
        y_pred_proba = model.predict_proba(X_scaled)[:,1]
        return y_pred_proba


def compute_pvalues(y_pred_proba):
        pvalues = np.zeros(len(y_pred_proba))
        params = stats.lognorm.fit(y_pred_proba)
        arg, loc, scale = params[:-2], params[-2], params[-1]
        for i in range(len(y_pred_proba)): pvalues[i] = 1.0 - stats.lognorm.cdf(y_pred_proba[i], loc=loc, scale=scale, *arg)
        return pvalues


def get_output(genes, X, structure_data, y_pred_proba, pvalues):
        cols = ['Target Length', 'upstream same strand', 'distance same strand', 'upstream any strand', 'distance any strand', 'has pair homologs', 'Homologs', 'Seed', 'Energy']
        names, products, mRNA_info, structures, previous = [], [], {}, [], {}
        UPSTREAM, DOWNSTREAM = 200, 100
        t_starts, t_stops, s_starts, s_stops = [], [], [], []
        idx = 0
        for i in range(len(genes)):
                g = genes[i]
                if (g[0] != 'CDS'): continue
                accession, name, synonym, product, sequence, ID = g[1], g[2], g[3], g[7], g[8], g[10]
                name2 = name
                if (name != synonym) and (len(synonym) > 0): name2 = name + ' (' + synonym + ')'
                if (name not in mRNA_info): mRNA_info[name] = {}
                mRNA_info[name][name2] = (accession, sequence, ID)
                parse_structure = structure_data[accession + '____' + name].split()
                if (len(parse_structure) > 3):
                        t_start, t_stop = parse_structure[1].split(',')
                        t_start, t_stop = int(t_start), int(t_stop)
                        if (t_start <= UPSTREAM): t_start = t_start - UPSTREAM - 1
                        else: t_start = t_start - UPSTREAM
                        if (t_stop <= UPSTREAM): t_stop = t_stop - UPSTREAM - 1
                        else: t_stop = t_stop - UPSTREAM
                        s_start, s_stop = parse_structure[3].split(',')
                        s_start, s_stop = int(s_start), int(s_stop)
                        if (s_start == 1): s_stop += 1
                if (name in previous):  # Handle duplicate gene names
                        old_idx = previous[name]
                        names[old_idx] = name2
                        products[old_idx] = product
                        if (len(parse_structure) > 3):
                                structures[old_idx] = parse_structure[0]
                                t_starts[old_idx] = t_start; t_stops[old_idx] = t_stop
                                s_starts[old_idx] = s_start; s_stops[old_idx] = s_stop
                        else:
                                structures[old_idx] = ''
                                t_starts[old_idx] = 0; t_stops[old_idx] = 0
                                s_starts[old_idx] = 0; s_stops[old_idx] = 0
                        continue
                names.append(name2)
                products.append(product)
                if (len(parse_structure) > 3):
                        structures.append(parse_structure[0])
                        t_starts.append(t_start); t_stops.append(t_stop)
                        s_starts.append(s_start); s_stops.append(s_stop)
                else:
                        structures.append('')
                        t_starts.append(0); t_stops.append(0)
                        s_starts.append(0); s_stops.append(0)
                previous[name] = idx
                idx += 1
        df = pd.DataFrame(X, columns=cols)
        df['P-value'] = pvalues
        df['Probability'] = y_pred_proba
        df.insert(0, 'Target', names)
        df.insert(len(df.columns), 'sRNA start', s_starts)
        df.insert(len(df.columns), 'sRNA stop', s_stops)
        df.insert(len(df.columns), 'Target start', t_starts)
        df.insert(len(df.columns), 'Target stop', t_stops)
        df.insert(len(df.columns), 'Annotation', products)
        df.insert(len(df.columns), 'Interaction', structures)
        df = df.astype({'Target Length':int , 'upstream same strand':int , 'distance same strand':int , 'upstream any strand':int , 'distance any strand':int , 'has pair homologs':int , 'Homologs':int , 'Seed':int })
        df_output = df[['Target', 'Energy', 'P-value', 'Probability', 'sRNA start', 'sRNA stop', 'Target start', 'Target stop', 'Annotation', 'Interaction']]

        df_output = df_output[(df_output['P-value']<=PVALUE) & (df_output['Probability']>=PROBABILITY)]
        df_output = df_output.sort_values(by=['Probability'], ignore_index=True, ascending=False)
        df_output.index += 1
        return df_output, mRNA_info


def run_TargetRNA3(sRNA_name, sRNA_sequence):
        genome, genes = read_in_genes(GENOME_DIR)
        seeds = getSeeds(sRNA_sequence, genes)
        is_upstream_same_strand, distances_same_strand, is_upstream_any_strand, distances_any_strand = getGeneDistances(genome, genes)
        sRNA_homologs, mRNA_homologs, pairs_homologs, sRNA_conservation = get_homologs(GENOME_DIR, SRNA_FILENAME, genome, genes)
        energies, structure_data = get_interaction_energies(GENOME_DIR, SRNA_FILENAME, sRNA_name, sRNA_sequence)
        X = generate_features(genes, len(sRNA_sequence), is_upstream_same_strand, distances_same_strand, is_upstream_any_strand, distances_any_strand, len(sRNA_homologs), mRNA_homologs, pairs_homologs, seeds, energies)
        y_pred_proba = predictions(X)
        pvalues = compute_pvalues(y_pred_proba)
        df_output, mRNA_info = get_output(genes, X, structure_data, y_pred_proba, pvalues)
        return df_output, mRNA_info, sRNA_conservation


##############################
##########   MAIN   ##########
##############################

if __name__ == "__main__":
        command()
        arguments()

        sRNA_name, sRNA_sequence = read_in_sRNA()
        df_output, mRNA_info, sRNA_conservation = run_TargetRNA3(sRNA_name, sRNA_sequence)

        if (OUTPUT_FILE != '') and (len(df_output) > 0):
                df_output.to_csv(OUTPUT_FILE, sep='\t', index_label="Index")
                sys.stderr.write('\nTargetRNA3 predicted ' + str(len(df_output)) + ' targets, which it wrote to the output file ' + OUTPUT_FILE + '\n\n')
        elif (OUTPUT_FILE == '') and (len(df_output) > 0): sys.stdout.write('\n' + df_output.to_string() + '\n\n')
        else: sys.stderr.write('\nTargetRNA3 did not find any significant targets. To increase the sensitivity and predict more targets, try lowering the probability threshold. To see all target predictions, use a probability of 0 and a p-value of 1.\n\n')

