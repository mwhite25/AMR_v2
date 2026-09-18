"""
This script is to be used on the CAMDA AMR Prediction Challenge 2025 train and test datasets
This script takes as input a path to a directory of assemblies (.fasta files) and returns a dataset for specified model

Key Options:
    - "model type": "sequence_based" or "matrix_based"
        - "sequence-based": e.g. DNABERT, model takes as input sequence (+ query_id, hit_count, and optionally species/antibiotic)
        - "matrix-based": e.g. Random Forest, model takes as input matrix of features (gene counts and optionally species/antibiotic)

    - "grouping": "full", "per_species", or "per_antibiotic"
        - grouping for DNABERT inference - depends on how we trained DNABERT models, downstream RF should always be per-species
        - "full": returns single dataset for all assemblies
        - "per-species": returns one dataset per species (9 datasets)
        - "per-antibiotic": returns one dataset per antibiotic (4 datasets)

    - "train": boolean, default is False
        - if True, will search for ground truth labels in metadata and return as final column in dataset
        - if False, will not return ground truth labels

    - "split": path to directory containing train/test/dev accession lists (train.txt, test.txt, dev.txt), default is None
        - TODO: NOT YET IMPLEMENTED
        - ***this is for split of the TRAIN SET for evaluation before submission*** to train on entire train set, dont specify
        - if not specified, will return single dataset

Paths to be specified:
    - "assemblies_dir": path to directory containing all assemblies (.fasta files)
    - "metadata_path": path to CSV containing accessions, genus, and species columns, and (if train is True) ground truth labels column
        - as downloaded from CAMDA website
    - "perspecies_dbgwas_dir": path to directory containing significant sequences from DBGWAS 
        - 9 fasta files: {species}_sig_sequences.fasta
    - "output_dir": base directory for output datasets, default is ./datasets

Outputs/Next Steps:
    - <output_dir>/<model_type>/<grouping>/<train/test (official CAMDA split)>/<[full]/[species]/[antibiotic].csv>
    - sequence-based datasets should be sent to dnabert/inference/inference.py for inference
        - if finetuning is required, pass full dataset (<output_dir>/<model_type>/<grouping>/<train/test (official CAMDA split)>/full_sequence_dataset.csv) to dnabert/finetune/run_finetune.sh to get a model
    - matrix-based datasets should be sent to rf or xgboost directory for finetuning/inference
    - if "split" is specified, <train>/<test>/... can be passed to (TODO) evaluation script for accession-level evaluation before submission



Additional parameters listed below, but should remain relatively constant and only be modified once at a time for consistency


"""


## parameters: keep these relatively constant for consistency

#universal parameters for blast
#BLAST_IDENTITY = 80
LEAKAGE = True #should only set to False to prevent data leakage when running an abalation study on spcecies/antibiotic features on full/per-antibiotic models
DBGWAS_SIG_LEVEL = str(0.05) #THIS IS JUST TO FIND APPROPRIATE FILES, we do NOT filter sequences in this script, this should be done in DBGWAS script
MAX_TARGET_SEQS = 10 #for BLAST


#sequence-based parameters
#SEQ_LENGTH = 1000

#matrix-based parameters
BINARY_FEATURES = False


## global mappings:
species_list = [
    'neisseria_gonorrhoeae', 
    'staphylococcus_aureus', 
    'streptococcus_pneumoniae', 
    'salmonella_enterica', 
    'klebsiella_pneumoniae', 
    'escherichia_coli', 
    'pseudomonas_aeruginosa', 
    'acinetobacter_baumannii', 
    'campylobacter_jejuni' 
]

antibiotic_list = ['GEN', 'ERY', 'CAZ', 'TET']

species_to_antibiotic = {
    'klebsiella_pneumoniae': 'GEN',
    'streptococcus_pneumoniae': 'ERY',
    'escherichia_coli': 'GEN',
    'campylobacter_jejuni': 'TET',
    'salmonella_enterica': 'GEN',
    'neisseria_gonorrhoeae': 'TET',
    'staphylococcus_aureus': 'ERY',
    'pseudomonas_aeruginosa': 'CAZ',
    'acinetobacter_baumannii': 'CAZ',
}

antibiotic_to_species = {
    'GEN': ['klebsiella_pneumoniae', 'escherichia_coli', 'salmonella_enterica'],
    'ERY': ['streptococcus_pneumoniae', 'staphylococcus_aureus'],
    'CAZ': ['pseudomonas_aeruginosa', 'acinetobacter_baumannii'],
    'TET': ['campylobacter_jejuni', 'neisseria_gonorrhoeae'],
}

# index mappings:
species_mapping = {
    'klebsiella_pneumoniae': 0,
    'streptococcus_pneumoniae': 1,
    'escherichia_coli': 2,
    'campylobacter_jejuni': 3,
    'salmonella_enterica': 4,
    'neisseria_gonorrhoeae': 5,
    'staphylococcus_aureus': 6,
    'pseudomonas_aeruginosa': 7,
    'acinetobacter_baumannii': 8
}

antibiotic_mapping = {'GEN': 0, 
    'ERY': 1,
    'CAZ': 2,
    'TET': 3,
    'tetracycline': 3
}

label_map = {'Resistant': 0, 'Intermediate': 0, 'Susceptible': 1}

# imports
import os
import argparse
import subprocess
import pandas as pd
import csv
from Bio import SeqIO
from tqdm import tqdm



# methods:


def group_sig_seqs(grouping, output_dir, perspecies_dir): #NOT YET TESTED
    #if grouping is full, concat all species files to full_sig_sequences.fasta
    if grouping == 'full':
        with open(os.path.join(output_dir, 'full_sig_sequences.fasta'), 'w') as output_file:
            for species in species_list:
                with open(os.path.join(perspecies_dir, f'{species}_sig_sequences.fasta'), 'r') as f:
                    for line in f:
                        if line.startswith('>'):
                            output_file.write(line)

    #if grouping is per_antibiotic, concat all antibiotic files to TET_sig_sequences.fasta, GEN_sig_sequences.fasta, etc. according to antibiotic_to_species
    elif grouping == 'per_antibiotic':
        for antibiotic in antibiotic_list:
            with open(os.path.join(output_dir, f'{antibiotic}_sig_sequences.fasta'), 'w') as output_file:
                for species in antibiotic_to_species[antibiotic]:
                    with open(os.path.join(perspecies_dir, f'{species}_sig_sequences.fasta'), 'r') as f:
                        for line in f:
                            if line.startswith('>'):
                                output_file.write(line)

    elif grouping == 'per_species':
        return

    else:
        raise ValueError(f'Invalid grouping: {grouping}')


def get_hits(accessions_list, assemblies_dir, species_to_accessions, dbgwas_dir, grouping, train_test, args):
    #runs BLAST to align sequences in dbgwas_dir to assemblies
    #if leakage is True (default), aligns assemblies to their species-specific features
    #if leakage is False, aligns to all sig seqs (if grouping is full), or to all sig seqs in antibiotic grouping (if grouping is per_antibiotic)

    def run_blast():
        for accession in tqdm(accessions, desc=f'Running BLAST for {species}'):
                assembly_path = os.path.join(assemblies_dir, f'{accession}.fasta')
                db_dir = os.path.join(db_base_dir, accession, accession) #make db in subdirectory of accession name
                os.makedirs(hits_base_dir, exist_ok=True)
                hits_file = os.path.join(hits_base_dir, f'{accession}_hits.tsv')

                if not os.path.exists(hits_file):
                    if not os.path.exists(db_dir):
                        cmd = f'makeblastdb -in {assembly_path} -dbtype nucl -out {db_dir}'
                        subprocess.run(cmd, shell=True)
                    
                    cmd = (f'blastn -query {sig_seqs_path} -db {db_dir} '
                        f'-max_target_seqs {MAX_TARGET_SEQS} -outfmt "6 qseqid sseqid pident length '
                        f'qstart qend sstart send sstrand bitscore" -perc_identity {BLAST_IDENTITY} '
                        f'-out {hits_file}')
                    print(hits_file)
                    subprocess.run(cmd, shell=True)

    db_base_dir = f'work/blast/dbs/{train_test}'
    hits_base_dir = f'work/blast/hits/{train_test}/{args.run_name}'
    os.makedirs(db_base_dir, exist_ok=True)
    os.makedirs(hits_base_dir, exist_ok=True)

    if LEAKAGE or grouping == 'per_species': 
        hits_base_dir = os.path.join(hits_base_dir, 'per_species') #should still use 'per_species' directory when doing full or per-antibiotic model unless LEAKAGE is False, when we're trying to do an abalation study
        for species in species_list:
            if species not in species_to_accessions:
                print(f'SKIP: no metadata rows for species {species}; not present in the current metadata subset')
                continue
            accessions = list(set(species_to_accessions[species]) & set(accessions_list))
            if not accessions:
                print(f'SKIP: species {species} has no accessions in the current metadata subset')
                continue
            sig_seqs_path = os.path.join(dbgwas_dir, f'{species}_sig_sequences.fasta')
            os.makedirs(hits_base_dir, exist_ok=True)

            run_blast()
                
    else:
        if grouping == 'full':
            sig_seqs_path = os.path.join(dbgwas_dir, 'full_sig_sequences.fasta')
            accessions = accessions_list
            hits_base_dir = os.path.join(hits_base_dir, 'full')
            os.makedirs(hits_base_dir, exist_ok=True)

            run_blast()

        elif grouping == 'per_antibiotic':
            hits_base_dir = os.path.join(hits_base_dir, 'per_antibiotic')
            os.makedirs(hits_base_dir, exist_ok=True)
            for antibiotic in antibiotic_list:
                antibiotic_species = antibiotic_to_species[antibiotic]
                accessions = []
                for species in antibiotic_species:
                    accessions.extend(species_to_accessions[species])
                sig_seqs_path = os.path.join(dbgwas_dir, f'{antibiotic}_sig_sequences.fasta')
                hits_base_dir = os.path.join(hits_base_dir, 'per_antibiotic')
                os.makedirs(hits_base_dir, exist_ok=True)

                run_blast()

def generate_sequence_dataset(accession_list, output_dir, grouping, train, accession_to_species, accession_to_antibiotic, accession_to_phenotype, assemblies_dir, train_test, args):
    #generates full csv dataset with the following columns:
    #sequence, accession, query_id, hit_count, species, antibiotic, and (if train =True) phenotype
    #labels will be mapped to numerical value according to universal mappings
    #hit_count is NOT normalized (as it will likely not be used as feature to DNABERT)

    def get_flanking_regions(accession, sseqids, sstarts, sends):
        #get universal SEQ_LENGTH flanking regions around each hit
        #read assembly file
        assembly_path = os.path.join(assemblies_dir, f'{accession}.fasta')
        #print(f'assembly path: {assembly_path}')
        contigs = {}
        #print(f'accession: {accession}')
        #print(f'sseqids: {set(sseqids)}')
        if not os.path.exists(assembly_path):
            print(f'Assembly file not found for accession: {accession}')
            return []
        for record in SeqIO.parse(assembly_path, 'fasta'):
            if record.id in set(sseqids):
                contig = str(record.seq)
                contigs[record.id] = contig
        #print(f'contig ids: {contigs.keys()}')
        sequences = []
        for sseqid, sstart, send in zip(sseqids, sstarts, sends):
            if sseqid not in contigs:
                print(f'WARNING: Contig not found for sseqid: {sseqid} for accession: {accession}')
                continue
            midpoint = (int(sstart) + int(send)) // 2
            contig = contigs[sseqid]
            flanking_region = contig[max(0, midpoint - SEQ_LENGTH//2) : min(len(contig), midpoint + SEQ_LENGTH//2)] #get region around midpoint of blast hit, make sure to not go out of bounds of contig length
            sequences.append(flanking_region)

        return sequences

    def get_random_seq(accession):
        #get random seq by finding first contig >seq_len, if none then just return empty string
        assembly_path = os.path.join(assemblies_dir, f'{accession}.fasta')
        if not os.path.exists(assembly_path):
            print(f'Assembly file not found for accession: {accession}')
            return ''
        for record in SeqIO.parse(assembly_path, 'fasta'):
            if len(str(record.seq)) > SEQ_LENGTH:
                sequence = str(record.seq)[0:SEQ_LENGTH]
                return sequence

        #if never returned, return empty string
        return ''





    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(os.path.join(output_dir, 'full_sequence_dataset.csv'), 'w') as dataset: #first generate full set, then split later based on grouping ( if grouping != 'full)
        writer = csv.writer(dataset)
        if train:
            writer.writerow(['sequence', 'accession', 'query_id', 'hit_count', 'species', 'antibiotic', 'phenotype'])
        else:
            writer.writerow(['sequence', 'accession', 'query_id', 'hit_count', 'species', 'antibiotic'])
        #iterate through accession list, find hit file, get query_id, hit_count, flanking regions (sequence), and map accession to species and antibiotic
        for accession in tqdm(accession_list, desc='Generating sequence dataset'):
            hit_file = os.path.join('work', 'blast', 'hits', train_test, args.run_name, 'per_species', f'{accession}_hits.tsv')
            hit_counts = {} #hit counts for each query_id
            sequences = []
            sseqids = []
            sstarts = []
            sends = []
            query_ids = []
            with open(hit_file, 'r') as f:
                for line in f:
                    qseqid, sseqid, pident, length, qstart, qend, sstart, send, sstrand, bitscore = line.split('\t')
                    query_ids.append(qseqid)
                    sseqids.append(sseqid)
                    sstarts.append(sstart)
                    sends.append(send)
                    if qseqid not in hit_counts.keys():
                        hit_counts[qseqid] = 0
                    hit_counts[qseqid] += 1

            species = species_mapping[accession_to_species[accession]]
            antibiotic = antibiotic_mapping[accession_to_antibiotic[accession]]
            if train:
                phenotype = label_map[accession_to_phenotype[accession]]

            #print(f'hit file: {hit_file}')
            sequences = get_flanking_regions(accession, sseqids, sstarts, sends)
            if not sequences: #if no sequences, we need to add something so that the accession stays in the set
                sequences.append(get_random_seq(accession))
                query_ids.append('random') #random is ID for no hits
                hit_counts['random'] = 1

            #write rows for accession by iterating through lists
            for sequence, query_id in zip(sequences, query_ids):
                if train:
                    writer.writerow([sequence, accession, query_id, hit_counts[query_id], species, antibiotic, phenotype])
                else:
                    writer.writerow([sequence, accession, query_id, hit_counts[query_id], species, antibiotic])
            #dataset.flush()

def group_dataset(grouping, output_dir, full_seq_path=None):
    #group if not full model
    #only for sequence-based grouping 

    if grouping == 'per_species':
        if full_seq_path is None:
            df = pd.read_csv(os.path.join(output_dir, 'full_sequence_dataset.csv'))
        else:
            df = pd.read_csv(full_seq_path)
        for species in species_list:
            df_species = df[df['species'] == species_mapping[species]]
            if not os.path.exists(os.path.join(output_dir, species)):
                os.makedirs(os.path.join(output_dir, species))
            df_species.to_csv(os.path.join(output_dir, species, f'{species}_sequence_dataset.csv'), index=False)
            print(f'Wrote {len(df_species)} rows to {species}_sequence_dataset.csv')

    elif grouping == 'per_antibiotic':
        if full_seq_path is None:
            df = pd.read_csv(os.path.join(output_dir, 'full_sequence_dataset.csv'))
        else:
            df = pd.read_csv(full_seq_path)
        for antibiotic in antibiotic_list:
            df_antibiotic = df[df['antibiotic'] == antibiotic_mapping[antibiotic]]
            if not os.path.exists(os.path.join(output_dir, antibiotic)):
                os.makedirs(os.path.join(output_dir, antibiotic))
            df_antibiotic.to_csv(os.path.join(output_dir, antibiotic, f'{antibiotic}_sequence_dataset.csv'), index=False)
            print(f'Wrote {len(df_antibiotic)} rows to {antibiotic}_sequence_dataset.csv')
                    

    
def main():


    # parse arguments
    parser = argparse.ArgumentParser(description='CAMDA AMR Prediction Challenge 2025 Data Pipeline')
    parser.add_argument('--assemblies_dir', type=str, help='Path to directory containing all assemblies (.fasta files)', default=None)
    parser.add_argument('--metadata_path', type=str, help='Path to CSV containing accessions, ground truth labels, genus, and species columns', default=None)
    parser.add_argument('--output_dir', type=str, help='Base directory for output datasets', default='./datasets')
    parser.add_argument('--model_type', type=str, help='Model type', choices=['sequence_based', 'matrix_based'], default='sequence_based')
    parser.add_argument('--grouping', type=str, help='Grouping', choices=['full', 'per_species', 'per_antibiotic'], default='per_species')
    parser.add_argument('--split', type=str, help='Path to directory containing train/test/dev accession lists (train.txt, test.txt, dev.txt)', default=None)
    parser.add_argument('--perspecies_dbgwas_dir', type=str, help='Path to directory containing significant sequences from DBGWAS', default=f'./data/dbgwas/p{DBGWAS_SIG_LEVEL}/per_species')
    parser.add_argument('--train', action='store_true')
    parser.add_argument('--seq_length', type=int, default=1000)
    parser.add_argument('--group_only', type=str, default=None, help='give path to full sequence dataset if already created to group to new grouping')
    parser.add_argument('--run_name', type=str, default=None)
    parser.add_argument('--blast_identity', type=float, default=80)
    args = parser.parse_args()
    print(f'Arguments: {args}')
    global SEQ_LENGTH
    SEQ_LENGTH = args.seq_length
    global BLAST_IDENTITY
    BLAST_IDENTITY = args.blast_identity

    # get base output directory
    train_test = 'train' if args.train else 'test'
    base_output_dir = os.path.join(args.output_dir, args.run_name, args.model_type, args.grouping, train_test)

    # parse metadata
    #read dataframes, create genus_species column from genus and species columns
    #create dictionaty mapping genus_species to list of accessions
    #create accession lists for train and test
    metadata = pd.read_csv(args.metadata_path)
    metadata['genus_species'] = metadata['genus'].str.lower() + '_' + metadata['species']
    accessions = metadata['accession'].tolist()
    accessions = list(set(accessions)) #remove duplicate accessions 
    species_to_accessions = metadata.groupby('genus_species')['accession'].apply(list).to_dict()
    accession_to_species = metadata.set_index('accession')['genus_species'].to_dict()
    accession_to_antibiotic = metadata.set_index('accession')['antibiotic'].to_dict()
    if args.train:
        accession_to_phenotype = metadata.set_index('accession')['phenotype'].to_dict()
    else:
        accession_to_phenotype = None





    ## start pipeline 

    if args.group_only is not None:
        group_dataset(args.grouping, base_output_dir, args.group_only)
        exit()

    # if running an abalation study on species/antibiotic features, prevent species/antibiotic data leakage by grouping significant sequences by species/antibiotic
    if not LEAKAGE:
        #group DBGWAS significant sequencesfasta files:
        #for full model, groups all to ./data/dbgwas/p{threshold}/full/full_sig_sequences.fasta
        #for per-antibiotic model, groups all to ./data/dbgwas/p{threshold}/per_antibiotic/TET_sig_sequences.fasta, GEN_sig_sequences.fasta, etc. according to antibiotic_to_species
        #for per-species model, does not group (input should be grouped by species, this is how DBGWAS is run), files stay in args.perspecies_dbgwas_dir
        dbgwas_dir = os.path.join('data', 'dbgwas', f'p{DBGWAS_SIG_LEVEL}', args.grouping)
        group_sig_seqs(args.grouping, dbgwas_dir, args.perspecies_dbgwas_dir)
    else:
        dbgwas_dir = args.perspecies_dbgwas_dir

    # run BLAST
    #generates {accession}_hits.tsv files
    #aligns sequences in dbgwas_dir to assemblies
    #normally aligns assemblies to their species-specific features, but if LEAKAGE is False, takes args.grouping and aligns to all sequences in dbgwas_dir
    #output {accession}_hits.tsv files in work/blast/hits/species_specific/{accession}_hits.tsv by default, unless LEAKAGE is False
    get_hits(accessions_list=accessions, assemblies_dir=args.assemblies_dir, species_to_accessions=species_to_accessions, dbgwas_dir=dbgwas_dir, grouping=args.grouping, train_test=train_test, args=args)
    
    # get datasets
    #run generate_matrix() for matrix-based models and generate_sequence_dataset() for sequence-based models
    if args.model_type == 'sequence_based':
        generate_sequence_dataset(accession_list=accessions, output_dir=base_output_dir, grouping=args.grouping, train=args.train, accession_to_species=accession_to_species, accession_to_antibiotic=accession_to_antibiotic, accession_to_phenotype=accession_to_phenotype, assemblies_dir=args.assemblies_dir, train_test=train_test, args=args)
    elif args.model_type == 'matrix_based':
        generate_matrix_dataset()

    #group dataset (if not full model):
    if args.model_type == 'sequence_based':
        group_dataset(args.grouping, base_output_dir)


    



if __name__ == "__main__":
    main()
