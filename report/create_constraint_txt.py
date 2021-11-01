import argparse
import yaml
from yaml.loader import SafeLoader
import sys
sys.path.append("../utilities")
sys.path.append("../predict-constraint")

from compute_mutation_counts import compute_mutation_counts
from process_gene_gff import gff_to_db,get_gene_feature

#def create_constraint_txt(gene_start, gene_end, chr_num, model, expected_out_filepath='expected_mut_counts.txt', observed_out_filepath='observed_mut_counts.txt', window_size=51, window_stride=25):

def create_constraint_txt(gene_start, gene_end, chr_num, model, expected_out_filepath='expected_mut_counts.txt', observed_out_filepath='observed_mut_counts.txt', window_size=51, window_stride=25):
    exp_fi = open(expected_out_filepath, 'w')
    obs_fi = open(observed_out_filepath, 'w')
    
    mutation_dict = compute_mutation_counts('{}:{}-{}'.format(chr_num,gene_start,gene_end), model, window_size, window_stride)
    xs = mutation_dict['windowPositions']
    y1s = mutation_dict['windowExpectedMutationCounts']
    y2s = mutation_dict['windowObservedMutationCounts']

    for (x,y) in list(zip(xs,y1s)): exp_fi.write('{}\t{}\t{}\n'.format(chr_num,x,y))
    for (x,y) in list(zip(xs,y2s)): obs_fi.write('{}\t{}\t{}\n'.format(chr_num,x,y))

    exp_fi.close()
    obs_fi.close()

# parser = argparse.ArgumentParser(description='create_constraint_txt')
# parser.add_argument('gff_path')
# parser.add_argument('model_path')
# parser.add_argument('gene_name')
# parser.add_argument('seqid')
parser = argparse.ArgumentParser(description='create_constraint_txt')
parser.add_argument('config_file')
args = parser.parse_args()

with open(args.config_file) as f:
    params = list(yaml.load_all(f, Loader=SafeLoader))
    plot_params = params[0]


# 'test.db' to make sure .db does not get overwritten during development    
#gff_db = gff_to_db(args.gff_path, args.gff_path + 'test.db')
gff_db = gff_to_db(plot_params['gff_path'], plot_params['gff_path'] + 'test.db')
#tmp = '/scratch/ucgd/lustre-work/quinlan/u6038618/constraint-tools/data/Homo_sapiens.GRCh38.104.gff3.db'

gene_feature = get_gene_feature(gff_db, plot_params['gene_name'])
#gene_feature = get_gene_feature(gff_db, args.gene_name)

#TODO maybe make this a config file
#create_constraint_txt(gene_feature.start, gene_feature.end, args.seqid, args.model_path, expected_out_filepath='expected_mut_counts.txt', observed_out_filepath='observed_mut_counts.txt', window_size=51, window_stride=25)
create_constraint_txt(gene_feature.start, gene_feature.end, plot_params['seqid'], plot_params['model_path'], expected_out_filepath='expected_mut_counts.txt', observed_out_filepath='observed_mut_counts.txt', window_size=51, window_stride=25)

#python create_constraint_txt.py '/scratch/ucgd/lustre-work/quinlan/u6038618/constraint-tools/data/Homo_sapiens.GRCh38.104.gff3.db' /scratch/ucgd/lustre-work/quinlan/data-shared/constraint-tools/model.ptm.json KCNQ2 chr20