#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
_stderr = sys.stderr
null = open(os.devnull,'wb')

# set the working directory (now set locally)
work_dir = '/g/scb/zaugg/deuner/SCENIC+/'
# set tmp directory
tmp_dir = '/g/scb/zaugg/deuner/SCENIC+/tmp/combined/'
# set the figures directory
fig_dir = '/g/scb/zaugg/deuner/SCENIC+/figures/'

import pickle
import sys
import argparse
import os
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *

#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
_stderr = sys.stderr
null = open(os.devnull,'wb')

# set the working directory (now set locally)
work_dir = '/g/scb/zaugg/deuner/SCENIC+/'
# set tmp directory
tmp_dir = '/g/scb/zaugg/deuner/SCENIC+/tmp/combined/'
# set the figures directory
fig_dir = '/g/scb/zaugg/deuner/SCENIC+/figures/'

def main():
    """
    The main executable function
    """
    # Load candidate enhancer regions identified in previous step.
    import pickle
    region_bin_topics_otsu = pickle.load(open(os.path.join(tmp_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
    region_bin_topics_top3k = pickle.load(open(os.path.join(tmp_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
    markers_dict = pickle.load(open(os.path.join(tmp_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))
    
    # Convert to dictionary of pyranges objects
    import pyranges as pr
    from pycistarget.utils import region_names_to_coordinates
    region_sets = {}
    region_sets['topics_otsu'] = {}
    region_sets['topics_top_3'] = {}
    region_sets['DARs'] = {}
    for topic in region_bin_topics_otsu.keys():
        regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    for topic in region_bin_topics_top3k.keys():
        regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
        region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

    for DAR in markers_dict.keys():
        regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
        if len(regions) > 0:
            region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

    # Define rankings, score and motif annotation database
    db_fpath = "/g/scb2/zaugg/deuner/SCENIC+/cistarget_databases/"
    motif_annot_fpath = "/g/scb2/zaugg/deuner/SCENIC+/cistarget_databases/"

    rankings_db = os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
    scores_db =  os.path.join(db_fpath, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
    motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

    # Run pycistarget using the run_pycistarget wrapper function
    if not os.path.exists(os.path.join(tmp_dir, 'motifs')):
        os.makedirs(os.path.join(tmp_dir, 'motifs'))

    from scenicplus.wrappers.run_pycistarget import run_pycistarget
    run_pycistarget(
        region_sets = region_sets,
        species = 'homo_sapiens',
        save_path = os.path.join(tmp_dir, 'motifs'),
        ctx_db_path = rankings_db,
        dem_db_path = scores_db,
        path_to_motif_annotations = motif_annotation,
        run_without_promoters = True,
        n_cpu = 24,
        _temp_dir = '/scratch/deuner/ray_spill',
        annotation_version = 'v10nr_clust',
        )

    # Explore some of the results
    import dill
    menr = dill.load(open(os.path.join(tmp_dir, 'motifs/menr.pkl'), 'rb'))

if __name__ == "__main__":
    main()

