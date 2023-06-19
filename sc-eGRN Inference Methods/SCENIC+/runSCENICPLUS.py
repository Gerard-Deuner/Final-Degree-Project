# Set up Environment
import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

# set working directory
work_dir = '/g/scb/zaugg/deuner/SCENIC+/'
# set tmp directory
tmp_dir = '/g/scb/zaugg/deuner/SCENIC+/tmp/combined/'
# set the figures directory
fig_dir = '/g/scb/zaugg/deuner/SCENIC+/figures/'
# set the output data directory
out_dir = '/g/scb/zaugg/deuner/SCENIC+/outputdata/'

def main():
    """
    The main executable function
    """
    
    # Load the AnnData object containing the scRNA-seq side of the analysis
    adata = sc.read_h5ad(os.path.join(tmp_dir, 'combined.nomicro.adata.h5ad'))

    # Load the cisTopic object containing the scATAC-seq side of the analysis.
    cistopic_obj = dill.load(open(os.path.join(tmp_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))

    # Load the motif enrichment dictionary containing the motif enrichment results.
    menr = dill.load(open(os.path.join(tmp_dir, 'motifs/menr.pkl'), 'rb'))

    
    # adapt barcodes of cistopic object
    new_bcs = []
    old_bcs = cistopic_obj.cell_names
    for i in range(len(old_bcs)):
        split_bc = str.split(old_bcs[i], "_")
        new_bc = split_bc[1] + "_" + split_bc[0]
        new_bcs.append(new_bc)
        cistopic_obj.selected_model.cell_topic.columns.values[i] = new_bc
    cistopic_obj.cell_names = new_bcs
    cistopic_obj.cell_data.index = new_bcs

    # maybe select the atac barcodes as defaults for adata
    list(set(adata.obs["barcode"]) & set(list(cistopic_obj.cell_names.copy())))
    l_bcs = []
    for i in range(len(adata.obs_names)):
        l_bc = adata.obs["orig.ident"][i] + "_" + adata.obs["barcode"][i]
        l_bcs.append(l_bc)
    adata.obs["long_barcode"] = l_bcs
    adata.obs_names = list(adata.obs["long_barcode"])

    # do the same for the cistopic object
    cistopic_obj.cell_data["long_barcode"] = cistopic_obj.cell_names

    print(len(adata.obs_names.copy(deep=True)) == len(set(adata.obs_names.copy(deep=True))))
    print(len(adata.obs_names.copy(deep=True).drop_duplicates(keep='first')))
    print(len(set(adata.obs_names.copy(deep=True))))
    print(len(list(cistopic_obj.cell_names.copy())) == len(set(list(cistopic_obj.cell_names.copy()))))
    print(len(list(cistopic_obj.cell_names.copy())))
    print(len(set(list(cistopic_obj.cell_names.copy()))))

    # Create the Scenic+ object
    from scenicplus.scenicplus_class import create_SCENICPLUS_object
    import numpy as np
    scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = adata,
        cisTopic_obj = cistopic_obj,
        menr = menr,
        gene_metadata = adata.var.copy(deep=True),
        bc_transform_func = None, #lambda x: f'{x}_timecourse' #None, #function to convert scATAC-seq barcodes to scRNA-seq ones
    )
    scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
    scplus_obj

    # Select the optimal gene names host
    ensembl_version_dict = {'105': 'http://www.ensembl.org',
                            '104': 'http://may2021.archive.ensembl.org/',
                            '103': 'http://feb2021.archive.ensembl.org/',
                            '102': 'http://nov2020.archive.ensembl.org/',
                            '101': 'http://aug2020.archive.ensembl.org/',
                            '100': 'http://apr2020.archive.ensembl.org/',
                            '99': 'http://jan2020.archive.ensembl.org/',
                            '98': 'http://sep2019.archive.ensembl.org/',
                            '97': 'http://jul2019.archive.ensembl.org/',
                            '96': 'http://apr2019.archive.ensembl.org/',
                            '95': 'http://jan2019.archive.ensembl.org/',
                            '94': 'http://oct2018.archive.ensembl.org/',
                            '93': 'http://jul2018.archive.ensembl.org/',
                            '92': 'http://apr2018.archive.ensembl.org/',
                            '91': 'http://dec2017.archive.ensembl.org/',
                            '90': 'http://aug2017.archive.ensembl.org/',
                            '89': 'http://may2017.archive.ensembl.org/',
                            '88': 'http://mar2017.archive.ensembl.org/',
                            '87': 'http://dec2016.archive.ensembl.org/',
                            '86': 'http://oct2016.archive.ensembl.org/',
                            '80': 'http://may2015.archive.ensembl.org/',
                            '77': 'http://oct2014.archive.ensembl.org/',
                            '75': 'http://feb2014.archive.ensembl.org/',
                            '54': 'http://may2009.archive.ensembl.org/'}

    import pybiomart as pbm
    def test_ensembl_host(scplus_obj, host, species):
        dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
        annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
        annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
        annot['Chromosome'] = annot['Chromosome'].astype('str')
        filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
        annot = annot[~filter]
        annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
        gene_names_release = set(annot['Gene'].tolist())
        print(len(list(set(gene_names_release) & set(scplus_obj.gene_names))) > 0)
        ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
        print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
        return ov

    n_overlap = {}
    for version in ensembl_version_dict.keys():
        print(f'host: {version}')
        try:
            n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
        except:
            print('Host not reachable')
    v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
    print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

    # Choose the best host
    biomart_host = "http://sep2019.archive.ensembl.org/"

    # Before running  also download a list of known human TFs from the human transcription factors database
    #!wget -O /g/scb/zaugg/deuner/SCENIC+/inputdata/utoronto_human_tfs_v_1.01.txt http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt

    # Also download a the program bedToBigBed this will be used to generate files which can be uploaded to the UCSC genome browser
    #!wget -O /g/scb/zaugg/deuner/SCENIC+/inputdata/bedToBigBed http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
    #!chmod +x /g/scb/zaugg/deuner/SCENIC+/inputdata/bedToBigBed

    # Run the analysis
    from scenicplus.wrappers.run_scenicplus import run_scenicplus
    try:
        run_scenicplus(
            scplus_obj = scplus_obj,
            variable = ['GEX_celltype'],
            species = 'hsapiens',
            assembly = 'hg38',
            tf_file = '/g/scb/zaugg/deuner/SCENIC+/inputdata/utoronto_human_tfs_v_1.01.txt',
            save_path = os.path.join(tmp_dir, 'scenicplus'),
            biomart_host = biomart_host,
            upstream = [1000, 150000],
            downstream = [1000, 150000],
            calculate_TF_eGRN_correlation = True,
            calculate_DEGs_DARs = True,
            export_to_loom_file = True,
            export_to_UCSC_file = True,
            path_bedToBigBed = '/g/scb/zaugg/deuner/SCENIC+/inputdata',
            n_cpu = 24,
            _temp_dir = None)#'/g/scb/zaugg/deuner/ray_spill')
    except Exception as e:
        #in case of failure, still save the object
        dill.dump(scplus_obj, open(os.path.join(out_dir, '/scplus_obj.pkl'), 'wb'), protocol=-1)
        raise(e)

if __name__ == "__main__":
    main()

