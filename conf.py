#Config file for the hiC_query pipeline

SNP_DATABASE_FP = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/snpIndex.db" #The database of SNPs to search for details on input dbSNP IDs
HIC_DATA_DIR = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/hic_data" #The directory containing the directories representing cell lines of HiC experiment tables
FRAGMENT_BED_FP = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/Homo_sapiens.ensembl.release-74.MboI.fragments.bed" #Bed file detailing the start and end points of HiC experiment fragments (hg19 fragmented by MboI by default.)
FRAGMENT_DATABASE_FP = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/fragmentIndex.db" #The database of fragments to search when assigning SNPs to fragments (GRCh37 digested by MboI by default.)
GENE_BED_FP = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/hg19_genes.bed" #Bed file detailing locations of genes in genome (the UCSC known gene list for hg19 by default).
GENE_DATABASE_FP = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/geneIndex.db" #Database constructed from the gene BED file.
EQTL_DATA_DIR = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/eQTLs" #The directory containing databases of eQTL data from various tissues (the GTEx Analysis V6 by default).
EXPRESSION_TABLE_FP = "/mnt/3dgenome/projects/cekb635/hiCquery_dev/lib/GTEx_Analysis_v6_RNA-seq_gene_median.txt" #The tab-separated table containing expression information about GTEx genes (the GTEx Analysis V6 RNA-seq gene median by default).
