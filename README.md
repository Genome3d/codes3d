The hiC_query script expects the following file structure:


	[base hiC_query directory]/
		eQTLs/
			[databases created from tissue-specific eQTL tables, provided by GTEx in the file GTEx_Analysis_V6_eQTLs.tar.gz and created using index_tables.py with the -e option]
		hic_data/
			[one directory for each cell line]/
				<replicate_name>.db (Databases made from HiC interaction table, e.g. tables by Rao et al. or tables in similar format, using index_tables.py with the -i option.)
		snps/ (Optional if snpIndex.db is present in base hiC_query directory)
			[BED files for dbSNP SNPs for each chromosome]
		snpIndex.db (Will be created automatically from the contents of snps/ if it does not exist)
		hiC_query.py
		index_tables.py (Not necessary for querying, but required if wishing to build databases from source tables, or use custom tables)

Some test data is available in the directories beginning with "test". To run the test data, use the programs as follows:

To index the HiC table, run:

		./index_tables.py -t test_hic_data/TEST/test_hic.txt -i (this will produce the file test_hic.db in the same directory by default.)

To query the HiC table against the eQTL databases using one of the supplied test SNPs, run:

		./hiC_query.py -i <any combination of [rs001 rs002 rs003 rs004 rs005 rs006 rs007 rs008 rs009 rs010 rs011 rs012 rs013 rs014 rs015 rs016 rs017 rs018 rs019 rs020 rs021 rs022 rs023 rs024 rs025]> -t test_snps -c test_hic_data -e test_eQTLs
