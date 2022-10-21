# edgeR4MOLAS User Manual

This script is used for RNA-Seq gene expression between-sample normalization via edgeR.</n> 

For normalization, users should provide sample metadata and read count table.</n>

For pairwise statistical test for differential expressed genes (DEGs), users should provide sample metadata, read count table, and compared groups.

## Input format
1. **FeatureCount.csv (read count table)**

|       | Sample1_Name | Sample2_Name | ... |
| ------|:------------:|:------------:|:---:|
| Gene1 | 10    | 20  | ... |
| Gene2 | 100   | 150 | ... |
| Gene3 | 0     | 15  | ... |

2. **metadata.csv**

| SampleID | Group |
| ---------|:--------:|
| Sample1_Name | GroupA |
| Sample2_Name | GroupA |
| Sample3_Name | GroupB |
| Sample4_Name | GroupB |
| Sample5_Name | GroupB |
| Sample6_Name | GroupC |
| ... | ... |

## Workflow

1. **Set environment**
```
# Establish default folders
mkdir edgeR4MOLAS edgeR4MOLAS/InputData edgeR4MOLAS/OutputData
## move FeatureCount.csv and metadata.csv to edgeR4MOLAS/InputData
mv edgeR4MOLAS
docker run -itd -v $(pwd):/data/ --name edgeR4MOLAS t050508/edger4molas:v3_stable
```

2. **Normalization**
```
docker exec -w /data edgeR4MOLAS edgeR4MOLAS.R -a normalize
```

3. **Statistical test**
```
docker exec -w /data edgeR4MOLAS edgeR4MOLAS.R -a statistica_test -c 'GroupA,GroupB'
```

4. **Turn off edgeR4MOLAS**
```
docker rm -f edgeR4MOLAS
```

## Option list
```
Options:
	-a CHARACTER, --analysis_type=CHARACTER
		select 'normalize' or 'statistical_test'

	-r CHARACTER, --feature_count_table=CHARACTER
		dataset file.csv path [default= InputData/FeatureCount.csv]

	-m CHARACTER, --metadata=CHARACTER
		metadata file.csv path [default= InputData/metadata.csv]

	-o CHARACTER, --output_folder=CHARACTER
		output folder name [default= OutputData]

	-c CHARACTER, --comparison=CHARACTER
		the name of compared groups, which should be 'GroupA,GroupB'

	-h, --help
		Show this help message and exit
```

