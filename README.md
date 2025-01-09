# scRNA embed

## Updates

### 01/08/2025

Data for NN search
- CZ: https://cellxgene.cziscience.com/collections/854c0855-23ad-4362-8b77-6b1639e7a9fc
- ENCODE: https://www.encodeproject.org/matrix/?type=Experiment&control_type!=*&assay_term_name=single-cell+RNA+sequencing+assay&status=released&assay_title=snRNA-seq

| Tissue       | CZ           | ENCODE        |
|--------------|--------------|---------------|
| Heart        | 931k cells   | 49 studies    |
| Hippocampus  | 676k cells   | 49 studies    |
| Liver        | 259k cells   | 8 studies     |
| Pancreas     | 121k cells   | 6 studies     |
| Kidney       | 194k cells   | 1 study       |
| Blood        | 335k cells   | 0 studies     |

### 01/07/2025

Nearest neighbor (NN) search for with PCA and JL embeddings

| Model | Dimensions | Recall at top k=5       | KNN runtime (s=200) |
|-------|------------|-------------------------|---------------------|
| None  | 29898      | (Ground truth)          | 316.45             |
| PCA   | 3537       | 0.386                   | 27.23              |
| JL    | 3537       | 0.713                   | 28.73              |
