# scRNA embed

## Updates 01/07/2025

Nearest neighbor (NN) search for with PCA and JL embeddings

| Model | Dimensions | Recall at top k=5       | KNN runtime (s=200) |
|-------|------------|-------------------------|---------------------|
| None  | 29898      | (Ground truth)          | 316.45             |
| PCA   | 3537       | 0.386                   | 27.23              |
| JL    | 3537       | 0.713                   | 28.73              |
