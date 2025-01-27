#!/usr/bin/env python3

# characterize mismatching tissues for queries and nearest neighbor

# goal: map cells to their types


import pandas as pd
import sys

# read
df_nn = pd.read_csv('25_01_27-knn_results.tsv',sep='\t')
print('Total number of rows:',df_nn.shape[0])
print(df_nn)


# set cell tag and tissues
df_nn['query_cell_tag'] = df_nn['query_index'].\
    apply(lambda x: x.split('_')[-1])
df_nn['hit_cell_tag1'] = df_nn['hit_index1'].\
    apply(lambda x: x.split('_')[-1])
df_nn['hit_cell_tag2'] = df_nn['hit_index2'].\
    apply(lambda x: x.split('_')[-1])
df_nn['hit_tissue1'] = df_nn['hit_index1'].\
    apply(lambda x: x.split('_')[0])
df_nn['hit_tissue2'] = df_nn['hit_index2'].\
    apply(lambda x: x.split('_')[0])
print(df_nn)
# order columns
df_nn = df_nn[['query_index','query_cell_tag','query_tissue','hit_index1','hit_cell_tag1','hit_tissue1','hit_index2','hit_cell_tag2','hit_tissue2']]

# read cell types
ct_brain_q = pd.read_csv('brain.4341_BA24_10x.ct',sep='\t')
ct_brain_d = pd.read_csv('brain.1823_BA24_10x.ct',sep='\t')
ct_eye_q = pd.read_csv('eye.GSM3745992.ct',sep='\t')
ct_eye_d = pd.read_csv('eye.GSM3745993.ct',sep='\t')
ct_lung_q = pd.read_csv('lung.GSM3926545.ct',sep='\t')
ct_lung_d = pd.read_csv('lung.GSM3926546.ct',sep='\t')
ct_heart_q = pd.read_csv('heart.BO_H46_LV0_R2_premrna.ct',sep='\t')
ct_heart_d = pd.read_csv('heart.BO_H51_LV0_premrna.ct',sep='\t')
ct_testes_q = pd.read_csv('testes.GSM3052919.ct',sep='\t')
ct_testes_d = pd.read_csv('testes.GSM3052918.ct',sep='\t')

# lookup cell types
def lookup_cell_type(cell_tag, tissue, type):
    if type == 'query':
        if tissue == 'brain':
            cell_type = ct_brain_q[ct_brain_q['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'eye':
            cell_type = ct_eye_q[ct_eye_q['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'lung':
            cell_type = ct_lung_q[ct_lung_q['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'heart':
            cell_type = ct_heart_q[ct_heart_q['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'testes':
            cell_type = ct_testes_q[ct_testes_q['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
    elif type == 'hit':
        if tissue == 'brain':
            cell_type = ct_brain_d[ct_brain_d['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'eye':
            cell_type = ct_eye_d[ct_eye_d['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'lung':
            cell_type = ct_lung_d[ct_lung_d['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'heart':
            cell_type = ct_heart_d[ct_heart_d['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
        elif tissue == 'testes':
            cell_type = ct_testes_d[ct_testes_d['cell_id'] == cell_tag]['cell_type'].values[0]
            return cell_type
    
# lookup cell types
cell_types = pd.DataFrame(columns=['query_cell_type','hit_cell_type1', 'hit_cell_type2'])
for i, row in df_nn.iterrows():
    # input
    t_q = row['query_tissue']
    tag_c = row['query_cell_tag']
    t_h1 = row['hit_tissue1']
    tag_h1 = row['hit_cell_tag1']
    t_h2 = row['hit_tissue2']
    tag_h2 = row['hit_cell_tag2']
    # output
    # query
    ct_q = lookup_cell_type(tag_c,t_q,'query')
    # hit
    ct_h1 = lookup_cell_type(tag_h1,t_h1,'hit')
    ct_h2 = lookup_cell_type(tag_h2,t_h2,'hit')
    cell_types = pd.concat([cell_types,pd.DataFrame({'query_cell_type':[ct_q],'hit_cell_type1':[ct_h1],'hit_cell_type2':[ct_h2]})])
# save
cell_types.to_csv('25_01_27-nn_cell_types.tsv',sep='\t',index=False)

# run
# paste mismatching_tissues.tsv mismatching_cell_types.tsv > tissue_mismatch_cell_types.tsv to hstack
