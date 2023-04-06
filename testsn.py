#!/usr/bin/env python
# coding: utf-8

# In[11]:


import pandas as pd
from Bio import Entrez, SeqIO
import ssl
ssl._create_default_https_context = ssl._create_unverified_context


# In[12]:


Entrez.email="XXXXXX"


# In[13]:


dataset = pd.read_csv('data/nextstrain_ncov_open_global_6m_metadata.tsv',sep="\t")
print(len(dataset))
dataset.dropna(subset=['genbank_accession'], inplace=True)
dataset.dropna(subset=['pango_lineage'], inplace=True)
dataset.sort_values(by='date', inplace=True)
print(len(dataset))


# In[14]:


accnums=""
for accnum in dataset['genbank_accession'].values:
    accnums = accnums + ", " + str(accnum)
accnums=accnums[2:]


# In[15]:


def accnum_to_seq(accnums):
    handle=Entrez.efetch(db='nucleotide',id=accnums,rettype='fasta',retmode='text')
    records=SeqIO.parse(handle,'fasta')
    seqs=[]
    for record in records:
        seqs.append(str(record.seq))
    return seqs


# In[16]:


from dna2vec.dna2vec.multi_k_model import MultiKModel
filepath = 'dna2vec/pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
mk_model = MultiKModel(filepath)


# In[17]:


import numpy as np
def seq_to_vec(seq):
    seq = seq.replace('R','')
    seq = seq.replace('N','')
    seq = seq.replace('Y','')
    seq = seq.replace('S','')
    seq = seq.replace('W','')
    seq = seq.replace('K','')
    seq = seq.replace('M','')
    seq = seq.replace('H','')
    total_vec = np.zeros((100,), dtype=np.float32)
    kmer = ""
    for i in range(len(seq[:-2])):
        kmer = seq[i] + seq[i+1] + seq[i+2]
        total_vec = total_vec + mk_model.vector(kmer)
    return(total_vec)


# In[18]:


seqs = accnum_to_seq(accnums)


# In[19]:


X_all = []
for seq in seqs:
    X_all.append(seq_to_vec(seq))


# In[20]:


print(len(X_all))


# In[21]:


pango = dataset['pango_lineage'].values


# In[22]:


pango = sorted(list(set(pango)))


# In[23]:


pango


# In[24]:


from sklearn.preprocessing import LabelEncoder
def variant_to_vec(pango):
    pango_arr = []
    for variant in pango:
        if '.' in variant:
            variant = variant.split('.') 
        else:
            variant = [variant]
        rem = 4-len(variant) 
        variant_arr=[]
        for v in variant:
            variant_arr.append(v)
        for r in range(rem):
            variant_arr.append('0')
        pango_arr.append(variant_arr)
    level_one = []
    level_two = []
    level_three = []
    level_four = []
    for i in range(len(pango_arr)):
        level_one.append(pango_arr[i][0])
        level_two.append(pango_arr[i][1])
        level_three.append(pango_arr[i][2])
        level_four.append(pango_arr[i][3])
    labels_one = LabelEncoder().fit_transform(level_one)
    labels_two = LabelEncoder().fit_transform(level_two)
    labels_three = LabelEncoder().fit_transform(level_three)
    labels_four = LabelEncoder().fit_transform(level_four)
    pango_arr_encoded = []
    for i in range(len(pango_arr)):
        variant_arr = []
        variant_arr.append(labels_one[i])
        variant_arr.append(labels_two[i])
        variant_arr.append(labels_three[i])
        variant_arr.append(labels_four[i])
        pango_arr_encoded.append(variant_arr)
    pango_encoding = {}
    for i in range(len(pango)):
        pango_encoding[pango[i]] = pango_arr_encoded[i]
    return pango_encoding


# In[25]:


pango_encoding = variant_to_vec(pango)


# In[26]:


y_all = []
for variant in dataset['pango_lineage'].values:
    y_all.append(pango_encoding[variant])


# In[27]:


X_train = np.array(X_all[:515])
X_test = np.array(X_all[515:])
y_train = np.array(y_all[:515])
y_test = np.array(y_all[515:])


# In[28]:


#arr1 = dataset['pango_lineage'].values


# In[29]:


#import numpy as np


# In[30]:


#np.where(arr1 == 'BA.1')


# In[1]:





# In[43]:





# In[45]:





# In[33]:





# In[34]:





# In[35]:





# In[36]:





# In[ ]:





# In[ ]:





# In[ ]:




