from model import *
from utils import *
import pickle
import numpy as np
import torch
import torch.nn as nn
import random
import os
from torch.utils.data import DataLoader
from sklearn import preprocessing

def seed_torch(RANDOM_SEED=12345):
    random.seed(RANDOM_SEED)
    os.environ['PYTHONHASHSEED'] = str(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)
    torch.manual_seed(RANDOM_SEED)
    torch.cuda.manual_seed(RANDOM_SEED)
    torch.cuda.manual_seed_all(RANDOM_SEED)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


seed_torch()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def load_data(filename):
    data = []
    gene_names = []
    data_labels = []  # ?
    lines = open(filename).readlines()
    data_names = lines[0].replace('\n', '').split('\t')[1:]  # name of columns
    dx = 1

    for line in lines[dx:]:
        values = line.replace('\n', '').split('\t')
        gene = str.upper(values[0])
        gene_names.append(gene)
        data.append(values[1:])
    data = np.array(data, dtype='float32')
    data = np.transpose(data)

    return data, data_labels, data_names, gene_names

# tumor expression data
data_exp_tcga, data_labels_exp_tcga, data_names_exp_tcga, gene_names_exp_tcga = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/tcga_exp0424.txt")

batch_size = 100

tcga_exp_loader = DataLoader(
    data_exp_tcga, batch_size=batch_size)
tcga_exp_loader = DataLoader(
    data_exp_tcga, batch_size=batch_size)
tcga_exp_loader = DataLoader(
    data_exp_tcga, batch_size=batch_size)

num_DepOI = 601

model = DepPred_multi(
    inplanes=601, input_dim_e=16082, prop=0.3).to(device)

y_pred = []
y_true = []
# weight_score=[]
test_loss = []
test_coral_loss = []
test_corr = []
#model.load_state_dict(torch.load('../model/DepPredMulti_broad_20240423.pth'))
model.load_state_dict(torch.load('../model/DepPredMulti_broad_20240608.pth'))
model.eval()


y_tcga_pred = []
y_tcga_true = []
# weight_score=[]
model.load_state_dict(torch.load('../model/DepPredMulti_broad_20240608.pth'))
model.eval()

tcga1_exp, tcga2_exp, tcga3_exp = list(enumerate(tcga_exp_loader)), list(enumerate(tcga_exp_loader)), list(
    enumerate(tcga_exp_loader))
test_steps = min(len(tcga1_exp), len(tcga2_exp),
                 len(tcga3_exp))

for batch_idx in range(test_steps):
    _, tcga1_exp_test = tcga1_exp[batch_idx]
    _, tcga2_exp_test = tcga2_exp[batch_idx]
    _, tcga3_exp_test = tcga3_exp[batch_idx]

    #_, source1_depscore_test = source1_depscore3[batch_idx]
    #_, source2_depscore_test = source2_depscore3[batch_idx]
    #_, target_depscore_test = target_depscore[batch_idx]

    tcga1_exp_test = tcga1_exp_test.to(device)
    tcga2_exp_test = tcga2_exp_test.to(device)
    tcga3_exp_test = tcga3_exp_test.to(device)
    #source1_depscore_test = source1_depscore_test.to(device)
    #source2_depscore_test = source2_depscore_test.to(device)
    #target_depscore_test = target_depscore_test.to(device)

    tcga1_exp_test = tcga1_exp_test.to(dtype=torch.float32)
    tcga2_exp_test = tcga2_exp_test.to(dtype=torch.float32)
    tcga3_exp_test = tcga3_exp_test.to(dtype=torch.float32)

    out_source1, out_source2, out_target, loss_f, l2 = model(
        tcga1_exp_test,tcga2_exp_test, tcga3_exp_test)

    y_tcga_pred += list(out_target.cpu().detach().numpy())
    #print(y_tcga_pred)

y_tcga_pred = np.array(y_tcga_pred)
np.savetxt("../data/tcga_broad_predict_20240608.csv", y_tcga_pred, delimiter=",")
print('TCGA prediction is saved.')