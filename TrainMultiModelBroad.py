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
from sklearn.metrics import r2_score

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

# broad expression data
data_exp_broad, data_labels_exp_broad, data_names_exp_broad, gene_names_exp_broad = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/broad_exp0424.txt")
# broad score data
data_depscore_broad, _, _, _ = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/broad_depscore0424.txt")

# sanger expression data
data_exp_sanger, data_labels_exp_sanger, data_names_exp_sanger, gene_names_exp_sanger = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/sanger_exp0424.txt")
# sanger score data
data_depscore_sanger, _, _, _ = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/sanger_depscore0424.txt")

# RNAi expression data
data_exp_rnai, data_labels_exp_rnai, data_names_exp_rnai, gene_names_exp_rnai = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/rnai_exp0424.txt")
# RNAi score data
data_depscore_rnai, _, _, _ = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/rnai_depscore0424.txt")

# tumor expression data
data_exp_tcga, data_labels_exp_tcga, data_names_exp_tcga, gene_names_exp_tcga = load_data(
    "/home/yushi/projects/def-hup-ab/yushi/data/tcga_exp0424.txt")

num_DepOI = 601
num_broad = data_exp_broad.shape[0]
num_sanger = data_exp_sanger.shape[0]
num_rnai = data_exp_rnai.shape[0]

# 80% Sanger CCLs for training, and 20% for validation,
id_rand = np.random.permutation(num_broad)  # randomly distribute CCL number
id_cell_train = id_rand[np.arange(0, round(num_broad*0.7))]
id_cell_val = id_rand[np.arange(round(num_broad*0.7), num_broad)]
#id_cell_test = id_rand[np.arange(round(num_sanger*0.85), num_sanger)]


exp_broad_train = data_exp_broad[id_cell_train]
exp_broad_val = data_exp_broad[id_cell_val]
#exp_sanger_test = data_exp_sanger[id_cell_test]
exp_broad_test = data_exp_broad

depscore_broad_train = data_depscore_broad[id_cell_train]
depscore_broad_val = data_depscore_broad[id_cell_val]
#depscore_sanger_test = data_depscore_sanger[id_cell_test]
depscore_broad_test = data_depscore_broad


learning_rate = 1e-3
momentum = 0.9

#sanger as target (sanger, rnai, broad)
#batch_size1 = [162, 20, 202]
#batch_size2 = [162, 20, 50]
#batch_size3 = [162, 20, 252]

batch_size1 = [162, 20, 176]
batch_size2 = [162, 20, 76]
batch_size3 = [200, 200, 200]

epoch = 3000
decay = 1e-6

target_exp_train_loader = DataLoader(
    exp_broad_train, batch_size=batch_size1[2])
target_exp_val_loader = DataLoader(
    exp_broad_val, batch_size=batch_size2[2])
target_exp_test_loader = DataLoader(
    exp_broad_test, batch_size=batch_size3[2])

target_depscore_train_loader = DataLoader(
    depscore_broad_train, batch_size=batch_size1[2])
target_depscore_val_loader = DataLoader(
    depscore_broad_val, batch_size=batch_size2[2])
target_depscore_test_loader = DataLoader(
    depscore_broad_test, batch_size=batch_size3[2])

source1_exp_train_loader = DataLoader(
    data_exp_sanger, batch_size=batch_size1[0])
source1_exp_val_loader = DataLoader(
    data_exp_sanger, batch_size=batch_size2[0])
source1_exp_test_loader = DataLoader(
    data_exp_sanger, batch_size=batch_size3[0])

source1_depscore_train_loader = DataLoader(
    data_depscore_sanger, batch_size=batch_size1[0])
source1_depscore_val_loader = DataLoader(
    data_depscore_sanger, batch_size=batch_size2[0])
source1_depscore_test_loader = DataLoader(
    data_depscore_sanger, batch_size=batch_size3[0])

source2_exp_train_loader = DataLoader(
    data_exp_rnai, batch_size=batch_size1[1])
source2_exp_val_loader = DataLoader(
    data_exp_rnai, batch_size=batch_size2[1])
source2_exp_test_loader = DataLoader(
    data_exp_rnai, batch_size=batch_size3[1])

source2_depscore_train_loader = DataLoader(
    data_depscore_rnai, batch_size=batch_size1[1])
source2_depscore_val_loader = DataLoader(
    data_depscore_rnai, batch_size=batch_size2[1])
source2_depscore_test_loader = DataLoader(
    data_depscore_rnai, batch_size=batch_size3[1])

model = DepPred_multi(
    inplanes=601, input_dim_e=16082, prop=0.3).to(device)

loss_fn = torch.nn.MSELoss()  # reduction=none 返回每个sample mse
# support different learning rate according to CORAL paper
# i.e. 10 times learning rate for the last two fc layers.
optimizer = torch.optim.Adam(model.parameters(
), lr=learning_rate)  # momentum parameters

source1_exp1, source2_exp1, target_exp1, source1_depscore1, source2_depscore1, target_depscore1 = list(
    enumerate(source1_exp_train_loader)), list(enumerate(source2_exp_train_loader)), list(enumerate(target_exp_train_loader)), list(enumerate(source1_depscore_train_loader)), list(enumerate(source2_depscore_train_loader)), list(enumerate(target_depscore_train_loader))

# print(len(source_exp_num), len(source_fprint_num))
train_steps = min(len(source1_exp1), len(source2_exp1), len(target_exp1))

min_loss = 1e99
max_corr = 0
min_epoch = 0
final_lambda = 0
output_loss = []
output_mse_loss = []
output_coral_loss = []
train_loss = []
train_mse_loss = []
train_coral_loss = []
for e in range(0, epoch):
    epoch_loss = []
    epoch_tot_loss = []
    #epoch_coral_loss = []
    epoch_corr = []
    val_loss = []
    val_tot_loss = []
    val_corr = []
    #val_coral_loss = []
    val_tot_corr = []
    model.train()
    _lambda = (e+1)/epoch
    #_lambda = 0.1
    for batch_idx in range(train_steps):
        _, source1_exp_train = source1_exp1[batch_idx]
        _, source2_exp_train = source2_exp1[batch_idx]

        _, target_exp_train = target_exp1[batch_idx]
        #np.savetxt("../data/train-data.csv", source_exp_train, delimiter=",")
        #np.savetxt("../data/test-data.csv", target_exp_train, delimiter=",")
        #print(CORAL(source_exp_train, target_exp_train))

        _, source1_depscore_train = source1_depscore1[batch_idx]
        _, source2_depscore_train = source2_depscore1[batch_idx]
        _, target_depscore_train = target_depscore1[batch_idx]
 
        source1_exp_train = source1_exp_train.to(device)
        source2_exp_train = source2_exp_train.to(device)
        target_exp_train = target_exp_train.to(device)
        source1_depscore_train = source1_depscore_train.to(device)
        source2_depscore_train = source2_depscore_train.to(device)
        target_depscore_train = target_depscore_train.to(device)

        source1_exp_train = source1_exp_train.to(dtype=torch.float32)
        source2_exp_train = source2_exp_train.to(dtype=torch.float32)
        target_exp_train = target_exp_train.to(dtype=torch.float32)

        optimizer.zero_grad()
        # print(batch_idx)
        out_source1, out_source2, out_target, loss_f, l2 = model(
            source1_exp_train, source2_exp_train, target_exp_train)

        mse_loss1 = loss_fn(out_source1.float(), source1_depscore_train.float())
        mse_loss2 = loss_fn(out_source2.float(), source2_depscore_train.float())
        corr = []
        for i in range(0, len(out_target)):
            corr.append(Corr(out_target[i] , target_depscore_train[i]))
        corr = torch.tensor(corr, dtype=torch.float32)
        corr_coeff = corr.mean()
        #mse_loss = mse_loss1+mse_loss2+mse_loss3+mse_loss4+mse_loss5
        loss = 0.9*(mse_loss1+mse_loss2)+_lambda*(loss_f+l2)
        #loss = _lambda*preprocessing.normalize(mse_loss.item()) + (1-_lambda)*preprocessing.normalize(coral_loss.item())
        loss.backward()
        optimizer.step()
        #epoch_loss.append(mse_loss.item())
        epoch_tot_loss.append(loss.item())
        epoch_loss.append(loss_f.item())
        epoch_corr.append(corr_coeff.item())

    source1_exp2, source2_exp2, target_exp2, source1_depscore2, source2_depscore2, target_depscore2= list(enumerate(source1_exp_val_loader)), list(enumerate(source2_exp_val_loader)), list(
        enumerate(target_exp_val_loader)), list(enumerate(source1_depscore_val_loader)), list(enumerate(source2_depscore_val_loader)), list(enumerate(target_depscore_val_loader))
    val_steps = min(len(source1_exp2), len(source2_exp2), len(target_exp2))

    # validation mode
    model.eval()
    y_pred = []
    y_true = []
    with torch.no_grad():
        for batch_idx in range(val_steps):
            _, source1_exp_val = source1_exp2[batch_idx]
            _, source2_exp_val = source2_exp2[batch_idx]
            _, target_exp_val = target_exp2[batch_idx]
            _, source1_depscore_val = source1_depscore2[batch_idx]
            _, source2_depscore_val = source2_depscore2[batch_idx]
            _, target_depscore_val = target_depscore2[batch_idx]

            source1_exp_val = source1_exp_val.to(device)
            source2_exp_val = source2_exp_val.to(device)
            target_exp_val = target_exp_val.to(device)
            source1_depscore_val = source1_depscore_val.to(device)
            source2_depscore_val = source2_depscore_val.to(device)
            target_depscore_val = target_depscore_val.to(device)
            
            source1_exp_val = source1_exp_val.to(dtype=torch.float32)
            source2_exp_val = source2_exp_val.to(dtype=torch.float32)
            target_exp_val = target_exp_val.to(dtype=torch.float32)
            
            out_source1, out_source2, out_target, loss_f, l2 = model(
                source1_exp_val, source2_exp_val, target_exp_val)

            mse_loss1 = loss_fn(out_source1.float(), source1_depscore_val.float())
            mse_loss2 = loss_fn(out_source2.float(), source2_depscore_val.float())

            corr = []
            for i in range(len(out_target)):
                corr.append(Corr(out_target[i] , target_depscore_val[i]))
            corr = torch.tensor(corr, dtype=torch.float32)
            corr_coeff = corr.mean() 

            loss = 0.9*(mse_loss1+mse_loss2)+final_lambda*(loss_f+l2)
            val_tot_loss.append(loss.item())
            val_corr.append(corr_coeff)
            val_tot_corr.append(Corr(out_source1.flatten(), source1_depscore_val.flatten()))
            #val_coral_loss.append(coral_loss.item())
            #y_pred += list(output.squeeze().cpu().detach().numpy())
            #y_true += list(source_depscore_val.squeeze().cpu().detach().numpy())
    #y_pred = np.array(y_pred)  # +=
    #y_true = np.array(y_true)
    val_tot_corr = torch.tensor(val_tot_corr, dtype=torch.float32)
    print('Epoch: %d,  Train Loss: %.4f, Train Part Loss: %.4f, Train Correlation: %.4f, Val Loss: %.4f, Val Correlation: %.4f' % (
        e, np.mean(epoch_tot_loss), np.mean(epoch_loss), np.mean(epoch_corr), np.mean(val_tot_loss), np.mean(val_corr)))
    train_loss.append(np.mean(epoch_tot_loss))
    train_mse_loss.append(np.mean(epoch_loss))
    #train_coral_loss.append(np.mean(epoch_coral_loss))
    output_loss.append(np.mean(val_tot_loss))
    output_mse_loss.append(np.mean(val_loss))
    #output_coral_loss.append(np.mean(val_coral_loss))
    if np.mean(val_corr) > max_corr:
        #min_loss = np.mean(val_tot_loss)
        max_corr = np.mean(val_corr)
        min_epoch = e
        final_lambda = _lambda
        torch.save(model.state_dict(), '../model/DepPredMulti_broad_20240608.pth')

y_pred = []
y_true = []
# weight_score=[]
test_loss = []
test_coral_loss = []
test_corr = []
model.load_state_dict(torch.load('../model/DepPredMulti_broad_20240608.pth'))
model.eval()


source1_exp3, source2_exp3, target_exp3, target_depscore = list(enumerate(target_exp_test_loader)), list(enumerate(target_exp_test_loader)), list(
    enumerate(target_exp_test_loader)), list(enumerate(target_depscore_test_loader))
test_steps = min(len(source1_exp3), len(source2_exp3),
                 len(target_exp3))

tot_corr = []
for batch_idx in range(test_steps):
    _, source1_exp_test = source1_exp3[batch_idx]
    _, source2_exp_test = source2_exp3[batch_idx]

    _, target_exp_test = target_exp3[batch_idx]

    #_, source1_depscore_test = source1_depscore3[batch_idx]
    #_, source2_depscore_test = source2_depscore3[batch_idx]
    _, target_depscore_test = target_depscore[batch_idx]

    source1_exp_test = source1_exp_test.to(device)
    source2_exp_test = source2_exp_test.to(device)
    target_exp_test = target_exp_test.to(device)
    #source1_depscore_test = source1_depscore_test.to(device)
    #source2_depscore_test = source2_depscore_test.to(device)
    target_depscore_test = target_depscore_test.to(device)

    source1_exp_test = source1_exp_test.to(dtype=torch.float32)
    source2_exp_test = source2_exp_test.to(dtype=torch.float32)
    target_exp_test = target_exp_test.to(dtype=torch.float32)

    out_source1, out_source2, out_target, loss_f, l2 = model(
        source1_exp_test, source2_exp_test, target_exp_test)
    mse_loss = loss_fn(out_target.float(), target_depscore_test.float())
    
    corr = []
    for i in range(len(out_target)):
        corr.append((Corr(out_target[i] , target_depscore_test[i])))
    corr = torch.tensor(corr, dtype=torch.float32)
    corr_coeff = corr.mean()
    #coral_loss = coral_loss.to(dtype=torch.float64)   
    tot_corr.append(Corr(out_target.flatten(), target_depscore_test.flatten()))

    loss = 0.9*mse_loss+final_lambda*(loss_f+l2)
    test_loss.append(mse_loss.item())
    #test_coral_loss.append(coral_loss.item())
    test_corr.append(corr_coeff.item())
    y_pred += list(out_target.cpu().detach().numpy())
    y_true += list(target_depscore_test.cpu().detach().numpy())
y_pred = np.array(y_pred)
y_true = np.array(y_true)
#tot_corr = np.corrcoef(y_pred, y_true)[0, 1]
#print(Corr(y_pred, y_true))
#per_corr_coeff = per_corr.mean()
tot_corr = torch.tensor(tot_corr, dtype=torch.float32)
print('Test Loss: %.4f,  Test Correlation: %.4f, Total Correlation: %.4f, Final Lambda: %.4f' %
      (np.mean(test_loss), np.mean(test_corr), tot_corr.mean(), final_lambda))
np.savetxt("../data/broad_true_20240608.csv", y_true, delimiter=",")
np.savetxt("../data/broad_predict_20240608.csv", y_pred, delimiter=",")
print(r2_score(y_true.flatten('F'), y_pred.flatten('F')))