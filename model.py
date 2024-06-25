# %%
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import os
import random
from utils import *


# %%
RANDOM_SEED = 123


def seed_torch():
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


# %%

class Feature_extractor_multi(nn.Module):
    def __init__(self, input_dim_e, prop):
        super(Feature_extractor_multi, self).__init__()
        # expression data feature extraction
        self.efeature_extraction1 = nn.Sequential(
            nn.Dropout(prop),
            nn.Linear(input_dim_e, 10000),
            nn.BatchNorm1d(10000, affine=True),
            nn.ReLU(inplace=True))
        self.efeature_extraction2 = nn.Sequential(
            nn.Dropout(prop),
            nn.Linear(10000, 5000),
            nn.BatchNorm1d(5000, affine=True),
            nn.ReLU(inplace=True))

        self.efeature_extraction3 = nn.Sequential(
            nn.Dropout(prop),
            nn.Linear(5000, 601),
            nn.BatchNorm1d(601, affine=True),
            nn.ReLU(inplace=True))
        
    def forward(self, x1s1, x1s2, x1t):
        #loss_f0 = CORAL(x1s, x1t)

        out1s1 = self.efeature_extraction1(x1s1)
        out1s2 = self.efeature_extraction1(x1s2)
        out1t = self.efeature_extraction1(x1t)
        loss_f1 = CORAL(out1s1, out1t) + CORAL(out1s2, out1t) + CORAL(out1s1, out1s2)

        out1s1 = self.efeature_extraction2(out1s1)
        out1s2 = self.efeature_extraction2(out1s2)
        out1t = self.efeature_extraction2(out1t)
        loss_f2 = CORAL(out1s1, out1t) + CORAL(out1s2, out1t) + CORAL(out1s1, out1s2)


        out1s1 = self.efeature_extraction3(out1s1)
        out1s2 = self.efeature_extraction3(out1s2)
        out1t = self.efeature_extraction3(out1t)
        loss_f5 = CORAL(out1s1, out1t) + CORAL(out1s2, out1t) + CORAL(out1s1, out1s2)

        # merge exp and fprint models
        #loss_f = loss_f1 + loss_f2 + loss_f3 + loss_f4 + loss_f5
        loss_f = loss_f5
        source1 = out1s1
        source2 = out1s2
        target = out1t

        return source1, source2, target, loss_f


# %%


class Prediction_network_multi(nn.Module):
    def __init__(self, inplanes):
        super(Prediction_network_multi, self).__init__()

        #self.predict = make_res_layer(Bottleneck, inplanes, blocks)
        self.predict1 = nn.Sequential(
            nn.Linear(inplanes, inplanes)
            #nn.ReLU(inplace=True),
            #nn.Linear(inplanes, inplanes)
            )

        self.predict2 = nn.Sequential(
            nn.Linear(inplanes, inplanes)
            #nn.ReLU(inplace=True),
            #nn.Linear(inplanes, inplanes)
        )

    def forward(self, source1, source2, target):
        source1 = self.predict1(source1)
        source2 = self.predict2(source2)
        target1 = self.predict1(target)
        target2 = self.predict2(target)

        return source1, source2, target1, target2
    
# %%

class DepPred_multi(nn.Module):
    def __init__(self, inplanes, input_dim_e, prop):
        super(DepPred_multi, self).__init__()
        self.f_extraction = Feature_extractor_multi(
            input_dim_e, prop)
        self.sharedNet = Prediction_network_multi(inplanes)

    # x1s source expression data, x2s source fingerprint data, x1t target expression data, x2t target fingerprint data
    def forward(self, x1s1, x1s2, x1t):
        #coral_loss = 0
        source1, source2, target, loss_f = self.f_extraction(
            x1s1, x1s2, x1t)
        source1, source2, target1, target2 = self.sharedNet(source1, source2, target)

        out_target = (target1 + target2) / 2.0
        #out_target = target1

        l2 = l2_distance(target1, target2)

        part_loss = loss_f + l2


        return source1, source2, out_target, loss_f, l2
    
class DepPred_multi_ablation(nn.Module):
    def __init__(self, inplanes, input_dim_e, prop):
        super(DepPred_multi_ablation, self).__init__()
        self.f_extraction = Feature_extractor_multi_ablation(
            input_dim_e, prop)
        self.sharedNet = Prediction_network(inplanes)

    # x1s source expression data, x2s source fingerprint data, x1t target expression data, x2t target fingerprint data
    def forward(self, x1s, x1t):
        #coral_loss = 0
        source, target, loss_f = self.f_extraction(
            x1s, x1t)
        source, target = self.sharedNet(source, target)

        #out_target = (target1 + target2) / 2.0
        #output = self.fc(source)
        #output = self.output(source)

        return source, target, loss_f
