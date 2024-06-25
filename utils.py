# %%
import numpy as np
import torch
import torch.nn as nn
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# %%


def CORAL(source, target):
    d = source.data.shape[1]

    # source covariance
    xm = torch.mean(source, 0, keepdim=True) - source
    xc = xm.t() @ xm

    # target covariance
    xmt = torch.mean(target, 0, keepdim=True) - target
    xct = xmt.t() @ xmt

    # frobenius norm between source and target
    loss = torch.mean(torch.mul((xc - xct), (xc - xct)))
    loss = loss/(4*d*d)

    return loss

def l2_distance(x, y):
    return torch.sqrt(torch.sum((x - y) ** 2))

# %%
def Corr(x, y):
    vx = x-torch.mean(x)
    vy = y-torch.mean(y)

    corr = torch.sum(vx * vy) / ((torch.sqrt(torch.sum(vx ** 2)))
                                 * torch.sqrt(torch.sum(vy ** 2)))

    return corr


def standardization(data):
    mu = np.mean(data, axis=0)
    sigma = np.std(data, axis=0)
    return (data-mu)/sigma


def add_noise(points, noise_type='gaussian'):
    if noise_type == 'gaussian':
        noise_factor = 0.1
        noisy_points = points + noise_factor * \
            torch.randn(*points.shape).to(device)
    else:
        raise ValueError(
            'Invalid noise type. Currently, only "gaussian" is supported for 2D data.')

    return torch.clamp(noisy_points, 0., 1.)


def pairwise_distance(x, y):

    if not len(x.shape) == len(y.shape) == 4:
        raise ValueError('Both inputs should be matrices.')

    if x.shape[1] != y.shape[1]:
        raise ValueError('The number of features should be the same.')

    #x = x.view(x.shape[0], x.shape[1], 1)
    #y = torch.transpose(y, 0, 1)
    output = torch.sum((x - y) ** 2, 1)
    output = torch.transpose(output, 0, 1)

    return output


def gaussian_kernel_matrix(x, y, sigmas):

    sigmas = sigmas.view(sigmas.shape[0], 1)
    beta = 1. / (2. * sigmas)
    dist = pairwise_distance(x, y).contiguous()
    dist_ = dist.view(1, -1)
    s = torch.matmul(beta, dist_)

    return torch.sum(torch.exp(-s), 0).view_as(dist)


def diverse_loss(feat1, feat2):
    loss = torch.sqrt(torch.sum((torch.mean(feat1 - feat2, 0)) ** 2))
    return loss


def agreement_loss(prob1, prob2):
    #return torch.abs(torch.mean(torch.sum(prob1 - prob2, -1)))
    return torch.mean(torch.sum(torch.abs(prob1 - prob2), -1))


def cross_entropy_loss(logits):
    prob = nn.functional.softmax(logits, 1)
    return torch.mean(torch.sum(- prob * torch.log(prob), 1))
