# ------------------------------------------------------------------------------
# --coding='utf-8'--
# Written by wissing
# ------------------------------------------------------------------------------
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

import numpy as np
import pandas as pd
import torch
import torch.optim as optim
import prettytable as pt

from networks import DeepSurv
from networks import NegativeLogLikelihood
from datasets import SurvivalDataset
from utils import read_config
from utils import c_index
from utils import adjust_learning_rate
from utils import create_logger

from lifelines import CoxPHFitter
from lifelines import NelsonAalenFitter
from lifelines import KaplanMeierFitter
from lifelines import WeibullFitter

import matplotlib.pyplot as plt


def test(ini_file):
    ''' Performs training according to .ini file

    :param ini_file: (String) the path of .ini file
    :return best_c_index: the best c-index
    '''
    # reads configuration from .ini file
    config = read_config(ini_file)
    # builds network|criterion|optimizer based on configuration
    model = DeepSurv(config['network']).to(device)
    criterion = NegativeLogLikelihood(config['network'], device).to(device)

    # cph = CoxPHFitter()
    # constructs data loaders based on configuration
    train_dataset = SurvivalDataset(config['train']['h5_file'], is_train=True, device=device)
    test_dataset = SurvivalDataset(config['train']['h5_file'], is_train=False, device=device)
    train_loader = torch.utils.data.DataLoader(
        train_dataset, batch_size=train_dataset.__len__())
    test_loader = torch.utils.data.DataLoader(
        test_dataset, batch_size=test_dataset.__len__())
    test_df = pd.read_csv(r'H:\project\DeepSurv\DeepSurv.pytorch-master\ours_test.csv', index_col=['PatientID'])
    train_df = pd.read_csv(r'H:\project\DeepSurv\DeepSurv.pytorch-master\ours_train.csv', index_col=['PatientID'])

    # train step
    best_c_index = 0
    # kmf = KaplanMeierFitter()
    naf = NelsonAalenFitter()
    # wf = WeibullFitter()
    naf.fit(test_df['Time_d'], event_observed=test_df['Event'])
    timeline = np.arange(0, 25000)
    base_risk = naf.predict(timeline)
    i = timeline[-1]
    while i > 0:
        base_risk[i] = base_risk[i] - base_risk[i-1]
        i -= 1
    np.savetxt('temp.txt', base_risk, '%.17f')
    # base_risk.to_csv('test_base_risk.csv', header=True)

    model.load_state_dict(torch.load(os.path.join(models_dir, ini_file.split('\\')[-1] + '.pth'))['model'])
    model.eval()

    for X, y, e in test_loader:
        with torch.no_grad():
            risk_pred = model(X)
            valid_loss = criterion(risk_pred, y, e, model)
            print(valid_loss)
            valid_c = c_index(-risk_pred, y, e)
            best_c_index = valid_c

            R = risk_pred.detach().cpu().numpy()[:, 0]
            for test_index in range(len(R)):
                # test_index = 120    # people
                _r = R[test_index]
                _y = y.detach().cpu().numpy()[test_index, 0]
                _e = e.detach().cpu().numpy()[test_index, 0]
                t0 = naf.predict(_y)

                risk = t0 * np.exp(_r)
                # print(np.exp(_r))
                print(risk, int(_e))

                # pre_y = 0.
                # m = np.min(np.where(p > 0.5))
                # print(int(_y), m, _e, p[int(_y)] >= 0.5)
                # if (p[int(_y)] >= 0.5) == bool(_e):
                #     ture += 1
            # print(ture/len(R))
            # if _e == pre_y:
            #     ture += 1
            # plt.plot(p)
            # plt.show()
    naf.plot()
    plt.show()
    return best_c_index


if __name__ == '__main__':
    # global settings
    logs_dir = 'logs'
    models_dir = os.path.join(logs_dir, 'models')
    if not os.path.exists(models_dir):
        os.makedirs(models_dir)
    logger = create_logger(logs_dir + '/test')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    configs_dir = 'configs'
    params = [
        ('Ours', 'ours.ini'),
        # ('Simulated Linear', 'linear.ini'),
        # ('Simulated Nonlinear', 'gaussian.ini'),
        # ('WHAS', 'whas.ini'),
        # ('SUPPORT', 'support.ini'),
        # ('METABRIC', 'metabric.ini'),
        # ('Simulated Treatment', 'treatment.ini'),
        # ('Rotterdam & GBSG', 'gbsg.ini')
    ]
    patience = 50
    # test
    headers = []
    values = []
    for name, ini_file in params:
        logger.info('Running {}({})...'.format(name, ini_file))
        fig = read_config(os.path.join(configs_dir, ini_file))
        t_fig = fig['train']
        n_fig = fig['network']
        logger.info("train: learning_rate = {:.6f}\nlr_decay_rate = {}\noptimizer = {}".format(t_fig['learning_rate'],
                                                                                               t_fig['lr_decay_rate'],
                                                                                               t_fig['optimizer']))
        logger.info("network: drop = {}\ndim = {}\nact = {}\nl2 = {}".format(n_fig['drop'], n_fig['dims'], n_fig['activation'],
                                                                             n_fig['l2_reg']))
        best_c_index = test(os.path.join(configs_dir, ini_file))
        headers.append(name)
        values.append('{:.6f}'.format(best_c_index))
        print('')
        logger.info("The best valid c-index: {}".format(best_c_index))
        logger.info('')
    # prints results
    tb = pt.PrettyTable()
    tb.field_names = headers
    tb.add_row(values)
    logger.info(tb)
