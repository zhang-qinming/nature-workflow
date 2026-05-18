import pickle
import numpy as np
import pandas as pd
import argparse
from functools import partial
import scipy

import xgboost as xgb
from ngboost import NGBRegressor
from ngboost.distns.distn import RegressionDistn
from ngboost.scores import LogScore

import torch
import torch.distributions as dist
import torch.nn.functional as F

from torchquad import Boole, set_up_backend

import torch.utils.data as data


torch.set_default_dtype(torch.float64)
set_up_backend("torch", data_type="float64")
print = partial(print, flush=True)


class PriorScore(LogScore):
    def score(self, y):
        """
        log score, -log P(y), for prior distribution P
        """
        params = torch.tensor(self.params_transf, device=DEVICE)
        params = params.unsqueeze(dim=-1)
        params.requires_grad = True

        score = 0
        grad = torch.zeros(self.n_params, y.shape[0], device=DEVICE)
        y = torch.tensor(y, device=DEVICE)

        for idx in torch.split(torch.randperm(y.shape[0]), split_size_or_sections=BATCH_SIZE):
            grid = torch.linspace(INTEGRATION_LB, INTEGRATION_UB, N_INTEGRATION_PTS).expand(len(idx), -1)
            eval_at_pts, _ = log_prob(Prior.distribution, idx, params, y, grid)
            result = torch.trapezoid(eval_at_pts, grid, dim=1)
            result = -torch.sum(torch.log(result))#+torch.max(y[idx],dim=1)[0])
            result.backward()
            score += result.item()
            grad += torch.squeeze(params.grad, dim=-1)
            params.grad = None

        params = params.squeeze(dim=-1)
        grad[Prior.positive, :] *= params[Prior.positive, :]
        self.gradient = grad.T.cpu().detach().numpy()

        return score

    def d_score(self, data):
        """
        derivative of the score
        """
        for i in range(Prior.n_params):
            p0 = np.min(self.params_transf[i])
            p1 = np.quantile(self.params_transf[i], 0.01)
            p99 = np.quantile(self.params_transf[i], 0.99)
            p100 = np.max(self.params_transf[i])
            print(f"param #{i} - min: {p0}, 1st percentile: {p1}, "
                  f"99th percentile: {p99}, max: {p100}")

        return self.gradient

    def metric(self):
        all_grad = []
        for i in range(N_METRIC):
            params = torch.tensor(self.params_transf, device=DEVICE)
            params.requires_grad = True

            gene_property = Prior.distribution(params).sample()
            loss = Prior.distribution(params).log_prob(gene_property)
            loss = -torch.sum(loss)
            loss.backward()

            params.grad[Prior.positive, :] *= params[Prior.positive, :]
            all_grad.append(params.grad.T.detach().cpu())
            params.grad = None

        grad = np.stack(all_grad)
        grad = np.mean(np.einsum("sik,sij->sijk", grad, grad), axis=0)

        return grad


def log_prob(prior_dist, idx, params, y, gene_property):
    prior_p = prior_dist([p[idx] for p in params]).log_prob(gene_property)
    if torch.trapezoid(torch.exp(prior_p), gene_property)[0] < 0.95:
        print("integration issue", torch.trapezoid(torch.exp(prior_p), gene_property)[0])

    lh = dist.Normal(gene_property, y[idx,1:]).log_prob(y[idx,:1])
    #y[idx]-torch.max(y[idx],dim=1,keepdim=True)[0]

    return (torch.exp(prior_p+lh),prior_p)


class MixGamma(dist.Distribution):
    def __init__(self, p, mu, sigma):
        self.bernoulli = dist.Bernoulli(logits=p)
        #self.gamma = dist.LogNormal(mu, sigma)
        #self.gamma = dist.TransformedDistribution(
        #                 dist.Normal(mu, sigma),
        #                 transforms=[dist.SigmoidTransform(), dist.AffineTransform(loc=0.,scale=2.)])
        #self.gamma = dist.TransformedDistribution(
        #                 dist.Beta(mu, sigma),
        #                 transforms=[dist.AffineTransform(loc=0.,scale=2.)])
        self.gamma = dist.Gamma(mu, sigma)

    def sample(self, sample_shape=()):
        sign = self.bernoulli.sample(sample_shape)
        sign = sign*2-1
        magnitude = self.gamma.sample(sample_shape)
        #return magnitude
        return sign*magnitude

    def log_prob(self, value):
        sign = (value>0).float()
        sign = self.bernoulli.log_prob(sign)
        magnitude = torch.abs(value)
        magnitude = self.gamma.log_prob(magnitude)

        #return magnitude
        return sign+magnitude


class Prior(RegressionDistn):
    n_params = 3
    positive = np.array([False, True, True])
    scores = [PriorScore]

    def __init__(self, params):
        params[0] = scipy.special.logit(scipy.special.expit(params[0])*0.99+0.005)
        self._params = params
        self.params_transf = np.copy(params)
        self.params_transf[Prior.positive] = np.exp(self.params_transf[Prior.positive])

    def distribution(params):
        return MixGamma(params[0],params[1],params[2])

    def fit(y):
        """
        fit initial prior distribution for all genes
        """
        print("fitting initial distribution...")
        params = []
        # 0, 0.1
        for param in [0.,0.5,5.]:
            params.append(torch.tensor(param, requires_grad=True, device=DEVICE))
        optimizer = torch.optim.AdamW(params, lr=LR_INIT)

        for i in range(Prior.n_params):
            params[i] = params[i].expand(len(y),1)

        lr_stage = 0

        for i in range(MAX_EPOCHS_INIT):
            loss = 0

            for idx in torch.split(torch.randperm(y.shape[0]), split_size_or_sections=BATCH_SIZE):
                optimizer.zero_grad()
                grid = torch.linspace(INTEGRATION_LB, INTEGRATION_UB, N_INTEGRATION_PTS).expand(len(idx), -1)
                eval_at_pts, prior_p = log_prob(Prior.distribution, idx, params, y, grid)
                result = torch.trapezoid(eval_at_pts, grid, dim=1)
                result = -torch.sum(torch.log(result))#+torch.max(y[idx],dim=1)[0])
                result.backward()
                optimizer.step()
                loss += result.item()

            print("loss", loss,
                  "params", torch.sigmoid(params[0][0]).item(), params[1][0].item(), params[2][0].item(),
                  "lr", optimizer.param_groups[0]['lr'])

            if i == 0 or loss < min_loss:
                min_loss = loss
                min_epoch = i
            if i - min_epoch >= PATIENCE_INIT:
                optimizer.param_groups[0]['lr'] = optimizer.param_groups[0]['lr'] / 10
                lr_stage += 1

            if lr_stage > 1:
                break
        
        params = [p[0].item() for p in params]
        print("initial params", params)
        params = np.array(params)
        params[Prior.positive] = np.log(params[Prior.positive])

        return params

    @property
    def params(self):
        return self.params_transf


def posterior(prob_y, idx, params, y, gene_property):
    '''
    return posterior pdf, p(y,gene_property)/p(y)
    '''
    prob_joint, _ = log_prob(Prior.distribution, idx, params, y, gene_property)
    return torch.exp(torch.log(prob_joint) - torch.log(prob_y))


def output_prior_posterior(model, feature_table, y, max_iter=None):
    print("calculating posterior distributions...")
    X = feature_table.to_numpy()
    params = torch.tensor(model.pred_dist(X, max_iter=max_iter)[:].params,
                          device=DEVICE)
    print(params.shape)
    # prior
    prior_mean = torch.mean(Prior.distribution(params).sample(torch.Size([10000])),
                            dim=0).detach().cpu().numpy()

    # posterior
    post_mean = []
    post2_mean = []
    lower_95, upper_95 = [], []

    for idx in range(params.shape[1]):
        idx = torch.tensor([idx], device=DEVICE)

        # compute p(y)
        grid = torch.linspace(INTEGRATION_LB, INTEGRATION_UB, N_INTEGRATION_PTS).expand(len(idx), -1)
        eval_at_pts, prior_p = log_prob(Prior.distribution, idx, params, y, grid)
        prob_y = torch.trapezoid(eval_at_pts, grid, dim=1)

        # compute posterior pdf + expected posterior gene_property
        post_pdf = posterior(prob_y, idx, params, y, grid)
        post = torch.trapezoid(post_pdf*grid, grid, dim=1)
        post2 = torch.trapezoid(post_pdf*grid**2, grid, dim=1)
        post_mean.append(post.item())
        post2_mean.append(post2.item())

        # compute lower/upper bounds
        grid_size = (INTEGRATION_UB - INTEGRATION_LB) / N_INTEGRATION_PTS
        cdf = torch.cumulative_trapezoid(post_pdf.flatten(), dx=grid_size, dim=0)
        print(cdf)
        lower = (torch.argsort(torch.abs(cdf - 0.05))[0] + 1).detach().cpu().numpy() * grid_size + INTEGRATION_LB
        upper = (torch.argsort(torch.abs(cdf - 0.95))[0] + 1).detach().cpu().numpy() * grid_size + INTEGRATION_LB
        lower_95.append(lower)
        upper_95.append(upper)

    params = params.detach().cpu().numpy()

    output = pd.DataFrame({"ensg":feature_table.index,
                           "param0":params[0],
                           "param1":params[1],
                           "param2":params[2],
                           "prior_mean":prior_mean,
                           "post_mean":post_mean,
                           "post2_mean":post2_mean,
                           "lower_95":lower_95,
                           "upper_95":upper_95})
    output.to_csv(args.out_prefix + ".per_gene_estimates.tsv", sep='\t', index=None)

    ### feature importance metrics ###
    #features = X.columns
    #importance = {"feature": features}
    #for i in range(model.feature_importances_.shape[0]):
    #    importance["param%s_importance" % i] = model.feature_importances_[i]
    #importance = pd.DataFrame(importance)
    #importance.to_csv(args.out_prefix + ".feature_importance.tsv", sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--train_genes", dest="train_genes", required=False)
    parser.add_argument("--val_genes", dest="val_genes", required=False)
    parser.add_argument("--gene_column", dest="gene_column", required=False, default="ensg",
                        help="Name of the column containing gene names/IDs.")
    #parser.add_argument("--chrom_column", dest="chrom_column", required=False, default="chrom",
    #                    help="Name of the column containing chromosome number.")
    parser.add_argument("--batch_size", dest="batch_size", type=int, default=1, required=False)
    parser.add_argument("--n_integration_pts", dest="n_integration_pts", type=int, default=1001, required=False,
                        help="Number of points for numerical integration. Larger values increase can improve accuracy but also increase training time and/or numerical instability.")
    parser.add_argument("--total_iterations", dest="total_iterations", type=int, default=500, required=False,
                        help="Maximum number of iterations. The actual number of iterations may be lower if using early stopping (see the 'train' option)")
    parser.add_argument("--early_stopping_iter", dest="early_stopping_iter", type=int, default=10, required=False,
                        help="If >0, chromosomes 2, 4, 6 will be held out for validation and training will end when the loss on these chromosomes stops decreasing for the specified number of iterations. Otherwise, the model will train on all genes for the number of iterations specified in 'total_iterations'.")
    parser.add_argument("--lr", dest="lr", type=float, default=0.05, required=False,
                        help="Learning rate for NGBoost. Smaller values may improve accuracy but also increase training time. Typical values: 0.01 to 0.2.")
    parser.add_argument("--max_depth", dest="max_depth", type=int, default=3, required=False,
                        help="XGBoost parameter. See https://xgboost.readthedocs.io/en/stable/python/python_api.html.")
    parser.add_argument("--n_trees_per_iteration", dest="n_estimators", type=int, default=1, required=False,
                        help="XGBoost parameter, n_estimators")
    parser.add_argument("--min_child_weight", dest="min_child_weight", type=float, required=False,
                        help="XGBoost parameter.")
    parser.add_argument("--reg_alpha", dest="reg_alpha", type=float, required=False,
                        help="XGBoost parameter.")
    parser.add_argument("--reg_lambda", dest="reg_lambda", type=float, required=False,
                        help="XGBoost parameter.")
    parser.add_argument("--subsample", dest="subsample", type=float, required=False,
                        help="XGBoost parameter.")

    required = parser.add_argument_group("required arguments")
    required.add_argument("--response", dest="response", required=True,
                          help="tsv containing data that can be related to the gene property of interest through a likelihood function. Also known as y or dependent variable.")
    required.add_argument("--features", dest="features", required=True,
                          help="tsv containing gene features. Also known as X or independent variables.")
    required.add_argument("--out", dest="out_prefix",
                          help="Prefix for the output files.")
    required.add_argument("--integration_lb", dest="integration_lb", type=float,
                          help="Lower bound for numerical integration - the smallest value that you expect for the gene property of interest.")
    required.add_argument("--integration_ub", dest="integration_ub", type=float,
                          help="Upper bound for numerical integration - the largest value you expect for the gene property of interest.")
    parser.add_argument("--model", dest="model", default=None, 
                        help="Provide a pretrained model to obtain predictions for that model.")
    #parser.add_argument("--lr_init", dest="lr_init", type=float, default=1e-2, )
    #parser.add_argument("--patience_init", dest="patience_init", type=int, default=1)
    #parser.add_argument("--max_epochs_init", dest="max_epochs_init", type=int, default=100)
    #parser.add_argument("--n_metric", dest="n_metric", type=int, default=1000)
    #parser.add_argument("--gpu_id", dest="gpu_id")

    args = parser.parse_args()
    print(vars(args))

    global LR_INIT, N_METRIC, MAX_EPOCHS_INIT, PATIENCE_INIT, INTEGRATION_LB, INTEGRATION_UB
    global GENES, N_CORES
    global DEVICE

    MAX_EPOCHS_INIT = 500 #args.max_epochs_init
    LR_INIT = 1e-2 #args.lr_init
    PATIENCE_INIT = 2 #args.patience_init
    N_METRIC = 1000 #args.n_metric
    INTEGRATION_LB = args.integration_lb
    INTEGRATION_UB = args.integration_ub
    N_INTEGRATION_PTS = args.n_integration_pts
    GENE_COLUMN = args.gene_column
    #CHROM_COLUMN = args.chrom_column
    BATCH_SIZE = args.batch_size
    N_CORES = 4
    #GPU_ID = args.gpu_id
   
    print(torch.cuda.is_available())
    DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(DEVICE)

    ### load likelihoods and training data ###
    X_all = pd.read_csv(args.features, sep='\t', index_col="ensg")
    y_all = pd.read_csv(args.response, sep=',', index_col="ensg")

    genes = list(set(X_all.index) & set(y_all.index))

    X_all = X_all.loc[genes]
    y_all = y_all.loc[genes]

    train_genes = list(set([x.strip() for x in open(args.train_genes)]) & set(genes))
    val_genes = list(set([x.strip() for x in open(args.val_genes)]) & set(genes))

    X_train = X_all.loc[train_genes]
    X_val = X_all.loc[val_genes]

    X_train = X_train.to_numpy()
    X_val = X_val.to_numpy()

    y_train = y_all.loc[train_genes]
    y_val = y_all.loc[val_genes]

    y_train = torch.tensor(y_train.to_numpy()).to(DEVICE)
    y_val = torch.tensor(y_val.to_numpy()).to(DEVICE)
    y_all = torch.tensor(y_all.to_numpy()).to(DEVICE)

    if args.model is not None:  # already have pretrained model, want to obtain predictions
        model = pickle.load(open(args.model, 'rb'))
    else:  # train model and make predictions
        xgb_params = {"max_depth": args.max_depth,
                      "reg_alpha": args.reg_alpha,
                      "reg_lambda": args.reg_lambda,
                      "min_child_weight": args.min_child_weight,
                      "eta": 1.0,
                      "subsample": args.subsample,
                      "n_estimators": args.n_estimators}
        xgb_params = {k: v for k, v in xgb_params.items() if v is not None}

        if not torch.cuda.is_available():
            learner = xgb.XGBRegressor(
                tree_method="hist",
                **xgb_params
            )
        else:
            learner = xgb.XGBRegressor(
                gpu_id=0,
                tree_method="gpu_hist",
                **xgb_params
            )

        if args.early_stopping_iter>0:
            model = NGBRegressor(n_estimators=args.total_iterations, Dist=Prior, Base=learner, Score=PriorScore,
                                 verbose_eval=1, learning_rate=args.lr, natural_gradient=False
                                 ).fit(X_train, y_train, X_val=X_val, Y_val=y_val, early_stopping_rounds=args.early_stopping_iter)
        else:
            model = NGBRegressor(n_estimators=args.total_iterations, Dist=Prior, Base=learner, Score=PriorScore,
                                 verbose_eval=1, learning_rate=args.lr, natural_gradient=False
                                 ).fit(X, y)

        f = open(args.out_prefix + ".model", 'wb')
        pickle.dump(model, f)
        f.close()

    if model.best_val_loss_itr is not None:
        output_prior_posterior(model, X_all, y_all, model.best_val_loss_itr)
    else:
        output_prior_posterior(model, X_all, y_all)
