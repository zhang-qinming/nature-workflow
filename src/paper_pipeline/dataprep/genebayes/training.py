from __future__ import annotations

import argparse
import pickle
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import torch
import torch.distributions as dist
import xgboost as xgb
from ngboost import NGBRegressor
from ngboost.distns.distn import RegressionDistn
from ngboost.scores import LogScore
from torchquad import set_up_backend


torch.set_default_dtype(torch.float64)
set_up_backend("torch", data_type="float64")
print = partial(print, flush=True)

LR_INIT = 1e-2
N_METRIC = 1000
MAX_EPOCHS_INIT = 500
PATIENCE_INIT = 2
INTEGRATION_LB = -3.99
INTEGRATION_UB = 3.99
N_INTEGRATION_PTS = 1001
BATCH_SIZE = 1
DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


class PriorScore(LogScore):
    def score(self, y):
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
            result = -torch.sum(torch.log(result))
            result.backward()
            score += result.item()
            grad += torch.squeeze(params.grad, dim=-1)
            params.grad = None

        params = params.squeeze(dim=-1)
        grad[Prior.positive, :] *= params[Prior.positive, :]
        self.gradient = grad.T.cpu().detach().numpy()
        return score

    def d_score(self, data):
        for i in range(Prior.n_params):
            p0 = np.min(self.params_transf[i])
            p1 = np.quantile(self.params_transf[i], 0.01)
            p99 = np.quantile(self.params_transf[i], 0.99)
            p100 = np.max(self.params_transf[i])
            print(
                f"param #{i} - min: {p0}, 1st percentile: {p1}, "
                f"99th percentile: {p99}, max: {p100}"
            )

        return self.gradient

    def metric(self):
        all_grad = []
        for _ in range(N_METRIC):
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

    likelihood = dist.Normal(gene_property, y[idx, 1:]).log_prob(y[idx, :1])
    return (torch.exp(prior_p + likelihood), prior_p)


class MixGamma(dist.Distribution):
    def __init__(self, p, mu, sigma):
        self.bernoulli = dist.Bernoulli(logits=p)
        self.gamma = dist.Gamma(mu, sigma)

    def sample(self, sample_shape=()):
        sign = self.bernoulli.sample(sample_shape)
        sign = sign * 2 - 1
        magnitude = self.gamma.sample(sample_shape)
        return sign * magnitude

    def log_prob(self, value):
        sign = (value > 0).float()
        sign = self.bernoulli.log_prob(sign)
        magnitude = torch.abs(value)
        magnitude = self.gamma.log_prob(magnitude)
        return sign + magnitude


class Prior(RegressionDistn):
    n_params = 3
    positive = np.array([False, True, True])
    scores = [PriorScore]

    def __init__(self, params):
        params[0] = scipy.special.logit(scipy.special.expit(params[0]) * 0.99 + 0.005)
        self._params = params
        self.params_transf = np.copy(params)
        self.params_transf[Prior.positive] = np.exp(self.params_transf[Prior.positive])

    def distribution(params):
        return MixGamma(params[0], params[1], params[2])

    def fit(y):
        print("fitting initial distribution...")
        params = []
        for param in [0.0, 0.5, 5.0]:
            params.append(torch.tensor(param, requires_grad=True, device=DEVICE))
        optimizer = torch.optim.AdamW(params, lr=LR_INIT)

        for i in range(Prior.n_params):
            params[i] = params[i].expand(len(y), 1)

        lr_stage = 0
        min_loss = None
        min_epoch = 0

        for epoch in range(MAX_EPOCHS_INIT):
            loss = 0

            for idx in torch.split(torch.randperm(y.shape[0]), split_size_or_sections=BATCH_SIZE):
                optimizer.zero_grad()
                grid = torch.linspace(INTEGRATION_LB, INTEGRATION_UB, N_INTEGRATION_PTS).expand(len(idx), -1)
                eval_at_pts, _ = log_prob(Prior.distribution, idx, params, y, grid)
                result = torch.trapezoid(eval_at_pts, grid, dim=1)
                result = -torch.sum(torch.log(result))
                result.backward()
                optimizer.step()
                loss += result.item()

            print(
                "loss",
                loss,
                "params",
                torch.sigmoid(params[0][0]).item(),
                params[1][0].item(),
                params[2][0].item(),
                "lr",
                optimizer.param_groups[0]["lr"],
            )

            if min_loss is None or loss < min_loss:
                min_loss = loss
                min_epoch = epoch
            if epoch - min_epoch >= PATIENCE_INIT:
                optimizer.param_groups[0]["lr"] = optimizer.param_groups[0]["lr"] / 10
                lr_stage += 1

            if lr_stage > 1:
                break

        final_params = np.array([param[0].item() for param in params])
        print("initial params", final_params)
        final_params[Prior.positive] = np.log(final_params[Prior.positive])
        return final_params

    @property
    def params(self):
        return self.params_transf


def posterior(prob_y, idx, params, y, gene_property):
    prob_joint, _ = log_prob(Prior.distribution, idx, params, y, gene_property)
    return torch.exp(torch.log(prob_joint) - torch.log(prob_y))


def output_prior_posterior(model, feature_table, y, out_prefix: str, max_iter=None):
    print("calculating posterior distributions...")
    x_all = feature_table.to_numpy()
    params = torch.tensor(model.pred_dist(x_all, max_iter=max_iter)[:].params, device=DEVICE)
    print(params.shape)
    prior_mean = torch.mean(
        Prior.distribution(params).sample(torch.Size([10000])),
        dim=0,
    ).detach().cpu().numpy()

    post_mean = []
    post2_mean = []
    lower_95 = []
    upper_95 = []

    for idx in range(params.shape[1]):
        idx_tensor = torch.tensor([idx], device=DEVICE)
        grid = torch.linspace(INTEGRATION_LB, INTEGRATION_UB, N_INTEGRATION_PTS).expand(len(idx_tensor), -1)
        eval_at_pts, _ = log_prob(Prior.distribution, idx_tensor, params, y, grid)
        prob_y = torch.trapezoid(eval_at_pts, grid, dim=1)

        post_pdf = posterior(prob_y, idx_tensor, params, y, grid)
        post = torch.trapezoid(post_pdf * grid, grid, dim=1)
        post2 = torch.trapezoid(post_pdf * grid**2, grid, dim=1)
        post_mean.append(post.item())
        post2_mean.append(post2.item())

        grid_size = (INTEGRATION_UB - INTEGRATION_LB) / N_INTEGRATION_PTS
        cdf = torch.cumulative_trapezoid(post_pdf.flatten(), dx=grid_size, dim=0)
        print(cdf)
        lower = (
            (torch.argsort(torch.abs(cdf - 0.05))[0] + 1).detach().cpu().numpy() * grid_size
            + INTEGRATION_LB
        )
        upper = (
            (torch.argsort(torch.abs(cdf - 0.95))[0] + 1).detach().cpu().numpy() * grid_size
            + INTEGRATION_LB
        )
        lower_95.append(lower)
        upper_95.append(upper)

    params = params.detach().cpu().numpy()
    pd.DataFrame(
        {
            "ensg": feature_table.index,
            "param0": params[0],
            "param1": params[1],
            "param2": params[2],
            "prior_mean": prior_mean,
            "post_mean": post_mean,
            "post2_mean": post2_mean,
            "lower_95": lower_95,
            "upper_95": upper_95,
        }
    ).to_csv(f"{out_prefix}.per_gene_estimates.tsv", sep="\t", index=None)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--train_genes", dest="train_genes", required=False)
    parser.add_argument("--val_genes", dest="val_genes", required=False)
    parser.add_argument(
        "--gene_column",
        dest="gene_column",
        required=False,
        default="ensg",
        help="Name of the column containing gene names/IDs.",
    )
    parser.add_argument("--batch_size", dest="batch_size", type=int, default=1, required=False)
    parser.add_argument(
        "--n_integration_pts",
        dest="n_integration_pts",
        type=int,
        default=1001,
        required=False,
        help="Number of points for numerical integration.",
    )
    parser.add_argument(
        "--total_iterations",
        dest="total_iterations",
        type=int,
        default=500,
        required=False,
        help="Maximum number of iterations.",
    )
    parser.add_argument(
        "--early_stopping_iter",
        dest="early_stopping_iter",
        type=int,
        default=10,
        required=False,
        help="Validation patience for early stopping. Use 0 to disable.",
    )
    parser.add_argument("--lr", dest="lr", type=float, default=0.05, required=False)
    parser.add_argument("--max_depth", dest="max_depth", type=int, default=3, required=False)
    parser.add_argument(
        "--n_trees_per_iteration",
        dest="n_estimators",
        type=int,
        default=1,
        required=False,
    )
    parser.add_argument("--min_child_weight", dest="min_child_weight", type=float, required=False)
    parser.add_argument("--reg_alpha", dest="reg_alpha", type=float, required=False)
    parser.add_argument("--reg_lambda", dest="reg_lambda", type=float, required=False)
    parser.add_argument("--subsample", dest="subsample", type=float, required=False)

    required = parser.add_argument_group("required arguments")
    required.add_argument("--response", dest="response", required=True)
    required.add_argument("--features", dest="features", required=True)
    required.add_argument("--out", dest="out_prefix", required=True)
    required.add_argument("--integration_lb", dest="integration_lb", type=float, required=True)
    required.add_argument("--integration_ub", dest="integration_ub", type=float, required=True)
    parser.add_argument(
        "--model",
        dest="model",
        default=None,
        help="Provide a pretrained model to obtain predictions for that model.",
    )
    return parser


def _load_gene_subset(path: str | None, genes: list[str]) -> list[str]:
    if path is None:
        return genes
    with open(path, encoding="utf-8") as handle:
        return list(set([line.strip() for line in handle]) & set(genes))


def train_gene_bayes_model(args: argparse.Namespace) -> None:
    global LR_INIT, N_METRIC, MAX_EPOCHS_INIT, PATIENCE_INIT
    global INTEGRATION_LB, INTEGRATION_UB, N_INTEGRATION_PTS, BATCH_SIZE, DEVICE

    MAX_EPOCHS_INIT = 500
    LR_INIT = 1e-2
    PATIENCE_INIT = 2
    N_METRIC = 1000
    INTEGRATION_LB = args.integration_lb
    INTEGRATION_UB = args.integration_ub
    N_INTEGRATION_PTS = args.n_integration_pts
    BATCH_SIZE = args.batch_size
    DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    print(torch.cuda.is_available())
    print(DEVICE)

    x_all = pd.read_csv(args.features, sep="\t", index_col="ensg")
    y_all = pd.read_csv(args.response, sep=",", index_col="ensg")

    genes = list(set(x_all.index) & set(y_all.index))
    x_all = x_all.loc[genes]
    y_all = y_all.loc[genes]

    if args.train_genes is None:
        if args.early_stopping_iter > 0:
            raise ValueError("--train_genes is required when early stopping is enabled")
        train_genes = genes
    else:
        train_genes = _load_gene_subset(args.train_genes, genes)

    if args.val_genes is None:
        if args.early_stopping_iter > 0:
            raise ValueError("--val_genes is required when early stopping is enabled")
        val_genes: list[str] = []
    else:
        val_genes = _load_gene_subset(args.val_genes, genes)

    x_train = x_all.loc[train_genes].to_numpy()
    x_val = x_all.loc[val_genes].to_numpy()
    y_train = torch.tensor(y_all.loc[train_genes].to_numpy()).to(DEVICE)
    y_val = torch.tensor(y_all.loc[val_genes].to_numpy()).to(DEVICE)
    y_all_tensor = torch.tensor(y_all.to_numpy()).to(DEVICE)

    if args.model is not None:
        with open(args.model, "rb") as handle:
            model = pickle.load(handle)
    else:
        xgb_params = {
            "max_depth": args.max_depth,
            "reg_alpha": args.reg_alpha,
            "reg_lambda": args.reg_lambda,
            "min_child_weight": args.min_child_weight,
            "eta": 1.0,
            "subsample": args.subsample,
            "n_estimators": args.n_estimators,
        }
        xgb_params = {key: value for key, value in xgb_params.items() if value is not None}

        if torch.cuda.is_available():
            learner = xgb.XGBRegressor(gpu_id=0, tree_method="gpu_hist", **xgb_params)
        else:
            learner = xgb.XGBRegressor(tree_method="hist", **xgb_params)

        regressor = NGBRegressor(
            n_estimators=args.total_iterations,
            Dist=Prior,
            Base=learner,
            Score=PriorScore,
            verbose_eval=1,
            learning_rate=args.lr,
            natural_gradient=False,
        )
        if args.early_stopping_iter > 0:
            model = regressor.fit(
                x_train,
                y_train,
                X_val=x_val,
                Y_val=y_val,
                early_stopping_rounds=args.early_stopping_iter,
            )
        else:
            model = regressor.fit(x_all.to_numpy(), y_all_tensor)

        out_path = Path(f"{args.out_prefix}.model")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "wb") as handle:
            pickle.dump(model, handle)

    if model.best_val_loss_itr is not None:
        output_prior_posterior(model, x_all, y_all_tensor, args.out_prefix, model.best_val_loss_itr)
    else:
        output_prior_posterior(model, x_all, y_all_tensor, args.out_prefix)


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    print(vars(args))
    train_gene_bayes_model(args)
