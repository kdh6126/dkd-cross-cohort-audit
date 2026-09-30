#!/usr/bin/env python
"""Generative robustness: does the signature survive on a synthetic population?

The parametric perturbations in stress_test.py (noise, masking, batch shift, subsampling) each
probe one axis. A generative model probes the joint distribution: sample a synthetic patient
population and ask whether the biomarker still discriminates, and whether the same genes would
still be chosen.

Two generators, both class-conditional and both trained on TRAINING cohorts only:
  VAE        conditional variational autoencoder, gaussian decoder
  DIFFUSION  small DDPM with an MLP denoiser, class-conditioned

Honest limitation, stated up front: these are trained on ~120 samples. A generative model at
that sample size cannot be trusted to extrapolate; it mostly re-expresses the training
covariance. So this is a CONSISTENCY check - "does the signature survive a resampling of the
learned joint distribution" - not evidence about real unseen patients. Three guards are
reported alongside so the reader can judge how much the synthetic data is worth:
  * nearest-neighbour distance from synthetic to real (memorisation check)
  * correlation of the per-gene CLASS-MEAN DIFFERENCE, real vs synthetic (fidelity check --
    NOT per-gene means, which are all ~0 because the data is per-cohort z-scored, making that
    correlation structurally meaningless)
  * a random-signature comparison on the same synthetic data (the null still applies)

Result at this sample size: the VAE works (class-difference correlation ~0.95); the diffusion
model does NOT. Its denoiser plateaus at eps-MSE ~0.57 (predicting zero scores 1.0), which is
not accurate enough for a 1000-step reverse process -- the error compounds and the sampler
diverges. That is a sample-size limit, not a bug, and it is reported rather than hidden.

Feature space is reduced to the signature genes plus high-variance context genes; a VAE over
9,900 dimensions from 120 samples would be pure memorisation.
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dkd_data import load_many, GLOM

torch.set_num_threads(4)


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- conditional VAE
class CVAE(nn.Module):
    def __init__(self, p, latent=16, hidden=128):
        super().__init__()
        self.enc = nn.Sequential(nn.Linear(p + 1, hidden), nn.ReLU(), nn.Linear(hidden, hidden), nn.ReLU())
        self.mu = nn.Linear(hidden, latent)
        self.lv = nn.Linear(hidden, latent)
        self.dec = nn.Sequential(nn.Linear(latent + 1, hidden), nn.ReLU(),
                                 nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, p))
        self.latent = latent

    def forward(self, x, y):
        h = self.enc(torch.cat([x, y], 1))
        mu, lv = self.mu(h), self.lv(h).clamp(-6, 2)
        z = mu + torch.randn_like(mu) * (0.5 * lv).exp()
        return self.dec(torch.cat([z, y], 1)), mu, lv


def train_vae(X, y, epochs=6000, latent=16, beta=1.0, seed=0):
    torch.manual_seed(seed)
    p = X.shape[1]
    m = CVAE(p, latent)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3, weight_decay=1e-5)
    xt = torch.tensor(X, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32).view(-1, 1)
    for e in range(epochs):
        opt.zero_grad()
        xr, mu, lv = m(xt, yt)
        rec = ((xr - xt) ** 2).sum(1).mean()
        kl = (-0.5 * (1 + lv - mu ** 2 - lv.exp()).sum(1)).mean()
        (rec + beta * kl).backward()
        opt.step()
    return m


@torch.no_grad()
def sample_vae(m, n, label, seed=0):
    torch.manual_seed(seed + 1)
    z = torch.randn(n, m.latent)
    y = torch.full((n, 1), float(label))
    return m.dec(torch.cat([z, y], 1)).numpy()


# ---------------------------------------------------------------- conditional DDPM
class Denoiser(nn.Module):
    def __init__(self, p, hidden=256):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(p + 2, hidden), nn.SiLU(),
                                 nn.Linear(hidden, hidden), nn.SiLU(), nn.Linear(hidden, p))

    def forward(self, x, t, y):
        return self.net(torch.cat([x, t, y], 1))


def train_ddpm(X, y, steps=1000, epochs=6000, seed=0):
    torch.manual_seed(seed)
    p = X.shape[1]
    betas = torch.linspace(1e-4, 0.02, steps)
    ab = torch.cumprod(1 - betas, 0)
    m = Denoiser(p)
    opt = torch.optim.Adam(m.parameters(), lr=1e-3)
    xt0 = torch.tensor(X, dtype=torch.float32)
    yt = torch.tensor(y, dtype=torch.float32).view(-1, 1)
    n = len(X)
    last = None
    for e in range(epochs):
        t = torch.randint(0, steps, (n,))
        a = ab[t].view(-1, 1)
        eps = torch.randn_like(xt0)
        xt = a.sqrt() * xt0 + (1 - a).sqrt() * eps
        opt.zero_grad()
        pred = m(xt, (t.float() / steps).view(-1, 1), yt)
        loss = ((pred - eps) ** 2).mean()
        loss.backward()
        opt.step()
        last = float(loss)
    # eps-MSE of 1.0 is what predicting zero achieves; anything near 1.0 means the denoiser
    # never learned the score and the reverse process will accumulate error and diverge.
    return m, betas, ab, last


@torch.no_grad()
def sample_ddpm(m, betas, ab, n, p, label, seed=0):
    torch.manual_seed(seed + 2)
    steps = len(betas)
    x = torch.randn(n, p)
    y = torch.full((n, 1), float(label))
    for i in reversed(range(steps)):
        t = torch.full((n, 1), i / steps)
        eps = m(x, t, y)
        a, b = ab[i], betas[i]
        mean = (x - b / (1 - a).sqrt() * eps) / (1 - b).sqrt()
        x = mean + (b.sqrt() * torch.randn_like(x) if i > 0 else 0)
        x = x.clamp(-20, 20)      # guard: without it a poorly-trained denoiser diverges to 1e2+
    return x.numpy()


# ---------------------------------------------------------------- evaluation
def auc(model, X, y, cols):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, model.predict_proba(X[:, cols])[:, 1])


def nn_distance(syn, real):
    d = ((syn[:, None, :] - real[None, :, :]) ** 2).sum(-1)
    return np.sqrt(d.min(1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--candidates', default='results/final_candidates_gated.tsv')
    ap.add_argument('--topn', type=int, default=30)
    ap.add_argument('--context', type=int, default=470, help='high-variance context genes')
    ap.add_argument('--n-syn', type=int, default=300)
    ap.add_argument('--n-random', type=int, default=30)
    ap.add_argument('--cohorts', default=','.join(GLOM) + ',GSE142025')
    ap.add_argument('--seed', type=int, default=20260824)
    ap.add_argument('--out', default='results/generative_robustness.tsv')
    args = ap.parse_args()

    cohorts = [c.strip() for c in args.cohorts.split(',') if c.strip()]
    data, genes = load_many(cohorts)
    genes = np.array([str(g) for g in genes])
    gidx = {g: i for i, g in enumerate(genes)}

    sym2id = pd.read_csv('data/processed/gene_space_all6.tsv', sep='\t', dtype=str)
    sym2id = sym2id.dropna(subset=['symbol']).set_index('symbol')['entrez_id'].to_dict()
    cand = pd.read_csv(args.candidates, sep='\t').head(args.topn)
    sig_full = [gidx[sym2id[s]] for s in cand['symbol'].dropna().astype(str)
                if s in sym2id and sym2id[s] in gidx]

    Xall = np.vstack([data[c][0].values for c in cohorts])
    var = Xall.var(0)
    ctx = [i for i in np.argsort(-var) if i not in set(sig_full)][:args.context]
    sub = np.array(sorted(set(sig_full) | set(ctx)))
    pos = {g: j for j, g in enumerate(sub)}
    sig = [pos[g] for g in sig_full]
    log('feature space for the generators: %d genes (%d signature + %d context)\n'
        % (len(sub), len(sig), len(sub) - len(sig)))

    rng = np.random.default_rng(args.seed)
    rows = []
    for held in cohorts:
        tr = [c for c in cohorts if c != held]
        Xtr = np.vstack([data[c][0].values for c in tr])[:, sub]
        ytr = np.concatenate([data[c][1] for c in tr])
        log('=== held out %s  (train n=%d, %d case) ===' % (held, len(ytr), int(ytr.sum())))

        clf = LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                                 class_weight='balanced').fit(Xtr[:, sig], ytr)
        rand_sets = [rng.choice(len(sub), len(sig), replace=False) for _ in range(args.n_random)]
        rand_clfs = [LogisticRegression(solver='lbfgs', C=1.0, max_iter=5000,
                                        class_weight='balanced').fit(Xtr[:, r], ytr)
                     for r in rand_sets]

        n1 = args.n_syn // 2
        for name in ('vae', 'diffusion'):
            if name == 'vae':
                m = train_vae(Xtr, ytr, seed=args.seed)
                s1 = sample_vae(m, n1, 1, seed=args.seed)
                s0 = sample_vae(m, n1, 0, seed=args.seed + 7)
            else:
                m, betas, ab, eps_mse = train_ddpm(Xtr, ytr, seed=args.seed)
                s1 = sample_ddpm(m, betas, ab, n1, len(sub), 1, seed=args.seed)
                s0 = sample_ddpm(m, betas, ab, n1, len(sub), 0, seed=args.seed + 7)

            Xs = np.vstack([s1, s0])
            ys = np.concatenate([np.ones(n1), np.zeros(n1)]).astype(int)

            a_sig = auc(clf, Xs, ys, sig)
            a_rnd = float(np.nanmean([auc(c_, Xs, ys, r) for c_, r in zip(rand_clfs, rand_sets)]))

            # ---- fidelity. NOT per-gene means: the data is per-cohort z-scored, so every
            # gene's mean is ~0 and SD ~1 and correlating those is structurally meaningless
            # (we measured ~0 for a generator that was in fact working). The quantity that
            # carries the biology is the per-gene CLASS-MEAN DIFFERENCE.
            real_diff = Xtr[ytr == 1].mean(0) - Xtr[ytr == 0].mean(0)
            syn_diff = s1.mean(0) - s0.mean(0)
            r_mu = float(np.corrcoef(real_diff, syn_diff)[0, 1])
            r_sd = float(Xs.std() / Xtr.std())      # overall dispersion ratio, 1.0 = matched
            d_syn = nn_distance(Xs, Xtr)
            d_real = nn_distance(Xtr, Xtr + 1e-9)  # self-distance baseline is ~0, use real-real 2nd NN
            dd = ((Xtr[:, None, :] - Xtr[None, :, :]) ** 2).sum(-1)
            np.fill_diagonal(dd, np.inf)
            d_real = np.sqrt(dd.min(1))
            mem = float(np.median(d_syn) / np.median(d_real))

            usable = (r_mu > 0.5) and (0.4 < r_sd < 2.0)
            log('  %-10s signature AUROC=%.3f random=%.3f margin=%+.3f | class-diff r=%+.2f '
                'dispersion=%.2f NNdist=%.2f  %s'
                % (name, a_sig, a_rnd, a_sig - a_rnd, r_mu, r_sd, mem,
                   'OK' if usable else 'GENERATOR FAILED - result not interpretable'))
            rows.append(dict(held_out=held, generator=name, signature_auroc=a_sig,
                             random_auroc=a_rnd, margin=a_sig - a_rnd,
                             class_diff_corr=r_mu, dispersion_ratio=r_sd,
                             nn_dist_ratio=mem, usable=bool(usable)))
        log('')

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep='\t', index=False)
    log('=== averaged over folds ===')
    for g, d in df.groupby('generator'):
        log('  %-10s signature=%.3f random=%.3f margin=%+.3f | class-diff r=%+.2f '
            'dispersion=%.2f NN-dist=%.2f  usable in %d/%d folds'
            % (g, d['signature_auroc'].mean(), d['random_auroc'].mean(), d['margin'].mean(),
               d['class_diff_corr'].mean(), d['dispersion_ratio'].mean(),
               d['nn_dist_ratio'].mean(), int(d['usable'].sum()), len(d)))
    log('\nNN-dist ratio near 1.0 means synthetic points sit as far from real data as real points')
    log('sit from each other - i.e. no memorisation. Much below 1.0 means the generator is')
    log('copying training samples and the test is uninformative.')
    log('\nwrote %s' % args.out)


if __name__ == '__main__':
    main()
