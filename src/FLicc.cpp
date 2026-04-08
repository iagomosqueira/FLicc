#include <TMB.hpp>

// ---------------------------------------------------------
// Selectivity
// sel_type:
//   1 = logistic
//   2 = dsnormal
//
// packed pars:
//   logistic : c(L50, log_steep)
//   dsnormal : c(mode, log_lsd, log_rsd)
// ---------------------------------------------------------
template<class Type>
vector<Type> calc_sel(const vector<Type>& L,
                      int sel_type,
                      const vector<Type>& pars,
                      Type Linf)
{
  int n = L.size();
  vector<Type> sel(n);

  if(sel_type == 1) {
    // logistic, LBSPR-style
    // pars = c(log_SL50_rel, log_dSL_rel)
    // where:
    //   SL50 = exp(pars[0]) * Linf
    //   dSL  = exp(pars[1]) * Linf
    //   SL95 = SL50 + dSL

    Type SL50 = exp(pars(0)) * Linf;
    Type dSL  = exp(pars(1)) * Linf;
    Type Delta = dSL;

    for(int i = 0; i < n; i++) {
      sel(i) = Type(1) / (Type(1) + exp(-log(Type(19.0)) * (L(i) - SL50) / Delta));
    }

  } else if(sel_type == 2) {
    // dsnormal, relative to Linf
    // pars = c(log_mode_rel, log_lsd_rel, log_rsd_rel)

    Type mode = exp(pars(0)) * Linf;
    Type lsd  = exp(pars(1)) * Linf;
    Type rsd  = exp(pars(2)) * Linf;
    Type two  = Type(2);

    for(int i = 0; i < n; i++) {
      Type d = L(i) - mode;
      if(L(i) < mode) {
        sel(i) = exp(-(d * d) / (two * lsd * lsd));
      } else {
        sel(i) = exp(-(d * d) / (two * rsd * rsd));
      }
    }

  } else if(sel_type == 3) {
    // symmetric normal, relative to Linf
    // pars = c(log_mode_rel, log_sd_rel)

    Type mode = exp(pars(0)) * Linf;
    Type sd   = exp(pars(1)) * Linf;
    Type two  = Type(2);

    for(int i = 0; i < n; i++) {
      Type d = L(i) - mode;
      sel(i) = exp(-(d * d) / (two * sd * sd));
    }

  } else {
    error("Unknown sel_type");
  }

  return sel;
}

// ---------------------------------------------------------
// Negative binomial parameterization:
// mean = mu
// variance = mu + mu^2 / phi
// implemented via size=phi, prob=phi/(phi+mu)
// ---------------------------------------------------------
template<class Type>
Type dnbinom_phi(Type x, Type mu, Type phi, int give_log = 0)
{
  Type logres =
    lgamma(x + phi) - lgamma(phi) - lgamma(x + Type(1)) +
    phi * log(phi / (phi + mu)) +
    x   * log(mu  / (phi + mu));

  if(give_log) return logres;
  return exp(logres);
}

// ---------------------------------------------------------
// fishblicc-style population at length using
// Gauss-Laguerre quadrature over gamma-distributed Linf
//
// Len  : lower bounds of length bins
// Zki  : total mortality by bin
// Galpha, Gbeta : gamma parameters for Linf variability
//
// Returns interval numbers-at-length NI, where
//   NI(i) = (S(i) - S(i+1)) / Z(i)
// and last bin is treated as a plus group:
//   NI(last) = S(last) / Z(last)
// ---------------------------------------------------------
template<class Type>
vector<Type> pop_len_glq(const vector<Type>& node,
                         const vector<Type>& quad_wt,
                         const vector<Type>& Len,
                         const vector<Type>& Zki,
                         Type Galpha,
                         Type Gbeta)
{
  int nv  = node.size();
  int LN  = Len.size();
  int LN1 = LN - 1;

  if(quad_wt.size() != nv) error("node and quad_wt must have same length");
  if(Zki.size() != LN) error("Zki and Len must have same length");

  Type lgamma_Galpha = lgamma(Galpha);
  Type Galpha_1      = Galpha - Type(1);

  vector<Type> surv(LN);
  surv.setZero();

  vector<Type> x_beta(nv);
  vector<Type> log_x_beta(nv);

  vector<Type> Zii(LN);
  Zii.setZero();

  vector<Type> NI(LN);
  NI.setZero();

  // differences in Z between adjacent bins
  Zii(0) = -Zki(0);
  for(int i = 1; i < (LN - 1); ++i) {
    Zii(i) = Zki(i - 1) - Zki(i);
  }

  // first two bins
  for(int i = 0; i < nv; ++i) {
    x_beta(i)     = node(i) / Gbeta;
    log_x_beta(i) = log(x_beta(i));

    surv(0) += exp(
      log(node(i) + Gbeta * Len(0)) * Galpha_1
    - Gbeta * Len(0)
      - lgamma_Galpha
    ) * quad_wt(i);

    if(LN > 1) {
      surv(1) += exp(
        log(x_beta(i) + Len(1) - Len(0)) * Zii(0)
      + log_x_beta(i) * Zki(0)
      + log(node(i) + Gbeta * Len(1)) * Galpha_1
      - Gbeta * Len(1)
      - lgamma_Galpha
      ) * quad_wt(i);
    }
  }

  // remaining bins
  for(int Li = 2; Li < LN; ++Li) {
    Type Ln = Len(Li);
    Type v1 = Gbeta * Ln + lgamma_Galpha;

    for(int i = 0; i < nv; ++i) {
      Type lim = Type(0);

      for(int Lii = 0; Lii < Li; ++Lii) {
        Type Lrange = Ln - Len(Lii);
        lim += log(x_beta(i) + Lrange) * Zii(Lii);
      }

      surv(Li) += exp(
        lim
        + log_x_beta(i) * Zki(Li - 1)
        + log(node(i) + Gbeta * Ln) * Galpha_1
        - v1
      ) * quad_wt(i);
    }
  }

  // interval abundance
  for(int i = 0; i < LN1; ++i) {
    NI(i) = (surv(i) - surv(i + 1)) / Zki(i);
  }

  // plus group
  NI(LN1) = surv(LN1) / Zki(LN1);

  return NI;
}

// ---------------------------------------------------------
// Objective
// ---------------------------------------------------------
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER(nlen);
  DATA_INTEGER(ngear);
  DATA_VECTOR(Lmid);           // length midpoints
  DATA_VECTOR(LLB);            // lower bin bounds
  DATA_MATRIX(obs);            // nlen x ngear counts
  DATA_VECTOR(catch_wt);       // retained for structure; not used in mu scaling here
  DATA_IVECTOR(sel_type);      // 1=logistic, 2=dsnormal
  DATA_IVECTOR(sm_start);      // 0-based start index in Sm
  DATA_IVECTOR(sm_n);          // number of pars for each gear
  DATA_VECTOR(Mscaler)       // scales MK(L) according to M model
  DATA_VECTOR(wt);            // vector weight-at-length
  DATA_VECTOR(mat);           // vector length-at-maturity

  DATA_SCALAR(catch_sd);

  // quadrature inputs for fishblicc recursion
  DATA_VECTOR(node);
  DATA_VECTOR(quad_wt);

  // Prior controls (simple normal penalties)
  DATA_VECTOR(prior_mu);
  DATA_VECTOR(prior_sd);
  DATA_IVECTOR(prior_code);

  // prior_code:
  // 1 log_Linf
  // 2 log_Galpha
  // 3 log_Mk
  // 4 log_phi
  // 100 + g => log_Fk[g]
  // 200 + j => Sm[j]

  // Parameters
  PARAMETER(log_Linf);
  PARAMETER(log_Galpha);
  PARAMETER(log_Mk);
  PARAMETER_VECTOR(log_Fk);    // length ngear
  PARAMETER_VECTOR(Sm);        // packed selectivity pars
  PARAMETER(log_phi);

  Type nll = 0.0;

  // Natural-scale parameters
  Type Linf   = exp(log_Linf);
  Type Galpha = exp(log_Galpha);
  Type Mk     = exp(log_Mk);
  Type phi    = exp(log_phi);
  Type Gbeta  = Galpha / Linf;

  vector<Type> Fk(ngear);
  for(int g = 0; g < ngear; g++) Fk(g) = exp(log_Fk(g));

  // Selectivity matrix
  matrix<Type> Sel(nlen, ngear);
  for(int g = 0; g < ngear; g++) {
    if(sm_start(g) < 0) error("sm_start < 0");
    if(sm_start(g) + sm_n(g) > Sm.size()) error("sm_start + sm_n exceeds Sm size");

    vector<Type> pars(sm_n(g));
    for(int j = 0; j < sm_n(g); j++) {
      pars(j) = Sm(sm_start(g) + j);
    }
    Sel.col(g) = calc_sel(Lmid, sel_type(g), pars, Linf);
  }

  // Natural mortality at length
  // kept on Lmid for now
  vector<Type> M(nlen);
  for(int l = 0; l < nlen; l++) {

      M(l) = Mk * Mscaler(l);

  }


  // Fishing mortality at length and total Z
  matrix<Type> F_len(nlen, ngear);
  vector<Type> Z(nlen);
  vector<Type> Z0(nlen);

  for(int l = 0; l < nlen; l++) {
    Z0(l) = M(l);
    Z(l)  = M(l);
    for(int g = 0; g < ngear; g++) {
      F_len(l, g) = Fk(g) * Sel(l, g);
      Z(l) += F_len(l, g);
    }
  }

  // ---------------------------------------------------------
  // fishblicc-style population recursion
  // uses lower bin bounds (LLB), not Lmid
  // ---------------------------------------------------------
  vector<Type> N0 = pop_len_glq(node, quad_wt, LLB, Z0, Galpha, Gbeta);
  vector<Type> N  = pop_len_glq(node, quad_wt, LLB, Z,  Galpha, Gbeta);

  // SPR = SBPR_F / SBPR_0
  Type spr0 = Type(0);
  Type sprf = Type(0);

  for(int l = 0; l < nlen; l++) {
    spr0 += N0(l) * mat(l) * wt(l);
    sprf += N(l)  * mat(l) * wt(l);
  }

  Type spr = Type(0);
  if(spr0 > Type(0)) spr = sprf / spr0;

  // Expected catch-at-length by gear:
  // Cpred(l,g) = N(l) * F_len(l,g)
  matrix<Type> Cpred(nlen, ngear);
  for(int l = 0; l < nlen; l++) {
    for(int g = 0; g < ngear; g++) {
      Cpred(l, g) = N(l) * F_len(l, g);
    }
  }

  // ---------------------------------------------------------
  // Global expected-count construction across all gears
  // This lets Fk influence allocation among gears as well
  // ---------------------------------------------------------
  vector<Type> pred_gear(ngear);
  vector<Type> obs_gear(ngear);
  vector<Type> pred_pgear(ngear);

  pred_gear.setZero();
  obs_gear.setZero();

  Type pred_total = Type(0);
  Type obs_total  = Type(0);

  for(int g = 0; g < ngear; g++) {
    for(int l = 0; l < nlen; l++) {
      pred_gear(g) += Cpred(l, g);
      obs_gear(g)  += obs(l, g);
      pred_total   += Cpred(l, g);
      obs_total    += obs(l, g);
    }
  }

  for(int g = 0; g < ngear; g++) {
    if(pred_total > Type(0)) {
      pred_pgear(g) = pred_gear(g) / pred_total;
    } else {
      pred_pgear(g) = Type(0);
    }
  }

  // ---------------------------------------------------------
  // Observed and predicted proportions at length
  // ---------------------------------------------------------
  matrix<Type> obs_p_lg(nlen, ngear);
  matrix<Type> pred_p_lg(nlen, ngear);
  vector<Type> obs_p_l(nlen);
  vector<Type> pred_p_l(nlen);

  obs_p_l.setZero();
  pred_p_l.setZero();

  for(int g = 0; g < ngear; g++) {
    for(int l = 0; l < nlen; l++) {
      // observed within-gear proportions at length
      if(obs_gear(g) > Type(0)) {
        obs_p_lg(l, g) = obs(l, g) / obs_gear(g);
      } else {
        obs_p_lg(l, g) = Type(0);
      }

      // predicted within-gear proportions at length
      if(pred_gear(g) > Type(0)) {
        pred_p_lg(l, g) = Cpred(l, g) / pred_gear(g);
      } else {
        pred_p_lg(l, g) = Type(0);
      }
    }
  }

  // catch-weighted overall proportions across gears
  for(int l = 0; l < nlen; l++) {
    for(int g = 0; g < ngear; g++) {
      obs_p_l(l)  += catch_wt(g) * obs_p_lg(l, g);
      pred_p_l(l) += catch_wt(g) * pred_p_lg(l, g);
    }
  }

  // build mu for the likelihood


  matrix<Type> mu(nlen, ngear);

  for(int g = 0; g < ngear; g++) {
    for(int l = 0; l < nlen; l++) {
      if(pred_gear(g) > Type(0)) {
        mu(l, g) = obs_gear(g) * Cpred(l, g) / pred_gear(g) + Type(1e-12);
      } else {
        mu(l, g) = Type(1e-12);
      }
    }
  }


  // Likelihood
  for(int g = 0; g < ngear; g++) {
    for(int l = 0; l < nlen; l++) {
      nll -= dnbinom_phi(obs(l, g), mu(l, g), phi, 1);
    }
  }

  // Normal prior penalties
  for(int i = 0; i < prior_mu.size(); i++) {
    Type x;
    int code = prior_code(i);

    if(code == 1) {
      x = log_Linf;
    } else if(code == 2) {
      x = log_Galpha;
    } else if(code == 3) {
      x = log_Mk;
    } else if(code == 4) {
      x = log_phi;
    } else if(code >= 100 && code < 200) {
      int g = code - 100;
      x = log_Fk(g);
    } else if(code >= 200) {
      int j = code - 200;
      x = Sm(j);
    } else {
      error("Unknown prior_code");
    }

    nll -= dnorm(x, prior_mu(i), prior_sd(i), true);
  }

  //Penalty
  // soft constraint on predicted gear shares using relative catch
  // smaller = stronger anchoring

  for(int g = 0; g < (ngear - 1); g++) {
    if(pred_pgear(g) > Type(0) &&
       pred_pgear(ngear - 1) > Type(0) &&
       catch_wt(g) > Type(0) &&
       catch_wt(ngear - 1) > Type(0)) {

      Type logratio_pred  = log(pred_pgear(g) / pred_pgear(ngear - 1));
      Type logratio_catch = log(catch_wt(g)   / catch_wt(ngear - 1));

      nll -= dnorm(logratio_pred, logratio_catch, catch_sd, true);
    }
  }


  REPORT(Linf);
  REPORT(Galpha);
  REPORT(Gbeta);
  REPORT(Mk);
  REPORT(phi);
  REPORT(Fk);
  REPORT(M);
  REPORT(Sel);
  REPORT(F_len);
  REPORT(Z0);
  REPORT(Z);
  REPORT(N0);
  REPORT(N);
  REPORT(Cpred);
  REPORT(mu);
  REPORT(mat);
  REPORT(wt);
  REPORT(spr0);
  REPORT(sprf);
  REPORT(spr);
  REPORT(pred_gear);
  REPORT(obs_gear);
  REPORT(pred_pgear);
  REPORT(obs_gear);
  REPORT(pred_gear);
  REPORT(pred_pgear);
  REPORT(obs_p_lg);
  REPORT(pred_p_lg);
  REPORT(obs_p_l);
  REPORT(pred_p_l);

  ADREPORT(Linf);
  ADREPORT(Galpha);
  ADREPORT(Mk);
  ADREPORT(phi);
  ADREPORT(Fk);
  ADREPORT(spr);

  return nll;
}
