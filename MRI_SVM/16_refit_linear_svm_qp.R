# Linear C=1 SVM with numerical optimality checks.
project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
if (!requireNamespace('quadprog', quietly=TRUE)) stop("Run install.packages('quadprog') first.")
base = file.path(project,'results','evaluation')
p = read.csv(file.path(base,'cv_folds.csv'), stringsAsFactors=FALSE)
p = p[order(p$feature_row),]
K = readRDS(file.path(project,'feature-prep','gm_1p5mm_fwhm8','linear_kernel_raw.rds'))
stopifnot(identical(as.integer(p$feature_row),1:175),
  identical(dim(K),c(175L,175L)),all(is.finite(K)),isSymmetric(K),
  all(p$diagnosis %in% c('CONTROL','SCHZ')),setequal(p$fold,1:5))
counts = table(p$fold,p$diagnosis)
stopifnot(all(counts[,'CONTROL']==25),all(counts[,'SCHZ']==10))
C = 1
pred = rep(NA_character_,175)
scores = rep(NA_real_,175)
models = diagnostics = vector('list',5)
for (f in 1:5) {
  tr = which(p$fold!=f); te = which(p$fold==f)
  y = ifelse(p$diagnosis[tr]=='SCHZ',1,-1)
  n = length(tr)
  Q = K[tr,tr]*outer(y,y)
  fit = quadprog::solve.QP(Q,rep(1,n),
    cbind(y,diag(n),-diag(n)),c(0,rep(0,n),rep(-C,n)),meq=1)
  alpha = fit$solution
  beta = alpha*y
  g = as.vector(K[tr,tr]%*%beta)
  bound_tol = 1e-8
  low = alpha <= bound_tol
  high = alpha >= C-bound_tol
  free = !low & !high
  if (any(free)) {
    intercept = median(y[free]-g[free])
  } else {
    lower = c((1-g)[y==1 & low],(-1-g)[y==-1 & high])
    upper = c((1-g)[y==1 & high],(-1-g)[y==-1 & low])
    if (!length(lower)||!length(upper)) stop('Cannot determine finite intercept.')
    intercept = (max(lower)+min(upper))/2
  }
  margin = y*(g+intercept)
  violation = ifelse(low,pmax(0,1-margin),
    ifelse(high,pmax(0,margin-1),abs(margin-1)))
  norm_sq = as.numeric(crossprod(beta,K[tr,tr]%*%beta))
  primal = .5*norm_sq+C*sum(pmax(0,1-margin))
  dual = sum(alpha)-.5*norm_sq
  gap = primal-dual
  box = max(c(0,-alpha,alpha-C))
  equality = abs(sum(beta))
  passed = box<1e-7 && equality<1e-7 && max(violation)<.01 &&
    abs(gap)<1e-4*max(1,abs(primal),abs(dual))
  diagnostics[[f]] = data.frame(fold=f,box_violation=box,
    equality_residual=equality,max_margin_violation=max(violation),
    primal=primal,dual=dual,objective_gap=gap,passed=passed)
  print(diagnostics[[f]],row.names=FALSE)
  if (!passed) stop('Optimisation check failed. No performance reported; send this output.')
  scores[te] = as.vector(K[te,tr,drop=FALSE]%*%beta+intercept)
  pred[te] = ifelse(scores[te]>0,'SCHZ','CONTROL')
  models[[f]] = list(training_rows=tr,test_rows=te,alpha=alpha,beta=beta,
    intercept=intercept,C=C,solver_iterations=fit$iterations)
}
stopifnot(!anyNA(pred))
metric = function(truth,guess) {
  sensitivity = mean(guess[truth=='SCHZ']=='SCHZ')
  specificity = mean(guess[truth=='CONTROL']=='CONTROL')
  data.frame(accuracy=mean(truth==guess),sensitivity=sensitivity,
    specificity=specificity,balanced_accuracy=(sensitivity+specificity)/2)
}
summary = rbind(cbind(model='Linear SVM C=1 (QP)',metric(p$diagnosis,pred)),
  cbind(model='Always CONTROL',metric(p$diagnosis,rep('CONTROL',175))))
fold_metrics = do.call(rbind,lapply(1:5,function(f) {
  i = p$fold==f
  cbind(fold=f,metric(p$diagnosis[i],pred[i]))
}))
out = file.path(base,paste0('linear_svm_C1_QP_',format(Sys.time(),'%Y%m%d_%H%M%S')))
if (dir.exists(out)) stop('Output already exists; retry in a second.')
dir.create(out,recursive=TRUE)
p$predicted_diagnosis = pred
p$correct = pred==p$diagnosis
p$decision_score = scores
write.csv(p,file.path(out,'out_of_fold_predictions.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'summary_metrics.csv'),row.names=FALSE)
write.csv(fold_metrics,file.path(out,'fold_metrics.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'optimisation_diagnostics.csv'),row.names=FALSE)
saveRDS(models,file.path(out,'fold_models_qp.rds'))
writeLines(c('Corrective alternative solver for the same unscaled raw-kernel C=1 objective.',
  'Same subjects, labels and folds as the first run; no class weights or tuning.',
  'Uses saved folds directly; earlier failed fits are not inputs.',
  'All fits passed declared feasibility, margin and objective-gap checks.',
  'Existing whole-cohort mask is transductive; not wholly independent preprocessing.',
  'No external validation; permutation, uncertainty and confound analyses pending.',
  'QP model objects have a different format from kernlab objects.',
  paste('quadprog version:',packageVersion('quadprog'))),file.path(out,'methods.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
cat('\nAll five optimisation checks passed.\n')
print(table(Actual=p$diagnosis,Predicted=factor(pred,levels=c('CONTROL','SCHZ'))))
print(summary,row.names=FALSE)
writeLines(normalizePath(out), file.path(base, 'latest_qp_run.txt'))
cat('\nSaved results to:',out,'\n')
