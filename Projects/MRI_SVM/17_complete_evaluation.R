# Subgroup summaries, sensitivity fits and resumable permutation tests.
main = function() {
  project = Sys.getenv("UCLA_MRI_PROJECT", unset = path.expand("~/ucla-mri-r"))
  run_dir = Sys.getenv("UCLA_MRI_RUN", unset = "")
  if (!nzchar(run_dir)) {
    pointer = file.path(project, "results", "evaluation", "latest_qp_run.txt")
    if (!file.exists(pointer)) stop("Run script 16 first or set UCLA_MRI_RUN to an existing QP result folder.")
    run_dir = readLines(pointer, warn = FALSE)[1]
  }
  if (is.na(run_dir) || !dir.exists(run_dir)) stop("QP result folder not found.")
  B = 1000L
  bootstrap_B = 5000L
  if (!requireNamespace('quadprog',quietly=TRUE)) stop("Install quadprog first.")
  p = read.csv(file.path(run_dir,'out_of_fold_predictions.csv'),stringsAsFactors=FALSE,check.names=FALSE)
  kernel_file = file.path(project,'feature-prep','gm_1p5mm_fwhm8','linear_kernel_raw.rds')
  K = readRDS(kernel_file)
  required = c('feature_row','diagnosis','fold','predicted_diagnosis',
                'ScannerSerialNumber','ghost_NoGhost','age','gender')
  if (!all(required %in% names(p))) stop('Missing columns: ',paste(setdiff(required,names(p)),collapse=', '))
  p = p[order(p$feature_row),]
  id_col = if ('subject_id' %in% names(p)) 'subject_id' else 'participant_id'
  stopifnot(id_col %in% names(p),!anyNA(p[[id_col]]),!anyDuplicated(p[[id_col]]),
    identical(as.integer(p$feature_row),1:175),identical(dim(K),c(175L,175L)),
    all(is.finite(K)),isSymmetric(K),all(p$diagnosis %in% c('CONTROL','SCHZ')),
    !anyNA(p$fold),setequal(p$fold,1:5))
  y = ifelse(p$diagnosis=='SCHZ',1,-1)
  scanner = as.character(p$ScannerSerialNumber)
  ghost = as.character(p$ghost_NoGhost)
  stopifnot(!anyNA(scanner),all(nzchar(scanner)),length(unique(scanner))==2L,
            !anyNA(ghost),all(ghost %in% c('ghost','No_ghost')))
  counts = table(p$fold,p$diagnosis)
  stopifnot(all(counts[,'CONTROL']==25),all(counts[,'SCHZ']==10))
  out = file.path(run_dir,'extended_evaluation')
  dir.create(out,showWarnings=FALSE)
  writeLines(c('Exploratory, post-initial-result evaluation; not a preregistered confirmatory analysis.',
    'Fixed raw-kernel C=1 SVM; no performance-based parameter selection.',
    'Primary permutation null: labels exchangeable within the existing diagnosis-stratified folds.',
    'Secondary null: labels exchangeable within fold x scanner x ghost-status blocks.',
    'Both tests refit all five models for each permutation; fold class counts remain fixed.',
    'Secondary test conditions on observed scanner/ghost association; does not control all confounding.',
    'Exchangeability assumes independent subjects within blocks; related subjects would require different blocks.',
    'P-values use (1 + null scores >= observed)/(B + 1); no failed fits silently omitted.',
    'Both raw and Holm-adjusted p-values reported for the two tests.',
    'Bootstrap intervals resample fixed out-of-fold predictions within diagnosis.',
    'These are descriptive conditional intervals: they omit retraining/split variability and CV dependence.',
    'Scanner-held-out tests assess transfer between these two scanners only.',
    'Sensitivity fits retain the original all-subject mask: no fully independent preprocessing claim.',
    'Ghost exclusion and scanner separation alter sample composition; differences are not causal effects.',
    'Metadata-only logistic baseline uses scanner, ghost flag, age and gender; no feature selection.',
    'No external validation; results are exploratory.'),file.path(out,'analysis_plan_and_limits.txt'))
  metric = function(truth,guess) {
    sen = if(any(truth==1)) mean(guess[truth==1]==1) else NA_real_
    spe = if(any(truth==-1)) mean(guess[truth==-1]==-1) else NA_real_
    c(accuracy=mean(truth==guess),sensitivity=sen,specificity=spe,balanced_accuracy=(sen+spe)/2)
  }
  fit_predict = function(tr,te,labels) {
    z = labels[tr]; n = length(tr)
    if(length(unique(z))!=2) stop('Training split lacks one diagnosis.')
    Q = K[tr,tr,drop=FALSE]*outer(z,z)
    fit = quadprog::solve.QP(Q,rep(1,n),cbind(z,diag(n),-diag(n)),
                             c(0,rep(0,n),rep(-1,n)),meq=1)
    a = fit$solution; beta = a*z
    g = as.vector(K[tr,tr,drop=FALSE]%*%beta)
    low = a<=1e-8; high = a>=1-1e-8; free = !low & !high
    if(any(free)) b = median(z[free]-g[free]) else {
      lo = c((1-g)[z==1 & low],(-1-g)[z==-1 & high])
      hi = c((1-g)[z==1 & high],(-1-g)[z==-1 & low])
      if(!length(lo)||!length(hi)) stop('Intercept interval unavailable.')
      b = (max(lo)+min(hi))/2
    }
    margin = z*(g+b)
    violation = ifelse(low,pmax(0,1-margin),ifelse(high,pmax(0,margin-1),abs(margin-1)))
    norm_sq = as.numeric(crossprod(beta,K[tr,tr,drop=FALSE]%*%beta))
    primal = .5*norm_sq+sum(pmax(0,1-margin)); dual = sum(a)-.5*norm_sq
    if(max(c(0,-a,a-1))>=1e-7 || abs(sum(beta))>=1e-7 || max(violation)>=.01 ||
       abs(primal-dual)>=1e-4*max(1,abs(primal),abs(dual))) stop('Optimisation failed; no result accepted.')
    ifelse(as.vector(K[te,tr,drop=FALSE]%*%beta+b)>0,1,-1)
  }
  cv = function(labels,indices=seq_len(nrow(p))) {
    guess = rep(NA_real_,nrow(p))
    for(f in sort(unique(p$fold[indices]))) {
      te = indices[p$fold[indices]==f]; tr = indices[p$fold[indices]!=f]
      if(length(te)) guess[te] = fit_predict(tr,te,labels)
    }
    stopifnot(!anyNA(guess[indices]))
    guess
  }
  cat('Reproducing the corrected fit with optimisation checks...\n')
  original = cv(y)
  saved = ifelse(p$predicted_diagnosis=='SCHZ',1,-1)
  stopifnot(all(original==saved))
  observed = metric(y,original)
  cat('Observed balanced accuracy:',observed['balanced_accuracy'],'\n')
  groups = list(fold=as.character(p$fold),scanner=scanner,ghost=ghost)
  subgroup = do.call(rbind,lapply(names(groups),function(v) {
    do.call(rbind,lapply(sort(unique(groups[[v]])),function(g) {
      i = which(groups[[v]]==g)
      cbind(data.frame(variable=v,group=g,n=length(i),n_CONTROL=sum(y[i]==-1),n_SCHZ=sum(y[i]==1)),
            as.data.frame(as.list(metric(y[i],original[i]))))
    }))
  }))
  write.csv(subgroup,file.path(out,'corrected_subgroup_metrics.csv'),row.names=FALSE)
  set.seed(202610171)
  strata = split(seq_along(y),y)
  boot = replicate(bootstrap_B,{
    i = unlist(lapply(strata,function(ix) ix[sample.int(length(ix),length(ix),replace=TRUE)]),use.names=FALSE)
    metric(y[i],original[i])
  })
  intervals = data.frame(metric=names(observed),estimate=as.numeric(observed),
    lower=apply(boot,1,quantile,.025),upper=apply(boot,1,quantile,.975),
    interpretation='Descriptive fixed-prediction bootstrap; not full CV uncertainty')
  write.csv(intervals,file.path(out,'descriptive_bootstrap_intervals.csv'),row.names=FALSE)
  collected = new.env(parent = emptyenv())
  collected$sensitivity = list(); collected$sensitivity_predictions = list()
  add_result = function(name,i,guess) {
    collected$sensitivity[[name]] = cbind(data.frame(analysis=name,n=length(i),
      n_CONTROL=sum(y[i]==-1),n_SCHZ=sum(y[i]==1)),as.data.frame(as.list(metric(y[i],guess))))
    collected$sensitivity_predictions[[name]] = data.frame(analysis=name,feature_row=i,
      actual=y[i],predicted=guess)
  }
  i = which(ghost=='No_ghost')
  ng = cv(y,i); add_result('No_ghost only, retrained five-fold',i,ng[i])
  for(s in sort(unique(scanner))) {
    i = which(scanner==s)
    within = cv(y,i); add_result(paste('Within scanner',s),i,within[i])
    tr = which(scanner!=s)
    cross = fit_predict(tr,i,y)
    add_result(paste('Held-out scanner',s),i,cross)
  }
  age = suppressWarnings(as.numeric(p$age))
  if(any(!is.finite(age)) || anyNA(p$gender)) stop('Missing age/gender requires an explicit baseline handling decision.')
  d = data.frame(response=as.integer(y==1),age=age,gender=factor(p$gender),
                  scanner=factor(scanner),ghost=factor(ghost))
  baseline = rep(NA_real_,length(y)); collected$baseline_warnings = character()
  for(f in 1:5) {
    tr = which(p$fold!=f); te = which(p$fold==f)
    for(v in c('gender','scanner','ghost')) {
      if(!all(as.character(d[[v]][te]) %in% as.character(d[[v]][tr])))
        stop('Unseen metadata category in a baseline test fold.')
    }
    m = withCallingHandlers(glm(response~age+gender+scanner+ghost,data=d[tr,],family=binomial()),
      warning=function(w) {collected$baseline_warnings = c(collected$baseline_warnings,paste('fold',f,conditionMessage(w)));invokeRestart('muffleWarning')})
    if(!isTRUE(m$converged)) collected$baseline_warnings = c(collected$baseline_warnings,paste('fold',f,'did not converge'))
    probability = predict(m,newdata=d[te,],type='response')
    if(anyNA(probability)) stop('Metadata baseline predictions missing.')
    baseline[te] = ifelse(probability>.5,1,-1)
  }
  add_result('Metadata-only logistic baseline',seq_along(y),baseline)
  writeLines(if(length(collected$baseline_warnings)) collected$baseline_warnings else 'No glm fitting warnings.',
             file.path(out,'metadata_baseline_warnings.txt'))
  sensitivity_table = do.call(rbind,collected$sensitivity)
  write.csv(sensitivity_table,file.path(out,'sensitivity_metrics.csv'),row.names=FALSE)
  write.csv(do.call(rbind,collected$sensitivity_predictions),file.path(out,'sensitivity_predictions_PRIVATE.csv'),row.names=FALSE)
  cat('Corrected subgroups, descriptive intervals and sensitivity fits saved.\n')
  schemes = list(within_fold=interaction(p$fold,drop=TRUE),
                  within_fold_scanner_ghost=interaction(p$fold,scanner,ghost,drop=TRUE))
  fingerprint = unname(tools::md5sum(c(kernel_file,file.path(run_dir,'out_of_fold_predictions.csv'))))
  permutation_results = list()
  for(scheme in names(schemes)) {
    blocks = split(seq_along(y),schemes[[scheme]])
    write.csv(do.call(rbind,lapply(names(blocks),function(g) {
      i = blocks[[g]];data.frame(block=g,n=length(i),CONTROL=sum(y[i]==-1),SCHZ=sum(y[i]==1))
    })),file.path(out,paste0('permutation_blocks_',scheme,'.csv')),row.names=FALSE)
    checkpoint = file.path(out,paste0('permutation_',scheme,'.rds'))
    state = list(version=1L,fingerprint=fingerprint,B=B,scheme=scheme,values=numeric())
    if(file.exists(checkpoint)) {
      state = readRDS(checkpoint)
      stopifnot(identical(state$version,1L),identical(state$fingerprint,fingerprint),
                identical(state$B,B),identical(state$scheme,scheme))
    }
    first = length(state$values)+1L
    if(first<=B) for(j in seq.int(first,B)) {
      set.seed(20261017L+j+if(scheme=='within_fold') 0L else 100000L)
      yp = y
      for(i in blocks) yp[i] = y[i][sample.int(length(i))]
      guess = cv(yp)
      state$values[j] = unname(metric(yp,guess)['balanced_accuracy'])
      tmp = paste0(checkpoint,'.tmp');saveRDS(state,tmp)
      if(!file.rename(tmp,checkpoint)) stop('Could not save permutation checkpoint.')
      if(j%%25L==0L || j==1L) cat(scheme,':',j,'/',B,'permutations complete\n')
    }
    exceed = sum(state$values>=observed['balanced_accuracy']-1e-12)
    mc = binom.test(exceed,B)$conf.int
    permutation_results[[scheme]] = data.frame(test=scheme,B=B,observed_BA=observed['balanced_accuracy'],
      null_mean=mean(state$values),exceedances=exceed,p_value=(exceed+1)/(B+1),
      monte_carlo_tail_lower=mc[1],monte_carlo_tail_upper=mc[2])
    pdf(file.path(out,paste0('null_distribution_',scheme,'.pdf')))
    hist(state$values,breaks=25,main=scheme,xlab='Null balanced accuracy',col='grey80',
      xlim=range(c(state$values,observed['balanced_accuracy'])))
    abline(v=observed['balanced_accuracy'],col='red',lwd=2)
    dev.off()
  }
  perm = do.call(rbind,permutation_results)
  perm$p_holm_two_tests = p.adjust(perm$p_value,method='holm')
  write.csv(perm,file.path(out,'permutation_results.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
  cat('\nDESCRIPTIVE INTERVALS (conditional on fixed predictions):\n');print(intervals,row.names=FALSE)
  cat('\nRETRAINED SENSITIVITIES AND METADATA BASELINE:\n');print(sensitivity_table,row.names=FALSE)
  cat('\nPERMUTATION RESULTS:\n');print(perm,row.names=FALSE)
  if(length(collected$baseline_warnings)) cat('\nMetadata baseline has warnings: inspect metadata_baseline_warnings.txt.\n')
  cat('\nEvaluation completed. Results:',out,'\n')
  cat('Aggregate results are ready for review.\n')
}
main()
