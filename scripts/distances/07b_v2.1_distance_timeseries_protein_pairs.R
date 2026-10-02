###############################################################################
# 07b_v2.1_multiPairs_fast.R
# Author: DenyCB
###############################################################################

suppressPackageStartupMessages({
  library(bio3d)
  library(parallel)
  library(ggplot2)
  library(dplyr)
  library(mclust)
})

###############################################################################
# USER SETTINGS
###############################################################################

INPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/input/md"
OUT_BASE  <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/MD_distances/07b_v2.0"

PAIRS <- list(
#  c("510","555","+0"),
#  c("501","510","+0"),
#  c("553","512","+0"),
#  c("501","505","+0"),
#  c("496","501","+0"),
#  c("450","470","+0"),
#  c("473","450","+0"),
  c("470","454","+0")
)

N_CORES <- 6

FRAME_PS <- 100
FRAME_RANGE <- NULL

SHOW_SD <- TRUE

PNG_W <- 3000
PNG_H <- 2200
PNG_DPI <- 600

###############################################################################
# COLORS
###############################################################################

COL_MUT <- c(
  WT    = "black",
  W426A = "orange",
  W697A = "forestgreen",
  Y441A = "purple"
)

COL_COND <- c(
  Apo  = "black",
  Cap  = "red",
  Heat = "blue"
)

###############################################################################
# FOLDERS
###############################################################################

dir.create(OUT_BASE, recursive=TRUE, showWarnings=FALSE)

# NUEVO: carpeta cache
CACHE_DIR <- file.path(OUT_BASE, "cache_xyz")
dir.create(CACHE_DIR, recursive=TRUE, showWarnings=FALSE)

LOG_FILE <- file.path(
  OUT_BASE,
  paste0("07b_v2.1_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

###############################################################################
# LOG
###############################################################################

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse=""))
  cat(txt, "\n")
  cat(txt, "\n", file=LOG_FILE, append=TRUE)
}

log_msg("===========================================================")
log_msg("START 07b_v2.1_multiPairs_fast")
log_msg("Pairs requested: ", length(PAIRS))
log_msg("Cores: ", N_CORES)

###############################################################################
# HELPERS
###############################################################################

shift_chain <- function(ch, sh) {
  ord <- c("A","B","C","D")
  i <- match(ch, ord)
  ord[((i - 1 + sh) %% 4) + 1]
}

dist3 <- function(a,b) sqrt(sum((a-b)^2))

get_atom_index <- function(atom, chain_id, resno) {
  
  x <- atom[
    atom$chain == chain_id &
      atom$resno == as.integer(resno) &
      atom$type == "ATOM",
  ]
  
  if (nrow(x) == 0) return(NA_integer_)
  
  cb <- which(x$elety == "CB")
  if (length(cb) > 0) return(as.integer(rownames(x)[cb[1]]))
  
  ca <- which(x$elety == "CA")
  if (length(ca) > 0) return(as.integer(rownames(x)[ca[1]]))
  
  return(NA_integer_)
}

classify_system <- function(sysname) {
  
  pdb <- substr(sysname,1,4)
  
  cond <- ifelse(pdb=="7LP9","Apo",
                 ifelse(pdb=="7LPB","Cap","Heat"))
  
  mut <- strsplit(sysname,"_")[[1]][2]
  
  list(cond=cond, mut=mut)
}

# NUEVO ------------------------------------------------------------
cache_file <- function(sysname, resid) {
  file.path(CACHE_DIR, paste0(sysname, "__", resid, ".rds"))
}

###############################################################################
# DISCOVER SYSTEMS
###############################################################################

systems <- list.files(INPUT_DIR, pattern="\\.dcd$")
systems <- sub("\\.dcd$", "", systems)

log_msg("Systems found: ", length(systems))

###############################################################################
# CORE WORKER
###############################################################################

process_system <- function(sysname, pair_list, input_dir, frame_ps, frame_range) {
  
  # residuos únicos solicitados
  needed_res <- unique(unlist(lapply(pair_list, function(x) x[1:2])))
  
  # revisar cache
  cache_exists <- sapply(needed_res, function(r)
    file.exists(cache_file(sysname, r)))
  
  need_read <- !all(cache_exists)
  
  # ------------------------------------------------------------
  # si falta algo -> leer PDB/DCD y construir faltantes
  # ------------------------------------------------------------
  if (need_read) {
    
    pdbfile <- file.path(input_dir, paste0(sysname, ".pdb"))
    dcdfile <- file.path(input_dir, paste0(sysname, ".dcd"))
    
    if (!file.exists(pdbfile) || !file.exists(dcdfile))
      return(NULL)
    
    pdb <- read.pdb(pdbfile)
    trj <- read.dcd(dcdfile)
    
    if (!is.null(frame_range))
      trj <- trj[frame_range,,drop=FALSE]
    
    atom <- pdb$atom
    rownames(atom) <- seq_len(nrow(atom))
    
    nf <- nrow(trj)
    
    # guardar caches faltantes
    for (rr in needed_res[!cache_exists]) {
      
      mat <- matrix(NA_real_, nrow=nf, ncol=12)
      
      colnames(mat) <- c(
        "A_x","A_y","A_z",
        "B_x","B_y","B_z",
        "C_x","C_y","C_z",
        "D_x","D_y","D_z"
      )
      
      for (ic in seq_along(c("A","B","C","D"))) {
        
        ch <- c("A","B","C","D")[ic]
        
        idx <- get_atom_index(atom, ch, rr)
        
        if (is.na(idx)) next
        
        xyz <- atom2xyz(idx)
        
        cols <- ((ic-1)*3+1):((ic-1)*3+3)
        
        mat[, cols] <- trj[, xyz, drop=FALSE]
      }
      
      saveRDS(mat, cache_file(sysname, rr))
    }
    
  } else {
    
    # para conocer nf desde cache
    tmp <- readRDS(cache_file(sysname, needed_res[1]))
    nf <- nrow(tmp)
  }
  
  # ------------------------------------------------------------
  # calcular pares usando SOLO cache
  # ------------------------------------------------------------
  pair_results <- list()
  raw_results <- list()###
  
  for (pp in seq_along(pair_list)) {
    
    pair <- pair_list[[pp]]
    
    r1 <- pair[1]
    r2 <- pair[2]
    sh <- as.integer(pair[3])
    
    xyz_res1 <- readRDS(cache_file(sysname, r1))
    xyz_res2 <- readRDS(cache_file(sysname, r2))
    
    meanv <- numeric(nf)
    sdv   <- numeric(nf)
    
    raw_mat <- matrix(NA_real_, nrow=nf, ncol=4)####
    colnames(raw_mat) <- c("chain1","chain2","chain3","chain4")####
    
    for (f in seq_len(nf)) {
      
      dvec <- c()
      
      for (ch in 1:4) {
        
        ch2 <- ((ch - 1 + sh) %% 4) + 1
        
        c1 <- ((ch-1)*3+1):((ch-1)*3+3)
        c2 <- ((ch2-1)*3+1):((ch2-1)*3+3)
        
        p1 <- as.numeric(xyz_res1[f, c1])
        p2 <- as.numeric(xyz_res2[f, c2])
        
        if (any(is.na(p1)) || any(is.na(p2))) next
        
        dvec <- c(dvec, dist3(p1,p2))
      }
      
      
      
      if(length(dvec) > 0){###
        raw_mat[f, seq_along(dvec)] <- dvec
      }
      
      meanv[f] <- mean(dvec)
      sdv[f]   <- sd(dvec)###
      
    }
    
    time_ns <- ((seq_len(nf)-1) * frame_ps) / 1000
    
    tag <- paste0(pair[1], "_", pair[2], "_p", sh)
    
    raw_results[[tag]] <- raw_mat###
    
    pair_results[[tag]] <- list(
      
      summary = data.frame(
        time_ns = time_ns,
        mean    = meanv,
        sd      = sdv,
        system  = sysname,
        stringsAsFactors = FALSE
      ),
      
      raw = raw_results[[tag]]
    )
  }
  
  pair_results
}

###############################################################################
# PARALLEL EXECUTION
###############################################################################

cl <- makeCluster(N_CORES)

clusterEvalQ(cl, library(bio3d))

clusterExport(
  cl,
  varlist = c(
    "PAIRS","INPUT_DIR","FRAME_PS","FRAME_RANGE",
    "process_system","shift_chain","dist3",
    "get_atom_index","CACHE_DIR","cache_file"
  ),
  envir=environment()
)

res_list <- parLapply(
  cl,
  systems,
  function(s)
    process_system(
      sysname=s,
      pair_list=PAIRS,
      input_dir=INPUT_DIR,
      frame_ps=FRAME_PS,
      frame_range=FRAME_RANGE
    )
)
str(res_list[[1]][[1]], max.level=2)

stopCluster(cl)

names(res_list) <- systems

log_msg("Parallel processing completed.")

###############################################################################
# REORGANIZE BY PAIR
###############################################################################

pair_db     <- list()
pair_raw_db <- list()

for (i in seq_along(res_list)) {
  
  sys <- names(res_list)[i]
  rr  <- res_list[[i]]
  
  if (is.null(rr)) next
  
  for (nm in names(rr)) {
    
    if (is.null(pair_db[[nm]]))
      pair_db[[nm]] <- list()
    
    if (is.null(pair_raw_db[[nm]]))
      pair_raw_db[[nm]] <- list()
    
    pair_db[[nm]][[sys]]     <- rr[[nm]]$summary
    pair_raw_db[[nm]][[sys]] <- rr[[nm]]$raw
  }
}

###############################################################################
# PLOT FUNCTION
###############################################################################

plot_group <- function(db, picks, cols, labels, title_txt, outfile, ylim_use) {
  
  if (length(picks) == 0) return(NULL)
  
  tt <- db[[picks[1]]]$time_ns
  
  mat_mean <- sapply(picks, function(x) db[[x]]$mean)
  
  avg <- rowMeans(mat_mean, na.rm=TRUE)
  sdd <- apply(mat_mean, 1, sd, na.rm=TRUE)
  
  ymin <- ylim_use[1]
  ymax <- ylim_use[2]
  
  png(outfile, width=PNG_W, height=PNG_H, res=PNG_DPI)
  
  plot(0,0,type="n",
       xlim=c(min(tt), max(tt)),
       ylim=c(ymin,ymax),
       xlab="Time (ns)",
       ylab="Distance (Å)",
       main=title_txt)
  
  groups <- split(seq_along(picks), factor(cols, levels = unique(cols)))
  
  leg_txt <- c()
  leg_col <- c()
  
  
  for (cc in names(groups)) {
    
    idx <- groups[[cc]]
    
    g <- picks[idx]
    
    mat <- sapply(g, function(x) db[[x]]$mean)
    
    if (is.vector(mat)) mat <- matrix(mat, ncol=1)
    
    mm <- rowMeans(mat, na.rm=TRUE)
    ss <- apply(mat, 1, sd, na.rm=TRUE)
    
    if (SHOW_SD) {
      polygon(
        c(tt, rev(tt)),
        c(mm-ss, rev(mm+ss)),
        border=NA,
        col=adjustcolor(cc, alpha.f=0.20)
      )
    }
    
    lines(tt, mm, col=cc, lwd=2)
    
    leg_txt <- c(leg_txt, labels[length(leg_txt)+1])
    leg_col <- c(leg_col, cc)
    
  }
  
  par(xpd=TRUE)
  
  legend("topright",
         inset=c(.97,1),
         legend=leg_txt,
         col=leg_col,
         lwd=3,
         bty="n")
  
  dev.off()
}

###############################################################################
# BUILD FIGURES FOR EACH PAIR
###############################################################################

for (pair_name in names(pair_db)) {
  
  log_msg("Plotting pair: ", pair_name)
  
  db <- pair_db[[pair_name]]
  
  # ------------------------------------------------------------
  # Global Y scale for all figures of this pair
  # ------------------------------------------------------------
  
  vals_low  <- unlist(lapply(db, function(x) x$mean - x$sd))
  vals_high <- unlist(lapply(db, function(x) x$mean + x$sd))
  
  all_min <- quantile(vals_low,  0.02, na.rm=TRUE)
  all_max <- quantile(vals_high, 0.98, na.rm=TRUE)
  
  yr <- all_max - all_min
  
  ylim_pair <- c(
    all_min - 0.2 * yr,
    all_max + 0.2 * yr
  )
  
  out_pair <- file.path(OUT_BASE, pair_name)
  dir.create(out_pair, showWarnings=FALSE)
  
  meta <- lapply(names(db), classify_system)
  names(meta) <- names(db)
  
  meta_mut  <- vapply(meta, function(x) x$mut,  character(1))
  meta_cond <- vapply(meta, function(x) x$cond, character(1))
  
  sp <- strsplit(pair_name, "_")[[1]]
  r1 <- sp[1]
  r2 <- sp[2]
  
  for (mut in c("WT","W426A","W697A","Y441A")) {
    
    picks <- names(db)[meta_mut == mut]
    
    if (length(picks) == 0) next
    
    ord <- c("Apo","Cap","Heat")
    picks <- picks[order(match(meta_cond[picks], ord))]
    
    cols <- COL_COND[meta_cond[picks]]
    
    plot_group(
      db=db,
      picks=picks,
      cols=cols,
      labels=c("Apo","Cap","Heat")[seq_along(picks)],
      title_txt=paste("Conditions -", mut, "-", r1, "vs", r2),
      outfile=file.path(
        out_pair,
        paste0("Conditions_", mut, "_", r1, "_", r2, ".png")
      ),
      ylim_use = ylim_pair
    )
  }
  
  for (cond in c("Apo","Cap","Heat")) {
    
    picks <- names(db)[meta_cond == cond]
    
    if (length(picks) == 0) next
    
    ord <- c("WT","W426A","W697A","Y441A")
    picks <- picks[order(match(meta_mut[picks], ord))]
    
    cols <- COL_MUT[meta_mut[picks]]
    
    plot_group(
      db=db,
      picks=picks,
      cols=cols,
      title_txt=paste("Mutants -", cond, "-", r1, "vs", r2),
      labels = c("WT","W426A","W697A","Y441A")[seq_along(picks)],
      outfile=file.path(
        out_pair,
        paste0("Mutants_", cond, "_", r1, "_", r2, ".png")
      ),
      ylim_use = ylim_pair
    )
  }
}
###############################################################################
# DISTRIBUTION PLOTS
###############################################################################

dist_out <- file.path(OUT_BASE, "distributions")
dir.create(dist_out, showWarnings=FALSE)

dir.create(file.path(dist_out,"violins"), showWarnings=FALSE)
dir.create(file.path(dist_out,"histograms"), showWarnings=FALSE)
dir.create(file.path(dist_out,"variance"), showWarnings=FALSE)

mutant_order <- c("WT","W426A","W697A","Y441A")
cond_order   <- c("Apo","Cap","Heat")

for(pair_name in names(pair_raw_db)){
  
  raw_db <- pair_raw_db[[pair_name]]
  
  all_df <- data.frame()
  
  for(sys in names(raw_db)){
    
    cls <- classify_system(sys)
    
    mat <- raw_db[[sys]]
    
    nf <- nrow(mat)
    time_ns <- ((seq_len(nf)-1) * FRAME_PS)/1000
    
    keep <- which(time_ns >= 10)
    
    for(i in keep){
      
      vals <- as.numeric(mat[i,])
      
      vals <- vals[is.finite(vals)]
      
      if(length(vals)==0) next
      
      tmp <- data.frame(
        time_ns = time_ns[i],
        angle = vals,
        system = sys,
        mutant = cls$mut,
        condition = cls$cond,
        label = paste(cls$mut, cls$cond)
      )
      
      all_df <- rbind(all_df, tmp)
    }
  }
  
  if(nrow(all_df)==0) next
  
  all_df$label <- factor(
    all_df$label,
    levels=c(
      "WT Apo","WT Cap","WT Heat",
      "W426A Apo","W426A Cap","W426A Heat",
      "W697A Apo","W697A Cap","W697A Heat",
      "Y441A Apo","Y441A Cap","Y441A Heat"
    )
  )
  
  g1 <- ggplot(all_df,aes(x=label,y=angle,fill=condition)) +
    geom_violin(trim=FALSE,color="black",alpha=0.7) +
    scale_fill_manual(values=c(Apo="white",Cap="red",Heat="blue")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    labs(
      title=paste("Distribution:",pair_name),
      x="System",
      y="Distance (A)"
    )
  
  ggsave(
    file.path(dist_out,"violins",paste0("violin_",pair_name,".png")),
    g1,width=14,height=6,dpi=300
  )
  
  g2 <- ggplot(all_df,aes(x=angle,fill=condition)) +
    geom_histogram(
      bins=35,
      alpha=0.45,
      color="black",
      position="dodge"
    ) +
    scale_fill_manual(values=c(Apo="white",Cap="red",Heat="blue")) +
    facet_wrap(~mutant,ncol=2) +
    theme_bw() +
    labs(
      title=paste("Histograms:",pair_name),
      x="Distance (A)",
      y="Count"
    )
  
  binwidth_global <- diff(range(all_df$angle, na.rm=TRUE)) / 35
  
  for(m in mutant_order){
    
    for(cond in cond_order){
      
      sub_df <- all_df %>%
        filter(mutant == m,
               condition == cond) %>%
        filter(is.finite(angle))
      
      if(nrow(sub_df) > 10){
        
        xmin <- min(sub_df$angle, na.rm=TRUE)
        xmax <- max(sub_df$angle, na.rm=TRUE)
        
        if(is.finite(xmin) &&
           is.finite(xmax) &&
           xmin < xmax){
          
          fit <- try(
            Mclust(sub_df$angle, G=1:2, verbose=FALSE),
            silent=TRUE
          )
          
          xs <- seq(xmin, xmax, length.out=300)
          
          use_two <- FALSE
          
          if(!inherits(fit, "try-error")){
            
            if(length(fit$parameters$pro) == 2){
              
              w1 <- fit$parameters$pro[1]
              w2 <- fit$parameters$pro[2]
              
              mu1 <- fit$parameters$mean[1]
              mu2 <- fit$parameters$mean[2]
              
              sd1 <- sqrt(fit$parameters$variance$sigmasq[1])
              sd2 <- sqrt(fit$parameters$variance$sigmasq[2])
              
              sep_test <- abs(mu1 - mu2) / mean(c(sd1, sd2))
              
              if(is.finite(w1) &&
                 is.finite(w2) &&
                 is.finite(sd1) &&
                 is.finite(sd2) &&
                 is.finite(mu1) &&
                 is.finite(mu2) &&
                 is.finite(sep_test)){
                
                if(min(w1, w2) >= 0.10 &&
                   sep_test >= 0.50){
                  use_two <- TRUE
                }
              }
            }
          }
          
          if(use_two){
            
            ys <- fit$parameters$pro[1] *
              dnorm(xs,
                    mean=fit$parameters$mean[1],
                    sd=sqrt(fit$parameters$variance$sigmasq[1])) +
              
              fit$parameters$pro[2] *
              dnorm(xs,
                    mean=fit$parameters$mean[2],
                    sd=sqrt(fit$parameters$variance$sigmasq[2]))
            
          } else {
            
            mu <- mean(sub_df$angle, na.rm=TRUE)
            sdv <- sd(sub_df$angle, na.rm=TRUE)
            
            if(!is.finite(sdv) || sdv <= 0){
              sdv <- 1
            }
            
            ys <- dnorm(xs, mean=mu, sd=sdv)
          }
          
          dens_df <- data.frame(
            angle = xs,
            density = ys,
            mutant = m
          )
          
          line_col <- if(cond == "Apo"){
            "black"
          } else if(cond == "Cap"){
            "red"
          } else {
            "blue"
          }
          
          g2 <- g2 +
            geom_line(
              data = dens_df,
              aes(
                x = angle,
                y = density *
                  nrow(sub_df) *
                  binwidth_global
              ),
              inherit.aes = FALSE,
              color = line_col,
              linewidth = 1
            )
        }
      }
    }
  }
  ggsave(
    file.path(dist_out,"histograms",paste0("hist_",pair_name,".png")),
    g2,width=10,height=8,dpi=300
  )
  
  ###########################################################
  # VARIANCE DISPERSION PLOT
  ###########################################################
  
  var_df <- all_df %>%
    group_by(system, mutant, condition, label) %>%
    summarise(
      var_angle = var(angle, na.rm=TRUE),
      .groups="drop"
    )
  
  var_sum <- var_df %>%
    group_by(mutant, condition, label) %>%
    summarise(
      mean_var = mean(var_angle, na.rm=TRUE),
      se_var   = sd(var_angle, na.rm=TRUE) / sqrt(n()),
      .groups="drop"
    )
  
  var_sum$label <- factor(
    var_sum$label,
    levels=c(
      "WT Apo","WT Cap","WT Heat",
      "gap1",
      "W426A Apo","W426A Cap","W426A Heat",
      "gap2",
      "W697A Apo","W697A Cap","W697A Heat",
      "gap3",
      "Y441A Apo","Y441A Cap","Y441A Heat"
    )
  )
  
  g3 <- ggplot(
    var_sum,
    aes(
      x = label,
      y = mean_var,
      color = condition,
      group = label
    )
  ) +
    
    geom_point(size=3) +
    
    geom_errorbar(
      aes(
        ymin = mean_var - se_var,
        ymax = mean_var + se_var
      ),
      width=0.25,
      linewidth=0.7
    ) +
    
    scale_x_discrete(
      drop=FALSE,
      labels=c(
        "WT Apo","WT Cap","WT Heat",
        "",
        "W426A Apo","W426A Cap","W426A Heat",
        "",
        "W697A Apo","W697A Cap","W697A Heat",
        "",
        "Y441A Apo","Y441A Cap","Y441A Heat"
      )
    ) +
    
    scale_color_manual(
      values=c(
        "Apo"="black",
        "Cap"="red",
        "Heat"="blue"
      )
    ) +
    
    labs(
      title=paste("Dispersion (variance):", pair_name),
      x="System",
      y="Mean variance ± SE"
    ) +
    
    theme_bw() +
    
    theme(
      axis.text.x=element_text(
        angle=45,
        hjust=1
      ),
      legend.position="none"
    )
  
  ggsave(
    filename=file.path(
      dist_out,
      "variance",
      paste0("variance_", pair_name, ".png")
    ),
    plot=g3,
    width=14,
    height=6,
    dpi=300
  )
}
###############################################################################
# END
###############################################################################

log_msg("Finished successfully.")
log_msg("===========================================================")