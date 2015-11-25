#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
plotSurvival <- function(sfit, sdiff, main,returns = FALSE,
                 xlabs = "Time", ylabs = "Probability",
                 ystratalabs = NULL, ystrataname = NULL,
                 timeby = 400,
                 pval = TRUE, ...) {
  
  if(is.null(ystratalabs)) {
    ystratalabs <- as.character(levels(summary(sfit)$strata))
  }
  m <- max(nchar(ystratalabs))
  if(is.null(ystrataname)) ystrataname <- "Treatment Pathway"
  times <- seq(0, max(sfit$time), by = timeby)
  .df <- data.frame(time = sfit$time, n.risk = sfit$n.risk,
                    n.event = sfit$n.event, surv = sfit$surv, strata = summary(sfit, censored = T)$strata,
                    upper = sfit$upper, lower = sfit$lower)
  levels(.df$strata) <- ystratalabs
  zeros <- data.frame(time = 0, surv = 1, strata = factor(ystratalabs, levels=levels(.df$strata)),
                      upper = 1, lower = 1)
  .df <- rbind.fill(zeros, .df)
  d <- length(levels(.df$strata))
  p <- ggplot(.df, aes(time, surv, group = strata)) +
    geom_step(aes(linetype = strata), size = 0.7) +
    theme_bw() +
    theme(axis.title.x = element_text(vjust = 0.5)) +
    scale_x_continuous(xlabs, breaks = times, limits = c(0, max(sfit$time))) +
    scale_y_continuous(ylabs, limits = c(0, 1)) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position = c(ifelse(m < 10, .28, .2), ifelse(d < 4, .25, .35))) +
    theme(legend.key = element_rect(colour = "grey")) +
    labs(linetype = ystrataname) +
    theme(plot.margin = unit(c(0, 1, .5, ifelse(m < 10, 1.5, 2.5)), "lines")) +
    theme(text = element_text(size=12))+
    ggtitle(main) 
  if(pval) {
    pval = as.numeric(1 - pchisq(sdiff$chisq, length(sdiff$n) - 1))
    pval <- format(round(as.numeric(as.character(pval)), 3), nsmall = 3)
    HR = (sdiff$obs[2]/sdiff$exp[2])/(sdiff$obs[1]/sdiff$exp[1])
    HRR <- format(round(as.numeric(as.character(HR)), 2), nsmall = 2)
    up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdiff$exp[2]+1/sdiff$exp[1]))
    up95 <- format(round(as.numeric(as.character(up95)), 2), nsmall = 2)
    low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdiff$exp[2]+1/sdiff$exp[1]))
    low95 <- format(round(as.numeric(as.character(low95)), 2), nsmall = 2)
    pvaltxt <- paste("p =", pval)
    HRText <- paste("HR = ",HRR," CI[",low95,",",up95,"]",sep="")
    p <- p + annotate("text", x = 0.11 * max(sfit$time), y = 0.1, label = pvaltxt)
    p <- p + annotate("text", x = 0.2 * max(sfit$time), y = 0.02, label = HRText)
  }
  print(p)
  if(returns) return(p) 
}
#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
ggsurv <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                   cens.col = 'red', lty.est = 1, lty.ci = 2,
                   cens.shape = 3, back.white = F, xlab = 'Time',
                   ylab = 'Survival', main = ''){
  
  library(ggplot2)
  strata <- ifelse(is.null(s$strata) ==T, 1, length(s$strata))
  stopifnot(length(surv.col) == 1 | length(surv.col) == strata)
  stopifnot(length(lty.est) == 1 | length(lty.est) == strata)
  
  ggsurv.s <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = ''){
    
    dat <- data.frame(time = c(0, s$time),
                      surv = c(1, s$surv),
                      up = c(1, s$upper),
                      low = c(1, s$lower),
                      cens = c(0, s$n.censor))
    dat.cens <- subset(dat, cens != 0)
    
    col <- ifelse(surv.col == 'gg.def', 'black', surv.col)
    
    pl <- ggplot(dat, aes(x = time, y = surv)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(col = col, lty = lty.est)
    
    pl <- if(CI == T | CI == 'def') {
      pl + geom_step(aes(y = up), color = col, lty = lty.ci) +
        geom_step(aes(y = low), color = col, lty = lty.ci)
    } else (pl)
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  
  ggsurv.m <- function(s, CI = 'def', plot.cens = T, surv.col = 'gg.def',
                       cens.col = 'red', lty.est = 1, lty.ci = 2,
                       cens.shape = 3, back.white = F, xlab = 'Time',
                       ylab = 'Survival', main = '') {
    n <- s$strata
    
    groups <- factor(unlist(strsplit(names
                                     (s$strata), '='))[seq(2, 2*strata, by = 2)])
    gr.name <-  unlist(strsplit(names(s$strata), '='))[1]
    gr.df <- vector('list', strata)
    ind <- vector('list', strata)
    n.ind <- c(0,n); n.ind <- cumsum(n.ind)
    for(i in 1:strata) ind[[i]] <- (n.ind[i]+1):n.ind[i+1]
    
    for(i in 1:strata){
      gr.df[[i]] <- data.frame(
        time = c(0, s$time[ ind[[i]] ]),
        surv = c(1, s$surv[ ind[[i]] ]),
        up = c(1, s$upper[ ind[[i]] ]),
        low = c(1, s$lower[ ind[[i]] ]),
        cens = c(0, s$n.censor[ ind[[i]] ]),
        group = rep(groups[i], n[i] + 1))
    }
    
    dat <- do.call(rbind, gr.df)
    dat.cens <- subset(dat, cens != 0)
    
    pl <- ggplot(dat, aes(x = time, y = surv, group = group)) +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      geom_step(aes(col = group, lty = group))
    
    col <- if(length(surv.col == 1)){
      scale_colour_manual(name = gr.name, values = rep(surv.col, strata))
    } else{
      scale_colour_manual(name = gr.name, values = surv.col)
    }
    
    pl <- if(surv.col[1] != 'gg.def'){
      pl + col
    } else {pl + scale_colour_discrete(name = gr.name)}
    
    line <- if(length(lty.est) == 1){
      scale_linetype_manual(name = gr.name, values = rep(lty.est, strata))
    } else {scale_linetype_manual(name = gr.name, values = lty.est)}
    
    pl <- pl + line
    
    pl <- if(CI == T) {
      if(length(surv.col) > 1 && length(lty.est) > 1){
        stop('Either surv.col or lty.est should be of length 1 in order
             to plot 95% CI with multiple strata')
      }else if((length(surv.col) > 1 | surv.col == 'gg.def')[1]){
        pl + geom_step(aes(y = up, color = group), lty = lty.ci) +
          geom_step(aes(y = low, color = group), lty = lty.ci)
      } else{pl +  geom_step(aes(y = up, lty = group), col = surv.col) +
               geom_step(aes(y = low,lty = group), col = surv.col)}
    } else {pl}
    
    
    pl <- if(plot.cens == T & length(dat.cens) > 0){
      pl + geom_point(data = dat.cens, aes(y = surv), shape = cens.shape,
                      col = cens.col)
    } else if (plot.cens == T & length(dat.cens) == 0){
      stop ('There are no censored observations')
    } else(pl)
    
    pl <- if(back.white == T) {pl + theme_bw()
    } else (pl)
    pl
  }
  pl <- if(strata == 1) {ggsurv.s(s, CI , plot.cens, surv.col ,
                                  cens.col, lty.est, lty.ci,
                                  cens.shape, back.white, xlab,
                                  ylab, main)
  } else {ggsurv.m(s, CI, plot.cens, surv.col ,
                   cens.col, lty.est, lty.ci,
                   cens.shape, back.white, xlab,
                   ylab, main)}
  pl
}
#'This function generates keyword and ignore lists based on the expansion of
#'concepts.
#'
#'@description Given any given concept_id or string of text this function
#'generates keyword and ignore lists based on the expansion of concepts (looking
#'at their synonyms).
#'
#'@param connection    The connection to the database server.
#'@param aphroditeConceptName  The string of text / concept name to use.
#'@param schema        The database schema being used.
#'@param dbms          The target DBMS for SQL to be rendered in.
#'
#'@details Takes the aphroditeConceptName looks for synonyms and builds a list
#'of related concepts using the vocabulary hierarchies
#'
#'@return A list with two elements: a list of positive keywords found
#'(keywordlist_ALL), and a list of ignore keywords (ignorelist_ALL)
#'
#' @examples \dontrun{
#'
#'wordLists <- buildKeywordList(conn, aphrodite_concept_name, cdmSchema, dbms)
#'
#' }
#'
#'@export
executeSQL_ro <- function (connection, schema, query, targetDBMS) {
  ### JMB: Added nasty hack to remove semi colon for Oracle instances ##
  ### Will probably remove after more testing - with sending no semi-colons on PostGre ##
  if (tolower(c(targetDBMS))=="oracle") {
    query= substr(query,1,nchar(query)-1)
  } #Not oracle? do nothing
  renderedSql <- renderSql(query, cdmSchema=schema)$sql
  translatedSql <- translateSql(renderedSql, sourceDialect = "sql server", targetDialect = targetDBMS)$sql
#  
  queryResults <- querySql(connection,translatedSql)
#  
  names(queryResults) <- tolower(names(queryResults))  #Hack added to make the field names lowercase - should/might be removed later
#  
  return(queryResults)
}