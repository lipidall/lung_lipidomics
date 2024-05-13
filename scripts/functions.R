library(ggplot2)
library(ggbiplot)
library(plyr)
library(extrafont)
library(magrittr)
library(gridExtra)
library(Mfuzz)
library(showtext)
library(stringr)
library(RColorBrewer)
library(cowplot)
library(ComplexHeatmap)
library(patchwork)
library(dplyr)
library(factoextra)
library(extrafont)
library(ggpubr)
library(xlsx)
library(data.table)
library(yyplot)

#set conflicts 
conflicted::conflicts_prefer(dplyr::arrange)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(dplyr::mutate)
conflicted::conflicts_prefer(base::union)

# coln is colnames of data.frame
# return a data.frame with name and class
getAttrs <- function(coln) {
    # name, class, #carbon, #double bond, color of class
    out <- cbind(
        name = as.character(coln),
        class = sapply(c(
            '^PA(?![A-z])', '^LPA(?![A-z])',
            '^PC(?![A-z])', '^LPC(?![A-z])', '^LysoPC(?![A-z])', 'PC$',
            '^PE(?![A-z])', '^LPE(?![A-z])', 'HETE-PE$', '^OX-PE$',
            '^PG(?![A-z])', '^LPG(?![A-z])',
            '^PI(?![A-z])', '^LPI(?![A-z])', 
            '^PS(?![A-z])', '^LPS(?![A-z])',
            '^CL(?![A-z])',
            '^Cer(?!1P)', '^Cer1P', 
            '^GluCer(?![A-z])',
            '^GalCer(?![A-z])',
            '^LacCer(?![A-z])',
            '^PhytoCer(?![A-z])',
            '^SM(?![A-z])', '^Lyso-SM(?![A-z])', '^LSM$',
            '^Sph(?![A-z])', '^S1P(?![A-z])',
            '^Gb3(?![A-z])', '^SL(?![A-z])',
            '^Cho$', '^CE(?![A-z])',
            '^DAG(?![A-z])', '^TAG(?![A-z])',  
            '^FFA(?![A-z])', 
            '^PAHSA(?![A-z])', '^FAHFA$',
            'carnitine$','carnitine-aq',
            '^MGDG(?![A-z])', '^DGDG(?![A-z])',
            '^GM1(?![A-z])', '^GM2(?![A-z])', '^GM3(?![A-z])', 
            '^LBPA(?![A-z])', 
            '^WE(?![A-z])', 
            '^BMP(?![A-z])',
            'CoA$'
            ),
            function(x) {
                res <- stringr::str_extract(coln, x)
                ## assign Cho to CE
                if (x == '^Cho$')
                    res[which(res == 'Cho')] <- 'CE'
                if (x == '^LysoPC(?![A-z])')
                    res[which(res == 'LysoPC')] <- 'LPC'
                if (x == '^Lyso-SM(?![A-z])')
                    res[which(res == 'Lyso-SM')] <- 'LSM'
                if (x == 'carnitine$')
                    res[which(res == 'carnitine')] <- 'Acylcarnitine'
                if (x == 'carnitine-aq')
                  res[which(res == 'carnitine-aq')] <- 'Acylcarnitine'
                if (x == '^PE(?![A-z])')
                    res[res == 'PE' & stringr::str_detect(coln, '\\[O\\]')] <- 'OX-PE'
                if (x == 'HETE-PE$')
                    res[which(res == 'HETE-PE')] <- 'OX-PE'
                if (x == '^PAHSA(?![A-z])')
                    res[which(res == 'PAHSA')] <- 'FAHFA'
                if (x == 'PC$')
                    res[which(res == 'PC')] <- 'OXPC'
                res
            }) %>% apply(1, function(x) {
                x[!is.na(x)][1]
            })
    ) %>% data.frame %>% dplyr::mutate(
        ## sphingolipid and phospholipids are named differently
        carbon = ifelse(class %in% c('SL', 'Cer', 'GM1', 'GM2', 'GM3',
                                     'SM', 'GluCer', 'GalCer',
                                     'LacCer', 'Gb3'),
                        stringr::str_extract(
                            stringr::str_extract(name, '[0-9]{1,3}:[0-9][h)]?$'),
                            '[0-9]+(?=:)'),
                        stringr::str_extract(name, '[0-9]+(?=:)')) %>% as.numeric,
        dbond = ifelse(class %in% c('SL', 'Cer', 'GM1', 'GM2', 'GM3', 'SM', 'GluCer',
                                    'LacCer', 'Gb3'),
                       stringr::str_extract(
                           stringr::str_extract(name, '[0-9]{1,3}:[0-9][h)]?$'),
                           '(?<=:)[0-9]+'),
                       stringr::str_extract(coln, '(?<=:)[0-9]+')) %>% as.numeric
    )  %>% dplyr::mutate(
        carbon_cat = ifelse(
            class %in% c('PC', 'PE', 'PI', 'PS', 'PA', 'PG', 'DAG', 'LBPA'),
            ifelse(carbon >= 44, 'Very long', ifelse(carbon >= 32, 'Long', 'Short')),
            ifelse(class %in% c('TAG'),
                   ifelse(carbon >= 66, 'Very long',
                          ifelse(carbon >= 48, 'Long', 'Short')),
                   ifelse(carbon >= 22, 'Very long',
                          ifelse(carbon >= 16, 'Long', 'Short'))
            )),
        dbond_cat = ifelse(
            class %in% c('PC', 'PE', 'PI', 'PS', 'PA', 'PG', 'DAG', 'LBPA'),
            ifelse(dbond > 4, 'Poly-unsaturated',
                   ifelse(dbond > 0, 'Mono/di-unsaturated', 'Saturated')),
            ifelse(class %in% c('TAG'),
                   ifelse(dbond > 6, 'Poly-unsaturated',
                          ifelse(dbond > 1, 'Mono/di-unsaturated', 'Saturated')),
                   ifelse(dbond > 2, 'Poly-unsaturated',
                          ifelse(dbond > 0, 'Mono/di-unsaturated', 'Saturated')))
        )
    ) %>% dplyr::mutate(
        carbon_cat = factor(
            carbon_cat,
            levels = c('Very long', 'Long', 'Short')
        ),
        dbond_cat = factor(
            dbond_cat,
            levels = c('Poly-unsaturated', 'Mono/di-unsaturated', 'Saturated'))
    )
    # color
    # remove black from color sequence
    out$color <- WGCNA::labels2colors(
        as.numeric(as.factor(out$class)),
        colorSeq = WGCNA::standardColors(50)[-7])
    # rowname
    rownames(out) <- out$name
    out
}

# In lipid classes which contain 1 carbon chain
# Short: [0, 16)
# Long:  [16, 22)
# Very long: [22, inf)
#
# Saturated: 0
# Mono/di-unsaturated: 1-2
# Poly-unsaturated: > 2
#
# In lipid classes which contain 2 carbon chains: 'PC', 'PE', 'PI', 'PS', 'PA', 'PG', 'DAG', 'LBPA'
# Short: [0, 32)
# Long:  [32, 44)
# Very long: [44, inf)
#
# Saturated: 0
# Mono/di-unsaturated: [1, 4]
# Poly-unsaturated: [5, Inf)
#
# In lipid classes which contain 3 carbon chains: TAG
# Short: [0, 48)
# Long:  [48, 66)
# Very long: [66, inf)
#
# Saturated: 0
# Mono-unsaturated: [1, 6]
# Poly-unsaturated: [7, Inf)



#' Generate discrete color palette as in ggplot2
#'
#' @param n Number of colors
#'
#' @return Character vector (length n) of color hex code
#' @export
#'
#' @examples
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


font_family = 'Arial'
fontSize_big=7
fontSize_small=5
#using ggpubr package set font size.
theme_parameter<- font("xylab", size = fontSize_big,family=font_family)+ 
  font("xy.text", size = fontSize_small,family=font_family)+
  font("legend.text", size = fontSize_small,family=font_family)+
  font("legend.title", size = fontSize_big,family=font_family)


batch_barplot <- function(d, g, dh, 
                          posthoc, 
                          add_sig = T, add_point = T,
                          x_angle = 0, h_just = 0.5, v_just = 0, 
                          conf.int = 0.95, 
                          ggtheme = theme_bw(), font_family = 'Arial') {
  # d: data.frame
  # g: group, length is the same as row number of d
  # dh: hypothesis test resutls
  # conf.int: confidence interval
  d.colnames <- colnames(d)
  d$group = g
  mult <- qnorm((1 + conf.int)/2)
  if (abs(conf.int - 0.68) < 0.01)
    mult = 1
  
  p_list <- lapply(d.colnames, function(var.i) {
    
    p <- ggplot(
      d ,
      aes(x = group,
          y = eval(parse(         
            text = paste0('`', var.i, '`')   #use eval function to select column.
          )))) +
      stat_summary(              #statistical average.
        fun = mean,
        geom = 'bar',
        fill = '#404040', 
        color = '#404040'
      ) +
      stat_summary(
        fun.data = mean_se,
        fun.args = list(mult = mult),
        geom = 'errorbar',
        width = 0.2
      )
    
    if (add_point) {
      p <- p +
        geom_jitter(position = position_jitter(width = 0.1))
    }
    
    if (add_sig) {
      if (nlevels(g) == 1) {    # if g only has one group.
        
      } else if(nlevels(g) > 2) { #if g contains more than two samples.
        p_matrix <- matrix(
          1, nrow = nlevels(g), ncol = nlevels(g),
          dimnames = list(levels(g), levels(g))    
        )
        ## multcompLetters work on lower.tri
        p_matrix[lower.tri(p_matrix)] <- dh[var.i, stringr::str_detect(colnames(dh), posthoc)]
        
        tryCatch(
          {
            sig_letters <<- multcompView::multcompLetters(
              p_matrix, threshold = 0.05
            )$Letters
            ## compute y limits
            ymax <- vaggregate(d[, var.i], g, mean_se) %>%
              unlist %>% max
            d_sig <<- data.frame(group = names(sig_letters), 
                                 y = 1.2 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
                                 label = sig_letters)                    
            
            suppressWarnings(
              p <- p + 
                geom_text(
                  data = d_sig,
                  aes(x = group, y = y, label = label),
                  family = font_family) + 
                scale_y_continuous(
                  expand = c(0.1, 0)
                )
            )
          },
          error = function(e) {
            print(paste(var.i, ':', e$message))
          }
        )
        
      } else {
        ## 2 groups
        p_col <- ifelse(posthoc == 'Parametric',
                        'parametric pvalue', 'non-parametric pvalue')
        d_sig <- data.frame(
          start = levels(g)[1],
          end = levels(g)[2],
          y = 1.1 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
          label = formatC(dh[var.i, p_col], digits = 2)
        )
        
        suppressWarnings(
          p <- p + ggsignif::geom_signif(
            data = d_sig,
            aes(xmin = start, xmax = end, y_position = y,
                annotations = label),
            manual = T,
            tip_length = min(
              0.01 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
              0.01),
            family = font_family
          ) + scale_y_continuous(expand = c(0.15, 0))
        )
      }
      
    }
    
    
    p <- p + labs(x = '', title = var.i, y ="MFP") +
      ggtheme + 
      theme(text = element_text(size = 8,family =  font_family))

    p
  })
  
  return(p_list)
}


hypothesisTest <- function(d1, g, pair = NA, method = union('none', p.adjust.methods)) {
  ## d1: data
  ## g: group
  method <- match.arg(method)
  
  if (!(length(pair) == 1 && is.na(pair))) {
    pair <- as.factor(pair)
    testthat::expect_equal(
      length(g), length(pair)
    )
    testthat::expect(
      all(table(g, pair) == 1),
      "Each pair must be in two groups."
    )
  }
  if (!is.data.frame(d1)) d1 <- as.data.frame(d1)
  if (!is.factor(g)) g <- as.factor(g)
  
  testthat::expect_equal(
    nrow(d1), length(g)
  )
  
  if (nlevels(g) == 1) {
    message('There is only 1 group. Return NULL.')
    return(NULL)
  }
  if (is.na(pair)) {
    ## not paired design
    test_result <- sapply(d1, function(x) {
      ## combination of group pairs
      g_pair <- combn(levels(g), 2)
      
      g_fold <- apply(g_pair, 2, function(pair.i) {
        mean(x[g == pair.i[2]], na.rm = T) /
          mean(x[g == pair.i[1]], na.rm = T)
      })
      
      names(g_fold) <- apply(g_pair, 2, function(x) {
        paste0('Fold: ', paste0(rev(x), collapse = '/'))
      })
      
      ## welch's ANOVA
      aov.out <- tryCatch(
        oneway.test(x ~ g),
        error = function(e) {
          ## Number of replicate is 1
          print(e)
          list(p.value = NA)
        }
      )
      
      ## non-parametric
      ## wilcox when #groups == 2
      ## kruskal-wallis when #groups > 2
      if (nlevels(g) > 2)
        nonparam.out <- kruskal.test(x ~ g)
      else
        nonparam.out <- wilcox.test(x ~ g)
      
      if (nlevels(g) > 2) {
        ## posthoc.tgh in functions.R
        posthoc.out <- userfriendlyscience::posthocTGH(
          x, g, 
          p.adjust = ifelse(method %in% c('bh', 'by'), toupper(method), method), 
          formatPvalue = F
        )
        
        ## replace - by :
        posthoc.names <- apply(
          combn(levels(g),2), 2, 
          function(x) {paste0(rev(x), collapse=':')})
        
        posthoc.out$output$tukey <- as.matrix(posthoc.out$output$tukey)
        rownames(posthoc.out$output$tukey) <- posthoc.names
        posthoc.out$output$games.howell <- as.matrix(posthoc.out$output$games.howell)
        rownames(posthoc.out$output$games.howell) <- posthoc.names
        
        ## dunn's test
        ## p value adjust
        
        dunn.out <- FSA::dunnTest(x ~ g, method = method)
        dunn.out$res[, 'Comparison'] <- stringr::str_replace(
          dunn.out$res[, 'Comparison'],
          ' - ', ':'
        ) ## make it same as Tukey
        ## dunnTest does not honor level order
        dunn.out$res <- dunn.out$res[
          match(
            combn(levels(g), 2, function(x) { paste0(sort(x), collapse=":") }),
            dunn.out$res[, 'Comparison']
          ),
        ]
        ## reverse the order
        paste0Rev <- function(x) {paste0(rev(x), collapse = ':')}
        dunn.out$res$Comparison <- combn(levels(g), 2, paste0Rev)
        
        res <- c(`parametric pvalue` = aov.out$p.value,
                 posthoc.out$output$tukey[, ifelse(method == 'none', 'p', 'p.adjusted')],
                 posthoc.out$output$games.howell[, ifelse(method == 'none', 'p', 'p.adjusted')],
                 `non-parametric pvalue` = nonparam.out$p.value,
                 dunn.out$res[, 'P.adj'],
                 g_fold)
        
        ## Need to set names because column of data.frame
        ## is not named
        ## number of pair-wise comparison
        n_comb <- nrow(posthoc.out$output$tukey)
        names(res)[2:(1 + 2 * n_comb)] <- paste0(
          rep(c('TukeyHSD: ', 'Games-Howell: '), each = n_comb),
          names(res)[2:(1 + 2 * n_comb)]
        )
        
        names(res)[(3 + 2 * n_comb):(2 + 3 * n_comb)] <- paste0(
          rep('Dunn: ', n_comb),
          dunn.out$res$Comparison
        )
        
      } else {
        res <- c(`parametric pvalue` = aov.out$p.value, 
                 `non-parametric pvalue` = nonparam.out$p.value,
                 g_fold)
      }
      
      res 
    }) %>% t
  } else {
    ## paired design
    ## paired t-test
    test_result <- sapply(d1, function(x) {
      order_pair <- order(pair)
      x_order <- x[order_pair]
      g_order <- g[order_pair]
      g1 <- which(g_order == levels(g)[1])
      g2 <- which(g_order == levels(g)[2])
      
      message('Fold change is the median of all pair-wise fold changes for paired design.')
      g_fold <- median(x_order[g2] / x_order[g1], na.rm = T)
      
      paired.t.out <- tryCatch(
        t.test(x_order[g1], x_order[g2], paired = T),
        error = function(e) {
          print(e)
          list(p.value = NA)
        }
      )
      
      paired.wilcox.out <- wilcox.test(x_order[g1], x_order[g2], paired = T)
      c(paired.t.out$p.value, 
        paired.wilcox.out$p.value,
        g_fold)
    }) %>% t
    
    colnames(test_result) <- c('parametric pvalue', 'non-parametric pvalue', paste0('Fold: ', levels(g)[2], '/', levels(g)[1]))
    
  }
  return(test_result)
}
