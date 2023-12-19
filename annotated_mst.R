#pacman is a progress management tool
pacman::p_load(dplyr,
               tidyr,
               ggplot2,
               catR,
               mstR,
               progress)


N_REP <- 100 # no. of replications per theta
THETA <- seq(-3, 3, 0.5) # theta ranging from -3 to 3 by 0.5 inc
TEST_LEN <- 42 # total test length
set.seed(123)

# Item generation
n_items <- 1500
# catR function to generate an itembank 
item_bank <- genDichoMatrix(items = n_items,
                            model = "3PL",
                            #sample a's from log normal distribution
                            aPrior = c("lnorm", 0, 0.25), 
                            #sample b's from normal distribution
                            bPrior = c("norm", 0, 1),
                            #sample c's from uniform distribution
                            cPrior = c("unif", 0.1, 0.2))

# vector of items to use as an ID variable, from 1 to number of items
id <- 1:n_items
# add the variable ID vector onto the simulated item bank
item_bank <- cbind(id, item_bank)

# Probability and Information
  #create a probability function that reads in theta and a set of items

  prob <- function(theta, items) {
      #outer product (-) of each b and theta inc
      p <- outer(items$b, theta, "-") 
      #finish the probabilities by applying the 3pl binary response function
      p <- items$c + (1 - items$c) / (1 + exp(1.7 * items$a * p)) 
      #return matrix 'p' containing probabilities with
      # row = items, col = theta
      p 
  }

info <- function(theta, items) {
    #outer product (-) of each b and theta inc
    p <- prob(theta, items)
    #finish the information matrix by applying the 3pl binary response function
    i <- (1.7 * items$a * (p - items$c) / (1 - items$c))^2 * (1 - p) / p
    #return matrix 'i' containing probabilities with
    # row = items, col = theta
    i 
}

# MST Setup Function
#function incorporates design, items per stage & assembly direction (forward/backward)
# 'design' argument is character string "1-3-3" or ""1-3-4"
# 'items_per_stage' is a character string "increasing", "decreasing" or "equal"
# 'assembly_priority' is a character string 'forward' or 'backward'

mst_setup <- function(design, items_per_stage, assembly_priority) {
  
    #load 'theta' object with ranges based on value of design (1-3-3 or 1-3-4)
    theta <- switch(design,
                    "1-3-3" = c(0, -1, 0, 1, -1, 0, 1),
                    "1-3-4" = c(0, -1, 0, 1, -1.25, -0.75, 0.75, 1.25))
    
    #load 'len' object with vector of item lengths based on 'items_per_stage' and overall test length
    len <- switch(items_per_stage,
                  "increasing" = TEST_LEN * c(1/6, 1/3, 1/2),
                  "equal"= TEST_LEN * c(1/3, 1/3, 1/3),
                  "decreasing" = TEST_LEN * c(1/2, 1/3, 1/6))
    
    #load 'design' object as a numeric vector, converted from character 
    design <- strsplit(design, split = "[-]") %>% unlist() %>% as.integer()
    
    #load 'n_stages' object with length of 'design' numeric vector
    n_stages <- length(design)
    
    #load 'n_modules' object with overall number (i.e., sum) of 'design' numeric vector
    n_modules <- sum(design)

    #instantiate a new object, 'module'
    module <- NULL
    
    #across each stage, 's'....
    for (s in 1:n_stages)
        #and across each design, 'm'...
        for (m in 1:design[s])
            #add a row with stage, design and length
            module <- rbind(module, c(stage = s, module = m, len = len[s]))
    #create a data frame with our module, id's (1:nrow), and theta vector
    module <- data.frame(cbind(module, index = 1:nrow(module), theta = theta))

    #create and export object, 'x' that includes: 
    x <- list(design = design, #design numeric vector
              n_panels = 5,    #number of panels
              n_stages = n_stages, #number of stages
              n_modules = n_modules, #number of modules
              module = module, #model matrix
              items_per_stage = items_per_stage, # text description for plotting later
              assembly_priority = assembly_priority)
    x #output/return x object
}

# MST Assemble Function
#function incorporates 'x' (object output from 'mst_setup' function),
#and assembly priority, e.g., forward/backward

mst_assemble <- function(x, assembly_priority) {
    #add an 'items' object to the 'x' object that captures
    #module length * number of panels as number of items, 
    #and columns that indicate item characteristics, IRT parms/administration
    x$items <- data.frame(matrix(nrow = sum(x$module$len) * x$n_panels,
                                 ncol = 10,
                                 dimnames = list(c(),
                                                 c("id", "a", "b", "c", "d",
                                                   "stage", "module", "panel",
                                                   "index", "form"))))
    #calculate information from simulated item bank at module theta increments
    info <- info(theta = x$module$theta, items = item_bank)
    
    #logic to how tests are constructed based on forward or backward assembly
    #RC Addition: create random and spiral assembly methods
    if (x$assembly_priority == "Forward")
        stage_order <- 1:x$n_stages
    if (x$assembly_priority == "Spiral"){
      stage_order <- c(ceiling(x$n_stages/2),c(c(-1,1)%*%t(c(1:floor(x$n_stages/2))))+ceiling(x$n_stages/2))
      stage_order <- stage_order[stage_order %in% 1:x$n_stages]}
    if (x$assembly_priority == "Random")
      stage_order <- sample(1:x$n_stages, x$n_stages)
    else
        stage_order <- x$n_stages:1

    temp <- 1

    #for each stage_order, 'i'....
    for (i in stage_order) {
        # randomize which module and which panel to assemble
        n_modules <- x$design[i] # no. of modules in stage i
        len <- with(x$module, len[stage == i])[1] # no. of items in stage i
        order <- cbind(rep(1:n_modules, times = len * x$n_panels),
                       rep(1:x$n_panels, each = len * n_modules))
        order <- order[sample(1:nrow(order)), ]
        for (j in 1:nrow(order)) {
            index <- with(x$module, index[stage == i & module == order[j, 1]])
            # select the item that maximizes the information at the specified theta
            max_info <- which.max(info[item_bank$id, index])
            selected_item <- item_bank[max_info, ]
            x$items[temp, ] <- c(selected_item, # id, a, b, c, d
                                 stage = i,
                                 module = order[j, 1],
                                 panel = order[j, 2],
                                 index = index,
                                 form = x$n_modules * (order[j, 2] - 1) + index)
            info[selected_item$id, ] <- NA
            temp <- temp + 1
        }
    }

    x$items <- x$items %>% arrange(form)
    x
}

sim <- function(x, routing_method) {
    # creating a matrix to list item membership in each module
    modules <- matrix(0, nrow(x$items) / x$n_panels, x$n_modules)
    for (i in 1:x$n_modules) {
        len1 <- sum(subset(x$module, subset = index < i)$len)
        len2 <- subset(x$module, subset = index == i)$len
        modules[(1 + len1):(len1 + len2), i] <- 1
    }

    # creating the transition matrix to the MST design
    trans <- matrix(0, x$n_modules, x$n_modules)
    if (paste0(x$design, collapse = "-") == "1-3-3"){
    trans[1, 2:4] <- trans[2, 5:6] <- trans[3, 5:7] <- trans[4, 6:7] <- 1
    }else{
    trans[1, 2:4] <- trans[2, 5:6] <- trans[3, 6:7] <- trans[4, 7:8] <- 1}


    x$mst_admin <- data.frame(matrix(ncol = 8, #<- this was 7, adjusted to 8 because of runtime
                                     nrow = length(THETA) * N_REP,
                                     dimnames = list(c(),
                                                     c("true_theta", "rep",
                                                       "item_no", "item_id",
                                                       "interim_theta", "interim_sem",
                                                       "panel", "path"))))

    x$cat_admin <- data.frame(matrix(ncol = 6,
                                     nrow = length(THETA) * N_REP * TEST_LEN,
                                     dimnames = list(c(),
                                                     c("true_theta", "rep",
                                                       "item_no", "item_id",
                                                       "interim_theta", "interim_sem"))))

    all_items <- x$items %>% select(a, b, c, d)
    
    #create population level responses if routing by NC
    if(routing_method=="NC")
    {popresp <- sapply(rnorm(1e3,0,1), genPattern, all_items, D = 1.7)}
    
    

    pb <- progress_bar$new(total = length(THETA) * N_REP)
    temp <- 1
    for (i in 1:length(THETA)) {
        for (j in 1:N_REP) {
            # generate responses for MST and CAT
            resp <- genPattern(THETA[i], all_items, D = 1.7)

            # select a panel for MST randomly
            rand <- sample(1:x$n_panels, 1)
            
            
            moduleInfo=apply(modules, 2, function(c){sum(c*info(THETA[i],x$items[which(x$items$panel == rand),]))} )
            if(paste0(x$design, collapse = "-") == "1-3-3"){paths=expand.grid(c(1),c(2:4),c(5:7))}
            if(paste0(x$design, collapse = "-") == "1-3-4"){paths=expand.grid(c(1),c(2:4),c(5:7))}
            #expected_path=paste(paths[which.max(apply(paths,1,function(c){sum(moduleInfo[as.numeric(c)])}))],collapse="-")
            #remove if unused
            
            if(routing_method=="NC"){
              #create population level responses to the specific panel
              poprespsubset=popresp[which(x$items$panel == rand),]
              
              #scores at each module
              moduleNC=apply(modules, 2, function(c){sum(c*resp[which(x$items$panel == rand)])} )
              popmoduleNC=apply(modules, 2, function(c){colSums(c*poprespsubset)} )
            
              #paths - 
              first=1
              
              middle=2+sum(moduleNC[1]>=quantile(popmoduleNC[,first], probs=seq(0,1,length.out=4))[2:3])
              #print(popmoduleNC[,first])
              if(paste0(x$design, collapse = "-")=="1-3-3"){
                if(middle==2){
                    final=5+sum(moduleNC[middle]>=quantile(popmoduleNC[,2], probs=seq(0,1,length.out=4))[3])}
                if(middle==3){
                    final=5+sum(moduleNC[middle]>=quantile(popmoduleNC[,3], probs=seq(0,1,length.out=4))[2:3])}
                if(middle==4){
                  final=6+sum(moduleNC[middle]>=quantile(popmoduleNC[,4], probs=seq(0,1,length.out=4))[2])}}
              
              if(paste0(x$design, collapse = "-")=="1-3-4"){
                if(middle==2){
                  final=5+sum(moduleNC[middle]>=quantile(popmoduleNC[,2], probs=seq(0,1,length.out=5))[4])}
                if(middle==3){
                  final=6+sum(moduleNC[middle]>=quantile(popmoduleNC[,3], probs=seq(0,1,length.out=5))[3])}
                if(middle==4){
                  final=7+sum(moduleNC[middle]>=quantile(popmoduleNC[,4], probs=seq(0,1,length.out=5))[2])}}        

              NC=sum(moduleNC[c(first,middle,final)])
              
              items=unlist(sapply(c(first,middle,final), function(c){which(modules[,c]==1)}))
              
              interim_theta=sapply(1:length(items),function(i){
                thetaEst(all_items[1:i,], resp[1:i], D=1.7, method="ML", range=c(-3.5,3.5))})
              
              interim_sem=sapply(1:length(items),function(i){
                semTheta(interim_theta[1:i], all_items[1:i,], D=1.7, method="ML")})

              # mst=data.frame(
              #   true_theta = rep(THETA[i], times=TEST_LEN),
              #   rep = rep(j, times=TEST_LEN),
              #   item_no = 1:TEST_LEN,
              #   item_id = items,
              #   interim_theta = interim_theta,
              #   interim_sem = interim_sem,
              #   panel = rep(rand, times=TEST_LEN),
              #   path = as.numeric(paste0(c(first,middle,final), collapse = ""))
              #   )
              
              mst=list(
                "trueTheta" = THETA[i],
                "testItems" = items,
                "allTheta" = cbind("th"=interim_theta,"se"=interim_sem),
                "selected.modules" = as.numeric(paste0(c(first,middle,final), collapse = "")))

            }else{
            #theta cutoff mat
            cutoffmat = switch(paste0(x$design, collapse = "-"),
                   "1-3-3"= matrix(ncol=3, byrow=TRUE, c(
                           c(2,3,-0.43),
                           c(3,4,0.43),
                           c(5,6,-0.43),
                           c(6,7,0.43))),

                   "1-3-4"= matrix(ncol=3, byrow=TRUE, c(
                           c(2,3,-0.43),
                           c(3,4,0.43),
                           c(5,6,-0.67),
                           c(6,7,0),
                           c(7,8,0.67)))
                     )

            testopts=switch(routing_method, 
                             "information"=list(method = "ML", D = 1.7,
                                                range = c(-3.5, 3.5),
                                                moduleSelect = "MFI"),
                             "theta"= list(method = "ML", D = 1.7,
                                           range = c(-3.5, 3.5),
                                           cutoff = cutoffmat))
  

            mst <- randomMST(trueTheta = THETA[i],
                             itemBank = subset(x$items,
                                               subset = panel == rand,
                                               select = c(a, b, c, d)),
                             modules = modules,
                             transMatrix = trans,
                             responses = resp[which(x$items$panel == rand)],
                             start = list(fixModule = 1, D = 1.7),
                             test = testopts,
                             final = list(method = "ML", D = 1.7,
                                          range = c(-3.5, 3.5)),
                             allTheta = TRUE)}
            
            # save output
            range <- (TEST_LEN * (temp - 1) + 1):(TEST_LEN * temp)
            x$mst_admin[range,] <- cbind(true_theta = mst$trueTheta,
                                         rep = j,
                                         item_no = 1:TEST_LEN,
                                         item_id = mst$testItems,
                                         interim_theta = mst$allTheta[, "th"],
                                         interim_sem = mst$allTheta[, "se"],
                                         panel = rand,
                                         path = as.numeric(paste0(mst$selected.modules, collapse = "")))#,
                                         #expectedpath=expected_path)
            x

            # simulate CAT
            cat <- randomCAT(trueTheta = THETA[i],
                             itemBank = all_items,
                             responses = resp,
                             start = list(nrItems = 1, D = 1.7,
                                          theta = runif(1, -0.5, 0.5),
                                          startSelect = "MFI"),
                             test = list(method = "ML", D = 1.7,
                                         range = c(-3.5, 3.5),
                                         itemSelect = "MFI"),
                             stop = list(rule = "length",
                                         thr = TEST_LEN),
                             final = list(method = "ML", D = 1.7,
                                          range = c(-3.5, 3.5)),
                             allTheta = TRUE)


            x$cat_admin[range,] <- cbind(true_theta = cat$trueTheta,
                                         rep = j,
                                         item_no = 1:TEST_LEN,
                                         item_id = cat$testItems,
                                         interim_theta = cat$thetaProv,
                                         interim_sem = cat$seProv)
            temp <- temp + 1
            pb$tick()
        }
    }
    x$mst_admin$path <- sapply(x$mst_admin$path,
                               function(p)
                                   p %>% as.character() %>%
                                   strsplit(split = "") %>%
                                   unlist() %>% as.integer() %>%
                                   paste0(collapse = "-"))
      
    # if(routing_method=="NC"){
    #   x$moduleNC=moduleNC
    #   x$popmoduleNC=popmoduleNC
    #   x$modules=modules
    #   x$resp=resp
    #   }
    
    x
}

stat <- function(x) {
    # for MST
    mst_final_est <- subset(x$mst_admin, item_no == TEST_LEN)

    # routing error
    path_freq <- as.data.frame.matrix(table(mst_final_est$true_theta, mst_final_est$path))

    expected_path_lookup <- data.frame(true_theta = THETA,
                                       expected_path = names(path_freq)[apply(path_freq, 1, which.max)])
    
    
    
    mst_final_est <- mst_final_est %>% left_join(expected_path_lookup, by = "true_theta")
    mst_final_est$n_routing_errors <- apply(mst_final_est, 1,
                                        function(p) {
                                            observed <- strsplit(p[7], split = "[-]")[[1]]
                                            expected <- strsplit(p[8], split = "[-]")[[1]]
                                            n_routing_errors <- sum(observed != expected)
                                            return(n_routing_errors)
                                        })

    x$routing_error_percent <- mst_final_est %>%
        group_by(true_theta) %>%
        summarise("1" = sum(n_routing_errors == 1) / N_REP * 100,
                  "2" = sum(n_routing_errors == 2) / N_REP * 100) %>%
        round(., 4) %>%
        gather(key = n_routing_errors, value = routing_error_percent, 2:3)

    x$routing_error_percent$design <- paste("MST", paste(x$design, collapse = "-"))
    x$routing_error_percent$items_per_stage <- switch(x$items_per_stage,
                                   "increasing" = "(7, 14, 21)",
                                   "equal" = "(14, 14, 14)",
                                   "decreasing" = "(21, 14, 7)")
    x$routing_error_percent$assembly_priority <- x$assembly_priority
    x$routing_error_percent$design <- factor(x$routing_error_percent$design)
    x$routing_error_percent$items_per_stage <- factor(x$routing_error_percent$items_per_stage,
                                   levels = c("(7, 14, 21)",
                                              "(14, 14, 14)",
                                              "(21, 14, 7)"))
    x$routing_error_percent$assembly_priority <- factor(x$routing_error_percent$assembly_priority)


    # mean bias, RMSE, mean SEM
    routing_error <- function(mst_final_est, n = NULL) {
        if(is.null(n)) {
            stat <- mst_final_est %>%
                group_by(true_theta) %>%
                summarise(mean_bias = mean(true_theta - interim_theta),
                          rmse = sqrt(mean((true_theta - interim_theta) ^ 2)),
                          mean_sem = mean(interim_sem)) %>%
                round(., 4)
        } else {
            stat <- mst_final_est %>%
                filter(n_routing_errors == n) %>%
                group_by(true_theta) %>%
                summarise(mean_bias = mean(true_theta - interim_theta),
                          rmse = sqrt(mean((true_theta - interim_theta) ^ 2)),
                          mean_sem = mean(interim_sem)) %>%
                mutate(n_routing_errors = n) %>%
                round(., 4)
            stat$n_routing_errors <- factor(stat$n_routing_errors)
        }
        stat$design <- paste("MST", paste(x$design, collapse = "-"))
        stat$items_per_stage <- switch(x$items_per_stage,
                                       "increasing" = "(7, 14, 21)",
                                       "equal" = "(14, 14, 14)",
                                       "decreasing" = "(21, 14, 7)")
        stat$assembly_priority <- x$assembly_priority

        stat$design <- factor(stat$design)
        stat$items_per_stage <- factor(stat$items_per_stage,
                                       levels = c("(7, 14, 21)",
                                                  "(14, 14, 14)",
                                                  "(21, 14, 7)"))
        stat$assembly_priority <- factor(stat$assembly_priority)
        return(stat)
    }
    x$expectedpath <- expected_path_lookup
    x$mst_stat <- routing_error(mst_final_est)
    x$mst_stat0 <- routing_error(mst_final_est, n = 0)
    x$mst_stat1 <- routing_error(mst_final_est, n = 1)
    x$mst_stat2 <- routing_error(mst_final_est, n = 2)

    # correlation
    x$mst_cor <- round(cor(mst_final_est$true_theta, mst_final_est$interim_theta), 4)


    # for CAT

    # no. of items to match MST's SEM
    routing_error <- function(n = NULL) {
        if (is.null(n)) {
            match_mst <- x$cat_admin %>%
                mutate(mst_sem = rep(mst_final_est$interim_sem, each = TEST_LEN),
                       diff = as.numeric(interim_sem > mst_sem)) %>%
                group_by(true_theta, rep) %>%
                mutate(match_mst = max(max(which(diff == 1)) + 1, TEST_LEN)) %>%
                group_by(true_theta) %>%
                summarize(match_mst = mean(match_mst))
        } else {
            temp <- x$mst_admin %>%
                inner_join(mst_final_est %>% select(true_theta, rep, n_routing_errors),
                           by = c("true_theta", "rep"))
            indices <- which(temp$n_routing_errors == n)
            match_mst <- x$cat_admin[indices, ] %>%
                mutate(mst_sem = rep((mst_final_est %>% filter(n_routing_errors == n))$interim_sem, each = TEST_LEN),
                       diff = as.numeric(interim_sem > mst_sem)) %>%
                group_by(true_theta, rep) %>%
                mutate(match_mst = max(which(diff == 1)) + 1,
                       match_mst = ifelse(match_mst == 43, 42, match_mst)) %>%
                group_by(true_theta) %>%
                summarize(match_mst = mean(match_mst)) %>%
                mutate(n_routing_errors = n)
            match_mst$n_routing_errors <- factor(match_mst$n_routing_errors)
        }
        match_mst$design <- paste("CAT", paste(x$design, collapse = "-"))
        match_mst$items_per_stage <- switch(x$items_per_stage,
                                             "increasing" = "(7, 14, 21)",
                                             "equal"="(14, 14, 14)",
                                             "decreasing" = "(21, 14, 7)")
        match_mst$assembly_priority <- x$assembly_priority
        match_mst$design <- factor(match_mst$design)
        match_mst$items_per_stage <- factor(match_mst$items_per_stage,
                                             levels = c("(7, 14, 21)",
                                                        "(14, 14, 14)",
                                                        "(21, 14, 7)"))
        match_mst$assembly_priority <- factor(match_mst$assembly_priority)
        return(match_mst)
    }
    match_mst <- routing_error()
    x$match_mst0 <- routing_error(n = 0)
    x$match_mst1 <- routing_error(n = 1)
    x$match_mst2 <- routing_error(n = 2)

    # mean bias, RMSE, mean SEM
    cat_final_est <- subset(x$cat_admin, item_no == TEST_LEN)

    x$cat_stat <- cat_final_est %>%
        group_by(true_theta) %>%
        summarise(mean_bias = mean(true_theta - interim_theta),
                  rmse = sqrt(mean((true_theta - interim_theta) ^ 2)),
                  mean_sem = mean(interim_sem)) %>%
        round(., 4) %>%
        inner_join(match_mst, by = "true_theta")

    x$cat_cor <- round(cor(cat_final_est$true_theta, cat_final_est$interim_theta), 4)

    x
}

compare <- function(design, items_per_stage, assembly_priority, routing_method) {
    x <- mst_setup(design,
                   items_per_stage = items_per_stage,
                   assembly_priority = assembly_priority)
    x <- mst_assemble(x)
    x <- sim(x, routing_method)
    x <- stat(x)
    x
}

compare2 <- function(design, items_per_stage, assembly_priority, routing_method) {
  x <- mst_setup(design,
                 items_per_stage = items_per_stage,
                 assembly_priority = assembly_priority)
  x <- mst_assemble(x)
  x <- sim(x, routing_method)
  #x <- stat(x)
   x
}




plot_module_info <- function(...) {
    mst <- list(...)
    theta <- seq(-3, 3, 0.1)
    data <- NULL
    for(i in 1:length(mst)) {
        for (j in unique(mst[[i]]$items$index)) {
            items <- subset(mst[[i]]$items, mst[[i]]$items$index == j)
            data <- rbind(data, data.frame(t = theta,
                                           info = colSums(info(theta = theta, items = items)) / mst[[i]]$n_panels,
                                           index = j,
                                           stage = items$stage[1],
                                           assembly_priority = mst[[i]]$assembly_priority))
        }
    }

    data <- data %>%
        group_by(t, index, stage, assembly_priority) %>%
        summarise(info = mean(info))

    data$assembly_priority <- factor(data$assembly_priority)
    data$index <- factor(data$index)
    data$stage <- factor(paste("Stage", data$stage))

    g <- ggplot(data, aes_string(x = "t", y = "info", color = "index")) +
        geom_line() + xlab(expression(theta)) + ylab("Information") +
        theme_bw() + theme(legend.position = "none") +
        guides(color = guide_legend("Modules")) +
        facet_grid(assembly_priority ~ stage)
    g
}


mst_plot <- function(x) {
    theta <- seq(-3, 3, 0.1)
    data <- NULL
    for (i in unique(x$items$form)) {
        items <- subset(x$items, x$items$form == i)
        data <- rbind(data, data.frame(t = theta,
                                       info = colSums(info(theta = theta, items = items)),
                                       form = i,
                                       panel = items$panel[1],
                                       stage = items$stage[1]))
    }
    data$form <- factor(data$form)
    data$panel <- factor(paste("Panel", data$panel))
    data$stage <- factor(paste("Stage", data$stage))
    g <- ggplot(data, aes_string(x = "t", y = "info", color = "form")) +
        geom_line() + xlab(expression(theta)) + ylab("Information") +
        theme_bw() + theme(legend.position = "none") +
        guides(color = guide_legend("Modules")) + facet_grid(panel ~ stage)
    g
}

# Simulation

Design=c("1-3-3","1-3-4")
NI=c("increasing","equal","decreasing")
Assembly=c("Forward","Backward","Spiral","Random")
Routing=c("information","theta","NC")


Design_Conditions=expand.grid("Design"=Design,"Number of Items"=NI, "Assembly Priority"=Assembly, "Routing"=Routing, stringsAsFactors = FALSE)
# Design_Conditions=Design_Conditions[1:2,]

Out=mapply(compare, SIMPLIFY=FALSE,
           Design_Conditions$Design, 
           Design_Conditions$'Number of Items', 
           Design_Conditions$'Assembly Priority',
           Design_Conditions$Routing)
names(Out)=apply(Design_Conditions,1,paste, names(Design_Conditions), collapse=" | ")

save.image(paste0(Sys.Date(), ".RData"))

View(Out[[1]]$mst_admin)

#########################################################################
### Table X - Average Path Error Percentage #############################
#########################################################################

AveragePathErrorByDesign=apply(Design_Conditions,1, paste, collapse="")

AveragePathErrorByDesign=Out

AveragePathErrorByDesignCondition=sapply(Routing, simplify=FALSE, function(routing){
  routingindex=grep(routing,names(AveragePathErrorByDesign))  
  tableout=do.call(cbind,
    lapply(Assembly, function(assembly){
      assemblyindex=grep(assembly, names(AveragePathErrorByDesign))
      routingassemblyindex=assemblyindex[assemblyindex %in% routingindex]
        sapply(Design, function(design){
          designindex=grep(design,names(AveragePathErrorByDesign))
          routingassemblydesignindex=designindex[designindex %in% routingassemblyindex]
          out=sapply(Out[routingassemblydesignindex], function(dat){
            dat=merge(
              dat$mst_admin[which(dat$mst_admin$item_no==42),],
              dat$expectedpath, 
              by="true_theta", all.x = TRUE)
            mean(dat$expected_path!=dat$path)
          })
          names(out)=gsub('^.*Design \\| \\s*|\\s*Number of Items.*$', '', names(Out)[routingassemblydesignindex])
          out
        })
    }))
  dimnames(tableout)[2][[1]]=c(outer(Design,Assembly, paste))
  tableout
})

for(routing in names(AveragePathErrorByDesignCondition)){
  write.csv(AveragePathErrorByDesignCondition[[routing]], paste(routing,"AveragePathErrorByDesignCondition.csv",sep="_"))
}


#% of Path Error for 1-3-3 and 1-3-4
#to be calculated




############################################################################
### Table X - Average Estimation Error #####################################
############################################################################


AveragePathErrorByDesign=apply(Design_Conditions,1, paste, collapse="")

AveragePathErrorByDesign=Out


AverageBiasByDesignCondition=sapply(Routing, simplify=FALSE, function(routing){
  routingindex=grep(routing,names(AveragePathErrorByDesign))  
  tableout=do.call(cbind,
                   lapply(Assembly, function(assembly){
                     assemblyindex=grep(assembly, names(AveragePathErrorByDesign))
                     routingassemblyindex=assemblyindex[assemblyindex %in% routingindex]
                     sapply(Design, function(design){
                       designindex=grep(design,names(AveragePathErrorByDesign))
                       routingassemblydesignindex=designindex[designindex %in% routingassemblyindex]
                       out=sapply(Out[routingassemblydesignindex], function(dat){
                         dat=merge(
                           dat$mst_admin[which(dat$mst_admin$item_no==42),],
                           dat$expectedpath, 
                           by="true_theta", all.x = TRUE)
                         
                         mean(dat$true_theta-dat$interim_theta)
                       })
                       names(out)=gsub('^.*Design \\| \\s*|\\s*Number of Items.*$', '', names(Out)[routingassemblydesignindex])
                       out
                     })
                   }))
  dimnames(tableout)[2][[1]]=c(outer(Design,Assembly, paste))
  tableout
})

for(routing in names(AverageBiasByDesignCondition)){
  write.csv(AverageBiasByDesignCondition[[routing]], paste(routing,"AverageBiasByDesignCondition.csv",sep="_"))
}

############################################################################
### Table X - RMSE #####################################
############################################################################


AveragePathErrorByDesign=apply(Design_Conditions,1, paste, collapse="")

AveragePathErrorByDesign=Out


RMSEByDesignCondition=sapply(Routing, simplify=FALSE, function(routing){
  routingindex=grep(routing,names(AveragePathErrorByDesign))  
  tableout=do.call(cbind,
                   lapply(Assembly, function(assembly){
                     assemblyindex=grep(assembly, names(AveragePathErrorByDesign))
                     routingassemblyindex=assemblyindex[assemblyindex %in% routingindex]
                     sapply(Design, function(design){
                       designindex=grep(design,names(AveragePathErrorByDesign))
                       routingassemblydesignindex=designindex[designindex %in% routingassemblyindex]
                       out=sapply(Out[routingassemblydesignindex], function(dat){
                         dat=merge(
                           dat$mst_admin[which(dat$mst_admin$item_no==42),],
                           dat$expectedpath, 
                           by="true_theta", all.x = TRUE)
                         

                         sqrt(sum((dat$true_theta-dat$interim_theta)^2)/length(dat$true_theta))
                       })
                       names(out)=gsub('^.*Design \\| \\s*|\\s*Number of Items.*$', '', names(Out)[routingassemblydesignindex])
                       out
                     })
                   }))
  dimnames(tableout)[2][[1]]=c(outer(Design,Assembly, paste))
  tableout
})

for(routing in names(RMSEByDesignCondition)){
  write.csv(RMSEByDesignCondition[[routing]], paste(routing,"RMSEByDesignCondition.csv",sep="_"))
}

###########################################################################################
###### Average % of Errors by Theta #######################################################
###########################################################################################

OverallPathErrors=unlist(lapply(Out, function(x){
  x=merge(x$mst_admin, x$expectedpath, by="true_theta")
  nerrors=apply(x,1,function(y){
    sum(!(strsplit(y["path"],"-")[[1]] %in% strsplit(y["expected_path"],"-")[[1]]))})}))

OverallPathErrorsTable=
  cbind(
    "n"=table(OverallPathErrors),
    "Percent"=round(prop.table(table(OverallPathErrors)),3))

write.csv(OverallPathErrorsTable,paste0(getwd(),"/Results/Tables/OverallPathErrorsTable.csv"))


###########################################################################################
###### Average % of Errors by Theta #######################################################
###########################################################################################


for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    png(paste(routing,"routing",design,"PercErrorsbyTrueThetaPlots.png", sep="_"), 1200,800)
    par(mfcol=c(length(Assembly),length(NI)),
        oma = c(3.5,3.5,3.5,3.5) + 0.1,
        mar = c(1.25,1.25,1.25,1.25) + 0.1)

    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(as.character(c(0:2)), function(n_errors){
          sapply(unique(temp$true_theta), function(t){
            thetaindx=(temp$true_theta==t)
            errorindx=(temp$n_errors==as.numeric(n_errors))
            if(sum(thetaindx & errorindx)>0){length(temp$interim_theta[thetaindx & errorindx])/length(temp$interim_theta[thetaindx])}
            else{NA}
          })
        })
        
        plot(NA, xlim=c(-3,3), ylim=c(0,1),type="l", yaxt="n", xlab="True Theta", ylab="Error Percentage")#, main=paste(assembly,ni))
        axis(2, at=seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        lines(seq(-3,3,0.5), temp2[,"1"], lty=2)
        lines(seq(-3,3,0.5), temp2[,"2"], lty=3)
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    title(paste(routing,"Routing",design,"Error Percentage by True Theta"),
          xlab = "Theta",
          ylab = "Error Percentage",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c( "One Error", "Two Errors"), lty=c(2:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    
    dev.off()
  }
}


###########################################################################################
###### RMSE Conditional Bias by # of Errors ###############################################
###########################################################################################

  #All Errors one graph per design (1-3-3 or 1-3-4)

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
  
    png(paste(routing,"routing",design,"AllErrorsRMSEbyTrueThetaPlots.png", sep="_"), 1200,800)
      
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(unique(temp$true_theta), function(t){
                thetaindx=(temp$true_theta==t)
                if(sum(thetaindx)>0){sqrt(sum((t-temp$interim_theta[thetaindx])^2)/sum(thetaindx))}
                else{NA}
              })
            
        
        plot(seq(-3,3,0.5), temp2, ylim=c(0,1.5),type="l", , yaxt="n", xlab="True Theta", ylab="RMSE")#, main=paste(assembly,ni))
        axis(2, at=seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    title(main= paste(routing,"Routing",design,"All Errors RMSE by True Theta"),
          xlab = "Theta",
          ylab = "RMSE",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    
    dev.off()
  }
}

#All Errors both designs (1-3-3 AND 1-3-4) in one graph 

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 

    png(paste(routing,"routing","AllErrorsBothDesignsRMSEbyTrueThetaPlots.png", sep="_"), 1200,800)
    
    Routing133Names=grep("1-3-3", RoutingNames, value=TRUE)
    Routing134Names=grep("1-3-4", RoutingNames, value=TRUE)
    
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      Routing133NINames=grep(ni, Routing133Names, value=TRUE)
      Routing134NINames=grep(ni, Routing134Names, value=TRUE)
      for(assembly in Assembly){
        Routing133NIAssemblyName=grep(assembly, Routing133NINames, value=TRUE)
        Routing134NIAssemblyName=grep(assembly, Routing134NINames, value=TRUE)
        dat133=merge(Out[[Routing133NIAssemblyName]]$mst_admin[which(Out[[Routing133NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing133NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        dat134=merge(Out[[Routing134NIAssemblyName]]$mst_admin[which(Out[[Routing134NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing134NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        dat133$n_errors=apply(dat133,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        dat134$n_errors=apply(dat134,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        dat133=sapply(unique(dat133$true_theta), function(t){
          thetaindx=(dat133$true_theta==t)
          if(sum(thetaindx)>0){sqrt(sum((t-dat133$interim_theta[thetaindx])^2)/sum(thetaindx))}
          else{NA}
        })
        dat134=sapply(unique(dat134$true_theta), function(t){
          thetaindx=(dat134$true_theta==t)
          if(sum(thetaindx)>0){sqrt(sum((t-dat134$interim_theta[thetaindx])^2)/sum(thetaindx))}
          else{NA}
        })
        
        
        plot(seq(-3,3,0.5), dat133, ylim=c(0,1.5),type="l", , yaxt="n", xlab="True Theta", ylab="RMSE")#, main=paste(assembly,ni))
        axis(2, at=seq(0,1.5,length.out=7), labels = seq(0,1.5,length.out=7))
        lines(seq(-3,3,0.5), dat134, lty=2)
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}

    title(main = paste(routing,"Routing",design,"All Errors Both Designs RMSE by True Theta"),
          xlab = "Theta",
          ylab = "RMSE",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    
    dev.off()
  
}



    
  #Zero Routing Errors one graph per design (1-3-3 or 1-3-4)

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    
    png(paste(routing,"routing",design,"ZeroErrorsRMSEbyTrueThetaPlots.png", sep="_"), 1200,800)
    
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(unique(temp$true_theta), function(t){
          thetaindx=(temp$true_theta==t)
          errorindx=(temp$n_errors==0)
          if(sum(thetaindx & errorindx)>0){sqrt(sum((t-temp$interim_theta[thetaindx & errorindx])^2)/sum(thetaindx & errorindx))}
          else{NA}
        })
        
        
        plot(seq(-3,3,0.5), temp2, ylim=c(0,1.5),type="l", , yaxt="n", xlab="True Theta", ylab="RMSE")#, main=paste(assembly,ni))
        axis(2, at=seq(0,1,length.out=5), labels = seq(0,1,length.out=5))
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    
    title(main = paste(routing,"Routing",design,"Zero Errors RMSE by True Theta"),
          xlab = "Theta",
          ylab = "RMSE",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    
    dev.off()
  }
}

#Zero Routing Errors both designs (1-3-3 AND 1-3-4) in one graph 

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  png(paste(routing,"routing","ZeroErrorsBothDesignsRMSEbyTrueThetaPlots.png", sep="_"), 1200,800)
  
  Routing133Names=grep("1-3-3", RoutingNames, value=TRUE)
  Routing134Names=grep("1-3-4", RoutingNames, value=TRUE)
  
  par(mfcol=c(length(Assembly),length(NI)))
  par(cex.lab=1.25)
  par(cex.axis=1.25)
  for(ni in NI){
    Routing133NINames=grep(ni, Routing133Names, value=TRUE)
    Routing134NINames=grep(ni, Routing134Names, value=TRUE)
    for(assembly in Assembly){
      Routing133NIAssemblyName=grep(assembly, Routing133NINames, value=TRUE)
      Routing134NIAssemblyName=grep(assembly, Routing134NINames, value=TRUE)
      dat133=merge(Out[[Routing133NIAssemblyName]]$mst_admin[which(Out[[Routing133NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing133NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat134=merge(Out[[Routing134NIAssemblyName]]$mst_admin[which(Out[[Routing134NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing134NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat133$n_errors=apply(dat133,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat134$n_errors=apply(dat134,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat133=sapply(unique(dat133$true_theta), function(t){
        thetaindx=(dat133$true_theta==t)
        errorindx=(dat133$n_errors==0)
        if(sum(thetaindx & errorindx)>0){sqrt(sum((t-dat133$interim_theta[thetaindx & errorindx])^2)/sum(thetaindx & errorindx))}
        else{NA}
      })
      dat134=sapply(unique(dat134$true_theta), function(t){
        thetaindx=(dat134$true_theta==t)
        errorindx=(dat134$n_errors==0)
        if(sum(thetaindx & errorindx)>0){sqrt(sum((t-dat134$interim_theta[thetaindx & errorindx])^2)/sum(thetaindx & errorindx))}
        else{NA}
      })
      
      plot(seq(-3,3,0.5), dat133, ylim=c(0,1.5),type="l", , yaxt="n", xlab="True Theta", ylab="RMSE")#, main=paste(assembly,ni))
      axis(2, at=seq(0,1.5,length.out=7), labels = seq(0,1.5,length.out=7))
      lines(seq(-3,3,0.5), dat134, lty=2)
      if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
      if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
      
    }
  }
  if(routing=="NC"){routing="PI-NC"}
  if(routing=="theta"){routing="PI-\u0398"}
  if(routing=="information"){routing="MI"}
  title(main = paste(routing,"Routing",design, "Zero Errors Both Designs RMSE by True Theta"),
        xlab = "Theta",
        ylab = "RMSE",
        outer = TRUE, line = 1, cex.lab=3)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
  
  dev.off()
  
}

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    png(paste0(getwd(),"/Results/Figures/",paste(routing,"routing",design,"RMSEbyTrueThetaPlots.png", sep="_")), 1200,800)
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(as.character(c(0:2)), function(n_errors){
          sapply(unique(temp$true_theta), function(t){
            thetaindx=(temp$true_theta==t)
            errorindx=(temp$n_errors==as.numeric(n_errors))
            if(sum(thetaindx & errorindx)>0){sqrt(sum((t-temp$interim_theta[thetaindx & errorindx])^2)/sum(thetaindx & errorindx))}
            else{NA}
          })
        })
        
        plot(seq(-3,3,0.5), temp2[,"0"], ylim=c(0,1.5),type="l", yaxt="n", xlab="True Theta", ylab="RMSE")#, main=paste(assembly,ni))
        axis(2, at=seq(0,1.5,length.out=7), labels = seq(0,1.5,length.out=7))
        lines(seq(-3,3,0.5), temp2[,"1"], lty=2)
        lines(seq(-3,3,0.5), temp2[,"2"], lty=3)
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    title(main = paste(routing,"Routing",design, "RMSE by True Theta"),
          xlab = "Theta",
          ylab = "RMSE",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)    
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')    
    
    dev.off()
  }
}

###########################################################################################
###### Mean Conditional Bias by # of Errors ###############################################
###########################################################################################

#All Errors one graph per design (1-3-3 or 1-3-4)

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    
    png(paste(routing,"routing",design,"AllErrorsMeanBiasbyTrueThetaPlots.png", sep="_"), 1200,800)
    
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(unique(temp$true_theta), function(t){
          thetaindx=(temp$true_theta==t)
          if(sum(thetaindx)>0){mean(t-temp$interim_theta[thetaindx])}
          else{NA}
        })
        
        
        plot(seq(-3,3,0.5), temp2, ylim=c(-1,1),type="l", yaxt="n", xlab="True Theta", ylab="Theta Bias")#, main=paste(assembly,ni))
        axis(2, at=seq(-1,1,length.out=5), labels = seq(-1,1,length.out=5))
        abline(h=0, lty=3)
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    title(main = paste(routing,"routing",design, "All Errors Mean Bias by True Theta"),
          xlab = "Theta",
          ylab = "Theta Bias",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)
    
    # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    # plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    # legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')    
    # 
    dev.off()
  }
}

#All Errors both designs (1-3-3 AND 1-3-4) in one graph 

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  png(paste(routing,"routing","AllErrorsBothDesignsMeanBiasbyTrueThetaPlots.png", sep="_"), 1200,800)
  
  Routing133Names=grep("1-3-3", RoutingNames, value=TRUE)
  Routing134Names=grep("1-3-4", RoutingNames, value=TRUE)
  
  par(mfcol=c(length(Assembly),length(NI)))
  par(cex.lab=1.25)
  par(cex.axis=1.25)
  for(ni in NI){
    Routing133NINames=grep(ni, Routing133Names, value=TRUE)
    Routing134NINames=grep(ni, Routing134Names, value=TRUE)
    for(assembly in Assembly){
      Routing133NIAssemblyName=grep(assembly, Routing133NINames, value=TRUE)
      Routing134NIAssemblyName=grep(assembly, Routing134NINames, value=TRUE)
      dat133=merge(Out[[Routing133NIAssemblyName]]$mst_admin[which(Out[[Routing133NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing133NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat134=merge(Out[[Routing134NIAssemblyName]]$mst_admin[which(Out[[Routing134NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing134NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat133$n_errors=apply(dat133,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat134$n_errors=apply(dat134,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat133=sapply(unique(dat133$true_theta), function(t){
        thetaindx=(dat133$true_theta==t)
        if(sum(thetaindx)>0){mean(t-dat133$interim_theta[thetaindx])}
        else{NA}
      })
      dat134=sapply(unique(dat134$true_theta), function(t){
        thetaindx=(dat134$true_theta==t)
        if(sum(thetaindx)>0){mean(t-dat134$interim_theta[thetaindx])}
        else{NA}
      })
      
      
      plot(seq(-3,3,0.5), dat133, ylim=c(-1,1),type="l", , yaxt="n", xlab="True Theta", ylab="Theta Bias")#, main=paste(assembly,ni))
      axis(2, at=seq(-1,1,length.out=5), labels = seq(-1,1,length.out=5))
      abline(h=0, lty=3)
      lines(seq(-3,3,0.5), dat134, lty=2)
      if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
      if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
      
    }
  }
  if(routing=="NC"){routing="PI-NC"}
  if(routing=="theta"){routing="PI-\u0398"}
  if(routing=="information"){routing="MI"}  
  
  title(main = paste(routing,"Routing",design, "All Errors Both Designs Mean Bias by True Theta"),
        xlab = "Theta",
        ylab = "Mean Bias",
        outer = TRUE, line = 1, cex.lab=3)  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')  
  
  dev.off()
  
}


#Zero Routing Errors one graph per design (1-3-3 or 1-3-4)

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    
    png(paste(routing,"routing",design,"ZeroErrorsMeanBiasbyTrueThetaPlots.png", sep="_"), 1200,800)
    
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(unique(temp$true_theta), function(t){
          thetaindx=(temp$true_theta==t)
          errorindx=(temp$n_errors==0)
          if(sum(thetaindx & errorindx)>0){mean(t-temp$interim_theta[thetaindx & errorindx])}
          else{NA}
        })
        
        
        plot(seq(-3,3,0.5), temp2, ylim=c(-1,1),type="l", , yaxt="n", xlab="True Theta", ylab="Theta Bias")#, main=paste(assembly,ni))
        axis(2, at=seq(-1,1,length.out=5), labels = seq(-1,1,length.out=5))
        abline(h=0, lty=3)
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    title(main = paste(routing,"routing",design, "Zero Errors Mean Bias by True Theta"),
          xlab = "Theta",
          ylab = "Mean Bias",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)    
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    
    dev.off()
  }
}

#Zero Routing Errors both designs (1-3-3 AND 1-3-4) in one graph 

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  png(paste(routing,"routing","ZeroErrorsBothDesignsMeanBiasbyTrueThetaPlots.png", sep="_"), 1200,800)
  
  Routing133Names=grep("1-3-3", RoutingNames, value=TRUE)
  Routing134Names=grep("1-3-4", RoutingNames, value=TRUE)
  
  par(mfcol=c(length(Assembly),length(NI)))
  par(cex.lab=1.25)
  par(cex.axis=1.25)
  for(ni in NI){
    Routing133NINames=grep(ni, Routing133Names, value=TRUE)
    Routing134NINames=grep(ni, Routing134Names, value=TRUE)
    for(assembly in Assembly){
      Routing133NIAssemblyName=grep(assembly, Routing133NINames, value=TRUE)
      Routing134NIAssemblyName=grep(assembly, Routing134NINames, value=TRUE)
      dat133=merge(Out[[Routing133NIAssemblyName]]$mst_admin[which(Out[[Routing133NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing133NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat134=merge(Out[[Routing134NIAssemblyName]]$mst_admin[which(Out[[Routing134NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing134NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat133$n_errors=apply(dat133,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat134$n_errors=apply(dat134,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat133=sapply(unique(dat133$true_theta), function(t){
        thetaindx=(dat133$true_theta==t)
        errorindx=(dat133$n_errors==0)
        if(sum(thetaindx & errorindx)>0){
          mean(t-dat133$interim_theta[thetaindx & errorindx])}
        else{NA}
      })
      dat134=sapply(unique(dat134$true_theta), function(t){
        thetaindx=(dat134$true_theta==t)
        errorindx=(dat134$n_errors==0)
        if(sum(thetaindx & errorindx)>0){
          mean(t-dat134$interim_theta[thetaindx & errorindx])}
        else{NA}
      })
      
      plot(seq(-3,3,0.5), dat133, ylim=c(-1,1),type="l", , yaxt="n", xlab="True Theta", ylab="Theta Bias")#, main=paste(assembly,ni))
      axis(2, at=seq(-1,1,length.out=5), labels = seq(-1,1,length.out=5))
      lines(seq(-3,3,0.5), dat134, lty=2)
      abline(h=0, lty=3)
      if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
      if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
      
    }
  }
  if(routing=="NC"){routing="PI-NC"}
  if(routing=="theta"){routing="PI-\u0398"}
  if(routing=="information"){routing="MI"}
  
  title(main = paste(routing,"Routing",design, "Zero Errors Both Designs Mean Bias byTrue Theta"),
        xlab = "Theta",
        ylab = "Mean Bias",
        outer = TRUE, line = 1, cex.lab=3)  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
  
  dev.off()
  
}

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    png(paste(routing,"routing",design,"MeanBiasbyTrueThetaPlots.png", sep="_"), 1200,800)
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(as.character(c(0:2)), function(n_errors){
          sapply(unique(temp$true_theta), function(t){
            thetaindx=(temp$true_theta==t)
            errorindx=(temp$n_errors==as.numeric(n_errors))
            if(sum(thetaindx & errorindx)>0){mean(t-temp$interim_theta[thetaindx & errorindx])}
            else{NA}
          })
        })
        
        plot(seq(-3,3,0.5), temp2[,"0"], ylim=c(-1,1),type="l", , yaxt="n", xlab="True Theta")#, main=paste(assembly,ni))
        axis(2, at=seq(-1,1,length.out=5), labels = seq(-1,1,length.out=5))
        abline(h=0, lty=3)
        lines(seq(-3,3,0.5), temp2[,"1"], lty=2)
        lines(seq(-3,3,0.5), temp2[,"2"], lty=3)
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}

    title(main = paste(routing,"Routing",design, "Mean Bias by True Theta"),
          xlab = "Theta",
          ylab = "Mean Bias",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)

    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')    
    
    dev.off()
  }
}


###########################################################################################
###### Mean SEM Conditional by # of Errors ################################################
###########################################################################################

#All Errors one graph per design (1-3-3 or 1-3-4)

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    
    png(paste(routing,"routing",design,"AllErrorsMeanSEMbyTrueThetaPlots.png", sep="_"), 1200,800)
    
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(unique(temp$true_theta), function(t){
          thetaindx=(temp$true_theta==t)
          if(sum(thetaindx)>0){mean(temp$interim_sem[thetaindx])}
          else{NA}
        })
        
        
        plot(seq(-3,3,0.5), temp2, ylim=c(0,5),type="l", , yaxt="n", xlab="True Theta", ylab="SEM")#, main=paste(assembly,ni))
        axis(2, at=seq(0,5,length.out=6), labels = seq(0,5,length.out=6))
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    title(main = paste(routing,"Routing",design, "All Errors Mean SEM by True Theta"),
          xlab = "Theta",
          ylab = "SEM",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)    
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')    
    
    dev.off()
  }
}

#All Errors both designs (1-3-3 AND 1-3-4) in one graph 

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  png(paste(routing,"routing","AllErrorsBothDesignsMeanSEMbyTrueThetaPlots.png", sep="_"), 1200,800)
  
  Routing133Names=grep("1-3-3", RoutingNames, value=TRUE)
  Routing134Names=grep("1-3-4", RoutingNames, value=TRUE)
  
  par(mfcol=c(length(Assembly),length(NI)))
  par(cex.lab=1.25)
  par(cex.axis=1.25)
  for(ni in NI){
    Routing133NINames=grep(ni, Routing133Names, value=TRUE)
    Routing134NINames=grep(ni, Routing134Names, value=TRUE)
    for(assembly in Assembly){
      Routing133NIAssemblyName=grep(assembly, Routing133NINames, value=TRUE)
      Routing134NIAssemblyName=grep(assembly, Routing134NINames, value=TRUE)
      dat133=merge(Out[[Routing133NIAssemblyName]]$mst_admin[which(Out[[Routing133NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing133NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat134=merge(Out[[Routing134NIAssemblyName]]$mst_admin[which(Out[[Routing134NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing134NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat133$n_errors=apply(dat133,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat134$n_errors=apply(dat134,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat133=sapply(unique(dat133$true_theta), function(t){
        thetaindx=(dat133$true_theta==t)
        if(sum(thetaindx)>0){mean(dat133$interim_sem[thetaindx])}
        else{NA}
      })
      dat134=sapply(unique(dat134$true_theta), function(t){
        thetaindx=(dat134$true_theta==t)
        if(sum(thetaindx)>0){mean(dat134$interim_sem[thetaindx])}
        else{NA}
      })
      
      
      plot(seq(-3,3,0.5), dat133, ylim=c(0,5),type="l", , yaxt="n", xlab="True Theta", ylab="SEM")#, main=paste(assembly,ni))
      axis(2, at=seq(0,5,length.out=6), labels = seq(0,5,length.out=6))
      lines(seq(-3,3,0.5), dat134, lty=2)
      if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
      if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
      
    }
  }
  if(routing=="NC"){routing="PI-NC"}
  if(routing=="theta"){routing="PI-\u0398"}
  if(routing=="information"){routing="MI"}
  title(main = paste(routing,"Routing",design, "All Errors Both Designs Mean SEM by True Theta"),
        xlab = "Theta",
        ylab = "SEM",
        outer = TRUE, line = 1, cex.lab=3)  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
  
  dev.off()
  
}


#Zero Routing Errors one graph per design (1-3-3 or 1-3-4)

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    
    png(paste(routing,"routing",design,"ZeroErrorsMeanSEMbyTrueThetaPlots.png", sep="_"), 1200,800)
    
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(unique(temp$true_theta), function(t){
          thetaindx=(temp$true_theta==t)
          errorindx=(temp$n_errors==0)
          if(sum(thetaindx & errorindx)>0){mean(temp$interim_sem[thetaindx & errorindx])}
          else{NA}
        })
        
        
        plot(seq(-3,3,0.5), temp2, ylim=c(0,5),type="l", , yaxt="n", xlab="True Theta", ylab="SEM")#, main=paste(assembly,ni))
        axis(2, at=seq(0,5,length.out=6), labels = seq(0,5,length.out=6))
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    title(main = paste(routing,"Routing",design, "Zero Errors Mean SEM by True Theta"),
          xlab = "Theta",
          ylab = "SEM",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)    
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    
    dev.off()
  }
}

#Zero Routing Errors both designs (1-3-3 AND 1-3-4) in one graph 

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  
  png(paste(routing,"routing","ZeroErrorsBothDesignsMeanSEMbyTrueThetaPlots.png", sep="_"), 1200,800)
  
  Routing133Names=grep("1-3-3", RoutingNames, value=TRUE)
  Routing134Names=grep("1-3-4", RoutingNames, value=TRUE)
  
  par(mfcol=c(length(Assembly),length(NI)))
  par(cex.lab=1.25)
  par(cex.axis=1.25)
  for(ni in NI){
    Routing133NINames=grep(ni, Routing133Names, value=TRUE)
    Routing134NINames=grep(ni, Routing134Names, value=TRUE)
    for(assembly in Assembly){
      Routing133NIAssemblyName=grep(assembly, Routing133NINames, value=TRUE)
      Routing134NIAssemblyName=grep(assembly, Routing134NINames, value=TRUE)
      dat133=merge(Out[[Routing133NIAssemblyName]]$mst_admin[which(Out[[Routing133NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing133NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat134=merge(Out[[Routing134NIAssemblyName]]$mst_admin[which(Out[[Routing134NIAssemblyName]]$mst_admin$item_no==42),], Out[[Routing134NIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
      dat133$n_errors=apply(dat133,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat134$n_errors=apply(dat134,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
      dat133=sapply(unique(dat133$true_theta), function(t){
        thetaindx=(dat133$true_theta==t)
        errorindx=(dat133$n_errors==0)
        if(sum(thetaindx & errorindx)>0){
          mean(dat133$interim_sem[thetaindx & errorindx])}
        else{NA}
      })
      dat134=sapply(unique(dat134$true_theta), function(t){
        thetaindx=(dat134$true_theta==t)
        errorindx=(dat134$n_errors==0)
        if(sum(thetaindx & errorindx)>0){
          mean(dat134$interim_sem[thetaindx & errorindx])}
        else{NA}
      })
      
      plot(seq(-3,3,0.5), dat133, ylim=c(0,5),type="l", , yaxt="n", xlab="True Theta", ylab="SEM")#, main=paste(assembly,ni))
      axis(2, at=seq(0,5,length.out=6), labels = seq(0,5,length.out=6))
      lines(seq(-3,3,0.5), dat134, lty=2)
      if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
      if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
      
    }
  }
  if(routing=="NC"){routing="PI-NC"}
  if(routing=="theta"){routing="PI-\u0398"}
  if(routing=="information"){routing="MI"}
  title(main = paste(routing,"Routing",design, "Zero Errors Both Designs Mean SEM by True Theta"),
        xlab = "Theta",
        ylab = "SEM",
        outer = TRUE, line = 1, cex.lab=3)  
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
  
  dev.off()
  
}

for(routing in Routing){
  RoutingNames=grep(routing, names(Out), value = TRUE) 
  for(design in Design){
    RoutingDesignNames=grep(design, RoutingNames, value=TRUE)
    png(paste(routing,"routing",design,"MeanSEMbyTrueThetaPlots.png", sep="_"), 1200,800)
    par(mfcol=c(length(Assembly),length(NI)),         oma = c(3.5,3.5,3.5,3.5) + 0.1,         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    for(ni in NI){
      RoutingDesignNINames=grep(ni, RoutingDesignNames, value=TRUE)
      for(assembly in Assembly){
        RoutingDesignNIAssemblyName=grep(assembly, RoutingDesignNINames, value=TRUE ) 
        temp=merge(Out[[RoutingDesignNIAssemblyName]]$mst_admin[which(Out[[RoutingDesignNIAssemblyName]]$mst_admin$item_no==42),], Out[[RoutingDesignNIAssemblyName]]$expectedpath, by="true_theta", all.x=TRUE)
        temp$n_errors=apply(temp,1,function(r){sum(!strsplit(r["path"],"-")[[1]]%in%strsplit(r["expected_path"],"-")[[1]])})
        temp2=sapply(as.character(c(0:2)), function(n_errors){
          sapply(unique(temp$true_theta), function(t){
            thetaindx=(temp$true_theta==t)
            errorindx=(temp$n_errors==as.numeric(n_errors))
            if(sum(thetaindx & errorindx)>0){mean(temp$interim_sem[thetaindx & errorindx])}
            else{NA}
          })
        })
        
        plot(seq(-3,3,0.5), temp2[,"0"], ylim=c(0,5),type="l", , yaxt="n", xlab="True Theta", ylab="SEM")#, main=paste(assembly,ni))
        axis(2, at=seq(0,5,length.out=6), labels = seq(0,5,length.out=6))
        lines(seq(-3,3,0.5), temp2[,"1"], lty=2)
        lines(seq(-3,3,0.5), temp2[,"2"], lty=3)
        if(assembly=="Forward"){mtext(ni,side=3, line=0.5, cex=1.5)}
        if(ni=="decreasing"){mtext(assembly, side=4, line=1, cex=1.5)}
        
      }
    }
    if(routing=="NC"){routing="PI-NC"}
    if(routing=="theta"){routing="PI-\u0398"}
    if(routing=="information"){routing="MI"}
    
    title(main = paste(routing,"Routing",design, "Mean SEM by True Theta"),
          xlab = "Theta",
          ylab = "SEM",
          outer = TRUE, line = 1.5, cex.lab=2.5, cex.main=3)    
    
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    
    dev.off()
  }
}

###########################################################################################
###### Module Info by Design ################################################
###########################################################################################

#information looks identical for different assembly methods... check that assembly is producing/selecting different items

  for(design in Design){
    png(paste(design,"ModuleInfo.png", sep="_"), 1200,800)
    par(mfrow=c(length(Assembly),3),
         oma = c(3.5,3.5,3.5,3.5) + 0.1,
         mar = c(1.25,1.25,1.25,1.25) + 0.1)
    
    par(cex.lab=1.25)
    par(cex.axis=1.25)
    par(cex.main=1.5)
    ni="equal"
      for(assembly in Assembly)
        for(stage in 1:3){
          infoplotnames=paste(design,"Design | equal Number of Items |", assembly, "Assembly Priority |", Routing[3], "Routing")

          temp=matrix(sapply(unique(Out[[infoplotnames[1]]]$items[Out[[infoplotnames[1]]]$items$stage==stage,"module"]), function(module){
              temp=sapply(infoplotnames, function(infoplotname){
                print(infoplotname)
                colSums(
                    info(seq(-3,3,0.1),
                         Out[[infoplotname]]$items[Out[[infoplotname]]$items$stage==stage & Out[[infoplotname]]$items$module==module, c("a","b","c")]))/5})
              print(sum(temp))
              rowMeans(temp)
          }),nrow=length(seq(-3,3,0.1)))
          
          
          print(assembly)
          print(ncol(temp))
          for(col in 1:ncol(temp)){
            if(col==1 & ncol(temp)==1 & assembly == Assembly[1]){
              plot(seq(-3,3,0.1),temp[,1], type="l",
                    ylim=c(0,20),
                   #xlab="",
                   #ylab="information",
                   main=paste("Stage",stage))}
            if(col==1 & ncol(temp)!=1 & assembly == Assembly[1]){
              plot(seq(-3,3,0.1),temp[,1], type="l",
                   ylim=c(0,20),
                   #xlab="",
                   main=paste("Stage",ncol(temp)))}
            if(col==1 & ncol(temp)==1 & assembly %in% Assembly[c(2,3)]){
              plot(seq(-3,3,0.1),temp[,col], type="l",
                   ylim=c(0,20),
                   # xlab="",
                   # ylab="information"
                   )}
            if(col==1 & ncol(temp)!=1 & assembly %in% Assembly[c(2,3)]){
              plot(seq(-3,3,0.1),temp[,col], type="l",
                   ylim=c(0,20)#,
                   #xlab="",
                   )}
            if(col==1 & ncol(temp)==1 & assembly == Assembly[4]){
              plot(seq(-3,3,0.1),temp[,col], type="l",
                   ylim=c(0,20)#,
                   # xlab="theta",
                   # ylab="information"
                   )}
            if(col==1 & ncol(temp)!=1 & assembly == Assembly[4]){
              plot(seq(-3,3,0.1),temp[,col], type="l",
                   ylim=c(0,20)#,
                   #xlab="theta"
                   )}
            if(col!=1){lines(seq(-3,3,0.1),temp[,col], lty=col)}
            if(stage==3){mtext(assembly,side=4,line=0.625)}
          }
        }
    
    
    title(main = paste(design, "Module Information by Stage "),
          xlab = "Theta",
          ylab = "Information",
          outer = TRUE, line =1.5, cex.lab=2.5, cex.main=3)    
    
    # par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    # plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    # legend('bottomright', legend = c("Zero Errors", "One Error", "Two Errors"), lty=c(1:3), xpd = TRUE, horiz = TRUE, cex = 1.75, seg.len=3, bty = 'n')
    dev.off()
      }


###########################################################################################
###### Item Descriptives  ################################################
###########################################################################################

ItemDescriptives=
  round(t(sapply(c("a","b","c"), function(i){
      c("mean"=mean(Out[[1]]$items[,i]),
        "sd"=sd(Out[[1]]$items[,i]),
        "min"=min(Out[[1]]$items[,i]),
        "max"=max(Out[[1]]$items[,i]))
})),2)

###########################################################################################
###### Descriptives by Condition ##########################################################
###########################################################################################

###Bias, RMSE, SEM
#swap this in when ready for all calculation across all levels
StatsByCondition=sapply(unique(unlist(Design_Conditions)), USE.NAMES=TRUE, simplify=FALSE, function(condition){
  out=do.call(rbind,
                lapply(grep(condition,names(Out)), function(x){
                  x=merge(Out[[x]]$mst_admin, Out[[x]]$expectedpath, by="true_theta")
                  x=x[which(x$item_no==42),]
                  x$nerrors=apply(x,1,function(y){
                    sum(!(strsplit(y["path"],"-")[[1]] %in% strsplit(y["expected_path"],"-")[[1]]))})
                  data.frame(
                    "Deviance"=x$true_theta-x$interim_theta,
                    "TrueTheta"=x$true_theta,
                    "SEM"=x$interim_sem,
                    "NErrors"=x$nerrors)}))

    sapply(c("all","1","2"), USE.NAMES=TRUE, simplify=FALSE, function(nerror){
        if(nerror=="all"){nerror=c(0:2)}else{nerror=as.numeric(nerror)}
        list(
          "Overall"=list(
            "Bias"=mean(out$Deviance[out$NErrors %in% nerror]),
            "RMSE"=sqrt(sum(out$Deviance[out$NErrors %in% nerror]^2)/length(out$Deviance[out$NErrors %in% nerror])),
            "SEM"=mean(out$SEM[out$NErrors %in% nerror])),
          "Theta"=list(
            "Bias"=sapply(as.character(seq(-3,3,0.5)), function(t){mean(out$Deviance[out$TrueTheta==t & (out$NErrors %in% nerror)])}),
            "RMSE"=sapply(as.character(seq(-3,3,0.5)), function(t){sqrt(sum(out$Deviance[out$TrueTheta==t& (out$NErrors %in% nerror)])^2/length(out$Deviance[out$TrueTheta==t & (out$NErrors %in% nerror)]))}),
            "SEM"=sapply(as.character(seq(-3,3,0.5)), function(t){mean(out$SEM[out$TrueTheta==t & (out$NErrors %in% nerror)])})))
        })
      })
  

#tabulation: Overall
StatsByCondition_OverallTable=sapply(unique(unlist(Design_Conditions)), function(condition){
  unlist(lapply(c("all","1","2"),function(x){
  out=unlist(StatsByCondition[[condition]][[x]]$Overall)
  names(out)=paste(names(out),x,sep="_")
  out}))})

write.csv(round(StatsByCondition_OverallTable,3),paste0(getwd(),"/Results/Tables/StatsByCondition_OverallTable.csv"))

#############################################################################################

#swap this in when ready for all calculation across all levels
StatsByCondition=sapply(unique(unlist(Design_Conditions)), USE.NAMES=TRUE, simplify=FALSE, function(condition){
  out=do.call(rbind,
              lapply(grep(condition,names(Out)), function(x){
                x=merge(Out[[x]]$mst_admin, Out[[x]]$expectedpath, by="true_theta")
                x=x[which(x$item_no==42),]
                x$nerrors=apply(x,1,function(y){
                  sum(!(strsplit(y["path"],"-")[[1]] %in% strsplit(y["expected_path"],"-")[[1]]))})
                data.frame(
                  "Deviance"=x$true_theta-x$interim_theta,
                  "TrueTheta"=x$true_theta,
                  "SEM"=x$interim_sem,
                  "NErrors"=x$nerrors)}))
  
  sapply(c("0","1","2"), USE.NAMES=TRUE, simplify=FALSE, function(nerror){
    if(nerror=="all"){nerror=c(0:2)}else{nerror=as.numeric(nerror)}
    list(
      "Overall"=list(
        "Bias"=mean(out$Deviance[out$NErrors %in% nerror]),
        "RMSE"=sqrt(sum(out$Deviance[out$NErrors %in% nerror]^2)/length(out$Deviance[out$NErrors %in% nerror])),
        "SEM"=mean(out$SEM[out$NErrors %in% nerror])),
      "Theta"=list(
        "Bias"=sapply(as.character(seq(-3,3,0.5)), function(t){mean(out$Deviance[out$TrueTheta==t & (out$NErrors %in% nerror)])}),
        "RMSE"=sapply(as.character(seq(-3,3,0.5)), function(t){sqrt(sum(out$Deviance[out$TrueTheta==t& (out$NErrors %in% nerror)])^2/length(out$Deviance[out$TrueTheta==t & (out$NErrors %in% nerror)]))}),
        "SEM"=sapply(as.character(seq(-3,3,0.5)), function(t){mean(out$SEM[out$TrueTheta==t & (out$NErrors %in% nerror)])})))
  })
})


#tabulation: Overall
StatsByCondition_OverallTable=sapply(unique(unlist(Design_Conditions)), function(condition){
  unlist(lapply(c("0","1","2"),function(x){
    out=unlist(StatsByCondition[[condition]][[x]]$Overall)
    names(out)=paste(names(out),x,sep="_")
    out}))})

write.csv(round(StatsByCondition_OverallTable,3),paste0(getwd(),"/Results/Tables/StatsByCondition_OverallTable_Revised.csv"))




###########################################################################################
###### Descriptives by Path ###############################################################
###########################################################################################

###Bias, RMSE, SEM
#swap this in when ready for all calculation across all levels

StatsByPath=sapply(Design, USE.NAMES=TRUE, simplify=FALSE, function(design){

sapply(unique(unlist(Design_Conditions[,-1])), USE.NAMES=TRUE, simplify=FALSE, function(condition){
  out=do.call(rbind,
              lapply(grep(paste(c(design,condition),collapse=".*"),names(Out)), function(x){
                x=merge(Out[[x]]$mst_admin, Out[[x]]$expectedpath, by="true_theta")
                x=x[which(x$item_no==42),]
                x$nerrors=apply(x,1,function(y){
                  sum(!(strsplit(y["path"],"-")[[1]] %in% strsplit(y["expected_path"],"-")[[1]]))})
                data.frame(
                  "Deviance"=x$true_theta-x$interim_theta,
                  "TrueTheta"=x$true_theta,
                  "SEM"=x$interim_sem,
                  "NErrors"=x$nerrors,
                  "Path"=x$path)}))
  
  if(design=="1-3-3"){paths=apply(expand.grid("1",c(2:4),c(5:7)),1,paste, collapse="-")}
  if(design=="1-3-4"){paths=apply(expand.grid("1",c(2:4),c(5:8)),1,paste, collapse="-")}
  
  sapply(paths, USE.NAMES=TRUE, simplify=TRUE, function(path){
    list(
        "Bias"=mean(out$Deviance[out$Path == path]),
        "RMSE"=sqrt(sum(out$Deviance[out$Path == path]^2)/length(out$Deviance[out$Path == path])),
        "SEM"=mean(out$SEM[out$Path == path]),
        "PropError"=mean(out$NErrors[out$Path == path]!=0))
    })
  })
})


#tabulation: Overall
StatsByPath_OverallTable=sapply(Design, USE.NAMES=TRUE, simplify=FALSE, function(design){
  do.call(rbind,lapply(unique(unlist(Design_Conditions[,-1])), function(condition){

    out=StatsByPath[[design]][[condition]]
    rownames(out)=paste(condition,rownames(out),sep="_")
    out}))})

write.csv(t(StatsByPath_OverallTable$`1-3-3`),paste0(getwd(),"/Results/Tables/1-3-3_StatsByPath_OverallTable.csv"))
write.csv(t(StatsByPath_OverallTable$`1-3-4`),paste0(getwd(),"/Results/Tables/1-3-4_StatsByPath_OverallTable.csv"))


###########################################################################################
###### Statistic by Path  ################################################
###########################################################################################

PathThetaStats_All=lapply(Out, Path_Statistics, OrgByPathOrTheta="both", ErrorOnly=FALSE)
ThetaStats_All=lapply(Out, Path_Statistics, OrgByPathOrTheta="Theta", ErrorOnly=FALSE)
PathStats_All=lapply(Out, Path_Statistics, OrgByPathOrTheta="Path", ErrorOnly=FALSE)

PathThetaStats_ErrorOnly=lapply(Out, Path_Statistics, OrgByPathOrTheta="both", ErrorOnly=TRUE)
ThetaStats_ErrorOnly=lapply(Out, Path_Statistics, OrgByPathOrTheta="Theta", ErrorOnly=TRUE)
PathStats_ErrorOnly=lapply(Out, Path_Statistics, OrgByPathOrTheta="Path", ErrorOnly=TRUE)

PathStats_ErrorOnly_Out=matrix(NA,nrow=0,ncol=ncol(PathStats_ErrorOnly[[1]])+1)
for(writeobj in names(PathStats_ErrorOnly)){
  PathStats_ErrorOnly_Out=rbind(PathStats_ErrorOnly_Out,cbind(c("condition"=writeobj),PathStats_ErrorOnly[[writeobj]]))
}
for(col in names(PathStats_ErrorOnly_Out)){
  PathStats_ErrorOnly_Out[,col]=as.character(PathStats_ErrorOnly_Out[,col])}

write.csv(PathStats_ErrorOnly_Out,"PathStats_ErrorOnly_Out.csv")




##########


plot_module_info(
  Out[["1-3-3 Design | increasing Number of Items | Forward Assembly Priority | information Routing"]],
  Out[["1-3-3 Design | increasing Number of Items | Backward Assembly Priority | information Routing"]] )

plot_module_info(
  Out[["1-3-4 Design | equal Number of Items | Forward Assembly Priority | information Routing"]],
  Out[["1-3-4 Design | equal Number of Items | Backward Assembly Priority | information Routing"]] )


for_133_in <- stat(for_133_in)
for_133_eq <- stat(for_133_eq)
for_133_de <- stat(for_133_de)
for_134_in <- stat(for_134_in)
for_134_eq <- stat(for_134_eq)
for_134_de <- stat(for_134_de)
bac_133_in <- stat(bac_133_in)
bac_133_eq <- stat(bac_133_eq)
bac_133_de <- stat(bac_133_de)
bac_134_in <- stat(bac_134_in)
bac_134_eq <- stat(bac_134_eq)
bac_134_de <- stat(bac_134_de)


# frequency table
data <- rbind(for_133_in$routing_error_percent,
              for_133_eq$routing_error_percent,
              for_133_de$routing_error_percent,
              bac_133_in$routing_error_percent,
              bac_133_eq$routing_error_percent,
              bac_133_de$routing_error_percent,
              for_134_in$routing_error_percent,
              for_134_eq$routing_error_percent,
              for_134_de$routing_error_percent,
              bac_134_eq$routing_error_percent,
              bac_134_in$routing_error_percent,
              bac_134_de$routing_error_percent)

data = 

#Table 3 Coding:

###Information Routing
Table3InfoRouting <- do.call(rbind, lapply(Out[grep("information", names(Out))], function(x) x$routing_error_percent)) %>% 
    group_by(design, items_per_stage, assembly_priority) %>%
    summarize(freq = sum(routing_error_percent) / 100 * N_REP) %>%
    mutate(percent = round(freq / (length(THETA) * N_REP) * 100, 3))
write.csv(Table3InfoRouting,"Table3InfoRouting.csv")

Table3ThetaRouting <- do.call(rbind, lapply(Out[grep("theta", names(Out))], function(x) x$routing_error_percent)) %>% 
  group_by(design, items_per_stage, assembly_priority) %>%
  summarize(freq = sum(routing_error_percent) / 100 * N_REP) %>%
  mutate(percent = round(freq / (length(THETA) * N_REP) * 100, 3))
write.csv(Table3ThetaRouting,"Table3ThetaRouting.csv")

# match_mst conditional on no. of routing errors
list_of_mst <- list(for_133_in, for_133_eq, for_133_de,
                    bac_133_in, bac_133_eq, bac_133_de,
                    for_134_in, for_134_eq, for_134_de,
                    bac_134_in, bac_134_eq, bac_134_de)

for (i in 1:12) {
  per0 <- (100 - tapply(list_of_mst[[i]]$routing_error_percent$routing_error_percent, list_of_mst[[i]]$routing_error_percent$true_theta, sum)) / 100
  err0 <- sum(list_of_mst[[i]]$match_mst0$match_mst * per0) / sum(per0)
  per1 <- subset(list_of_mst[[i]]$routing_error_percent, n_routing_errors == 1, routing_error_percent) / 100
  err1 <- sum(list_of_mst[[i]]$match_mst1$match_mst * per1) / sum(per1)
  per2 <- subset(list_of_mst[[i]]$routing_error_percent, n_routing_errors == 2, routing_error_percent) / 100
  err2 <- sum(list_of_mst[[i]]$match_mst2$match_mst * per2) / sum(per2)
  print(c(err0, err1, err2))
}

sum(for_133_in$match_mst1$match_mst * subset(for_133_in$routing_error_percent, n_routing_errors == 1, routing_error_percent) / 100 * 500) / sum(subset(for_133_in$routing_error_percent, n_routing_errors == 1, routing_error_percent) / 100 * 500)
sum(for_133_in$match_mst2$match_mst * subset(for_133_in$routing_error_percent, n_routing_errors == 2, routing_error_percent) / 100 * 500) / sum(subset(for_133_in$routing_error_percent, n_routing_errors == 2, routing_error_percent) / 100 * 500)


#################################################################################################


Path_Statistics=function(input, OrgByPathOrTheta="both", ErrorOnly=FALSE){
  
  #install.packages("data.tree")
  library(data.tree)
  
  Path_Theta=list()
  if(OrgByPathOrTheta=="both"){Path_Theta=expand.grid(true_theta=unique(input$mst_admin$true_theta), path=sort(unique(input$mst_admin$path)))}
  if(OrgByPathOrTheta=="Path"){Path_Theta$path=sort(unique(input$mst_admin$path)); Path_Theta$true_theta=replicate(length(Path_Theta$path),NA)}
  if(OrgByPathOrTheta=="Theta"){Path_Theta$true_theta=unique(input$mst_admin$true_theta); Path_Theta$path=replicate(length(Path_Theta$true_theta),NA)}    
  
  #ErrorOnly=TRUE
  
  Path_Statistics=data.frame(t(mapply(function(path,true_theta){
    freq=NA 
    Bias=NA 
    RMSE=NA
    TrueThetaMSDRng=NA
    
    index=switch(OrgByPathOrTheta,
                 "both"=which(input$mst_admin$true_theta==true_theta & input$mst_admin$path==path),
                 "Path"=which(input$mst_admin$path==path),
                 "Theta"=which(input$mst_admin$true_theta==true_theta))
    
    
    pathdata=input$mst_admin[index,c("true_theta", "interim_theta", "path")]
    
    pathdata=merge(pathdata, input$expectedpath, by="true_theta")
    
    if(ErrorOnly!=FALSE){pathdata=pathdata[which(pathdata$path!=pathdata$expected_path),]}
    
    if(nrow(pathdata)>0){ 
      
      Theta_Error=pathdata$true_theta-pathdata$interim_theta
      
      Squared_Error=Theta_Error^2
      
      freq=nrow(pathdata)
      Bias=round(sum(Theta_Error)/length(Theta_Error),2)
      RMSE=round(sqrt(sum(Squared_Error)/length(Squared_Error)),2)
      TrueThetaMSDRng=paste0(
        round(mean(pathdata$true_theta),1),
        "(",round(sd(pathdata$true_theta),2),"), ",
        min(pathdata$true_theta, na.rm = TRUE)," - ",
        max(pathdata$true_theta, na.rm = TRUE))
      
    }
    
    list(
      # "path"=as.character(unique(pathdata$path)),
      # "expectedpath"=as.character(unique(pathdata$expected_path)),
      "path"=paste(as.character(unique(pathdata$path)), collapse=", "),
      "expectedpath"=paste(as.character(unique(pathdata$expected_path)), collapse = ", "),
      "TrueTheta"=true_theta,
      #"pathString"=pathString,
      "freq"=freq,
      "Bias"=Bias,
      "RMSE"=RMSE,
      "TrueThetaMSDRng"=TrueThetaMSDRng)
    
  },Path_Theta$path,Path_Theta$true_theta)))
  
  Path_Statistics=Path_Statistics[!is.na(Path_Statistics$freq),]
  
  ########################
  
}

PathThetaStats_All=lapply(Out, Path_Statistics, OrgByPathOrTheta="both", ErrorOnly=FALSE)
ThetaStats_All=lapply(Out, Path_Statistics, OrgByPathOrTheta="Theta", ErrorOnly=FALSE)
PathStats_All=lapply(Out, Path_Statistics, OrgByPathOrTheta="Path", ErrorOnly=FALSE)

PathThetaStats_ErrorOnly=lapply(Out, Path_Statistics, OrgByPathOrTheta="both", ErrorOnly=TRUE)
ThetaStats_ErrorOnly=lapply(Out, Path_Statistics, OrgByPathOrTheta="Theta", ErrorOnly=TRUE)
PathStats_ErrorOnly=lapply(Out, Path_Statistics, OrgByPathOrTheta="Path", ErrorOnly=TRUE)

###############################################################################################
##
###############################################################################################




#### Notes/ To-Do's
#1) apply to all output

output=list(
  "for_133_in"=for_133_in,
  "for_133_eq"=for_133_eq,
  "for_133_de"=for_133_de,
  "for_134_in"=for_134_in,
  "for_134_eq"=for_134_eq,
  "for_134_de"=for_134_de,
  "bac_133_in"=bac_133_in,
  "bac_133_eq"=bac_133_eq,
  "bac_133_de"=bac_133_de,
  "bac_134_in"=bac_134_in,
  "bac_134_eq"=bac_134_eq,
  "bac_134_de"=for_134_de)

formatted=lapply(output, function(element){print(Path_Statistics(element),"freq", "Bias", "RMSE", "TrueThetaMSDRng")})
for(element in names(formatted)){formatted[[element]][1:nrow(formatted[[element]]),"MSTadmin"]=element}

outout=do.call(rbind, formatted)
write.csv(outout, "AllPathsByError.csv")





