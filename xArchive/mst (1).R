pacman::p_load(dplyr,
               tidyr,
               ggplot2,
               catR,
               mstR,
               progress)
N_REP <- 15 # no. of replications per theta
THETA <- seq(-3, 3, 0.5)
TEST_LEN <- 42 # total test length
set.seed(123)

# Item generation
n_items <- 1500
item_bank <- genDichoMatrix(items = n_items,
                            model = "3PL",
                            aPrior = c("lnorm", 0, 0.25),
                            bPrior = c("norm", 0, 1),
                            cPrior = c("unif", 0.1, 0.2))
id <- 1:n_items
item_bank <- cbind(id, item_bank)

# Probability and Information
prob <- function(theta, items) {
    p <- outer(items$b, theta, "-")
    p <- items$c + (1 - items$c) / (1 + exp(1.7 * items$a * p))
    p # row = items, col = theta
}

info <- function(theta, items) {
    p <- prob(theta, items)
    i <- (1.7 * items$a * (p - items$c) / (1 - items$c))^2 * (1 - p) / p
    i # row = items, col = theta
}

# MST
mst_setup <- function(design, items_per_stage, assembly_priority) {
    theta <- switch(design,
                    "1-3-3" = c(0, -1, 0, 1, -1, 0, 1),
                    "1-3-4" = c(0, -1, 0, 1, -1.25, -0.75, 0.75, 1.25))
    len <- switch(items_per_stage,
                  "increasing" = TEST_LEN * c(1/6, 1/3, 1/2),
                  "equal"= TEST_LEN * c(1/3, 1/3, 1/3),
                  "decreasing" = TEST_LEN * c(1/2, 1/3, 1/6))
    design <- strsplit(design, split = "[-]") %>% unlist() %>% as.integer()
    n_stages <- length(design)
    n_modules <- sum(design)

    module <- NULL
    for (s in 1:n_stages)
        for (m in 1:design[s])
            module <- rbind(module, c(stage = s, module = m, len = len[s]))
    module <- data.frame(cbind(module, index = 1:nrow(module), theta = theta))

    x <- list(design = design,
              n_panels = 5,
              n_stages = n_stages,
              n_modules = n_modules,
              module = module,
              items_per_stage = items_per_stage, # text description for plotting later
              assembly_priority = assembly_priority)
    x
}

mst_assemble <- function(x, assembly_priority) {
    x$items <- data.frame(matrix(nrow = sum(x$module$len) * x$n_panels,
                                 ncol = 10,
                                 dimnames = list(c(),
                                                 c("id", "a", "b", "c", "d",
                                                   "stage", "module", "panel",
                                                   "index", "form"))))
    info <- info(theta = x$module$theta, items = item_bank)

    if (x$assembly_priority == "Forward")
        stage_order <- 1:x$n_stages
    else
        stage_order <- x$n_stages:1

    temp <- 1

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

sim <- function(x) {
    # creating a matrix to list item membership in each module
    modules <- matrix(0, nrow(x$items) / x$n_panels, x$n_modules)
    for (i in 1:x$n_modules) {
        len1 <- sum(subset(x$module, subset = index < i)$len)
        len2 <- subset(x$module, subset = index == i)$len
        modules[(1 + len1):(len1 + len2), i] <- 1
    }

    # creating the transition matrix to the MST design
    trans <- matrix(0, x$n_modules, x$n_modules)
    if (paste0(x$design, collapse = "-") == "1-3-3")
        trans[1, 2:4] <- trans[2, 5:6] <- trans[3, 5:7] <- trans[4, 6:7] <- 1
    else
        trans[1, 2:4] <- trans[2, 5:6] <- trans[3, 6:7] <- trans[4, 7:8] <- 1


    x$mst_admin <- data.frame(matrix(ncol = 7,
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

    pb <- progress_bar$new(total = length(THETA) * N_REP)
    temp <- 1
    for (i in 1:length(THETA)) {
        for (j in 1:N_REP) {
            # generate responses for MST and CAT
            resp <- genPattern(THETA[i], all_items, D = 1.7)

            # select a panel for MST randomly
            rand <- sample(1:x$n_panels, 1)

            # simulate MST
            mst <- randomMST(trueTheta = THETA[i],
                             itemBank = subset(x$items,
                                               subset = panel == rand,
                                               select = c(a, b, c, d)),
                             modules = modules,
                             transMatrix = trans,
                             responses = resp[which(x$items$panel == rand)],
                             start = list(fixModule = 1, D = 1.7),
                             test = list(method = "ML", D = 1.7,
                                         range = c(-3.5, 3.5),
                                         moduleSelect = "MFI"),
                             final = list(method = "ML", D = 1.7,
                                          range = c(-3.5, 3.5)),
                             allTheta = TRUE)

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

            # save output
            range <- (TEST_LEN * (temp - 1) + 1):(TEST_LEN * temp)
            x$mst_admin[range,] <- cbind(true_theta = mst$trueTheta,
                                         rep = j,
                                         item_no = 1:TEST_LEN,
                                         item_id = mst$testItems,
                                         interim_theta = mst$allTheta[, "th"],
                                         interim_sem = mst$allTheta[, "se"],
                                         panel = rand,
                                         path = as.numeric(paste0(mst$selected.modules, collapse = "")))

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

compare <- function(design, items_per_stage, assembly_priority) {
    x <- mst_setup(design,
                   items_per_stage = items_per_stage,
                   assembly_priority = assembly_priority)
    x <- mst_assemble(x)
    x <- sim(x)
    x <- stat(x)
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
for_133_in <- compare("1-3-3", "increasing", "Forward")
for_133_eq <- compare("1-3-3", "equal", "Forward")
for_133_de <- compare("1-3-3", "decreasing", "Forward")
for_134_in <- compare("1-3-4", "increasing", "Forward")
for_134_eq <- compare("1-3-4", "equal", "Forward")
for_134_de <- compare("1-3-4", "decreasing", "Forward")
bac_133_in <- compare("1-3-3", "increasing", "Backward")
bac_133_eq <- compare("1-3-3", "equal", "Backward")
bac_133_de <- compare("1-3-3", "decreasing", "Backward")
bac_134_in <- compare("1-3-4", "increasing", "Backward")
bac_134_eq <- compare("1-3-4", "equal", "Backward")
bac_134_de <- compare("1-3-4", "decreasing", "Backward")

save.image(paste0(Sys.Date(), ".RData"))

plot_module_info(for_133_eq, bac_133_eq)
plot_module_info(for_134_eq, bac_134_eq)


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
data <- data %>%
    group_by(n_routing_errors, design, items_per_stage, assembly_priority) %>%
    summarize(freq = sum(routing_error_percent) / 100 * N_REP) %>%
    mutate(percent = round(freq / (length(THETA) * N_REP) * 100, 3))

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
