# Copyright (C) 2016  Aghiles Salah. All rights reserved.
# License: Apache 2.0 License.


# Normalized Mutual Information
NMI <- function(row_cluster, true_labels, n){
    N_kl = as.matrix(table(row_cluster, true_labels))
    N_k = rowSums(N_kl)
    N_l = colSums(N_kl)
    Num = 0
    denum_k = 0
    for (i in 1 : nrow(N_kl)) {
        denum_k = denum_k + (N_k[i] / n) * log((N_k[i] / n))
        for (j in 1 : ncol(N_kl)) {
            if (N_kl[i, j] != 0) {
                num = log(n) + log(N_kl[i, j])
                den = log(N_k[i]) + log(N_l[j])
                Num = Num + (N_kl[i, j] / n) * (num - den)
            }
        }
    }
    denum_l = 0
    for (j in 1 : ncol(N_kl)) {
        denum_l = denum_l + (N_l[j] / n) * log((N_l[j] / n))
    }
    resnmi = Num / (sqrt(denum_k * denum_l))
    resnmi
}

# Partition generator
par_gen <-
function(n, k, nb_par=1){

    par = as.integer(sample(as.numeric(1 : k), n, replace=TRUE))
    if (nb_par > 1) {
        for (i in 2 : nb_par) {
            par = rbind(par, as.integer(sample(as.numeric(1 : k), n, replace=TRUE)))
        }
        par = as.matrix(par)
    }
    par
}
