# Normalized Mutual Information
NMI <- function(rowcluster, truelabels, n){
    N_kl = as.matrix(table(rowcluster, truelabels))
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