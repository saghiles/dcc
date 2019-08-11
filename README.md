# Directional Co-clustering with a Conscience (DCC)

The code in this repository implements the DCC co-clustering algorithm presented in the paper: 
- **[Model-based von Mises-Fisher Co-clustering with a Conscience](https://epubs.siam.org/doi/pdf/10.1137/1.9781611974973.28)** <br/>Aghiles Salah, Mohamed Nadif<br/>*SIAM International Conference on Data Mining*. 2017

## Usage example
The following code is an example of how we can fit DCC to a real-world dataset and assess the quality of the obtained clustering.
```R
# Load NG2 datatset
ng2 <- readMat("./data/NG2.mat")
ng2_mat <- ng2$mat
ng2_class_labels <- as.vector(ng2$class.labels)


# TF-IDF representation
ng2_tfidf <- tf_idf(ng2_mat)

# Fit DCC to NG2 data
res= dcc(X=ng2_tfidf, k=2, iter.max=100, n_init=5, stoch_iter.max=70)

# Compare DCC clutering with the ground truth
NMI(res$rowcluster,ng2_class_labels,dim(ng2_tfidf)[1])
adjustedRandIndex(res$rowcluster,ng2_class_labels)
```
**Output:**

    NMI = 0.714    ARI = 0.810
Results may vary slightly from one run to another due to random initialization as well as stochastic assignments in early iterations. For more details, please refer to section 6.3 in DCC paper.   
