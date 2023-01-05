import ot
import numpy as np
import scipy.io as sio
import time
import sys

path = sys.argv[1] 
print(path)

ot_cost_matrix = np.loadtxt(open(path + "/ot_cost.csv", "rb"), delimiter=",")
gene_expr_matrix = sio.mmread(path + "/gene_expression.mtx").todense()
emd_mat = np.zeros((gene_expr_matrix.shape[1], gene_expr_matrix.shape[1]))

N_cal = (np.square(emd_mat.shape[0])- emd_mat.shape[0])/2
start = time.time()
COUNTER = 0
FLAG = 1
for i in range(0, emd_mat.shape[0]):
  for j in range(i+1, emd_mat.shape[0]):
    COUNTER += 1
    gene_i = np.array(gene_expr_matrix[:, i]).reshape(-1)
    gene_j = np.array(gene_expr_matrix[:, j]).reshape(-1)
    gene_i = gene_i/sum(gene_i)
    gene_j = gene_j/sum(gene_j)

    emd_mat[i, j] = ot.emd2(gene_i[np.nonzero(gene_i)],
                        gene_j[np.nonzero(gene_j)], 
                        np.array(ot_cost_matrix)[np.nonzero(gene_i)[0],:][:,np.nonzero(gene_j)[0]],
                        numItermax=5000)
    if COUNTER >= FLAG/20*N_cal:
        print(str(np.round(FLAG/20*100,1)) + "% finished.")
        FLAG += 1

end = time.time()
print('elapsed time = %.4f s' % (end-start))


emd_mat2 = emd_mat + np.transpose(emd_mat)
np.savetxt(path + "/emd.csv", emd_mat2, delimiter=",")





