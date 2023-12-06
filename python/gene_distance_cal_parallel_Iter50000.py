import ot
import numpy as np
import scipy.io as sio
import time
import sys
import multiprocessing
from multiprocessing.pool import Pool
import time

path = sys.argv[1]
print(path)

ot_cost_matrix = np.loadtxt(open(path + "/ot_cost.csv", "rb"), delimiter=",")
gene_expr_matrix = sio.mmread(path + "/gene_expression.mtx").todense()


def cal_ot(i, j):
    # report a message
    #print(f'Compute the distance between gene {i} and gene {j}', flush=True)
    # compute the distance
    gene_i = np.array(gene_expr_matrix[:, i]).reshape(-1)
    gene_j = np.array(gene_expr_matrix[:, j]).reshape(-1)
    gene_i = gene_i/sum(gene_i)
    gene_j = gene_j/sum(gene_j)

    emd_dist = ot.emd2(gene_i[np.nonzero(gene_i)],
                        gene_j[np.nonzero(gene_j)], 
                        np.array(ot_cost_matrix)[np.nonzero(gene_i)[0],:][:,np.nonzero(gene_j)[0]],
                        numItermax=50000)
    # return the generated value
    return (emd_dist)


if __name__ == '__main__':
    start_time = time.perf_counter()
    # create and configure the process pool
    with Pool(processes=8) as pool:
        # prepare arguments
        items = [(i, j) for i in range(0,gene_expr_matrix.shape[1]-1) for j in range(i+1,gene_expr_matrix.shape[1])]
        # execute tasks and process results in order
        result = pool.starmap(cal_ot, items)
        #for result in pool.starmap(cal_ot, items):
            #print(f'Got result: {result}', flush=True)
    finish_time = time.perf_counter()
    print("Program finished in {} seconds - using multiprocessing".format(finish_time-start_time))
    print("---")
    
ind = 0
emd_mat = np.zeros((gene_expr_matrix.shape[1], gene_expr_matrix.shape[1]))
for i in range(0,gene_expr_matrix.shape[1]-1): 
    for j in range(i+1,gene_expr_matrix.shape[1]):
        emd_mat[i, j] = result[ind]
        ind += 1

emd_mat2 = emd_mat + np.transpose(emd_mat)
np.savetxt(path + "/emd.csv", emd_mat2, delimiter=",")
