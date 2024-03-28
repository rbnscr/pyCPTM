from scipy.sparse import csc_matrix

def tracer_model_sparse(t, y, Q: csc_matrix):
    # t: time
    # y: Conservative tracer concentration
    # Q: Exchange flow
    # Delafosse (2014)
    return scalar_conv(Q, y)

def scalar_conv(Q: csc_matrix, y):
    return Q.dot(y)