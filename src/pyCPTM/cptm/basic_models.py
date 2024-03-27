def tracer_model_sparse(t, y, Q):
    # t: time
    # y: Tracer concentration
    # Q: Exchange flow
    # Delafosse (2014)
    return scalar_conv(Q, y)

def scalar_conv(Q, y):
    return Q.dot(y)