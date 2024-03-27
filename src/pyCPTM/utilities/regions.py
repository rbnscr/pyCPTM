import numpy as np
# import matplotlib.pyplot as plt

def region_cylinder(compartmentModel, origin, end, radius):
    # https://math.stackexchange.com/questions/3518495/check-if-a-general-point-is-inside-a-given-cylinder
    cellCoordinates = np.copy(np.asarray(compartmentModel.cellCoordinates))
    origin = np.asarray(origin)
    end = np.asarray(end)
    region = np.zeros(compartmentModel.nCell,dtype=np.int64)

    e = end - origin
    m = np.cross(origin,end)

    # d must be less than or equal to radius
    d = np.linalg.norm(m + np.cross(e,cellCoordinates),axis = 1)/np.linalg.norm(e)

    # return d
    rq = cellCoordinates + np.cross(e,(m+np.cross(e,cellCoordinates)))/np.linalg.norm(e)**2
    # return rq
    # setup solver
    a = np.array([[1 + np.dot(origin, origin),1 + np.dot(origin, end)],[1 + np.dot(origin, end),1 + np.dot(end, end)]])
    b = np.array([1+np.dot(rq,origin),1+np.dot(rq,end)])
    w = np.linalg.solve(a,b).T


    for i in range(len(w)):
        if (d[i] <= radius and w[i,0] <= 1 and 0 <= w[i,0] and w[i,1] <= 1 and 0 <= w[i,1]):
            region[i] = 1

    print(f"{np.sum(region)} cells found in region.")
    return region


# test = region_cylinder(cm,[0.005,0.05,0.01],[0.005,0.05,0.09],0.01)

# fig = plt.figure()
# ax = fig.add_subplot(projection="3d")
# im = ax.scatter(
#     cm.cellCoordinates[:, 0],
#     cm.cellCoordinates[:, 1],
#     cm.cellCoordinates[:, 2],
#     c=test,
# )