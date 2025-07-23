import open3d as o3d
import numpy as np
import random
import sys
import copy
from sklearn.decomposition import PCA
import time
from sklearn.neighbors import KDTree
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.sparse import lil_matrix
from scipy.sparse.csgraph import connected_components

volSize = 0.05
disThreshold = 1
lengthThreshold = 0.05
fitThreshold = 2.5 * volSize

def Distance(p1, p2):
    return np.linalg.norm(p1 - p2)

def notTooClose2(p1, p2):
    return np.linalg.norm(p1 - p2) > disThreshold

def notTooClose3(p):
    return notTooClose2(p[0], p[1]) and notTooClose2(p[0], p[2]) and notTooClose2(p[1], p[2])

def farAway(l1, l2):
    ss = np.linalg.norm(l1)
    tt = np.linalg.norm(l2)
    return abs(ss - tt) > lengthThreshold * (ss + tt)

def display_inlier_outlier(pcd, ind):
    inlier_cloud = pcd.select_by_index(ind)
    outlier_cloud = pcd.select_by_index(ind, invert=True)
    outlier_cloud.paint_uniform_color([1, 0, 0])
    inlier_cloud.paint_uniform_color([0.8, 0.8, 0.8])
    o3d.visualization.draw_geometries([inlier_cloud, outlier_cloud])

def estimateAvgDis(points):
    sample = random.sample(list(points), 10)
    dis = [Distance(p1, p2) for p1 in sample for p2 in sample if (p1 != p2).all()]
    global disThreshold
    disThreshold = np.mean(dis)/2

def topological_noise_filter(pcd, epsilon=0.1, tau=5, sigma=0.1):
    points = np.asarray(pcd.points)
    n = len(points)
    A = lil_matrix((n, n))
    tree = cKDTree(points)
    
    for i in range(n):
        neighbors = tree.query_ball_point(points[i], epsilon)
        for j in neighbors:
            if i != j:
                A[i,j] = 1
    
    degrees = np.array(A.sum(axis=1)).flatten()
    densities = np.zeros(n)
    
    for i in range(n):
        neighbor_indices = A[i].nonzero()[1]
        if len(neighbor_indices) > 0:
            distances = np.linalg.norm(points[neighbor_indices] - points[i], axis=1)
            densities[i] = np.sum(np.exp(-(distances**2)/(2*sigma**2)))
    
    filtered_indices = [i for i in range(n) if degrees[i] >= tau and densities[i] > np.mean(densities)/2]
    return pcd.select_by_index(filtered_indices)

def adaptive_downsample(pcd, initial_voxel_size=0.05, alpha=1.1, target_threshold=10000):
    points = np.asarray(pcd.points)
    N = len(points)
    V = initial_voxel_size
    
    while N > target_threshold:
        pcd = pcd.voxel_down_sample(voxel_size=V)
        N = len(pcd.points)
        V *= alpha
    
    return pcd

def pca_fpfh(fpfh_data, k=32):
    X = fpfh_data.T
    mu = np.mean(X, axis=0)
    X_centered = X - mu
    cov_matrix = np.cov(X_centered, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
    idx = np.argsort(eigenvalues)[::-1]
    W = eigenvectors[:, idx[:k]]
    Y = X_centered @ W
    return Y.T

def TCSAC(src_points, tgt_points, src_fpfh, tgt_fpfh, lmin=0.1, lmax=1.0, theta_min=30, theta_max=120):
    fmean = np.mean(src_fpfh, axis=0)
    D = np.array([np.linalg.norm(f - fmean) for f in src_fpfh])
    candidate_indices = np.argsort(D)[-int(0.3*len(D)):]
    
    best_R = None
    best_T = None
    max_inliers = 0
    
    for _ in range(1000):
        idx1, idx2, idx3 = np.random.choice(candidate_indices, 3, replace=False)
        P1, P2, P3 = src_points[idx1], src_points[idx2], src_points[idx3]
        
        a = Distance(P2, P3)
        b = Distance(P1, P3)
        c = Distance(P1, P2)
        
        if not (lmin < a < lmax and lmin < b < lmax and lmin < c < lmax):
            continue
            
        alpha = np.arccos((b**2 + c**2 - a**2)/(2*b*c)) * 180/np.pi
        beta = np.arccos((a**2 + c**2 - b**2)/(2*a*c)) * 180/np.pi
        gamma = np.arccos((a**2 + b**2 - c**2)/(2*a*b)) * 180/np.pi
        
        if not (theta_min < alpha < theta_max and theta_min < beta < theta_max and theta_min < gamma < theta_max):
            continue
            
        src_tri = np.array([P1, P2, P3])
        tgt_tri = []
        
        for idx in [idx1, idx2, idx3]:
            f = src_fpfh[idx]
            dists = np.linalg.norm(tgt_fpfh - f, axis=1)
            tgt_idx = np.argmin(dists)
            tgt_tri.append(tgt_points[tgt_idx])
            
        tgt_tri = np.array(tgt_tri)
        
        R, T = calculateTrans(src_tri, tgt_tri)
        transformed = (R @ src_points.T + T).T
        dists = np.min(np.linalg.norm(transformed[:, np.newaxis] - tgt_points, axis=2), axis=1)
        inliers = np.sum(dists < fitThreshold)
        
        if inliers > max_inliers:
            max_inliers = inliers
            best_R = R
            best_T = T
            
    return best_R, best_T

def prepare(path, color, downSave=False, outlier=False, draw=False, pcaTag=False):
    pcd = o3d.io.read_point_cloud(path)
    pcd.paint_uniform_color(color)
    oldPcd = copy.deepcopy(pcd)
    
    pcd = topological_noise_filter(pcd)
    
    if downSave:
        pcd = adaptive_downsample(pcd)
    else:
        pcd = pcd.voxel_down_sample(voxel_size=volSize)
            
    if outlier:
        pcd, ind = pcd.remove_statistical_outlier(nb_neighbors=20, std_ratio=0.95)
        if draw:
            display_inlier_outlier(oldPcd, ind)
            
    pcd.estimate_normals(o3d.geometry.KDTreeSearchParamKNN(knn=30))
    KDT = o3d.geometry.KDTreeFlann(pcd)
    fpfh = o3d.pipelines.registration.compute_fpfh_feature(pcd, o3d.geometry.KDTreeSearchParamKNN(knn=200))
    
    if pcaTag:
        fpfh.data = pca_fpfh(fpfh.data, pcaTag)
        
    fpfhKDT = o3d.geometry.KDTreeFlann(fpfh)
    global fitThreshold
    fitThreshold = 2.5 * volSize

    return KDT, fpfhKDT, oldPcd, pcd, fpfh.data.T

def calculateTrans(src, tgt):
    assert src.shape == tgt.shape
    src = np.array(src)
    tgt = np.array(tgt)
    num = src.shape[0]
    srcAvg = np.mean(src, axis=0).reshape(1,3)
    tgtAvg = np.mean(tgt, axis=0).reshape(1,3)
    src -= np.tile(srcAvg, (num, 1))
    tgt -= np.tile(tgtAvg, (num, 1))
    H = np.transpose(src) @ tgt
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T @ U.T
    T = -R @ srcAvg.T + tgtAvg.T
    return R, T

def ICP(src, tgt):
    limit = fitThreshold
    retR = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    retT = np.array([[0], [0], [0]])
    trace = []
    for _ in range(400):
        tgtCorr = []
        srcCorr = []
        for point in src:
            k, idx, dis2 = tgtKDT.search_knn_vector_3d(point, knn=1)
            if dis2[0] < (limit)**2:
                srcCorr.append(point)
                tgtCorr.append(tgt[idx[0]])
        trace.append([limit, len(srcCorr)])
        R, T = calculateTrans(np.array(srcCorr), np.array(tgtCorr))
        retR = R @ retR
        retT = R @ retT + T
        src = np.transpose((R @ src.T) + np.tile(T, (1, srcNum)))
        limit = (limit - fitThreshold/1.5) * 0.95 + fitThreshold/1.5
        if len(trace) > 50 and len(set([x[1] for x in trace[-20:]])) == 1:
            break
    return retR, retT, srcCorr, tgtCorr

def RANSAC():
    maxCount = 0
    jisuan=0
    j = 0
    while True:
        j += 1
        srcCorr = random.sample(range(srcNum), 3)
        if not notTooClose3([srcPoints[x] for x in srcCorr]):
            continue
        tgtCorr = []
        for id in srcCorr:
            k, idx, dis2 = tgtFpfhKDT.search_knn_vector_xd(srcFpfhq[id], knn=1)
            tgtCorr.append(idx[0])
            
        if True in [farAway(srcPoints[i[0]] - srcPoints[j[0]], 
                    tgtPoints[i[1]] - tgtPoints[j[1]]) 
                for i in zip(srcCorr, tgtCorr) 
                for j in zip(srcCorr, tgtCorr)]:
            continue
        jisuan += 1
        R, T = calculateTrans(np.array([srcPoints[i] for i in srcCorr]), 
                          np.array([tgtPoints[i] for i in tgtCorr]))
        A = np.transpose((R @ srcPoints.T) + np.tile(T, (1, srcNum)))
        
        count = 0
        for point in range(0, srcNum, 1):           
            k, idx, dis2 = tgtKDT.search_hybrid_vector_3d(A[point], 
                                                      radius=fitThreshold, max_nn=1)
            count += k
        if count > maxCount:
            maxCount = count
            bestR, bestT = R, T
        if jisuan > 50 and j > 1000:
            break
        
    return bestR, bestT

if __name__ == "__main__":
    start = time.time()
    srcPath = "data/1/2plys/1-1.ply"
    tgtPath = "data/1/2plys/1-2.ply"

    srcKDT, srcFpfhKDT, oldSrc, src, srcFpfhq = prepare(srcPath, [1, 0.706, 0], downSave=True, outlier=False)
    tgtKDT, tgtFpfhKDT, oldTgt, tgt, tgtFpfhq = prepare(tgtPath, [0, 0.651, 0.929], outlier=False)

    srcPoints = np.array(src.points)
    tgtPoints = np.array(tgt.points)
    srcNum = np.asarray(srcPoints).shape[0]
    tgtNum = np.asarray(tgtPoints).shape[0]

    estimateAvgDis(srcPoints)
    
    R1, T1 = TCSAC(srcPoints, tgtPoints, srcFpfhq, tgtFpfhq)
    
    A = np.transpose((R1 @ srcPoints.T) + np.tile(T1, (1, srcNum)))
    A=o3d.utility.Vector3dVector(A)
    src.points = A
    
    R2, T2, srcCorr, tgtCorr = ICP(np.array(A), tgtPoints)

    tgtCorr = np.array(tgtCorr).reshape(-1, 3)
    R = R2 @ R1
    T = R2 @ T1 + T2

    A = np.array(oldSrc.points)
    A = np.transpose((R @ np.array(A).T) + np.tile(T, (1, A.shape[0])))
    A = o3d.utility.Vector3dVector(A)
    oldSrc.points = A

    oldSrc.paint_uniform_color([1, 0.706, 0])
    oldTgt.paint_uniform_color([0, 0.651, 0.929])
    o3d.visualization.draw_geometries([oldSrc, oldTgt])

# if __name__ == "__main__":
#     # srcPath = sys.argv[1]
#     # tgtPath = sys.argv[2]
#     # savePath = sys.argv[3]
#     start = time.time()
#     # srcPath = "data/3DMatch/0ply/2.ply"
#     # tgtPath = "data/3DMatch/0ply/3.ply"

#     srcPath = "data/1/2plys/1-1.ply"
#     tgtPath = "data/1/2plys/1-2.ply"

#     srcKDT, srcFpfhKDT, oldSrc, src, srcFpfhq = prepare(srcPath, [1, 0.706, 0], downSave=True, outlier=False)
#     print("srcFpfh:",srcFpfhq.shape[0])
#     print("src:",len(src.points))
#     tgtKDT, tgtFpfhKDT, oldTgt, tgt, tgtFpfhq = prepare(tgtPath, [0, 0.651, 0.929], outlier=False)

#     print("tgtFpfh:",tgtFpfhq.shape[0])
#     print("tgt:",len(tgt.points))
   
#     ###########################################################
#     # KDT：点云的 KD 树。
#     # fpfhKDT：FPFH 特征的 KD 树。
#     # oldPcd：未经处理的原始点云数据的深拷贝。
#     # pcd：经过预处理的点云数据。
#     # fpfh.data.T：FPFH 特征的转置，是一个 numpy 数组，存储了点云的 FPFH 特征向量。



#     # print("fpfh wd:",srcFpfh)

#     srcPoints = np.array(src.points)
#     tgtPoints = np.array(tgt.points)
#     ###########################################################
#     # src2 = np.load("data/workpiece/cloud_bin_1_fpfh.npz")
#     # tgt2 = np.load("data/workpiece/cloud_bin_2_fpfh.npz")
#     # srcPoints = src2["xyz"]
#     # tgtPoints = tgt2["xyz"]
#     # srcFpfh2 = src2["features"]
#     # tgtFpfh2 = tgt2["features"]
#     # print("tgtKDT",len(tgtKDT.points))
#     # print("tgtFpfhKDT",tgtFpfhKDT.shape[0])
#     ###########################################################


#     srcNum = np.asarray(srcPoints).shape[0]
#     tgtNum = np.asarray(tgtPoints).shape[0]


#     estimateAvgDis(srcPoints)

#     # o3d.visualization.draw_geometries([src, tgt])
    
#     print("srcNum: %d\ntgtNum: %d" % (srcNum, tgtNum))
#     R1, T1 = RANSAC()
#     ###########################################################
#     # RANSAC()后的旋转矩阵和平移向量

#     A = np.transpose((R1 @ srcPoints.T) + np.tile(T1, (1, srcNum)))
#     A=o3d.utility.Vector3dVector(A)
#     src.points = A
#     # o3d.visualization.draw_geometries([src, tgt])
    
    
#     R2, T2, srcCorr, tgtCorr = ICP(np.array(A), tgtPoints)
#     ###########################################################
#     # ICP()后的旋转矩阵和平移向量，以及对应源点云和对应目标点云

#     tgtCorr = np.array(tgtCorr).reshape(-1, 3)
#     # print("123123:",tgtCorr.size)

#     R = R2 @ R1
#     T = R2 @ T1 + T2


# ################################################################################################################# 

#     A = np.array(oldSrc.points)
#     # print("A:\n",A)
#     A = np.transpose((R @ np.array(A).T) + np.tile(T, (1, A.shape[0])))
#     A = o3d.utility.Vector3dVector(A)
#     oldSrc.points = A
#     print("全局配准花费了： %.3f 秒.\n" % (time.time() - start))

#     oldSrc.paint_uniform_color([1, 0.706, 0])
#     oldTgt.paint_uniform_color([0, 0.651, 0.929])
#     o3d.visualization.draw_geometries([oldSrc, oldTgt])
#     # o3d.io.write_point_cloud("results/" + savePath, oldSrc + oldTgt)


# #################################################################################################################
#     # TT = compose_transform_matrix(R, T)
#     # print("TT:\n",TT)
#     # # 计算变换矩阵 T 的逆矩阵
#     # TT_inv = np.linalg.inv(TT)
#     # print("TT_inv:\n",TT_inv)
#     # np.save('data/3DMatch/temp/2tc.npy',TT)
#     # R_inv = TT_inv[:3, :3]
#     # T_inv = TT_inv[:3, 3]
#     # # srcCorr = np.transpose((R_inv @ np.array(tgtCorr).T) + np.tile(T_inv, (1, tgtCorr.shape[0])))
#     # A_minus_T = tgtCorr - T.T
#     # R_inv = np.linalg.inv(R)
#     # A_rotated = R_inv @ A_minus_T.T
#     # srcCorr= np.transpose(A_rotated)
#     # np.save('data/3DMatch/2corr/2tc.npy', srcCorr)
#     # np.save('data/3DMatch/2corr/3tc.npy', tgtCorr)

# #################################################################################################################
