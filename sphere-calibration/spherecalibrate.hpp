#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fstream>
#include <vector>
#include <cmath>
#include <Eigen/SVD>
#include <math.h>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/LU>
#include <Eigen/QR>
#include <stdlib.h>

#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/core/persistence.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d.hpp>

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix4d Pose;
typedef std::vector<Pose> Poses;

Eigen::Matrix3d skew(Eigen::Vector3d u);
Pose SolveX(Poses &AS, Poses &BS);
Pose SolveX_(Poses &AS, Poses &BS);
int GetPose(std::vector<std::vector<double>> &robot_points, std::vector<std::vector<double>> &cam_points, Poses &AS, Poses &BS);
cv::Mat EMat2RT(Eigen::Matrix4d &Rs);
int GetPose_(std::vector<std::vector<double>> &robot_points, std::vector<std::vector<double>> &cam_points, Poses &robots, Poses &cams);
cv::Mat Solve_LS(Poses &robot_mats, Poses &cam_mats);

bool isRotationMatrix(const cv::Mat & R);

int ComputeRotMat(std::vector<std::vector<double>> &robot, std::vector<std::vector<double>> &sphere, Eigen::Matrix3d &Rs);
int ComputeTrans(std::vector<std::vector<double>> &poses_rob, std::vector<std::vector<double>> &poses_cam, Eigen::Vector3d &Ts, Eigen::Matrix3d &Rs);
int Test(std::vector<std::vector<double>> &robot, std::vector<std::vector<double>> &sphere, cv::Mat &hcg);

cv::Mat R_T2RT(Eigen::Matrix3d &Rs, Eigen::Vector3d &Ts);
void csvRead(const std::string& inputfile, std::vector<std::vector<double>>& _robot_points, const int _linecount);
void outputEyeHandParams(std::string& yaml_file, cv::Mat& hcg);