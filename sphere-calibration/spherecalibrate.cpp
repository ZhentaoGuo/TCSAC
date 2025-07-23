#include "spherecalibrate.hpp"

static int NUM_PointCloud = 9; //9;

Eigen::Matrix3d GetTFromPose(double rx, double ry, double rz) {
    double rx_ = rx * 180 / M_PI;
    double ry_ = ry * 180 / M_PI;
    double rz_ = rz * 180 / M_PI;
    Eigen::MatrixXd Rx(3, 3);
    Eigen::MatrixXd Ry(3, 3);
    Eigen::MatrixXd Rz(3, 3);
    Rx << 1, 0, 0, 0, cos(rx_), -sin(rx_), 0, sin(rx_), cos(rx_);
    Ry << cos(ry_), 0, sin(ry_), 0, 1, 0, -sin(ry_), 0, cos(ry_);
    Rz << cos(rz_), -sin(rz_), 0, sin(rz_), cos(rz_), 0, 0, 0, 1;
    Eigen::Matrix3d rot_mat(3, 3);
    rot_mat = Rz * Ry * Rx;
    return rot_mat;
}
Eigen::Matrix3d rotZYZ(double rz_1, double ry, double rz_2)
{
	double rz_1_ = rz_1 * 180 / M_PI;
    double ry_ = ry * 180 / M_PI;
    double rz_2_ = rz_2 * 180 / M_PI;
    auto _1 = rz_1_, _2 = ry_, _3 = rz_2_;
    auto c1 = std::cos(_1); auto c2 = std::cos(_2); auto c3 = std::cos(_3);
    auto s1 = std::sin(_1); auto s2 = std::sin(_2); auto s3 = std::sin(_3);
    cv::Mat_<double> rot_mat = (cv::Mat_<double>(3, 3) << c1 * c2 * c3 - s1 * s3, -1 * c3 * s1 - c1 * c2 * s3, c1 * s2,
                      c1 * s3 + c2 * c3 * s1, c1 * c3 - c2 * s1 * s3, s1 * s2,
                      -1* c3 * s2, s2 * s3, c2);
    if (!isRotationMatrix(rot_mat)) {
        cv::error(cv::Error::StsAssert, "Euler angle can not convert to rotated matrix",
            __FUNCTION__, __FILE__, __LINE__);
    }
	Eigen::Matrix3d rot_mat_;
	rot_mat_ << rot_mat[0][0], rot_mat[0][1], rot_mat[0][2],
			    rot_mat[1][0], rot_mat[1][1], rot_mat[1][2],
				rot_mat[2][0], rot_mat[2][1], rot_mat[2][2];
    //std::cout << rot_mat_ << std::endl;
    return rot_mat_;
}

bool isRotationMatrix(const cv::Mat & R) {
	cv::Mat tmp33 = R({ 0,0,3,3 });
	cv::Mat shouldBeIdentity;
	shouldBeIdentity = tmp33.t()*tmp33;
	cv::Mat I = cv::Mat::eye(3, 3, shouldBeIdentity.type());
	//std::cout << cv::norm(I, shouldBeIdentity) << std::endl;
	return  cv::norm(I, shouldBeIdentity) < 1e-6;
}
int ComputeRotMat(std::vector<std::vector<double>> &robot, std::vector<std::vector<double>> &sphere, Eigen::Matrix3d &Rs) {
    Eigen::Matrix3d Ro;
	std::vector<Eigen::Vector3d> Pcs(4), Tos(4);// 四组点
    Ro = GetTFromPose(robot[0][3], robot[0][4], robot[0][5]);
    for (size_t i = 0; i < 4; i++) {
        Pcs[i](0) = robot[i][0]; Pcs[i](1) = robot[i][1]; Pcs[i](2) = robot[i][2];
        Tos[i](0) = sphere[i][0]; Tos[i](1) = sphere[i][1]; Tos[i](2) = sphere[i][2];
    }
	Eigen::Matrix3d P;
	Eigen::Matrix3d T;
	P << Pcs[0](0) - Pcs[1](0), Pcs[0](1) - Pcs[1](1), Pcs[0](2) - Pcs[1](2),
	     Pcs[0](0) - Pcs[2](0), Pcs[0](1) - Pcs[2](1), Pcs[0](2) - Pcs[2](2),
		 Pcs[0](0) - Pcs[3](0), Pcs[0](1) - Pcs[3](1), Pcs[0](2) - Pcs[3](2);
	T << Tos[1](0) - Tos[0](0), Tos[1](1) - Tos[0](1), Tos[1](2) - Tos[0](2),
	     Tos[2](0) - Tos[0](0), Tos[2](1) - Tos[0](1), Tos[2](2) - Tos[0](2),
		 Tos[3](0) - Tos[0](0), Tos[3](1) - Tos[0](1), Tos[3](2) - Tos[0](2);
	std::cout << "P: \n" << P << "\nT: \n" << T << std::endl;
	Eigen::MatrixXd B;
	B = T * Ro;
	Rs = ((P).colPivHouseholderQr().solve(B)).transpose();
	std::cout << " QR RS: \n" << Rs << std::endl;
	cv::Mat_<double> R1 = (cv::Mat_<double>(3, 3) << Rs(0, 0), Rs(0, 1), Rs(0, 2),
		Rs(1, 0), Rs(1, 1), Rs(1, 2),
		Rs(2, 0), Rs(2, 1), Rs(2, 2));
	if (isRotationMatrix(R1)) { std::cout << "Rs为旋转矩阵" << std::endl; }
	return 0;
}

int ComputeTrans(std::vector<std::vector<double>> &poses_rob, std::vector<std::vector<double>> &poses_cam, Eigen::Vector3d &Ts, Eigen::Matrix3d &Rs) {
	//std::vector<std::vector<double>> poses_rob; // inputData
	//std::vector<std::vector<double>> poses_cam; // inputData

	Eigen::Matrix3d Ro;
	Eigen::Vector3d To, To_, Pc, Pc_;
	std::vector<Eigen::Matrix3d> Ros;
	std::vector<Eigen::Vector3d> Ds;
	if (poses_cam.size() != poses_rob.size()) {
		return -1;
	}
	int n = poses_rob.size();
	std::cout << "n: " << n << std::endl;
	for (int i = 1; i < n; i++) {
		Ro = GetTFromPose(poses_rob[i][3], poses_rob[i][4], poses_rob[i][5]) - GetTFromPose(poses_rob[i - 1][3], poses_rob[i - 1][4], poses_rob[i - 1][5]);
		Ros.push_back(Ro);
		Pc << poses_cam[i][0], poses_cam[i][1], poses_cam[i][2];
		To << poses_rob[i][0], poses_rob[i][1], poses_rob[i][2];
		Pc_ << poses_cam[i-1][0], poses_cam[i-1][1], poses_cam[i-1][2];
		To_ << poses_rob[i-1][0], poses_rob[i-1][1], poses_rob[i-1][2];
		Eigen::Vector3d dd = GetTFromPose(poses_rob[i][3], poses_rob[i][4], poses_rob[i][5]) * Rs * Pc - GetTFromPose(poses_rob[i - 1][3], poses_rob[i - 1][4], poses_rob[i - 1][5]) * Rs * Pc_ + To - To_;
		Ds.push_back(dd);
	}
	Eigen::MatrixXd C(15, 3), D(15, 1);
	C << Ros[0], Ros[1], Ros[2], Ros[3], Ros[4];
	D << Ds[0], Ds[1], Ds[2], Ds[3], Ds[4];
	
	Ts = C.colPivHouseholderQr().solve(D);
	std::cout << "Ts: " << Ts << std::endl;
	return 0;
}

void readAlgParas(const std::string &path, cv::Mat &trans_mat) {
    std::string filename = path + "/calibdata.yml";
    cv::FileStorage fs(filename, cv::FileStorage::READ);
    cv::FileNode n;
    if (!fs.isOpened()) {
        std::cerr << "failed to open " << filename << std::endl;
        return;
    }
    else {
        n = fs["translate matrix"];
        if (n.type() != cv::FileNode::NONE) {
			n["data"] >> trans_mat;
        }
    }
    fs.release();
}

void csvRead(const std::string& inputfile, std::vector<std::vector<double>>& _robot_points, const int _linecount) {
	const int MAX_LINE = 3000;
	FILE *fp = nullptr;
	char line[MAX_LINE];
	if ((fp = fopen(inputfile.c_str(), "r")) != nullptr)
	{
		char delims[] = ",";
		while (fgets(line, MAX_LINE, fp))
		{
			char* _string;
			std::vector<double> _line_point;
			_line_point.push_back(std::strtod(strtok(line, delims), &_string));
			std::cout << _line_point.back() << ',';
			for (size_t _i = 0; _i < _linecount-1; _i++)
			{
				_line_point.push_back(std::strtod(strtok(nullptr, delims), &_string));
				std::cout << _line_point.back() << ',';
			}
			_robot_points.push_back(_line_point);
			std::cout << std::endl;
		}
        std::fclose(fp);
	}
}

cv::Mat R_T2RT(Eigen::Matrix3d &Rs, Eigen::Vector3d &Ts) {
    cv::Mat RT;
    cv::Mat_<double> R1 = (cv::Mat_<double>(4, 3) << Rs(0, 0), Rs(0, 1), Rs(0, 2),
		Rs(1, 0), Rs(1, 1), Rs(1, 2),
		Rs(2, 0), Rs(2, 1), Rs(2, 2),
		0.0, 0.0, 0.0);
    cv::Mat_<double> T1 = (cv::Mat_<double>(4, 1) << Ts(0), Ts(1), Ts(2), 1.0);
    cv::hconcat(R1, T1, RT);//C=A+B左右拼接
	return RT;
}
cv::Mat EMat2RT(Eigen::Matrix4d &Rs) {
    cv::Mat_<double> RT = (cv::Mat_<double>(4, 4) << Rs(0, 0), Rs(0, 1), Rs(0, 2), Rs(0, 3),
		Rs(1, 0), Rs(1, 1), Rs(1, 2), Rs(1, 3),
		Rs(2, 0), Rs(2, 1), Rs(2, 2), Rs(2, 3),
		0.0, 0.0, 0.0, 1.0);
	return RT;
}
void outputEyeHandParams(std::string& yaml_file, cv::Mat& hcg)
{
    cv::FileStorage fs(yaml_file, cv::FileStorage::WRITE);

    if(!fs.isOpened())
    {
        std::cout << "open file error!" << std::endl;
        return;
    }

    fs << "hcg" << hcg;
    std::cout << "hcg: " << hcg << std::endl;
    fs.release();
}

int Test(std::vector<std::vector<double>> &robot, std::vector<std::vector<double>> &sphere, cv::Mat &hcg) {
    Eigen::Matrix3d Ro;
    cv::Mat RT_rob;
    cv::Mat Pose_sphere;
    std::cout << "----手眼标定测试----" <<std::endl;
    for (size_t i = 0; i < 9; i++) {
		Ro = GetTFromPose(robot[i][3], robot[i][4], robot[i][5]);
        RT_rob = (cv::Mat_<double>(4, 4) << Ro(0, 0), Ro(0, 1), Ro(0, 2), robot[i][0],
		Ro(1, 0), Ro(1, 1), Ro(1, 2), robot[i][1],
		Ro(2, 0), Ro(2, 1), Ro(2, 2), robot[i][2],
        0.0, 0.0, 0.0, 1.0);
        Pose_sphere = (cv::Mat_<double>(4, 1) << 
                sphere[i][0], sphere[i][1], sphere[i][2], 1.0);
       cv::Mat worldPos = RT_rob * hcg * Pose_sphere;
       std::cout << i << ": " << worldPos.t() << std::endl;
    }
    return 0;
}

Eigen::Matrix3d skew(Eigen::Vector3d u)
{
  Eigen::Matrix3d u_hat = Eigen::MatrixXd::Zero(3,3);
  u_hat(0,1) = u(2);
  u_hat(1,0) = -u(2);
  u_hat(0,2) = -u(1);
  u_hat(2,0) = u(1);
  u_hat(1,2) = u(0);
  u_hat(2,1) = -u(0);

  return u_hat;
}

Pose SolveX(Poses &AS, Poses &BS) {
	Poses A_ = AS;
	Poses B_ = BS;
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(12*A_.size(),12);
  for(int i=0;i<A_.size();i++) {
    //extract R,t from homogophy matrix
    Eigen::Matrix3d Ra = A_[i].topLeftCorner(3,3);
    Eigen::Vector3d Ta = A_[i].topRightCorner(3,1);
    Eigen::Matrix3d Rb = B_[i].topLeftCorner(3,3);
    Eigen::Vector3d Tb = B_[i].topRightCorner(3,1);

    m.block<9,9>(12*i,0) = Eigen::MatrixXd::Identity(9,9) - Eigen::kroneckerProduct(Ra,Rb);
    Eigen::Matrix3d Ta_skew = skew(Ta);
    m.block<3,9>(12*i+9,0) = Eigen::kroneckerProduct(Ta_skew,Tb.transpose());
    m.block<3,3>(12*i+9,9) = Ta_skew - Ta_skew*Ra;
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> svd( m, Eigen::ComputeFullV | Eigen::ComputeFullU );
  //CHECK(svd.computeV())<<"fail to compute V";

  Eigen::Matrix3d R_alpha;
  R_alpha.row(0) = svd.matrixV().block<3,1>(0,11).transpose();
  R_alpha.row(1) = svd.matrixV().block<3,1>(3,11).transpose();
  R_alpha.row(2) = svd.matrixV().block<3,1>(6,11).transpose();
  //double a = std::fabs(R_alpha.determinant());
  //double alpha = R_alpha.determinant()/(pow(std::fabs(R_alpha.determinant()),4./3.));
  double det = R_alpha.determinant();
  double alpha = std::pow(std::abs(det),4./3.)/det;
  Eigen::HouseholderQR<Eigen::Matrix3d> qr(R_alpha/alpha);

  Pose handeyetransformation = Pose::Identity(4,4);
  Eigen::Matrix3d Q = qr.householderQ();
  Eigen::Matrix3d Rwithscale = alpha*Q.transpose()*R_alpha;
  Eigen::Vector3d R_diagonal = Rwithscale.diagonal();
  for(int i=0;i<3;i++)
  {
    handeyetransformation.block<3,1>(0,i) = int(R_diagonal(i)>=0?1:-1)*Q.col(i);
  }

  handeyetransformation.topRightCorner(3,1) = svd.matrixV().block<3,1>(9,11)/alpha;
  return handeyetransformation;
}

Pose SolveX_(Poses &AS, Poses &BS) {
	Poses A_ = AS;
	Poses B_ = BS;
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(12*A_.size(),12);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(12*A_.size());
  for(int i=0;i<A_.size();i++)
  {
    //extract R,t from homogophy matrix
    Eigen::Matrix3d Ra = A_[i].topLeftCorner(3,3);
    Eigen::Vector3d Ta = A_[i].topRightCorner(3,1);
    Eigen::Matrix3d Rb = B_[i].topLeftCorner(3,3);
    Eigen::Vector3d Tb = B_[i].topRightCorner(3,1);

    m.block<9,9>(12*i,0) = Eigen::MatrixXd::Identity(9,9) - Eigen::kroneckerProduct(Ra,Rb);
    Eigen::Matrix3d Ta_skew = skew(Ta);
    m.block<3,9>(12*i+9,0) = Eigen::kroneckerProduct(Eigen::MatrixXd::Identity(3,3),Tb.transpose());
    m.block<3,3>(12*i+9,9) = Eigen::MatrixXd::Identity(3,3) - Ra;
    b.block<3,1>(12*i+9,0) = Ta;
  }

  Eigen::Matrix<double, 12, 1> x = m.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  Eigen::Matrix3d R = Eigen::Map< Eigen::Matrix<double, 3, 3, Eigen::RowMajor> >(x.data()); //row major

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(R, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Pose handeyetransformation = Pose::Identity(4,4);
  handeyetransformation.topLeftCorner(3,3) = svd.matrixU() * svd.matrixV().transpose();
  handeyetransformation.topRightCorner(3,1) = x.block<3,1>(9,0);
  return handeyetransformation;
}

int GetPose(std::vector<std::vector<double>> &robot_points, std::vector<std::vector<double>> &cam_points, Poses &AS, Poses &BS) {
	Pose robot, cam;
	Poses robots, cams;

	for (int i = 0; i < NUM_PointCloud; i++) {
		cam << 1, 0, 0, cam_points[i][0],
		       0, 1, 0, cam_points[i][1],
			   0, 0, 1, cam_points[i][2],
			   0, 0, 0, 1;
		std::cout << "camera_pose: " << std::endl;
		std::cout << cam << std::endl;
		cams.push_back(cam);
	}
	
	for (int i = 0; i < NUM_PointCloud; i++) {
		Eigen::Matrix3d Ro = GetTFromPose(robot_points[i][3], robot_points[i][4], robot_points[i][5]);
		// rotZYZ();
		robot << Ro(0, 0), Ro(0, 1), Ro(0, 2), robot_points[i][0],
		Ro(1, 0), Ro(1, 1), Ro(1, 2), robot_points[i][1], 
		Ro(2, 0), Ro(2, 1), Ro(2, 2), robot_points[i][2], 
		0, 0, 0, 1;
		std::cout << "robot pose: " << std::endl;
		std::cout << robot << std::endl;
		robots.push_back(robot);
	}
	for (int i = 0; i < NUM_PointCloud-2; i++) {
		Pose A, B;
		A = (robots[i+1]).inverse()*robots[i];
		B = cams[i+1]*(cams[i].inverse());
		AS.push_back(A);
		BS.push_back(B);
	}	

	return 0;
}
int GetPose_(std::vector<std::vector<double>> &robot_points, std::vector<std::vector<double>> &cam_points, Poses &robots, Poses &cams) {
	Pose robot, cam;

	for (int i = 0; i < NUM_PointCloud; i++) {
		cam << 1, 0, 0, cam_points[i][0],
		       0, 1, 0, cam_points[i][1],
			   0, 0, 1, cam_points[i][2],
			   0, 0, 0, 1;
		std::cout << "camera_pose: " << std::endl;
		std::cout << cam << std::endl;
		cams.push_back(cam);
	}
	
	for (int i = 0; i < NUM_PointCloud; i++) {
		Eigen::Matrix3d Ro = GetTFromPose(robot_points[i][3], robot_points[i][4], robot_points[i][5]);
		// rotZYZ();
		robot << Ro(0, 0), Ro(0, 1), Ro(0, 2), robot_points[i][0],
		Ro(1, 0), Ro(1, 1), Ro(1, 2), robot_points[i][1], 
		Ro(2, 0), Ro(2, 1), Ro(2, 2), robot_points[i][2], 
		0, 0, 0, 1;
		std::cout << "robot pose: " << std::endl;
		std::cout << robot << std::endl;
		robots.push_back(robot);
	}

	return 0;
}
cv::Mat Solve_LS(Poses &robot_mats, Poses &cam_mats) {
	if (robot_mats.size() != cam_mats.size() ) {
		printf("input data error!");
		//return ;
	}
	std::vector<Eigen::MatrixXd> Avs;
	std::vector<Eigen::MatrixXd> Bvs;
	for (size_t i = 0; i < robot_mats.size()-1; i++) {
		Eigen::MatrixXd temp(3, 12);
		Eigen::MatrixXd temp_b(3, 1);
		temp << (robot_mats[i](0,0)*cam_mats[i](0,3)-robot_mats[i+1](0,0)*cam_mats[i+1](0,3)), 
				(robot_mats[i](0,0)*cam_mats[i](1,3)-robot_mats[i+1](0,0)*cam_mats[i+1](1,3)),
				(robot_mats[i](0,0)*cam_mats[i](2,3)-robot_mats[i+1](0,0)*cam_mats[i+1](2,3)),
				(robot_mats[i](0,1)*cam_mats[i](0,3)-robot_mats[i+1](0,1)*cam_mats[i+1](0,3)),
				(robot_mats[i](0,1)*cam_mats[i](1,3)-robot_mats[i+1](0,1)*cam_mats[i+1](1,3)),
				(robot_mats[i](0,1)*cam_mats[i](2,3)-robot_mats[i+1](0,1)*cam_mats[i+1](2,3)),
				(robot_mats[i](0,2)*cam_mats[i](0,3)-robot_mats[i+1](0,2)*cam_mats[i+1](0,3)),
				(robot_mats[i](0,2)*cam_mats[i](1,3)-robot_mats[i+1](0,2)*cam_mats[i+1](1,3)),
				(robot_mats[i](0,2)*cam_mats[i](2,3)-robot_mats[i+1](0,2)*cam_mats[i+1](2,3)),
				(robot_mats[i](0,0)-robot_mats[i+1](0,0)),
				(robot_mats[i](0,1)-robot_mats[i+1](0,1)),
				(robot_mats[i](0,2)-robot_mats[i+1](0,2)),
				(robot_mats[i](1,0)*cam_mats[i](0,3)-robot_mats[i+1](1,0)*cam_mats[i+1](0,3)), 
				(robot_mats[i](1,0)*cam_mats[i](1,3)-robot_mats[i+1](1,0)*cam_mats[i+1](1,3)),
				(robot_mats[i](1,0)*cam_mats[i](2,3)-robot_mats[i+1](1,0)*cam_mats[i+1](2,3)),
				(robot_mats[i](1,1)*cam_mats[i](0,3)-robot_mats[i+1](1,1)*cam_mats[i+1](0,3)),
				(robot_mats[i](1,1)*cam_mats[i](1,3)-robot_mats[i+1](1,1)*cam_mats[i+1](1,3)),
				(robot_mats[i](1,1)*cam_mats[i](2,3)-robot_mats[i+1](1,1)*cam_mats[i+1](2,3)),
				(robot_mats[i](1,2)*cam_mats[i](0,3)-robot_mats[i+1](1,2)*cam_mats[i+1](0,3)),
				(robot_mats[i](1,2)*cam_mats[i](1,3)-robot_mats[i+1](1,2)*cam_mats[i+1](1,3)),
				(robot_mats[i](1,2)*cam_mats[i](2,3)-robot_mats[i+1](1,2)*cam_mats[i+1](2,3)),
				(robot_mats[i](1,0)-robot_mats[i+1](1,0)),
				(robot_mats[i](1,1)-robot_mats[i+1](1,1)),
				(robot_mats[i](1,2)-robot_mats[i+1](1,2)),
				(robot_mats[i](2,0)*cam_mats[i](0,3)-robot_mats[i+1](2,0)*cam_mats[i+1](0,3)), 
				(robot_mats[i](2,0)*cam_mats[i](1,3)-robot_mats[i+1](2,0)*cam_mats[i+1](1,3)),
				(robot_mats[i](2,0)*cam_mats[i](2,3)-robot_mats[i+1](2,0)*cam_mats[i+1](2,3)),
				(robot_mats[i](2,1)*cam_mats[i](0,3)-robot_mats[i+1](2,1)*cam_mats[i+1](0,3)),
				(robot_mats[i](2,1)*cam_mats[i](1,3)-robot_mats[i+1](2,1)*cam_mats[i+1](1,3)),
				(robot_mats[i](2,1)*cam_mats[i](2,3)-robot_mats[i+1](2,1)*cam_mats[i+1](2,3)),
				(robot_mats[i](2,2)*cam_mats[i](0,3)-robot_mats[i+1](2,2)*cam_mats[i+1](0,3)),
				(robot_mats[i](2,2)*cam_mats[i](1,3)-robot_mats[i+1](2,2)*cam_mats[i+1](1,3)),
				(robot_mats[i](2,2)*cam_mats[i](2,3)-robot_mats[i+1](2,2)*cam_mats[i+1](2,3)),
				(robot_mats[i](2,0)-robot_mats[i+1](2,0)),
				(robot_mats[i](2,1)-robot_mats[i+1](2,1)),
				(robot_mats[i](2,2)-robot_mats[i+1](2,2));
		temp_b << (robot_mats[i+1](0,3) - robot_mats[i](0,3)),
				  (robot_mats[i+1](1,3) - robot_mats[i](1,3)),
				  (robot_mats[i+1](2,3) - robot_mats[i](2,3));
		Avs.push_back(temp);
		Bvs.push_back(temp_b);
	}
	
	int num = Avs.size();
	Eigen::MatrixXd A(num*3, 12), B(num*3, 1);
	for (int i = 0; i < num; i++) {
		A.block((i*3),0,3,12) = Avs[i];
		B.block((i*3),0, 3, 1) = Bvs[i];
	} 
	
	std::cout << "A: " << A << "\n" << " B: " << B << std::endl;
	Eigen::MatrixXd H;
	H = A.householderQr().solve(B);
	std::cout << "H_qr:" << H.transpose() << std::endl;
	//H = A.ldlt().solve(B);
	//std::cout << "H_ldlt:" << H.transpose() << std::endl;
	
	cv::Mat_<double> RT = (cv::Mat_<double>(4, 4) << H(0, 0), H(1, 0), H(2, 0), H(9, 0),
		H(3, 0), H(4, 0), H(5, 0), H(10, 0),
		H(6, 0), H(7, 0), H(8, 0), H(11, 0),
		0.0, 0.0, 0.0, 1.0);
	std::cout << RT << std::endl;
	return RT;
}