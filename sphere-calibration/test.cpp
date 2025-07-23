#include "spherecalibrate.hpp"


int main() {
    // std::string robot_r_file = "robot_r_.csv";
    // std::string cam_r_file = "cam_r_.csv";
    // std::string robot_t_file = "robot_t_.csv";
    // std::string cam_t_file = "cam_t_.csv";
    // std::vector<std::vector<double>> robot_r_points, cam_r_points, robot_t_points, cam_t_points;
    // std::cout << "rotate cam data:" << std::endl;
    // csvRead(cam_r_file, cam_r_points, 3);
    // std::cout << "rotate robot data:" << std::endl;
    // csvRead(robot_r_file, robot_r_points, 6);
    // std::cout << "translate robot data:" << std::endl;
    // csvRead(robot_t_file, robot_t_points, 6);
    // std::cout << "translate cam data:" << std::endl;
    // csvRead(cam_t_file, cam_t_points, 3);
    
    // Eigen::Matrix3d Rs;
    // ComputeRotMat(robot_r_points, cam_r_points, Rs);
    
    // Eigen::Vector3d Ts;
    // ComputeTrans(robot_t_points, cam_t_points, Ts, Rs);
    
    // std::string yaml_file = "calib.yaml"; 
    // cv::Mat hcg = R_T2RT(Rs, Ts);
    // std::cout << "hcg是否为旋转矩阵: " << isRotationMatrix(hcg) << std::endl;
    // outputEyeHandParams(yaml_file, hcg);
    
    // Test(robot_t_points, cam_t_points, hcg);
    
    std::string cam_com_file = "cam_com.csv";
    std::string robot_com_file = "robot_com.csv";
    std::vector<std::vector<double>> cams, robots;
    Poses A_rs, B_cs;
    csvRead(cam_com_file, cams, 3);
    csvRead(robot_com_file, robots, 6);
    //GetPose(robots, cams, A_rs, B_cs);
    //Eigen::Matrix4d hcg_ = SolveX(A_rs, B_cs);
    //Eigen::Matrix4d hcg_test = SolveX_(A_rs, B_cs);
    
    Poses robot_mats, cam_mats;
    GetPose_(robots, cams, robot_mats, cam_mats);
    cv::Mat ls_mat = Solve_LS(robot_mats, cam_mats);
    //cv::Mat hcg_cv = EMat2RT(hcg_);
    //cv::Mat hcg_cv_test = EMat2RT(hcg_test);
    std::cout << "hcg: " << ls_mat << std::endl;
    //std::cout << "hcg: " << hcg_cv_test << std::endl;
    //std::cout << "hcg_cv是否为旋转矩阵: " << isRotationMatrix(hcg_cv) << std::endl;
    //std::cout << "hcg_cv_test是否为旋转矩阵: " << isRotationMatrix(hcg_cv_test) << std::endl;
    Test(robots, cams, ls_mat);
    //Test(robots, cams, hcg_cv);
    //Test(robots, cams, hcg_cv_test);
    return 0;
}