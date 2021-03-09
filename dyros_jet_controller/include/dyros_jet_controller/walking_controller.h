#ifndef WALKING_CONTROLLER_H
#define WALKING_CONTROLLER_H


#include "dyros_jet_controller/dyros_jet_model.h"
#include "math_type_define.h"
#include <vector>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <thread>
#include <mutex>

#define ZERO_LIBRARY_MODE


const int FILE_CNT = 15;

const std::string FILE_NAMES[FILE_CNT] =
{
  ///change this directory when you use this code on the other computer///
 

};

using namespace std;
namespace dyros_jet_controller
{

class WalkingController
{
public:
  fstream file[FILE_CNT];


  static constexpr unsigned int PRIORITY = 8;


  WalkingController(DyrosJetModel& model, const VectorQd& current_q, const VectorQd& current_qdot, const double hz, const double& control_time) :
    total_dof_(DyrosJetModel::HW_TOTAL_DOF), model_(model), current_q_(current_q), current_qdot_(current_qdot), hz_(hz), current_time_(control_time), start_time_{}, end_time_{}, slowcalc_thread_(&WalkingController::slowCalc, this), calc_update_flag_(false), calc_start_flag_(false), ready_for_thread_flag_(false), ready_for_compute_flag_(false), foot_step_planner_mode_(false), walking_end_foot_side_ (false), foot_plan_walking_last_(false), foot_last_walking_end_(false)
  {
    walking_state_send = false;
    walking_end_ = false;

  } 

  void compute();
  void setTarget(int walk_mode, bool hip_compensation, bool lqr, int ik_mode, bool heel_toe,
                 bool is_right_foot_swing, double x, double y, double z, double height, double theta,
                 double step_length, double step_length_y, bool walking_pattern);
  void setEnable(bool enable);
  void setFootPlan(int footnum, int startfoot, Eigen::MatrixXd footpose);
  void updateControlMask(unsigned int *mask);
  void writeDesired(const unsigned int *mask, VectorQd& desired_q);

  void parameterSetting();
  //functions in compute
  //
  double CP_X_ref = 0, CP_Y_ref = 0;
  double CP_X_act = 0, CP_Y_act = 0;
  double del_zmp_x = 0, del_zmp_y = 0;
  double del_ux_LPF = 0, del_uy_LPF = 0;
  void getCPTrajectory_M();
  Eigen::Vector2d p;
  Eigen::Vector2d cp_eos; 
  Eigen::Vector2d cp;
  Eigen::Vector2d cp_pre;
  Eigen::Vector2d com;
  Eigen::Vector2d com_dot;
  Eigen::Vector2d com_pre;
  Eigen::Vector2d dT, b;
  Eigen::Vector2d cp_measured_;
  Eigen::Vector2d p_d;
  Eigen::Vector2d p_d_pre;
  double zmp_x = 0;
  double zmp_y = 0; 
  double del_t = 0.005;
  //

  Eigen::Vector2d zmp_measured_;
  Eigen::Vector2d zmp_measured_LPF_;
  Eigen::Vector2d zmp_err_;

  void getRobotState();
  void getComTrajectory();
  void getZmpTrajectory();
  void getPelvTrajectory();
  void getFootTrajectory();
  void computeIkControl_MJ(Eigen::Isometry3d float_trunk_transform, Eigen::Isometry3d float_lleg_transform, Eigen::Isometry3d float_rleg_transform, Eigen::Vector12d& desired_leg_q);
  void compensator();

  void supportToFloatPattern();
  void updateNextStepTime();
  void updateInitialState(); 
  //functions for getFootStep()
  void floatToSupportFootstep();
  void addZmpOffset();
  void zmpGenerator(const unsigned int norm_size, const unsigned planning_step_num);
  void comGenerator(const unsigned int norm_size, const unsigned planning_step_num);
  void onestepZmp(unsigned int current_step_number, Eigen::VectorXd& temp_px, Eigen::VectorXd& temp_py);
  void onestepCom(unsigned int current_step_number, Eigen::VectorXd& temp_px, Eigen::VectorXd& temp_py); 
  void OfflineCoM_MJ(unsigned int current_step_number, Eigen::VectorXd& temp_cx, Eigen::VectorXd& temp_cy);
  void calculateFootStepSeparate();
  void calculateFootStepTotal();
  void calculateFootStepTotal_MJ();
  void usingFootStepPlanner();
  void hip_compensator(); // HW6
  void Compliant_control(Eigen::Vector12d desired_leg_q);

  //PreviewController    
  void previewcontroller(double dt, int NL, int tick, double x_i, double y_i, Eigen::Vector3d xs, Eigen::Vector3d ys, double& UX, double& UY, 
       Eigen::MatrixXd Gi, Eigen::VectorXd Gd, Eigen::MatrixXd Gx, Eigen::MatrixXd A, Eigen::VectorXd B, Eigen::MatrixXd C, Eigen::Vector3d &XD, Eigen::Vector3d &YD);  
  void preview_Parameter(double dt, int NL, Eigen::MatrixXd& Gi, Eigen::VectorXd& Gd, Eigen::MatrixXd& Gx, Eigen::MatrixXd& A, Eigen::VectorXd& B, Eigen::MatrixXd& C);
  void preview_Parameter_CPM(double dt, int NL, Eigen::Matrix3d& K, Eigen::Vector3d com_support_init_, Eigen::MatrixXd& Gi, Eigen::VectorXd& Gd, Eigen::MatrixXd& Gx, 
  Eigen::MatrixXd& A, Eigen::VectorXd& B, Eigen::MatrixXd& C, Eigen::MatrixXd& D, Eigen::MatrixXd& A_bar, Eigen::VectorXd& B_bar);
  void previewcontroller_CPM(double dt, int NL, int tick, double x_i, double y_i, Eigen::Vector3d xs, Eigen::Vector3d ys, double& UX, double& UY, 
       Eigen::MatrixXd Gi, Eigen::VectorXd Gd, Eigen::MatrixXd Gx, Eigen::MatrixXd A, Eigen::VectorXd B, Eigen::MatrixXd A_bar, Eigen::VectorXd B_bar, Eigen::Vector2d &XD, Eigen::Vector2d &YD, Eigen::VectorXd& X_bar_p, Eigen::VectorXd& Y_bar_p);
  void slowCalc();  

  VectorQd desired_q_not_compensated_;

  bool walking_end_foot_side_;
  bool walking_end_;
  bool foot_plan_walking_last_;
  bool foot_last_walking_end_;
  bool walking_state_send;

  //CapturePoint
  void getCapturePointTrajectory();
  void getCapturePoint_init_ref();
  void CapturePointModify();
  void zmptoInitFloat();
  Eigen::VectorXd capturePoint_refx, capturePoint_refy;
  Eigen::VectorXd zmp_refx, zmp_refy;
  Eigen::VectorXd capturePoint_ox, capturePoint_oy, zmp_dx, zmp_dy;
  Eigen::Vector3d capturePoint_measured_;
  double last_time_;
  int capturePoint_current_num_;
  Eigen::Vector3d com_float_prev_;
  Eigen::Vector4d com_float_prev_dot_;
  Eigen::Vector4d com_float_prev;
  Eigen::Vector3d com_support_prev;
  double ux_1, uy_1;
  Eigen::Vector3d xs, ys;
  int currentstep;
  bool firsttime = false;
  Eigen::Vector2d capturePoint_offset_;
  Eigen::Isometry3d float_support_init;
  Eigen::Isometry3d current_step_float_support_;

  Eigen::Isometry3d support_float_init;
  Eigen::Isometry3d current_step_support_float_;

  Eigen::Vector6d q_sim_virtual_;
  Eigen::Vector6d q_sim_dot_virtual_;
  Eigen::VectorXd com_refx;
  Eigen::VectorXd com_refy;
  Eigen::VectorXd com_dot_refx;
  Eigen::VectorXd com_dot_refy;
  Eigen::Vector2d com_initx;
  Eigen::Vector2d com_inity;


private:

  const double hz_;
  const double &current_time_; // updated by control_base
  unsigned int walking_tick_ = 0;
  unsigned int com_tick_ = 0;
  double walking_time_ = 0;
  unsigned int tick_d1 = 0, tick_d2 = 0;
  //sensorData
  Eigen::Vector6d r_ft_;
  Eigen::Vector6d l_ft_;
  Eigen::Vector3d imu_acc_;
  Eigen::Vector3d imu_ang_;
  Eigen::Vector3d imu_grav_rpy_;

  double total_mass = 0;

  //parameterSetting()
  double t_last_;
  double t_start_;
  double t_start_real_;
  double t_temp_;
  double t_imp_;
  double t_rest_init_;
  double t_rest_last_;
  double t_double1_;
  double t_double2_;
  double t_total_;
  double foot_height_;
  double com_height_;

  bool com_control_mode_;
  bool com_update_flag_; // frome A to B
  bool gyro_frame_flag_;
  bool ready_for_thread_flag_;
  bool ready_for_compute_flag_;
  bool estimator_flag_;

  int ik_mode_;
  int walk_mode_;
  bool hip_compensator_mode_;
  bool lqr_compensator_mode_; 
  int is_right_foot_swing_;
  bool foot_step_planner_mode_;

  bool walking_enable_;
  bool joint_enable_[DyrosJetModel::HW_TOTAL_DOF];
  double step_length_x_;
  double step_length_y_;

  //double step_angle_theta_;
  unsigned int print_flag = 0;
  unsigned int print_flag_1 = 0;
  double target_x_;
  double target_y_;
  double target_z_;
  double target_theta_;
  double total_step_num_;
  double current_step_num_;
  int foot_step_plan_num_;
  int foot_step_start_foot_;
  bool walkingPatternDCM_;
  Eigen::MatrixXd foot_pose_;

  Eigen::MatrixXd foot_step_;
  Eigen::MatrixXd foot_step_support_frame_;
  Eigen::MatrixXd foot_step_support_frame_offset_;

  Eigen::MatrixXd org_ref_zmp_;
  Eigen::MatrixXd ref_zmp_;
  Eigen::MatrixXd modified_ref_zmp_;
  Eigen::MatrixXd ref_com_;
  Eigen::MatrixXd ref_zmp_float_;

  VectorQd start_q_;
  VectorQd desired_q_;
  VectorQd target_q_;
  const VectorQd& current_q_;
  const VectorQd& current_qdot_;

  double prev_zmp_error_y = 0, prev_zmp_error_x = 0;


  //const double &current_time_;
  const unsigned int total_dof_;
  double start_time_[DyrosJetModel::HW_TOTAL_DOF];
  double end_time_[DyrosJetModel::HW_TOTAL_DOF];

  //Step initial state variable//
  Eigen::Isometry3d pelv_support_init_2;
  Eigen::Isometry3d pelv_support_init_;
  Eigen::Isometry3d lfoot_support_init_;
  Eigen::Isometry3d rfoot_support_init_;
  Eigen::Isometry3d pelv_float_init_;
  Eigen::Isometry3d lfoot_float_init_;
  Eigen::Isometry3d rfoot_float_init_;

  Eigen::Vector3d pelv_support_euler_init_;
  Eigen::Vector3d lfoot_support_euler_init_;
  Eigen::Vector3d rfoot_support_euler_init_;
  VectorQd q_init_;

  Eigen::Vector6d supportfoot_float_init_;
  Eigen::Vector6d supportfoot_support_init_;
  Eigen::Vector6d supportfoot_support_init_offset_;
  Eigen::Vector6d swingfoot_float_init_;
  Eigen::Vector6d swingfoot_support_init_;
  Eigen::Vector6d swingfoot_support_init_offset_;

  Eigen::Isometry3d pelv_support_start_;

  Eigen::Vector3d com_float_init_;
  Eigen::Vector3d com_support_init_;
  Eigen::Vector3d com_support_init_2;
  double lfoot_zmp_offset_;   //have to be initialized
  double rfoot_zmp_offset_;
  Eigen::Vector3d com_offset_;

  //Step current state variable//
  Eigen::Vector3d com_support_current_Euler;
  Eigen::Vector3d com_support_current_CLIPM_b;
  Eigen::Vector3d com_support_current_CLIPM;
  Eigen::Vector3d com_support_current_b;
  Eigen::Vector3d com_support_current_dot;
  Eigen::Vector3d com_support_current_; 
  Eigen::Vector3d com_middle_support_current_;
  Eigen::Vector3d com_support_dot_current_;//from support foot
  Eigen::Vector3d com_support_ddot_current_;//from support foot

  Eigen::Vector3d com_support_current_dot_b; 
  Eigen::Vector3d com_support_current_dot_LPF;
  
  double R_angle = 0, P_angle = 0;
  double R_angle_i = 0, P_angle_i = 0;
  double leg_length_L = 0, leg_length_R = 0;
  ///simulation
  Eigen::Vector3d com_sim_current_;
  Eigen::Vector3d com_sim_dot_current_;
  Eigen::Isometry3d lfoot_sim_global_current_;
  Eigen::Isometry3d rfoot_sim_global_current_;
  Eigen::Isometry3d base_sim_global_current_;
  Eigen::Isometry3d lfoot_sim_float_current_;
  Eigen::Isometry3d rfoot_sim_float_current_;
  Eigen::Isometry3d supportfoot_float_sim_current_;

  Eigen::Vector3d gyro_sim_current_;
  Eigen::Vector3d accel_sim_current_;
  
  Eigen::Isometry3d supportfoot_float_current_Euler;
  Eigen::Isometry3d supportfoot_float_current_;
  Eigen::Isometry3d pelv_support_current_Euler;
  Eigen::Isometry3d pelv_support_current_;
  Eigen::Isometry3d lfoot_support_current_;
  Eigen::Isometry3d rfoot_support_current_;
  Eigen::Isometry3d lfoot_support_current_ZMP;
  Eigen::Isometry3d rfoot_support_current_ZMP;

  Eigen::Vector3d com_float_current_;
  Eigen::Vector3d com_float_current_RPY;
  Eigen::Vector3d com_global_current_;
  Eigen::Vector3d com_float_current_dot_;
  Eigen::Isometry3d pelv_float_current_;
  Eigen::Isometry3d lfoot_float_current_;
  Eigen::Isometry3d rfoot_float_current_;
  Eigen::Isometry3d lfoot_float_current_Euler;
  Eigen::Isometry3d rfoot_float_current_Euler;
  Eigen::Isometry3d R_;
  Eigen::Matrix3d R;
  
  Eigen::Matrix6d current_leg_jacobian_l_;
  Eigen::Matrix6d current_leg_jacobian_r_;
  DyrosJetModel &model_;

  Eigen::Vector3d com_desired_float_;
  double final_ref_zmp_print = 0 ;
  double final_com_print = 0 ;
  //desired variables
  Eigen::Vector12d q_des;
  Eigen::Vector12d desired_leg_q_;
  Eigen::Vector12d desired_leg_q_dot_;
  Eigen::Vector3d com_desired_;
  Eigen::Vector3d com_dot_desired_;
  Eigen::Vector2d zmp_desired_;
  // 수업용
  Eigen::Vector3d com_desired_dot_, com_desired_ddot_ , com_desired_dot_b_;
  //
  Eigen::Isometry3d rfoot_trajectory_support_;  //local frame
  Eigen::Isometry3d lfoot_trajectory_support_;
  Eigen::Vector3d rfoot_trajectory_euler_support_;
  Eigen::Vector3d lfoot_trajectory_euler_support_;
  Eigen::Vector6d rfoot_trajectory_dot_support_; //x,y,z translation velocity & roll, pitch, yaw velocity
  Eigen::Vector6d lfoot_trajectory_dot_support_;
 
  Eigen::Isometry3d pelv_trajectory_support_; //local frame
  Eigen::Isometry3d pelv_trajectory_float_; //pelvis frame
 //
  Eigen::Isometry3d rfoot_trajectory_float_;  //pelvis frame
  Eigen::Isometry3d lfoot_trajectory_float_;
  Eigen::Vector3d rfoot_trajectory_euler_float_;
  Eigen::Vector3d lfoot_trajectory_euler_float_;
  Eigen::Vector3d rfoot_trajectory_dot_float_;
  Eigen::Vector3d lfoot_trajectory_dot_float_;

  //getComTrajectory() variables
  double xi_;
  double yi_;
  Eigen::Vector3d xs_;
  Eigen::Vector3d ys_;
  Eigen::Vector3d xd_;
  Eigen::Vector3d yd_; 

  //Preview Control
  Eigen::Vector3d preview_x, preview_y, preview_x_b, preview_y_b, preview_x_b2, preview_y_b2;
  double ux_, uy_, ux_1_, uy_1_;
  double zc_;
  double gi_;
  double zmp_start_time_; //원래 코드에서는 start_time, zmp_ref 시작되는 time같음
  Eigen::Matrix4d k_;
  Eigen::Matrix4d K_act_;
  Eigen::VectorXd gp_l_;
  Eigen::Matrix1x3d gx_;
  Eigen::Matrix3d a_;
  Eigen::Vector3d b_;
  Eigen::Matrix1x3d c_;

  //Preview CLIPM MJ
  Eigen::MatrixXd A_;
  Eigen::VectorXd B_;
  Eigen::MatrixXd C_;
  Eigen::MatrixXd D_;
  Eigen::Matrix3d K_;
  Eigen::MatrixXd Gi_;
  Eigen::MatrixXd Gx_;
  Eigen::VectorXd Gd_;
  Eigen::MatrixXd A_bar_;
  Eigen::VectorXd B_bar_;
  Eigen::Vector2d Preview_X, Preview_Y, Preview_X_b, Preview_Y_b;
  Eigen::VectorXd X_bar_p_, Y_bar_p_;
  Eigen::Vector2d XD_;
  Eigen::Vector2d YD_;
  double UX_, UY_; 

  //resolved momentum control
  Eigen::Vector3d p_ref_;
  Eigen::Vector3d l_ref_;

  //Gravitycompensate
  Eigen::Vector12d joint_offset_angle_;
  Eigen::Vector12d grav_ground_torque_;

  //vibrationCotrol
  std::mutex slowcalc_mutex_;
  std::thread slowcalc_thread_;

  Eigen::Vector12d current_motor_q_leg_;
  Eigen::Vector12d current_link_q_leg_;
  Eigen::Vector12d pre_motor_q_leg_;
  Eigen::Vector12d pre_link_q_leg_;
  Eigen::Vector12d lqr_output_;
  Eigen::Vector12d lqr_output_pre_;
  Eigen::Vector12d DOB_IK_output_;
  Eigen::Vector12d DOB_IK_output_b_;

  VectorQd thread_q_;
  unsigned int thread_tick_;
 
  Eigen::Matrix<double, 12, 48> kkk_copy_;
  Eigen::Matrix<double, 48, 48> ad_total_copy_;
  Eigen::Matrix<double, 48, 12> bd_total_copy_;
  Eigen::Matrix<double, 36, 36> ad_copy_;
  Eigen::Matrix<double, 36, 12> bd_copy_;

  Eigen::Matrix<double, 36, 36> ad_right_;

  Eigen::Matrix<double, 36, 12> bd_right_;
  Eigen::Matrix<double, 48, 48> ad_total_right_;
  Eigen::Matrix<double, 48, 12> bd_total_right_;
  Eigen::Matrix<double, 12, 48> kkk_motor_right_;

  Eigen::Vector12d d_hat_b;

  bool calc_update_flag_;
  bool calc_start_flag_;


  Eigen::Matrix<double, 18, 18> mass_matrix_;
  Eigen::Matrix<double, 18, 18> mass_matrix_pc_;
  Eigen::Matrix<double, 12, 12> mass_matrix_sel_;
  Eigen::Matrix<double, 36, 36> a_right_mat_;
  Eigen::Matrix<double, 36, 12> b_right_mat_;
  Eigen::Matrix<double, 12, 36> c_right_mat_;
  Eigen::Matrix<double, 36, 36> a_disc_;
  Eigen::Matrix<double, 36, 12> b_disc_;
  Eigen::Matrix<double, 48, 48> a_disc_total_;
  Eigen::Matrix<double, 48, 12> b_disc_total_;
  Eigen::Matrix<double, 48, 48> kkk_;



  //HW 1 ///////////////////////////////////////////////
  void circling_motion();
  //////////////////////////////////////////////////////
  //HW 2 ///////////////////////////////////////////////

  //////////////////////////////////////////////////////
};

}
#endif // WALKING_CONTROLLER_H