#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/walking_controller.h"
#include "cvxgen_6_8_0/cvxgen/solver.h"
#include <fstream>

Vars vars;
Params params;
Workspace work;
Settings settings;

namespace dyros_jet_controller
{

ofstream Zmp_tra_graph("/home/myeongju/Zmp_tra_graph.txt");
ofstream Com_tra_graph("/home/myeongju/Com_tra_graph.txt");
ofstream Swing_tra_graph("/home/myeongju/Swing_tra_graph.txt");
ofstream Des_q_tra_graph("/home/myeongju/Des_q_tra_graph.txt");
ofstream LQR_q_tra_graph("/home/myeongju/LQR_q_tra_graph.txt");

void WalkingController::compute()
{
  if((walking_enable_ == true))
  {

    updateInitialState(); // Step change 하기 1tick 
    getRobotState(); // 매 tick 돈다.
    floatToSupportFootstep();
  
    if(ready_for_thread_flag_ == false)
    {
      ready_for_thread_flag_ = true;
    }

    if(ready_for_compute_flag_ == true || lqr_compensator_mode_ == false)
    {

      if(current_step_num_< total_step_num_)
      {
        getZmpTrajectory();
        getComTrajectory_MJ();
        getPelvTrajectory();
        getFootTrajectory_MJ();
        supportToFloatPattern();       
          
        /////////////////////compute/////////////////////////
        if (ik_mode_ == 0)
        {
          //computeIkControl(pelv_trajectory_float_, lfoot_trajectory_float_, rfoot_trajectory_float_, desired_leg_q_);
          computeIkControl_MJ(pelv_trajectory_float_, lfoot_trajectory_float_, rfoot_trajectory_float_, q_des);
          
          for(int i=0; i<12; i++)
          {            
            desired_q_(i) = q_des(i);
          }
        }
        else if (ik_mode_ == 1)
        {
          computeJacobianControl(lfoot_trajectory_float_, rfoot_trajectory_float_, lfoot_trajectory_euler_float_, rfoot_trajectory_euler_float_, desired_leg_q_dot_);
          for(int i=0; i<12; i++)
          {
            if(walking_tick_ == 0)
            {
              desired_q_(i) = q_init_(i);
            }
            desired_q_(i) = desired_leg_q_dot_(i)/hz_ + current_q_(i);
          }
        }
        //////////////////////////////////////////////////////

        desired_q_not_compensated_ = desired_q_ ; // IK 풀어서 나온 Desired Joint angle

        compensator();
        updateNextStepTime();
      }
      else
      {
        desired_q_ = current_q_;
      }
    }
    else
    {
      desired_q_ = current_q_;
    }
  }
}

void WalkingController::setTarget(int walk_mode, bool hip_compensation, bool lqr, int ik_mode, bool heel_toe,
                                  bool is_right_foot_swing, double x, double y, double z, double height, double theta,
                                  double step_length_x, double step_length_y, bool walking_pattern)
{
  target_x_ = x;
  target_y_ = y;
  target_z_ = z;
  com_height_ = height;
  target_theta_ = theta;
  step_length_x_ = step_length_x;
  step_length_y_ = step_length_y;
  ik_mode_ = ik_mode;
  walk_mode_ = walk_mode;
  hip_compensator_mode_ = hip_compensation; //uint32 compensator_mode[0] : HIP_COMPENSTOR    uint32 compensator_mode[1] : EXTERNAL_ENCODER
  lqr_compensator_mode_ = lqr;
  
  heel_toe_mode_ = heel_toe;
  is_right_foot_swing_ = is_right_foot_swing;
  
  walkingPatternDCM_ = walking_pattern;
  calc_start_flag_ = lqr;

  parameterSetting();
}

void WalkingController::setEnable(bool enable)
{
  walking_enable_ = enable;
  desired_q_ = current_q_;
}

void WalkingController::updateControlMask(unsigned int *mask)
{
  if(walking_enable_)
  {
    for (int i=0; i<total_dof_-18; i++) //control only leg
    {
      mask[i] = (mask[i] | PRIORITY);
    }
    mask[total_dof_-1] = (mask[total_dof_-1] & ~PRIORITY); //Gripper
    mask[total_dof_-2] = (mask[total_dof_-2] & ~PRIORITY); //Gripper
    mask[total_dof_-3] = (mask[total_dof_-3] & ~PRIORITY); //Head
    mask[total_dof_-4] = (mask[total_dof_-4] & ~PRIORITY); //Head
  }
  else
  {
    for (int i=0; i<total_dof_; i++)
    {
      mask[i] = (mask[i] & ~PRIORITY);
    }
  }
}

void WalkingController::writeDesired(const unsigned int *mask, VectorQd& desired_q)
{
  for(unsigned int i=0; i<total_dof_; i++)
  { 
    
    if(lqr_compensator_mode_ == true) //MJ DOB
    {
      if( mask[i] >= PRIORITY && mask[i] < PRIORITY * 2 )
      { 
        if(walking_tick_ == 0)
        {
          desired_q(i) = desired_q_(i);
        }
        else
        {
          desired_q(i) = DOB_IK_output_b_(i);
          desired_q(12) = desired_q_(12); // 허리
          desired_q(13) = desired_q_(13); // 허리
        }
      }
    }
    else
    {
      if( mask[i] >= PRIORITY && mask[i] < PRIORITY * 2 )
      {
        desired_q(i) = desired_q_(i);
      }
    }     
  }
}

void WalkingController::parameterSetting()
{
  t_double1_ = 0.10*hz_;
  t_double2_ = 0.10*hz_;
  t_rest_init_ = 0.05*hz_;
  t_rest_last_ = 0.05*hz_;
  t_total_= 1.2*hz_;
  t_temp_ = 3.0*hz_;
  t_last_ = t_total_ + t_temp_; //4.2*hz_;
  t_start_ = t_temp_ + 1;

  //initialize (should revise)
  t_start_real_ = t_start_ + t_rest_init_;

  current_step_num_ = 0;
  walking_tick_ = 0;
  walking_time_ = 0;

  foot_height_ = 0.05;
  //com_update_flag_ = true; // frome A to B1
  gyro_frame_flag_ = false;
  com_control_mode_ = true;
  estimator_flag_ = false; 

}

void WalkingController::getRobotState()
{
  Eigen::Matrix<double, DyrosJetModel::MODEL_WITH_VIRTUAL_DOF, 1> q_temp, qdot_temp;
  q_temp.setZero();
  qdot_temp;
  q_temp.segment<28>(6) = current_q_.segment<28>(0); //segment<포함갯수>(시작점) 
  qdot_temp.segment<28>(6)= current_qdot_.segment<28>(0);
  
 /*if(walking_tick_ > 0) // Using desired joint angle for RBDL update
  {
    q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0); //segment<포함갯수>(시작점), desired_q_not_compensated_ 는 IK에서 구해지는 Desired Joint angle.     
  }
  */
  model_.updateKinematics(q_temp, qdot_temp);
  
  com_float_current_ = model_.getCurrentCom(); 
  com_float_current_dot_= model_.getCurrentComDot();
  
  // Pelvis에서 왼쪽, 오른쪽 발 FK
  lfoot_float_current_ = model_.getCurrentTrasmfrom((DyrosJetModel::EndEffector)0); 
  rfoot_float_current_ = model_.getCurrentTrasmfrom((DyrosJetModel::EndEffector)1);  
  r_ft_ = model_.getRightFootForce();
  l_ft_ = model_.getLeftFootForce();
  imu_acc_ = model_.getImuAccel();
  imu_ang_ = model_.getImuAngvel();
  imu_grav_rpy_ = model_.getImuGravityDirection();
  pelv_float_current_.setIdentity();
  
  if(foot_step_(current_step_num_, 6) == 0) //right foot support
  {
    supportfoot_float_current_ = rfoot_float_current_;
  }
  else if(foot_step_(current_step_num_, 6) == 1)
  {
    supportfoot_float_current_ = lfoot_float_current_;
  }

  Eigen::Vector3d com_support_prev_;
  Eigen::Vector3d com_support_dot_prev_;
  
  com_support_prev_ = com_support_current_;
  com_support_dot_prev_ = com_support_dot_current_;
  pelv_support_current_ = DyrosMath::inverseIsometry3d(supportfoot_float_current_);
  lfoot_support_current_ = DyrosMath::multiplyIsometry3d(pelv_support_current_,lfoot_float_current_);
  rfoot_support_current_ = DyrosMath::multiplyIsometry3d(pelv_support_current_,rfoot_float_current_);
  com_support_current_ =  DyrosMath::multiplyIsometry3dVector3d(pelv_support_current_, com_float_current_);
  current_leg_jacobian_l_ = model_.getLegJacobian((DyrosJetModel::EndEffector) 0);
  current_leg_jacobian_r_ = model_.getLegJacobian((DyrosJetModel::EndEffector) 1);

  com_support_dot_current_ = (com_support_current_ - com_support_prev_) / 0.005 ;
  com_support_ddot_current_ = (com_support_dot_current_ - com_support_dot_prev_) / 0.005 ;
  
  if(walking_tick_ == 0)
  {
    com_float_prev_ = com_float_current_;
  }  
  slowcalc_mutex_.lock();
  thread_q_ = current_q_;
  current_motor_q_leg_ = current_q_.segment<12>(0); // 시뮬레이션에서 모터 각도와 링크 각도에 같은 값이 들어가고 있음.
  current_link_q_leg_ = current_q_.segment<12>(0);
  slowcalc_mutex_.unlock();

}
void WalkingController::calculateFootStepTotal()
{
  /*
   * this function calculate foot steps which the robot should put on
   * algorith: set robot orientation to the destination -> go straight -> set target orientation on the destination
   *
   * foot_step_(current_step_num_, i) is the foot step where the robot will step right after
   * foot_step_(crrennt_step_num_, 6) = 0 means swingfoot is left(support foot is right)
   */
  // 목표 각도와 거리를 정함 -> 로봇을 회전 시키고, 직진 -> 다시 원상태로 회전 복구.
   
  double initial_rot = 0.0;
  double final_rot = 0.0;
  double initial_drot = 0.0;
  double final_drot = 0.0;

  initial_rot = atan2(target_y_, target_x_); // GUI 명령에서 받아오는 값으로 회전 각도 계산.
  
  if(initial_rot > 0.0) // 반시계 방향으로 회전, GUI 에서는 이것만 됨.
    initial_drot = 10*DEG2RAD;
  else // 시계 방향으로 회전.
    initial_drot = -10*DEG2RAD;

  unsigned int initial_total_step_number = initial_rot/initial_drot; // 처음 실제 회전해야하는 각도 / 10deg , 정수만 추려냄. 
  double initial_residual_angle = initial_rot - initial_total_step_number*initial_drot; // 정수만 추려낸 나머지 회전해야 하는 각도.

  final_rot = target_theta_ - initial_rot; // 원하는 회전 각도가 없으면, 처음각도에 음수만 붙이면 다시 원상태로 돌릴 수 있음.
  if(final_rot > 0.0)
    final_drot = 10*DEG2RAD; // 반시계 방향 회전.
  else
    final_drot = -10*DEG2RAD; // 시계 방향 회전.

  unsigned int final_total_step_number = final_rot/final_drot; // 마지막 실제 회전해야 하는 각도 / 10deg, 정수만 추려냄.
  double final_residual_angle = final_rot - final_total_step_number*final_drot; // 정수만 추려낸 나머지 회전해야 하는 각도.
  double length_to_target = sqrt(target_x_*target_x_+target_y_*target_y_); // 목표까지의 최종 거리.
  double dlength = step_length_x_; 
  unsigned int middle_total_step_number = length_to_target/(dlength); // 목표 최종 거리를 Step 보폭으로 나눈 정수로 표현한 스텝 갯수.
  double middle_residual_length = length_to_target - middle_total_step_number*(dlength); // 나머지 거리.
  
  if(length_to_target == 0) // 제자리 걸음
  {
    middle_total_step_number = 50; //walking on the spot 10 times
    dlength = 0;
  }
  ////// 목표 x, y 변위, 설정된 step length를 이용해서 middle_total_step number를 구함.
  
  unsigned int number_of_foot_step;

  int del_size;

  del_size = 1;
  number_of_foot_step = initial_total_step_number*del_size + middle_total_step_number*del_size + final_total_step_number*del_size;
  
  if(initial_total_step_number != 0 || abs(initial_residual_angle) >= 0.0001) // 회전이 있는 경우.
  {
    if(initial_total_step_number % 2 == 0) // 회전해야 하는 Step 수가 짝수 일 경우. 
      number_of_foot_step = number_of_foot_step + 2*del_size; // 예를 들어, 45도 회전일 경우, initial_total_step_number = 4, 4번째 발이 목표에 오고, 5번째 발이 4번째 발과 맞춰짐, 마지막으로 한 발자국은 그냥 제자리 걸음. 
    else // 회전해야 하는 Step 수가 홀수 일 경우
    {
      if(abs(initial_residual_angle)>= 0.0001) // 회전해야 하는 Step 수가 홀수면서 나머지 회전 해야하는 각도가 남았을 경우, 기존에서 3 Step 추가
        number_of_foot_step = number_of_foot_step + 3*del_size; // 예를 들어, 15도 회전일 경우, initial_total_step_number = 1, 2번째 발이 목표에 오고 3번째 발이 2번쨰 발과 맞춰짐, 마지막 발은 제자리.
      else
        number_of_foot_step = number_of_foot_step + del_size; // 회전해야 하는 Step 수가 홀수면서 나머지 회전 해야하는 각도가 남지 않았을 경우, 기존에서 1 Step (반대 발) 추가
    }
  }
  
  if(middle_total_step_number != 0 || abs(middle_residual_length)>=0.0001) //제자리 걸음이 아닌 경우.
  {
    if(middle_total_step_number % 2 == 0) // 목표거리가 Step 길이의 짝수배 일 경우. 예를 들어 목표 길이가 1.2m 인 경우 middle_total_step_number는 6임. 6번째 목표에 맞춰지고 7번째 발이 6번째 발 옆에, 8번째 발은 제자리걸음.
      number_of_foot_step = number_of_foot_step + 2*del_size; // 첫 발이 오른발인 경우 마지막 Swing이 왼발이고, 오른발은 마지막 Swing 발 바로 옆에 위치, 왼발은 다시한번 제자리 step
    else // 목표거리가 Step 길이의 홀수배 일 경우
    {
      if(abs(middle_residual_length)>= 0.0001) 
        number_of_foot_step = number_of_foot_step + 3*del_size; // 예를 들어 목표 길이가 1.1m 인 경우, middle_total_step_number = 5, 나머지 10cm 때문에 6번째 발이 절반 Swing, 7번재 Swing발이 바로 옆에 위치, 8번째 발이 마지막 제자리 step, 더해져 8이됨.   
      else
        number_of_foot_step = number_of_foot_step + del_size; // 예를 들어 목표 길이가 1m 인 경우 middle_total_step_number는 5임. 거기서 마지막 Swing 발이 이전 Swing발 바로 옆에 위치, 1이 더해져 6이 됨. (목표 변위에는 영향 없음.)
    }
  }
  
  if(final_total_step_number != 0 || abs(final_residual_angle) >= 0.0001) // 마지막으로 다시 로봇의 방향을 되돌리는 회전이 있는 경우.
  {
    if(abs(final_residual_angle) >= 0.0001) // 나머지 회전 각도가 남았을 경우.
      number_of_foot_step = number_of_foot_step + 2*del_size; // 최종적으로 2 step 더해줌.
    else
      number_of_foot_step = number_of_foot_step + del_size; // 나머지 회전 각도가 남지 않았을 경우, 한걸음만 추가.
  }
  
  // 여기까지 발걸음 개수 계산.

  foot_step_.resize(number_of_foot_step, 7);
  foot_step_.setZero();
  foot_step_support_frame_.resize(number_of_foot_step, 7);
  foot_step_support_frame_.setZero();
  foot_step_support_frame_offset_.resize(number_of_foot_step, 7);
  foot_step_support_frame_offset_.setZero();

  int index = 0;
  int temp, temp2, temp3, is_right;

  if(is_right_foot_swing_ == true)
    is_right = 1;
  else
    is_right = -1;


  temp = -is_right; //right foot will be first swingfoot
  temp2 = -is_right;
  temp3 = -is_right;


  if(initial_total_step_number != 0 || abs(initial_residual_angle) >= 0.0001) // 회전이 있는 경우.
  {
    for (int i =0 ; i < initial_total_step_number; i++) // 45도 회전하면, index가 0 ~ 3이므로 4발자국을 먼저 찍는 것.
    {
      temp *= -1;
      foot_step_(index,0) = temp*0.127794*sin((i+1)*initial_drot); // 10도씩 회전 시키는데, 회전 했을때의 X방향 좌표.
      foot_step_(index,1) = -temp*0.127794*cos((i+1)*initial_drot); // 10도씩 회전 시킬때, 회전 했을때의 Y방향 좌표.
      foot_step_(index,5) = (i+1)*initial_drot; // Z방향 10도 돌린 Foot step.
      foot_step_(index,6) = 0.5 + 0.5*temp; // 오른 발이 Swing 일때는, 1.
      index++;
    }

    if(temp == is_right)
    {
      if(abs(initial_residual_angle) >= 0.0001) // 예를 들어 15 deg 회전이면, 첫번째 오른 Swing 발이 10deg, 두번째 왼 Swing 발이 15deg, 세번째 오른 Swing발이 15도, 나머지 제자리.
      {
        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;

        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot+initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;

        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot+initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;

      }
      else
      {
        temp *= -1;

        foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
        foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
        foot_step_(index,6) = 0.5 + 0.5*temp;
        index++;
      }
    }
    else if(temp == -is_right) // 여분으로 움직여야 할 각도가 있고 시작을 오른 발 Swing으로 할 때.
    {
      temp *= -1;

      foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle); // 나머지 각도 채워줌.
      foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
      foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
      foot_step_(index,6) = 0.5 + 0.5*temp;
      index ++;

      temp *= -1;

      foot_step_(index,0) = temp*0.127794*sin((initial_total_step_number)*initial_drot + initial_residual_angle); // 반대 발도 나머지 각도 채워줌.
      foot_step_(index,1) = -temp*0.127794*cos((initial_total_step_number)*initial_drot + initial_residual_angle);
      foot_step_(index,5) = (initial_total_step_number)*initial_drot + initial_residual_angle;
      foot_step_(index,6) = 0.5 + 0.5*temp;
      index ++;
    }
  }

  if(middle_total_step_number != 0 || abs(middle_residual_length) >= 0.0001)
  {
    for (int i = 0 ; i < middle_total_step_number; i++) // 예를 들어 1m 면, middle_total_step_number = 5,
    {
      temp2 *= -1; // Swing발이 오른발이면 temp는 1.

      foot_step_(index,0) = cos(initial_rot)*(dlength*(i+1)) + temp2*sin(initial_rot)*(0.127794); // Z방향으로 initial_rot 만큼 회전 시켰을 때, X방향 좌표.
      foot_step_(index,1) = sin(initial_rot)*(dlength*(i+1)) - temp2*cos(initial_rot)*(0.127794); // Z방향으로 initial_rot 만큼 회전 시켰을 때, Y방향 좌표.
      foot_step_(index,5) = initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp2;
      index ++;
    }

    if(temp2 == is_right)
    {
      if(abs(middle_residual_length) >= 0.0001) // 예를 들어 목표가 1.1m이면 20cm 보폭으로 나눴을 때 10cm의 나머지 거리가 생김.
      { // 3 Step 모두 같은 발자국 위치를 가짐.
        temp2 *= -1; // temp를 -1로 바꿈 (왼발을 Swing 발로)

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2; // 0 -> 오른 발이 Support foot.
        //나머지 거리만큼 왼발로 감. 목표치를 만족 시킴.
        index++;

        temp2 *= -1; // temp를 1로 바꿈 (오른발을 Swing 발로)

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2;
        index++;

        temp2 *= -1; // 제자리 Step.

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2;
        index++;
      }
      else // 예를 들어 1m이면 20cm 보폭으로 나눴을때 나머지가 없기 때문에, 
      {
        temp2 *= -1; // 왼발을 Swing 발로 바꿔주고, 제자리에서 Step. (residual_length가 0이기 때문에)

        foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
        foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
        foot_step_(index,5) = initial_rot;
        foot_step_(index,6) = 0.5 + 0.5*temp2; // 왼발이 Swing 발이기 때문에 temp는 -1, 그러므로 foot_step(i,6) = 0
        index++;
      }
    }
    else if(temp2 == -is_right) // 끝이 왼발 swing으로 끝났을 때.
    {
      temp2 *= -1; // temp를 1로 바꿔주고 (오른쪽 발 Swing)

      foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
      foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
      foot_step_(index,5) = initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp2; 
      index++;

      temp2 *= -1;
      // 제자리 Step.
      foot_step_(index,0) = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) + temp2*sin(initial_rot)*(0.127794);
      foot_step_(index,1) = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length) - temp2*cos(initial_rot)*(0.127794);
      foot_step_(index,5) = initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp2;
      index++;
    }
  }

  double final_position_x = cos(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length); // 최종적으로 도착한 목적지의 X 좌표.
  double final_position_y = sin(initial_rot)*(dlength*(middle_total_step_number) + middle_residual_length); // 최종적으로 도착한 목적지의 Y 좌표.


  if(final_total_step_number != 0 || abs(final_residual_angle) >= 0.0001) // 최종적으로 로봇을 원래 방향대로 회전 해야 할 경우.
  {
    for(int i = 0 ; i < final_total_step_number; i++)
    {
      temp3 *= -1; // temp를 1로 만듬 (오른발 Swing), 예를 들어 45도 회전해야 한다고 하면, 4번 회전 먼저 하고

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin((i+1)*final_drot + initial_rot);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos((i+1)*final_drot + initial_rot);
      foot_step_(index,5) = (i+1)*final_drot + initial_rot;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;
    }

    if(abs(final_residual_angle) >= 0.0001) // 남은 각도는 10deg 미만이기 때문에 한번에 감.
    {
      temp3 *= -1; // 두 Step의 발자국 위치는 같고 첫 Step 먼저 목표에 위치하고 두번째 Step이 다음으로 위치.

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin(target_theta_);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos(target_theta_);
      foot_step_(index,5) = target_theta_;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;

      temp3 *= -1;

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin(target_theta_);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos(target_theta_);
      foot_step_(index,5) = target_theta_;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;
    }
    else
    { // 여분의 각도가 없다면, 반대발은 그냥 제자리 걸음 한번.
      temp3 *= -1;

      foot_step_(index,0) = final_position_x + temp3*0.127794*sin(target_theta_);
      foot_step_(index,1) = final_position_y - temp3*0.127794*cos(target_theta_);
      foot_step_(index,5) = target_theta_;
      foot_step_(index,6) = 0.5 + 0.5*temp3;
      index++;
    }
  }
}

void WalkingController::getZmpTrajectory()
{
  unsigned int planning_step_number = 3;

  unsigned int norm_size = 0;

  if(current_step_num_ >= total_step_num_ - planning_step_number)
    norm_size = (t_last_-t_start_+1)*(total_step_num_-current_step_num_)+20*hz_;
  else
    norm_size = (t_last_-t_start_+1)*(planning_step_number);
  if(current_step_num_ == 0)
    norm_size = norm_size + t_temp_+1;
  addZmpOffset();
  zmpGenerator(norm_size, planning_step_number);  
}

void WalkingController::floatToSupportFootstep()
{
  Eigen::Isometry3d reference;

  if(current_step_num_ == 0) 
  {
    if(foot_step_(0,6) == 0) //right support
    {
      reference.translation() = rfoot_float_init_.translation(); 
      reference.translation()(2) = 0.0;
      reference.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(rfoot_float_init_.linear())(2));
      reference.translation()(0) = 0.0;
    }
    else //left support
    {
      reference.translation() = lfoot_float_init_.translation();  
      reference.translation()(2) = 0.0; // Z는 Pelvis 기준으로 발까지 - 변위 였고,
      reference.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(lfoot_float_init_.linear())(2)); // Pelvis부터 왼쪽 발의 Rotation Matirx를 alpha, beta, gamma로 바꾸고, gamma(Z 축)만 빼서 Rotation matrix로 다시 바꿈.
      reference.translation()(0) = 0.0; // X는 Pevlis 기준으로 Offset 되어 있는 변위 였다.
    } 
  }
  else
  {
    reference.linear() = DyrosMath::rotateWithZ(foot_step_(current_step_num_-1,5)); // 현재 current_step_num에 맞는 계획되어 있는 Gloabal에서 Foot step의 Z방향 회전 각도를 Z회전 Rotation matrix의 변수에 넣는다.
    for(int i=0 ;i<3; i++)
      reference.translation()(i) = foot_step_(current_step_num_-1,i); // 현재 current_step_num에 맞는 계획되어 있는 Gloabal에서 Foot step의 X,Y,Z 좌표를 넣음.
  } 

  Eigen::Vector3d temp_local_position;
  Eigen::Vector3d temp_global_position;
  // temp_local,global_position, 여기서는 Global에서 본 Foot step 좌표들을 Current_step_num_에 따른 현재 지지발 좌표계에서 보기위해 사용됨. 
  if(current_step_num_ == 0)
  {
    for(int i = 0; i < total_step_num_; i++) // total_step_num_은 foot_step_.col(1).size()에 의해 결정된다.
    {
      for(int j = 0; j < 3; j ++)
        temp_global_position(j) = foot_step_(i,j); // current_step_num_과 상관 없이 계획된 Foot step의 x,y,z를 전부 넣는다.

      temp_local_position = reference.linear().transpose()*(temp_global_position - reference.translation()); 
      // local_position이 의미하는 것은, 계획된 Foot step이 모두 temp_global_position에 들어가있고, reference.translation()은 current_step_num에 따라 차례로 증가되서
      // 그 차이는 현재 지지발에서 본 Global의 Footstep 들이라고 할 수 있음.
      for(int j=0; j<3; j++)
        foot_step_support_frame_(i,j) = temp_local_position(j); // 지지발에서 본 계획된 Foot step들.

      foot_step_support_frame_(i,3) = foot_step_(i,3); // roll
      foot_step_support_frame_(i,4) = foot_step_(i,4); // pitch
      foot_step_support_frame_(i,5) = foot_step_(i,5) - supportfoot_float_init_(5);      
    }
  }
  else
  {
    for(int i = 0; i < total_step_num_; i ++) // total_step_num_은 foot_step_.col(1).size()에 의해 결정된다.
    {
      for(int j = 0; j < 3; j ++)
        temp_global_position(j) = foot_step_(i,j); // current_step_num_과 상관 없이 계획된 Foot step의 x,y,z를 전부 넣는다.

      temp_local_position = reference.linear().transpose()*(temp_global_position - reference.translation());
      
      for(int j = 0; j < 3; j ++)
        foot_step_support_frame_(i,j) = temp_local_position(j); 
       
      foot_step_support_frame_(i,3) = foot_step_(i,3);
      foot_step_support_frame_(i,4) = foot_step_(i,4);
      foot_step_support_frame_(i,5) = foot_step_(i,5) - foot_step_(current_step_num_-1,5); // 지지발 기준에서 본 거니까 계속 빼줘야함.
      
    }
  }
   
  for(int j=0;j<3;j++)
    temp_global_position(j) = swingfoot_float_init_(j); // swingfoot_float_init_은 walking_tick = 0 일때 Pelvis에서 본 Swing 발의 Position, orientation. 
  
  temp_local_position = reference.linear().transpose()*(temp_global_position - reference.translation());
  // 여기서의 temp_local_position은 현재 current_step_num_에 따라 계획된 global에서의 footstep 좌표계에서 맨 처음 Swing발 좌표계를 본 것.
 
  for(int j=0;j<3;j++)
    swingfoot_support_init_(j) = temp_local_position(j);

  swingfoot_support_init_(3) = swingfoot_float_init_(3);
  swingfoot_support_init_(4) = swingfoot_float_init_(4);

  if(current_step_num_ == 0)
    swingfoot_support_init_(5) = swingfoot_float_init_(5) - supportfoot_float_init_(5);
  else
    swingfoot_support_init_(5) = swingfoot_float_init_(5) - foot_step_(current_step_num_-1,5);

  for(int j=0;j<3;j++)
    temp_global_position(j) = supportfoot_float_init_(j);

  temp_local_position = reference.linear().transpose()*(temp_global_position-reference.translation());
  // 여기서의 temp_local_position은 현재 current_step_num_에 따라 계획된 global에서의 footstep 좌표계에서 맨 처음 지지발 좌표계를 본 것.
  for(int j=0;j<3;j++)
    supportfoot_support_init_(j) = temp_local_position(j);
  
  supportfoot_support_init_(3) = supportfoot_float_init_(3);
  supportfoot_support_init_(4) = supportfoot_float_init_(4);

  if(current_step_num_ == 0)
    supportfoot_support_init_(5) = 0;
  else
    supportfoot_support_init_(5) = supportfoot_float_init_(5) - foot_step_(current_step_num_-1,5);  
   
}

void WalkingController::updateInitialState()
{
  if(walking_tick_ == 0)
  {
    calculateFootStepTotal();
    thread_tick_ = 0;  

    q_init_ = current_q_;
    lfoot_float_init_ = model_.getCurrentTrasmfrom((DyrosJetModel::EndEffector)(0));
    rfoot_float_init_ = model_.getCurrentTrasmfrom((DyrosJetModel::EndEffector)(1));
    com_float_init_ = model_.getCurrentCom();
    pelv_float_init_.setIdentity();
    
    Eigen::Isometry3d ref_frame;

    if(foot_step_(0, 6) == 0)  //right foot support
    {
      ref_frame = rfoot_float_init_;
    }
    else if(foot_step_(0, 6) == 1)
    {
      ref_frame = lfoot_float_init_;
    }
    
    lfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),lfoot_float_init_);
    rfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),rfoot_float_init_);
    pelv_support_init_ = DyrosMath::inverseIsometry3d(ref_frame);
    
    com_support_init_ = pelv_support_init_.linear()*com_float_init_ + pelv_support_init_.translation();
    
    pelv_support_euler_init_ = DyrosMath::rot2Euler(pelv_support_init_.linear());
    rfoot_support_euler_init_ = DyrosMath::rot2Euler(rfoot_support_init_.linear());
    lfoot_support_euler_init_ = DyrosMath::rot2Euler(lfoot_support_init_.linear());

    supportfoot_float_init_.setZero();
    swingfoot_float_init_.setZero();


    if(foot_step_(0,6) == 1) //left suppport foot
    {
      for(int i=0; i<2; i++)
        supportfoot_float_init_(i) = lfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        supportfoot_float_init_(i+3) = DyrosMath::rot2Euler(lfoot_float_init_.linear())(i);

      for(int i=0; i<2; i++)
        swingfoot_float_init_(i) = rfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        swingfoot_float_init_(i+3) = DyrosMath::rot2Euler(rfoot_float_init_.linear())(i);

      supportfoot_float_init_(0) = 0.0;
      swingfoot_float_init_(0) = 0.0;
    }
    else
    {
      for(int i=0; i<2; i++)
        supportfoot_float_init_(i) = rfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        supportfoot_float_init_(i+3) = DyrosMath::rot2Euler(rfoot_float_init_.linear())(i);

      for(int i=0; i<2; i++)
        swingfoot_float_init_(i) = lfoot_float_init_.translation()(i);
      for(int i=0; i<3; i++)
        swingfoot_float_init_(i+3) = DyrosMath::rot2Euler(lfoot_float_init_.linear())(i);

      supportfoot_float_init_(0) = 0.0;
      swingfoot_float_init_(0) = 0.0;
    }

    zc_ = com_support_init_(2);
    
    pelv_suppprt_start_ = pelv_support_init_;

    total_step_num_ = foot_step_.col(1).size();
    
    xi_ = com_support_init_(0); // preview parameter
    yi_ = com_support_init_(1);
    
  }
  else if(current_step_num_ != 0 && walking_tick_ == t_start_) // step change 할 때.
  {  
    Eigen::Matrix<double, DyrosJetModel::MODEL_WITH_VIRTUAL_DOF, 1> q_temp, qdot_temp;
    q_temp.setZero();
    qdot_temp;
    q_temp.segment<28>(6) = current_q_.segment<28>(0);
    qdot_temp.segment<28>(6)= current_qdot_.segment<28>(0);  
    //q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);
     
    model_.updateKinematics(q_temp, qdot_temp);
     
    lfoot_float_init_ = model_.getCurrentTrasmfrom((DyrosJetModel::EndEffector)(0));
    rfoot_float_init_ = model_.getCurrentTrasmfrom((DyrosJetModel::EndEffector)(1));
    com_float_init_ = model_.getCurrentCom();  

    pelv_float_init_.setIdentity();

    Eigen::Isometry3d ref_frame;

    if(foot_step_(current_step_num_, 6) == 0)  //right foot support
    {
      ref_frame = rfoot_float_init_;
    }
    else if(foot_step_(current_step_num_, 6) == 1)
    {
      ref_frame = lfoot_float_init_;
    }

    lfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),lfoot_float_init_);
    rfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),rfoot_float_init_);
    pelv_support_init_ = DyrosMath::inverseIsometry3d(ref_frame);
    com_support_init_ = pelv_support_init_.linear()*com_float_init_ + pelv_support_init_.translation();
    //pelv_support_init -> 지지발에서 본 pelvis rot 
    pelv_support_euler_init_ = DyrosMath::rot2Euler(pelv_support_init_.linear());
    rfoot_support_euler_init_ = DyrosMath::rot2Euler(rfoot_support_init_.linear());
    lfoot_support_euler_init_ = DyrosMath::rot2Euler(lfoot_support_init_.linear());

  }  
}

void WalkingController::updateNextStepTime()
{
  if(walking_tick_ == t_last_) // 840 tick 부터 시작해서 240 tick 주기
  {
    if(current_step_num_ != total_step_num_-1)
    {
      t_start_ = t_last_ +1; //841 tick
      t_start_real_ = t_start_ + t_rest_init_; // 841 + 10 = 851 tick
      t_last_ = t_start_ + t_total_ -1; // 840 + 240 = 1080 tick 
      
      current_step_num_ ++; // 시작 단계에서, walking tick이 600, 한 step 걸었을 때, 즉 t_last의 시간과 같게 되면 current_step_num_ 이 증가, 매 step은 240tick, 1.2초가 걸림.
    }    
  }
   if(current_step_num_ == total_step_num_-1 && walking_tick_ >= t_last_ + t_total_)// + 160) // 여기에 tick을 늘려주는 만큼 보행 완료 대기 시간을 늘릴 수 있음.
   {
     walking_enable_ = false;
   }
   walking_tick_ ++;
}

void WalkingController::addZmpOffset()
{
  lfoot_zmp_offset_ = -0.02;
  rfoot_zmp_offset_ = 0.02;

  foot_step_support_frame_offset_ = foot_step_support_frame_;

  if(foot_step_(0,6) == 0) //right support foot
  {
    supportfoot_support_init_offset_(1) = supportfoot_support_init_(1) + rfoot_zmp_offset_;
    swingfoot_support_init_offset_(1) = swingfoot_support_init_(1) + lfoot_zmp_offset_;
  }
  else
  {
    supportfoot_support_init_offset_(1) = supportfoot_support_init_(1) + lfoot_zmp_offset_;
    swingfoot_support_init_offset_(1) = swingfoot_support_init_(1) + rfoot_zmp_offset_;
  }

  for(int i=0; i<total_step_num_; i++)
  {
    if(foot_step_(i,6) == 0)//right support, left swing
    {
      foot_step_support_frame_offset_(i,1) += lfoot_zmp_offset_;
    }
    else
    {
      foot_step_support_frame_offset_(i,1) += rfoot_zmp_offset_;
    }
  }
}

void WalkingController::zmpGenerator(const unsigned int norm_size, const unsigned planning_step_num)
{
  ref_zmp_.resize(norm_size, 2);
  ref_com_.resize(norm_size, 2);
  com_offset_.setZero();
  Eigen::VectorXd temp_px;
  Eigen::VectorXd temp_py;
  Eigen::VectorXd temp_cx;
  Eigen::VectorXd temp_cy;
  
  unsigned int index = 0;
  // 매 tick 마다 zmp가 3발 앞까지 계산 된다. 

  if(current_step_num_ == 0) // Walking을 수행 할 때, 정지 상태 일때 3초 동안 Ref X ZMP를 0으로 보냄. Y ZMP는 제자리 유지.  
  {
    for (int i = 0; i <= t_temp_; i++) //600 tick
    {
      if(i <= 0.5*hz_) 
      {
        ref_zmp_(i,0) = com_support_init_(0) + com_offset_(0);
        ref_zmp_(i,1) = com_support_init_(1) + com_offset_(1);
      }
      else if(i < 1.5*hz_) 
      {
        double del_x = i - 0.5*hz_;
        ref_zmp_(i,0) = com_support_init_(0) + com_offset_(0) - del_x * (com_support_init_(0) + com_offset_(0))/(1.0*hz_);
        ref_zmp_(i,1) = com_support_init_(1) + com_offset_(1);
      }
      else 
      {
        ref_zmp_(i,0) = 0.0;
        ref_zmp_(i,1) = com_support_init_(1) + com_offset_(1);
      }      
      index++;
      ref_com_(i,0) = DyrosMath::cubic(i, 0, 2*hz_,com_support_init_(0) + com_offset_(0), 0, 0, 0);
      
      if(i < 2*hz_)
      {
        ref_com_(i,1) = com_support_init_(1) + com_offset_(1);
      }
      else
      { 
        Eigen::Vector3d ref_y_com ;
        ref_y_com = (DyrosMath::QuinticSpline(i, 2*hz_, 3*hz_,com_support_init_(1)+com_offset_(1), 0, 0, com_support_init_(1)+com_offset_(1), 0.289384/hz_,0));
        ref_com_(i,1) = ref_y_com(0);
      }
    }    
  }
  if(current_step_num_ >= total_step_num_ - planning_step_num)
  {    
    for(unsigned int i = current_step_num_; i < total_step_num_; i++)
    {
      onestepZmp_MJ(i, temp_px, temp_py);
      OfflineCoM_MJ(i, temp_cx, temp_cy);
     
      for(unsigned int j = 0; j < t_total_; j++)
      {
        ref_zmp_(index + j, 0) = temp_px(j);
        ref_zmp_(index + j, 1) = temp_py(j);
        ref_com_(index + j, 0) = temp_cx(j);
        ref_com_(index + j, 1) = temp_cy(j);        
      }
      index = index + t_total_;
    }
    
    for(unsigned int j = 0; j < 20*hz_; j++)
    {
      ref_zmp_(index + j, 0) = ref_zmp_(index -1, 0);
      ref_zmp_(index + j, 1) = ref_zmp_(index -1, 1);
      ref_com_(index + j, 0) = ref_com_(index -1, 0);
      ref_com_(index + j, 1) = ref_com_(index -1, 1);
    }
    
    if((current_step_num_ == total_step_num_ - 1)) //마지막 Step에만 사용 됨. 오른발이 끝발 -2, 왼발이 끝발 -1
    { Eigen::Vector3d ref_y_com ;
            
      for(int i = 0; i < 240; i++)
      {
        ref_y_com = DyrosMath::QuinticSpline(i+239, 120, 479, 0.031081, -2.60209e-18, 1.05331e-05, 0.12779, 0, 0);
        ref_com_(index + i, 1) = ref_y_com(0) ;
      }
    } 
    index = index + 20*hz_;     
  }
  else // 보행 중 사용 하는 Ref ZMP
  {
    for(unsigned int i = current_step_num_; i < current_step_num_ + planning_step_num; i++)  
    {
      onestepZmp_MJ(i, temp_px, temp_py);
      OfflineCoM_MJ(i, temp_cx, temp_cy);
      for (unsigned int j = 0; j < t_total_; j++) // 1 step 보행은 1.2초, 240 tick
      {
        ref_zmp_(index+j,0) = temp_px(j);
        ref_zmp_(index+j,1) = temp_py(j);
        ref_com_(index+j, 0) = temp_cx(j);
        ref_com_(index+j, 1) = temp_cy(j);
      }      
      index = index + t_total_; // 참조 zmp가 이만큼 쌓였다.      
      // 결국 실제 로봇 1Hz마다 720개의 ref_zmp를 생성함. 3.6초
    }   
  }    
}

void WalkingController::OfflineCoM_MJ(unsigned int current_step_number, Eigen::VectorXd& temp_cx, Eigen::VectorXd& temp_cy)
{
  temp_cx.resize(t_total_); 
  temp_cy.resize(t_total_);
  temp_cx.setZero();
  temp_cy.setZero();

  double wn = 0;
  double A = 0, B = 0, Cx1 = 0, Cx2 = 0, Cy1 = 0, Cy2 = 0, Kx = 0, Ky = 0 ;
  if(current_step_number == 0)
  {
    wn = sqrt(GRAVITY / com_support_init_(2));
    A = -(foot_step_support_frame_(current_step_number, 1) )/2 ;
    B =  (supportfoot_support_init_(0) + foot_step_support_frame_(current_step_number, 0))/2;
    Kx = (B * 0.15 * wn) / ((0.15*wn) + tanh(wn*(0.45)));
    Ky = (A * 0.15 * wn * tanh(wn*0.45))/(1 + 0.15 * wn * tanh(wn*0.45)); 
    Cx1 = Kx - B;
    Cx2 = Kx/(wn*0.15);
    Cy1 = Ky - A;
    Cy2 = Ky/(wn*0.15);      
    for(int i = 0; i < t_total_; i++)
    {
      temp_cx(i) = DyrosMath::cubic(i, 0, t_total_-1,0 ,0.10154 , 0, Kx/(t_rest_init_ + t_double1_));
      if(i >= 0 && i < (t_rest_init_ + t_double1_))
      {  
        temp_cy(i) = com_offset_(1) + com_support_init_(1) + Ky / (t_rest_init_ + t_double1_) * (i+1);
      }
      else if(i >= (t_rest_init_ + t_double1_) && i < t_total_ - t_rest_last_ - t_double2_ )
      {
        temp_cy(i) = A + com_offset_(1) + com_support_init_(1) + Cy1 *cosh(wn*(i/hz_ - 0.15)) + Cy2*sinh(wn*(i/hz_-0.15)) ;
      }
      else if(i >= t_total_ - t_rest_last_ - t_double2_  && i < t_total_) //1.05 ~ 1.15초 , 210 ~ 230 tick 
      {
        temp_cy(i) = Ky + (supportfoot_support_init_(1) + foot_step_support_frame_(current_step_number, 1))/2 + Ky/(t_rest_last_ + t_double2_)*-(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }     
    } 
  }
  else if(current_step_number == 1)
  {
    wn = sqrt(GRAVITY / com_support_init_(2));
    A = (foot_step_support_frame_(current_step_number-1, 1) - supportfoot_support_init_(1))/2 ;
    B = foot_step_support_frame_(current_step_number-1, 0) - (supportfoot_support_init_(0) + foot_step_support_frame_(current_step_number-1, 0))/2;
    Kx = (B * 0.15 * wn) / ((0.15*wn) + tanh(wn*0.45));
    Ky = (A * 0.15 * wn * tanh(wn*0.45))/(1 + 0.15 * wn * tanh(wn*0.45)); 
    Cx1 = Kx - B;
    Cx2 = Kx/(wn*0.15);
    Cy1 = Ky - A;
    Cy2 = Ky/(wn*0.15);    
    for(int i = 0; i < t_total_; i++)
    {
      if(i >= 0 && i < (t_rest_init_ + t_double1_))
      { 
        temp_cx(i) = (foot_step_support_frame_(current_step_number-1, 0) + supportfoot_support_init_(0))/2 + Kx / (t_rest_init_+ t_double1_) * (i+1);
        temp_cy(i) = (foot_step_support_frame_(current_step_number-1, 1) + supportfoot_support_init_(1))/2 + Ky / (t_rest_init_+ t_double1_) * (i+1);
      }
      else if(i >= (t_rest_init_ + t_double1_) && i < (t_total_ - t_rest_last_ - t_double2_)) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_cx(i) = (supportfoot_support_init_(0) + foot_step_support_frame_(current_step_number-1, 0))/2 + Cx1 *cosh(wn*(i/hz_ - 0.15)) + Cx2*sinh(wn*(i/hz_-0.15)) + B;
        temp_cy(i) = A + (supportfoot_support_init_(1) + foot_step_support_frame_(current_step_number-1, 1))/2 + Cy1 *cosh(wn*(i/hz_ - 0.15)) + Cy2*sinh(wn*(i/hz_-0.15)) ;
      }
      else if(i >= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_) //1.05 ~ 1.2초 , 210 ~ 240 tick 
      { 
        temp_cx(i) = (foot_step_support_frame_(current_step_number, 0)+ foot_step_support_frame_(current_step_number-1, 0)) /2 -Kx + Kx/(t_rest_last_ + t_double2_)*(i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_cy(i) = Ky + (foot_step_support_frame_(current_step_number-1, 1) + foot_step_support_frame_(current_step_number, 1))/2 + Ky/(t_rest_last_ + t_double2_)*-(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
    }    
  }
  else
  {
    wn = sqrt(GRAVITY / com_support_init_(2));
    A = (foot_step_support_frame_(current_step_number-1, 1) - foot_step_support_frame_(current_step_number-2, 1))/2 ;
    B = foot_step_support_frame_(current_step_number-1, 0) - (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2;
    Kx = (B * 0.15 * wn) / ((0.15*wn) + tanh(wn*0.45));
    Ky = (A * 0.15 * wn * tanh(wn*0.45))/(1 + 0.15 * wn * tanh(wn*0.45)); 
    Cx1 = Kx - B;
    Cx2 = Kx/(wn*0.15);
    Cy1 = Ky - A;
    Cy2 = Ky/(wn*0.15);
    for(int i = 0; i < t_total_; i++)
    {
      if(i >= 0 && i < (t_rest_init_ + t_double1_))
      {
        temp_cx(i) = (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 + Kx/(t_rest_init_ + t_double1_)*(i+1);
        temp_cy(i) = (foot_step_support_frame_(current_step_number-2, 1) + foot_step_support_frame_(current_step_number-1, 1))/2 + Ky/(t_rest_init_ + t_double1_)*(i+1);
      }            
      else if(i >= (t_rest_init_ + t_double1_) && i < (t_total_ - t_rest_last_ - t_double2_)) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_cx(i) = (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 + Cx1 *cosh(wn*(i/hz_ - 0.15)) + Cx2*sinh(wn*(i/hz_-0.15)) + B;
        temp_cy(i) = A + (foot_step_support_frame_(current_step_number-2, 1) + foot_step_support_frame_(current_step_number-1, 1))/2 + Cy1 *cosh(wn*(i/hz_ - 0.15)) + Cy2*sinh(wn*(i/hz_-0.15)) ;
         
      }
      else if(i >= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_) //1.05 ~ 1.2초 , 210 ~ 240 tick 
      {
        temp_cx(i) = (foot_step_support_frame_(current_step_number, 0)+ foot_step_support_frame_(current_step_number-1, 0)) /2 -Kx + Kx/(t_rest_last_ + t_double2_)*(i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_cy(i) = Ky + (foot_step_support_frame_(current_step_number-1, 1) + foot_step_support_frame_(current_step_number, 1))/2 - Ky/(t_rest_last_ + t_double2_)*(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
      if(i >= (t_rest_init_ + t_double1_) && (current_step_num_ == total_step_num_ - 1) && i < t_total_ ) //X방향, 마지막 Step에만 사용 됨. 오른발이 끝발 -2, 왼발이 끝발 -1
      { 
        Eigen::Vector3d ref_x_com ;
        ref_x_com = DyrosMath::QuinticSpline(i, 30, 239, (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 + Kx, 0.289384/hz_, 0, (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 + B, 0, 0);//com_support_init_(1)+com_offset_(1), 0.289384/hz_,0));
        temp_cx(i) = ref_x_com(0);
      }
      if( i >= 120 && i < t_total_ && (current_step_num_ == total_step_num_ - 1)) //마지막 Step에만 사용 됨. 오른발이 끝발 -2, 왼발이 끝발 -1
      { 
        Eigen::Vector3d ref_y_com ;
        ref_y_com = DyrosMath::QuinticSpline(i, 120, 479, 0.031081, (Cy1 *sinh(wn*(120/hz_ - 0.15))*wn/hz_ + Cy2*cosh(wn*(120/hz_-0.15))*wn/hz_), (Cy1 *cosh(wn*(120/hz_ - 0.15))*wn/hz_*wn/hz_ + Cy2*sinh(wn*(120/hz_-0.15))*wn/hz_*wn/hz_), 0.12779, 0, 0);
        temp_cy(i) = ref_y_com(0);
      }      
    }
  }    
}

void WalkingController::onestepZmp_MJ(unsigned int current_step_number, Eigen::VectorXd& temp_px, Eigen::VectorXd& temp_py)
{
  temp_px.resize(t_total_); // 함수가 실행 될 때 마다, 240 tick의 참조 ZMP를 담는다. Realtime = 1.2초
  temp_py.resize(t_total_);
  temp_px.setZero();
  temp_py.setZero();

  double Kx = 0, Kx2 = 0, Ky = 0, Ky2 = 0, A = 0, B = 0, wn = 0;
  if(current_step_number == 0)
  { 
    wn = sqrt(GRAVITY / com_support_init_(2));
    A = -(foot_step_support_frame_(current_step_number, 1))/2 ;
    B =  (supportfoot_support_init_(0) + foot_step_support_frame_(current_step_number, 0))/2;
    Kx = (B * 0.15 * wn) / ((0.15*wn) + tanh(wn*(0.45)));
    Ky = (A * 0.15 * wn * tanh(wn*0.45))/(1 + 0.15 * wn * tanh(wn*0.45));
        
    for(int i = 0; i < t_total_; i++)
    {
      if(i >= 0 && i < t_rest_init_ + t_double1_) //0 ~ 0.15초 , 0 ~ 30 tick
      {
        temp_px(i) = 0;
        temp_py(i) = (com_offset_(1) + com_support_init_(1)) + Ky / (t_rest_init_ + t_double1_)* (i+1);
      }
      else if(i >= t_rest_init_ + t_double1_ && i < t_total_ - t_rest_last_ - t_double2_ ) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_px(i) = supportfoot_support_init_(0);
        temp_py(i) = supportfoot_support_init_(1);
      }
      else if(i >= t_total_ - t_rest_last_ - t_double2_  && i < t_total_) //1.05 ~ 1.2초 , 210 ~ 240 tick 
      {
        temp_px(i) = B - Kx + Kx / (t_rest_last_ + t_double2_) * (i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_py(i) = Ky + (supportfoot_support_init_(1) + foot_step_support_frame_(current_step_number, 1))/2 + Ky/(t_rest_last_ + t_double2_)*-(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }

    }
  
  }
  else if(current_step_number == 1)
  { 
    wn = sqrt(GRAVITY / com_support_init_(2));
    A = (foot_step_support_frame_(current_step_number-1, 1) - supportfoot_support_init_(1))/2 ;
    B = foot_step_support_frame_(current_step_number-1, 0) - (supportfoot_support_init_(0) + foot_step_support_frame_(current_step_number-1, 0))/2;
    Kx = (B * 0.15 * wn) / ((0.15*wn) + tanh(wn*(0.45)));
    Ky = (A * 0.15 * wn * tanh(wn*0.45))/(1 + 0.15 * wn * tanh(wn*0.45)); 
    
    for(int i = 0; i < t_total_; i++)
    {
      if(i >= 0 && i < t_rest_init_ + t_double1_) //0 ~ 0.15초 , 10 ~ 30 tick
      {
        temp_px(i) = (foot_step_support_frame_(current_step_number-1, 0) + supportfoot_support_init_(0))/2 + Kx / (t_rest_init_+ t_double1_) * (i+1);
        temp_py(i) = (foot_step_support_frame_(current_step_number-1, 1) + supportfoot_support_init_(1))/2 + Ky / (t_rest_init_+ t_double1_) * (i+1);
      }
      else if(i >= t_rest_init_ + t_double1_ && i < t_total_ - t_rest_last_ - t_double2_) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_px(i) = foot_step_support_frame_(current_step_number-1, 0);
        temp_py(i) = foot_step_support_frame_(current_step_number-1, 1);
      }
      else if(i >= t_total_ - t_rest_last_ - t_double2_ && i < t_total_) //1.05 ~ 1.2초 , 210 ~ 240 tick 
      {               
        temp_px(i) = (foot_step_support_frame_(current_step_number-1, 0) + foot_step_support_frame_(current_step_number, 0))/2 - Kx + Kx /(t_rest_last_ + t_double2_) * (i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_py(i) = Ky + (foot_step_support_frame_(current_step_number-1, 1) + foot_step_support_frame_(current_step_number, 1))/2 + Ky/(t_rest_last_ + t_double2_)*-(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
    }
  }
  else
  { 
    wn = sqrt(GRAVITY / com_support_init_(2));
    A = (foot_step_support_frame_(current_step_number-1, 1) - foot_step_support_frame_(current_step_number-2, 1))/2 ;
    B = foot_step_support_frame_(current_step_number-1, 0) - (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2;
    Kx = (B * 0.15 * wn) / ((0.15*wn) + tanh(wn*(0.45))) ;
    Ky = (A * 0.15 * wn * tanh(wn*0.45))/(1 + 0.15 * wn * tanh(wn*0.45)); 
    for(int i = 0; i < t_total_; i++)
    {      
      if(i >= 0 && i < t_rest_init_ + t_double1_) //0 ~ 0.15초 , 0 ~ 30 tick
      { 
        temp_px(i) = (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 + Kx/(t_rest_init_ + t_double1_)*(i+1);
        temp_py(i) = (foot_step_support_frame_(current_step_number-2, 1) + foot_step_support_frame_(current_step_number-1, 1))/2 + Ky/(t_rest_init_ + t_double1_)*(i+1);
      }
      else if(i >= (t_rest_init_ + t_double1_) && i < (t_total_ - t_rest_last_ - t_double2_)) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_px(i) = foot_step_support_frame_(current_step_number-1, 0) ;
        temp_py(i) = foot_step_support_frame_(current_step_number-1, 1) ;
      }
      else if( i >= (t_total_ - t_rest_last_ - t_double2_) && (i < t_total_) && (current_step_num_ == total_step_num_ - 1))
      {
        temp_px(i) = (foot_step_support_frame_(current_step_number, 0) + foot_step_support_frame_(current_step_number-1, 0))/2;
        temp_py(i) = Ky + (foot_step_support_frame_(current_step_number-1, 1) + foot_step_support_frame_(current_step_number, 1))/2 + Ky/(t_rest_last_ + t_double2_)*-(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }       
      else if(i >= (t_total_ - t_rest_last_ - t_double2_) && i < t_total_) //1.05 ~ 1.2초 , 210 ~ 240 tick 
      { 
        temp_px(i) = (foot_step_support_frame_(current_step_number, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 -Kx + Kx/(t_rest_last_ + t_double2_)*(i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_py(i) = Ky + (foot_step_support_frame_(current_step_number-1, 1) + foot_step_support_frame_(current_step_number, 1))/2 + Ky/(t_rest_last_ + t_double2_)*-(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
    } 
  }
  
}
/*
void WalkingController::onestepZmp_MJ(unsigned int current_step_number, Eigen::VectorXd& temp_px, Eigen::VectorXd& temp_py)
{
  temp_px.resize(t_total_); // 함수가 실행 될 때 마다, 240 tick의 참조 ZMP를 담는다. Realtime = 1.2초
  temp_py.resize(t_total_);
  temp_px.setZero();
  temp_py.setZero();

  double Kx = 0, Kx2 = 0, Ky = 0, Ky2 = 0;
  if(current_step_number == 0)
  {
    Kx = supportfoot_support_init_offset_(0); 
    Ky = supportfoot_support_init_offset_(1) - (com_offset_(1) + com_support_init_(1));
    Kx2 = (foot_step_support_frame_(current_step_number,0) - supportfoot_support_init_(0))/2 - supportfoot_support_init_offset_(0);
    Ky2 = (foot_step_support_frame_(current_step_number,1) - supportfoot_support_init_(1))/2 - supportfoot_support_init_offset_(1);
    
    for(int i = 0; i < t_total_; i++)
    {
      if(i < t_rest_init_) //0 ~ 0.05초 , 0 ~ 10 tick
      {
        temp_px(i) = 0;
        temp_py(i) = com_offset_(1) + com_support_init_(1);
      }
      else if(i >= t_rest_init_ && i < t_rest_init_ + t_double1_) //0.05 ~ 0.15초 , 10 ~ 30 tick
      {
        temp_px(i) = Kx / t_double1_ * (i+1 - t_rest_init_);
        temp_py(i) = com_offset_(1) + com_support_init_(1) + Ky / t_double1_ * (i+1 - t_rest_init_);
      }
      else if(i >= t_rest_init_ + t_double1_ && i < t_total_ - t_rest_last_ - t_double2_ ) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_px(i) = supportfoot_support_init_offset_(0);
        temp_py(i) = supportfoot_support_init_offset_(1);
      }
      else if(i >= t_total_ - t_rest_last_ - t_double2_  && i < t_total_ - t_rest_last_) //1.05 ~ 1.15초 , 210 ~ 230 tick 
      {
        temp_px(i) = supportfoot_support_init_offset_(0) + Kx2 / t_double2_ * (i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_py(i) = supportfoot_support_init_offset_(1) + Ky2 / t_double2_ * (i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
      else //1.15 ~ 1.2초 , 230 ~ 240 tick
      {
        temp_px(i) = temp_px(i-1);
        temp_py(i) = temp_py(i-1);
      }
    }
  
  }
  else if(current_step_number == 1)
  {
    Kx = foot_step_support_frame_offset_(current_step_number-1, 0) - (foot_step_support_frame_(current_step_number-1, 0) + supportfoot_support_init_(0))/2;
    Ky = foot_step_support_frame_offset_(current_step_number-1, 1) - (foot_step_support_frame_(current_step_number-1, 1) + supportfoot_support_init_(1))/2;
    Kx2 = (foot_step_support_frame_(current_step_number, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 - foot_step_support_frame_offset_(current_step_number-1, 0);
    Ky2 = (foot_step_support_frame_(current_step_number, 1) + foot_step_support_frame_(current_step_number-1, 1))/2 - foot_step_support_frame_offset_(current_step_number-1, 1);

    for(int i = 0; i < t_total_; i++)
    {
      if(i < t_rest_init_) //0 ~ 0.05초 , 0 ~ 10 tick
      { // current_step_num = 1 일때, supportfoot_support_init이 foot_step_support_frame_(current_step_number-2, 0) 의 역할을 해서 -0.2의 값을 갖는다.
        temp_px(i) = (foot_step_support_frame_(current_step_number-1, 0) + supportfoot_support_init_(0))/2 ;  
        temp_py(i) = (foot_step_support_frame_(current_step_number-1, 1) + supportfoot_support_init_(1))/2 ;
      }
      else if(i >= t_rest_init_ && i < t_rest_init_ + t_double1_) //0.05 ~ 0.15초 , 10 ~ 30 tick
      {
        temp_px(i) = (foot_step_support_frame_(current_step_number-1, 0) + supportfoot_support_init_(0))/2 + Kx / t_double1_ * (i+1 - t_rest_init_);
        temp_py(i) = (foot_step_support_frame_(current_step_number-1, 1) + supportfoot_support_init_(1))/2 + Ky / t_double1_ * (i+1 - t_rest_init_);
      }
      else if(i >= t_rest_init_ + t_double1_ && i < t_total_ - t_rest_last_ - t_double2_) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_px(i) = foot_step_support_frame_offset_(current_step_number-1, 0);
        temp_py(i) = foot_step_support_frame_offset_(current_step_number-1, 1);
      }
      else if(i >= t_total_ - t_rest_last_ - t_double2_ && i < t_total_ - t_rest_last_) //1.05 ~ 1.15초 , 210 ~ 230 tick 
      {
        temp_px(i) = foot_step_support_frame_offset_(current_step_number-1, 0) + Kx2 / t_double2_ * (i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_py(i) = foot_step_support_frame_offset_(current_step_number-1, 1) + Ky2 / t_double2_ * (i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
      else //1.15 ~ 1.2초 , 230 ~ 240 tick
      {
        temp_px(i) = temp_px(i-1);
        temp_py(i) = temp_py(i-1);
      }
    }
  }
  else
  {
    Kx = foot_step_support_frame_offset_(current_step_number-1, 0) - ((foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2);
    Ky = foot_step_support_frame_offset_(current_step_number-1, 1) - ((foot_step_support_frame_(current_step_number-2, 1) + foot_step_support_frame_(current_step_number-1, 1))/2);
    Kx2 = (foot_step_support_frame_(current_step_number, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 - foot_step_support_frame_offset_(current_step_number-1, 0);
    Ky2 = (foot_step_support_frame_(current_step_number, 1) + foot_step_support_frame_(current_step_number-1, 1))/2 - foot_step_support_frame_offset_(current_step_number-1, 1);

    for(int i = 0; i < t_total_; i++)
    {
      if(i < t_rest_init_) //0 ~ 0.05초 , 0 ~ 10 tick
      {
        temp_px(i) = (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2;
        temp_py(i) = (foot_step_support_frame_(current_step_number-2, 1) + foot_step_support_frame_(current_step_number-1, 1))/2;
      }
      else if(i >= t_rest_init_ && i < t_rest_init_ + t_double1_) //0.05 ~ 0.15초 , 10 ~ 30 tick
      {
        temp_px(i) = (foot_step_support_frame_(current_step_number-2, 0) + foot_step_support_frame_(current_step_number-1, 0))/2 + Kx/t_double1_*(i+1 - t_rest_init_);
        temp_py(i) = (foot_step_support_frame_(current_step_number-2, 1) + foot_step_support_frame_(current_step_number-1, 1))/2 + Ky/t_double1_*(i+1 - t_rest_init_);
      }
      else if(i >= t_rest_init_ + t_double1_ && i < t_total_ - t_rest_last_ - t_double2_) //0.15 ~ 1.05초 , 30 ~ 210 tick
      {
        temp_px(i) = foot_step_support_frame_offset_(current_step_number-1, 0);
        temp_py(i) = foot_step_support_frame_offset_(current_step_number-1, 1);
      }
      else if(i >= t_total_ - t_rest_last_ - t_double2_ && i < t_total_ - t_rest_last_) //1.05 ~ 1.15초 , 210 ~ 230 tick 
      {
        temp_px(i) = foot_step_support_frame_offset_(current_step_number-1, 0) + Kx2/t_double2_*(i+1 - (t_total_ - t_rest_last_ - t_double2_));
        temp_py(i) = foot_step_support_frame_offset_(current_step_number-1, 1) + Ky2/t_double2_*(i+1 - (t_total_ - t_rest_last_ - t_double2_));
      }
      else //1.15 ~ 1.2초 , 230 ~ 240 tick
      {
        temp_px(i) = temp_px(i-1);
        temp_py(i) = temp_py(i-1);
      }
      
    }
  }
  
}
*/
void WalkingController::getComTrajectory_MJ()
{ 
  xs_(0) = com_support_current_(0); ys_(0) = com_support_current_(1);
    
  modifiedPreviewControl_MJ();  
  
  xs_(1) = xd_(1); ys_(1) = yd_(1);  
  xs_(2) = xd_(2); ys_(2) = yd_(2);
  
  // CPM Preview
  //xs_(1) = X_bar_p_(1)/0.005; ys_(1) = Y_bar_p_(1)/0.005;
  //xs_(2) = X_bar_p_(2)/0.005; ys_(2) = Y_bar_p_(2)/0.005;

  if (walking_tick_ == t_start_ + t_total_-1 && current_step_num_ != total_step_num_-1) // 지지 발이 바뀌기 1tick 전 
  { // 지지 발이 바뀌기 1tick 전에 아래의 작업을 하는 이유는 지지 발이 바뀌는 순간에도 이전 tick의 CoM desired를 써야하는데 이전 tick의 CoM desired는 이전 지지 발에서 생성한 것이기 때문.
    Eigen::Vector3d com_pos_prev;
    Eigen::Vector3d com_pos;
    Eigen::Vector3d com_vel_prev;
    Eigen::Vector3d com_vel;
    Eigen::Vector3d com_acc_prev;
    Eigen::Vector3d com_acc;

    Eigen::Matrix3d temp_rot;
    Eigen::Vector3d temp_pos;
    
    temp_rot = DyrosMath::rotateWithZ(-foot_step_support_frame_(current_step_num_,5)); // 회전
    for(int i=0; i<3; i++)
      temp_pos(i) = foot_step_support_frame_(current_step_num_,i);

    com_pos_prev(0) = xs_(0);
    com_pos_prev(1) = ys_(0);
    com_pos = temp_rot*(com_pos_prev - temp_pos);
     
    com_vel_prev(0) = xs_(1);
    com_vel_prev(1) = ys_(1);
    com_vel_prev(2) = 0.0;
    com_vel = temp_rot*com_vel_prev;

    com_acc_prev(0) = xs_(2);
    com_acc_prev(1) = ys_(2);
    com_acc_prev(2) = 0.0;
    com_acc = temp_rot*com_acc_prev;
    // xs_ => com의 위치, 속도, 가속도
    xs_(0) = com_pos(0);
    ys_(0) = com_pos(1);
    xs_(1) = com_vel(0);
    ys_(1) = com_vel(1);
    xs_(2) = com_acc(0);
    ys_(2) = com_acc(1);

    //Preview Actual
    preview_x(0) = xs_(0); 
    preview_y(0) = ys_(0); 
  }

  double start_time;

  if(current_step_num_ == 0)
    start_time = 0;
  else
    start_time = t_start_;

  if(com_control_mode_ == true)
  { 
    /*   
    if(walking_tick_ <= 840)
    {
      com_desired_(0) = ref_com_(walking_tick_,0);
      com_desired_(1) = ref_com_(walking_tick_,1);
      Com_tra_graph << com_desired_(0) << "," << com_desired_(1)  << std::endl;
    }    
    else if(walking_tick_ > 840 )
    { 
      if(com_tick_ % 240 == 0 )//&& print_flag != 10)
      {
        com_tick_ = 0;
        print_flag ++;
        if(print_flag % (int)total_step_num_ == 0)
        {com_tick_ = 240; print_flag = 0;}
      }
      com_desired_(0) = ref_com_(com_tick_,0);
      com_desired_(1) = ref_com_(com_tick_,1);
      Com_tra_graph << com_desired_(0) << "," << com_desired_(1)  << std::endl;
      com_tick_ ++;      
    }*/
    /*
    if(walking_tick_ <= 840) 
    {
      com_desired_(0) = ref_com_(walking_tick_,0);
      if(walking_tick_ >= 600)
      {
        com_desired_(0) = 0;//ref_com_(walking_tick_,0);
      }
      com_desired_(1) = ref_com_(walking_tick_,1);
      
      Com_tra_graph << com_desired_(0) << "," << com_desired_(1) + 0.129510  << "," << com_support_current_(0) << "," << com_support_current_(1) + 0.12960 << std::endl;
    }    
    else if(walking_tick_ > 840 )
    { 
      if(com_tick_ % 240 == 0)//&& print_flag != 10)
      {
        com_tick_ = 0;
        //print_flag ++;
        //if(print_flag % (int)total_step_num_ == 0)
        //{com_tick_ = 240; print_flag = 0;}
      }
      com_desired_(0) = 0;//ref_com_(com_tick_,0);
      com_desired_(1) = ref_com_(com_tick_,1);
      if((int)current_step_num_ % 2 == 0)
      {
        Com_tra_graph << com_desired_(0) << "," << com_desired_(1) + 0.12707  << "," << com_support_current_(0) << "," << com_support_current_(1) + 0.12707 << std::endl;
      }
      else if((int)current_step_num_ % 2 == 1)
      {
        Com_tra_graph << com_desired_(0) << "," << com_desired_(1) - 0.12707  << "," << com_support_current_(0) << "," << com_support_current_(1) - 0.12707 << std::endl;
      }
      //Com_tra_graph << com_desired_(0) << "," << com_desired_(1) << "," << com_support_current_(0) << "," << com_support_current_(1) << std::endl;
      com_tick_ ++;      
    }*/
    //CPM Preview
    //com_desired_(0) = XD_(0); 
    //com_desired_(1) = YD_(0);
    com_desired_(0) = xd_(0); // xd_, yd_ 대신 xs_, ys_를 쓰면 지지 발이 바뀌는 순간 이전 tick의 CoM desired를 바라보는 시점이 바뀌지 않아서 안됨.
    com_desired_(1) = yd_(0);
    com_desired_(2) = DyrosMath::cubic(walking_tick_, t_start_, t_start_real_, pelv_support_init_.translation()(2), pelv_suppprt_start_.translation()(2), 0, 0);
    
  }
  else
  {
    com_desired_(0) = xd_(0);
    com_desired_(1) = yd_(0);
    com_desired_(2) = DyrosMath::cubic(walking_tick_, t_start_, t_start_real_, pelv_support_init_.translation()(2), pelv_suppprt_start_.translation()(2), 0, 0);
  }
}

void WalkingController::getPelvTrajectory()
{
  double z_rot = foot_step_support_frame_(current_step_num_,5);

  //Trunk Position
  if(com_control_mode_ == true)
  {
    double kp = 1.50;
    if(estimator_flag_ == false || l_ft_(2)+r_ft_(2)<0.7*51*9.81)
    {
      pelv_trajectory_support_.translation()(0) = pelv_support_current_.translation()(0) + kp*(com_desired_(0) - com_support_current_(0));
      pelv_trajectory_support_.translation()(1) = pelv_support_current_.translation()(1) + kp*(com_desired_(1) - com_support_current_(1));
    }
    else
    {
      pelv_trajectory_support_.translation()(0) = pelv_support_current_.translation()(0) + kp*(com_desired_(0) - X_hat_post_2_(0));
      pelv_trajectory_support_.translation()(1) = pelv_support_current_.translation()(1) + kp*(com_desired_(1) - X_hat_post_2_(1));
    }
    pelv_trajectory_support_.translation()(2) = com_desired_(2); //_T_Trunk_support.translation()(2) + kp*(_COM_desired(2) - _COM_real_support(2));
  }  

  //Trunk orientation
  Eigen::Vector3d Trunk_trajectory_euler;
  //Trunk_trajectory_euler.setZero();
  if(walking_tick_ < t_start_real_+t_double1_)
  {
    for(int i=0; i<2; i++)
      Trunk_trajectory_euler(i) = DyrosMath::cubic(walking_tick_,t_start_,t_start_real_+t_double1_,pelv_support_euler_init_(i),0.0,0.0,0.0);
    Trunk_trajectory_euler(2) = pelv_support_euler_init_(2);
  }
  else if(walking_tick_ >= t_start_real_+t_double1_ && walking_tick_ < t_start_+t_total_-t_double2_-t_rest_last_)
  {
    for(int i=0; i<2; i++)
      Trunk_trajectory_euler(i) = 0.0;

    if(foot_step_(current_step_num_,6) == 2)
      Trunk_trajectory_euler(2) = pelv_support_euler_init_(2);
    else
      Trunk_trajectory_euler(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_double2_-t_rest_last_, pelv_support_euler_init_(2),z_rot/2.0,0.0,0.0);
  }
  else
  {
    for(int i=0; i<2; i++)
      Trunk_trajectory_euler(i) = 0.0;

    if(foot_step_(current_step_num_,6) == 2)
      Trunk_trajectory_euler(2) = pelv_support_euler_init_(2);
    else
      Trunk_trajectory_euler(2) = z_rot/2.0;
  }  

  pelv_trajectory_support_.linear() = DyrosMath::rotateWithZ(Trunk_trajectory_euler(2))*DyrosMath::rotateWithY(Trunk_trajectory_euler(1))*DyrosMath::rotateWithX(Trunk_trajectory_euler(0));

}

void WalkingController::getFootTrajectory_MJ()
{
  Eigen::Vector6d target_swing_foot;

  for(int i=0; i<6; i++)
  { // 지지 발 기준에서 Foot step들을 봤을 때 current_step_num_은 다음 발의 좌표 뜻 함. Swing 발이 가야하는 좌표를 뜻 함.
    target_swing_foot(i) = foot_step_support_frame_(current_step_num_,i); 
  }
    
  if(walking_tick_ < t_start_real_ + t_double1_) // 0~630, 841~870 ... 30 tick , 0.15초 // DSP 구간
  { 
    lfoot_trajectory_support_.translation() = lfoot_support_init_.translation(); // t_start일때 Desired trajectory의 시작을 FK풀어서 나온값의 초기값으로 바꿔줌.
    lfoot_trajectory_dot_support_.setZero(); // 지지발에서 본 왼발 속도의 시작은 0.
    rfoot_trajectory_support_.translation() = rfoot_support_init_.translation();
    rfoot_trajectory_dot_support_.setZero();
    rfoot_trajectory_euler_support_ = rfoot_support_euler_init_;
    lfoot_trajectory_euler_support_ = lfoot_support_euler_init_;

    if(foot_step_(current_step_num_,6) == 1) // foot_step_(, 6) == 1은 왼발이 지지발.
    { 
      lfoot_trajectory_support_.translation().setZero();
      lfoot_trajectory_euler_support_.setZero();
    }      
    else // 오른 발이 지지 발인 경우
    {
      if(current_step_num_ == 0)
      {
        lfoot_trajectory_support_.translation()(2) = lfoot_support_init_.translation()(2);
      }
      else
      {        
        if(current_step_num_ >= 2)
        {
          lfoot_trajectory_support_.translation()(0) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,0),lfoot_support_init_.translation()(0),0.0,0.0);
          lfoot_trajectory_support_.translation()(1) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,1),lfoot_support_init_.translation()(1),0.0,0.0);
          lfoot_trajectory_euler_support_(2) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,5),lfoot_support_euler_init_(2),0.0,0.0);
          lfoot_trajectory_dot_support_(5) = DyrosMath::cubicDot(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,5),lfoot_support_euler_init_(2),0.0,0.0,hz_);        
        }
        else if(current_step_num_ >= 1)
        {
          lfoot_trajectory_support_.translation()(0) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,supportfoot_support_init_(0),lfoot_support_init_.translation()(0),0.0,0.0);
          lfoot_trajectory_support_.translation()(1) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,supportfoot_support_init_(1),lfoot_support_init_.translation()(1),0.0,0.0);
          lfoot_trajectory_euler_support_(2) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,supportfoot_support_init_(5),lfoot_support_euler_init_(2),0.0,0.0);
          lfoot_trajectory_dot_support_(5) = DyrosMath::cubicDot(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,supportfoot_support_init_(5),lfoot_support_euler_init_(2),0.0,0.0,hz_);
        }
        lfoot_trajectory_support_.translation()(2) = 0;
        lfoot_trajectory_euler_support_(0) = 0;
        lfoot_trajectory_euler_support_(1) = 0;
      }
    }   

    if(foot_step_(current_step_num_,6) == 0) // foot_step_(, 6) == 0은 오른발이 지지발.
    { 
      rfoot_trajectory_support_.translation().setZero();
      rfoot_trajectory_euler_support_.setZero(); 
    }
    else // 왼발이 지지 발인 경우, 왼발 기준으로 오른발을 보면 FK를 풀어서 나온 값을 초기값으로 바꿔주고 0으로 바꿔줌. 스윙발을 지지발 기준의 높이로 맞춰줌.
    {
      if(current_step_num_ == 0)
      {
        rfoot_trajectory_support_.translation()(2) = rfoot_support_init_.translation()(2);
      }        
      else
      {
        if(current_step_num_ >=2 ) // 오른발이 먼저 나가니까 current_step_num_ 이 1일때 고려 안했음.
        {
          rfoot_trajectory_support_.translation()(0) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,0),rfoot_support_init_.translation()(0),0.0,0.0);
          rfoot_trajectory_support_.translation()(1) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,1),rfoot_support_init_.translation()(1),0.0,0.0);
          rfoot_trajectory_euler_support_(2) = DyrosMath::cubic(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,5),rfoot_support_euler_init_(2),0.0,0.0);
          rfoot_trajectory_dot_support_(5) = DyrosMath::cubicDot(walking_tick_,t_start_,t_start_ + t_double1_ + t_rest_init_,foot_step_support_frame_(current_step_num_-2,5),rfoot_support_euler_init_(2),0.0,0.0,hz_);        
        }
      
       rfoot_trajectory_support_.translation()(2) = 0;
       rfoot_trajectory_euler_support_(0) = 0;
       rfoot_trajectory_euler_support_(1) = 0;
      }
    }

    lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
    rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
  
  }
  
  else if(walking_tick_ >= t_start_real_+t_double1_ && walking_tick_ < t_start_+t_total_-t_double2_-t_rest_last_) // 30 ~ 209 tick , 0.9초
  {  
    double t_rest_temp = 0.05*hz_;

    if(foot_step_(current_step_num_,6) == 1) // 왼발이 지지발일때, 지지발은 고정, 오른발은 목표 위치로 스윙
    {
      lfoot_trajectory_support_.translation() = lfoot_support_init_.translation(); // 왼발이 지지 발이니 전부 0임.            
      lfoot_trajectory_euler_support_.setZero(); // 지지발의 발의 Orientation (Roll, Pitch, Yaw)를 0으로 계속 보내고 있음.(자세 유지, 회전할때도 마찬가지)
      // 지지발좌표계에서 지지발의 발끝의 Orientation이든 XYZ든 전부 0으로 보내는게 맞다. 발은 고정하되 Pelvis가 움직이는것.
      lfoot_trajectory_dot_support_.setZero(); // 지지 발 속도 0.
      lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
      // 왼 발(지지발)의 Roll, Pitch, Yaw를 Rotation matrix로 바꿔서 lfoot_trajectory_support의 Homogeneous t matrix를 정의함. (Rotation 부분은 I3 matrix)
      ///여기까지 지지발 정의.
      if(walking_tick_ < t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0) // 스윙발 시작
      { // 90 tick 동안 발 높이 만큼 든다.
        rfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0,foot_height_,0.0,0.0);
        rfoot_trajectory_dot_support_(2) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0,foot_height_,0.0,0.0,hz_);
        
        rfoot_trajectory_euler_support_(1) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0.0,0.0,0.0,0.0);
        rfoot_trajectory_dot_support_(4) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0.0,0.0,0.0,0.0,hz_);
        // Ankle Pitch 궤적을 0도로.     
      }  
      else
      { // 0.6 ~ 1.05
        rfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-t_rest_temp,foot_height_,target_swing_foot(2),0.0,0.0);
        rfoot_trajectory_dot_support_(2) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-t_rest_temp,foot_height_,target_swing_foot(2),0.0,0.0,hz_);
        //target_swing_foot(2) = 0 임.
        rfoot_trajectory_euler_support_(1) = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-t_rest_temp,t_start_+t_total_-t_rest_last_,0.0,0.0,0.0,0.0);
        rfoot_trajectory_dot_support_(4) = DyrosMath::cubicDot(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-t_rest_temp,t_start_+t_total_-t_rest_last_,0.0,0.0,0.0,0.0,hz_);
      } 
      // X방향인 Roll 부분도 필요하면 시간에 나눠서 궤적안에 넣어줘도 됨. Ankle Pitch 처럼.
      rfoot_trajectory_euler_support_(0) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,0.0,target_swing_foot(0+3),0.0,0.0);
      rfoot_trajectory_dot_support_(0+3) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,0.0,target_swing_foot(0+3),0.0,0.0,hz_);
      // target_swing_foot(3) = 0 임.
      for(int i=0; i<2; i++) // X, Y 방향 궤적 생성, 여기서 불연속 문제를 회복해야함. 스윙발이 완벽히 따라가지 못하고 지지발이 되었을때, 현재 지지발 기준으로 봤을때 스윙발이 오차가 생긴것이기 때문에 3차곡선으로 이어줘야함. 
      {
        rfoot_trajectory_support_.translation()(i) = DyrosMath::cubic(walking_tick_,t_start_real_ + t_double1_ + 2*t_rest_temp, t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-2*t_rest_temp,rfoot_support_init_.translation()(i),target_swing_foot(i),0.0,0.0);
        rfoot_trajectory_dot_support_(i) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+2*t_rest_temp,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-2*t_rest_temp,rfoot_support_init_.translation()(i),target_swing_foot(i),0.0,0.0,hz_);
      } 

      rfoot_trajectory_euler_support_(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,rfoot_support_euler_init_(2),target_swing_foot(5),0.0,0.0);
      rfoot_trajectory_dot_support_(5) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,rfoot_support_euler_init_(2),target_swing_foot(5),0.0,0.0,hz_);
      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
    }
    else if(foot_step_(current_step_num_,6) == 0) // 오른발이 지지발일때, 지지발은 고정, 왼발은 목표 위치로 스윙
    { // 왼쪽 발도 똑같이.

      rfoot_trajectory_support_.translation() = rfoot_support_init_.translation(); 
      rfoot_trajectory_euler_support_.setZero(); 
      rfoot_trajectory_dot_support_.setZero();
      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
 
      if(walking_tick_ < t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0)
      {
        lfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0,foot_height_,0.0,0.0);
        lfoot_trajectory_dot_support_(2) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0,foot_height_,0.0,0.0,hz_);

        lfoot_trajectory_euler_support_(1) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0.0,0.0,0.0,0.0);
        lfoot_trajectory_dot_support_(4) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+t_rest_temp,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,0.0,0.0,0.0,0.0,hz_);
      }
      else
      {
        lfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-t_rest_temp,foot_height_,target_swing_foot(2),0.0,0.0);
        lfoot_trajectory_dot_support_(2) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_-t_imp_)/2.0,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-t_rest_temp,foot_height_,target_swing_foot(2),0.0,0.0,hz_);
        
        lfoot_trajectory_euler_support_(1) = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-t_rest_temp,t_start_+t_total_-t_rest_last_,0.0,0.0,0.0,0.0);
        lfoot_trajectory_dot_support_(4) = DyrosMath::cubicDot(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-t_rest_temp,t_start_+t_total_-t_rest_last_,0.0,0.0,0.0,0.0,hz_);
      }
      lfoot_trajectory_euler_support_(0) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,0.0,target_swing_foot(0+3),0.0,0.0);
      lfoot_trajectory_dot_support_(0+3) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,0.0,target_swing_foot(0+3),0.0,0.0,hz_);

      for(int i=0; i<2; i++)
      {
        lfoot_trajectory_support_.translation()(i) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+2*t_rest_temp,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-2*t_rest_temp,lfoot_support_init_.translation()(i),target_swing_foot(i),0.0,0.0);
        lfoot_trajectory_dot_support_(i) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_+2*t_rest_temp,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_-2*t_rest_temp,lfoot_support_init_.translation()(i),target_swing_foot(i),0.0,0.0,hz_);
      }  

      lfoot_trajectory_euler_support_(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,lfoot_support_euler_init_(2),target_swing_foot(5),0.0,0.0);
      lfoot_trajectory_dot_support_(5) = DyrosMath::cubicDot(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_-t_imp_,lfoot_support_euler_init_(2),target_swing_foot(5),0.0,0.0,hz_);
      
      lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
    } 
  }
  else // 210 ~ 239 tick , 0.15초 
  { 
    if(foot_step_(current_step_num_,6) == 1) // 왼쪽 발 지지
    {
      //lfoot_trajectory_support_.translation() = lfoot_support_init_.translation();
      lfoot_trajectory_euler_support_.setZero();
      lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
      lfoot_trajectory_dot_support_.setZero();

      for(int i=0; i<3; i++)
      {
        rfoot_trajectory_support_.translation()(i) = target_swing_foot(i);
        rfoot_trajectory_euler_support_(i) = target_swing_foot(i+3);
      }
      rfoot_trajectory_dot_support_.setZero();
      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
    }
    else if (foot_step_(current_step_num_,6) == 0) 
    {
      //rfoot_trajectory_support_.translation() = rfoot_support_init_.translation();
      rfoot_trajectory_euler_support_.setZero();
      rfoot_trajectory_dot_support_.setZero();
      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));

      for(int i=0; i<3; i++)
      {
        lfoot_trajectory_support_.translation()(i) = target_swing_foot(i);
        lfoot_trajectory_euler_support_(i) = target_swing_foot(i+3);
      }
      lfoot_trajectory_dot_support_.setZero();
      lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
    }
  }
  
  Swing_tra_graph << rfoot_trajectory_support_.translation()(0) << "," << rfoot_trajectory_support_.translation()(1) << "," << rfoot_trajectory_support_.translation()(2) 
  << "," << lfoot_trajectory_support_.translation()(0) << "," << lfoot_trajectory_support_.translation()(1) << "," << lfoot_trajectory_support_.translation()(2) << 
  "," << lfoot_trajectory_euler_support_(0) << "," << lfoot_trajectory_euler_support_(1) << "," << lfoot_trajectory_euler_support_(2) << "," << rfoot_trajectory_euler_support_(0) << "," << rfoot_trajectory_euler_support_(1) << "," << rfoot_trajectory_euler_support_(2) 
  << std::endl;  
}

void WalkingController::supportToFloatPattern()
{
  if(gyro_frame_flag_ == true)
  {
    Eigen::Isometry3d reference = pelv_trajectory_float_;
    DyrosMath::floatGyroframe(pelv_trajectory_support_,reference,pelv_trajectory_float_);
    DyrosMath::floatGyroframe(lfoot_trajectory_support_,reference,lfoot_trajectory_float_);
    DyrosMath::floatGyroframe(rfoot_trajectory_support_,reference,rfoot_trajectory_float_);
    lfoot_trajectory_euler_float_ = DyrosMath::rot2Euler(lfoot_trajectory_float_.linear());
    rfoot_trajectory_euler_float_ = DyrosMath::rot2Euler(rfoot_trajectory_float_.linear());
  }
  else
  {
    pelv_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*pelv_trajectory_support_;
    lfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*lfoot_trajectory_support_;
    rfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*rfoot_trajectory_support_;
    lfoot_trajectory_euler_float_ = DyrosMath::rot2Euler(lfoot_trajectory_float_.linear());
    rfoot_trajectory_euler_float_ = DyrosMath::rot2Euler(rfoot_trajectory_float_.linear());
   }
}

void WalkingController::computeIkControl_MJ(Eigen::Isometry3d float_trunk_transform, Eigen::Isometry3d float_lleg_transform, Eigen::Isometry3d float_rleg_transform, Eigen::Vector12d& q_des)
{
  //float = World/ trunk = pelvis
  // 명주 정리 (KAJITA 책 <-> Code 구성)
  // float_trunk_transform.rotation() : float에서 바라본 pelvis rotation -> R1
  // float_trunk_transform.translation() : float에서 바라본 pelvis 좌표 -> P1
  // float_rleg_transform.rotation() : float에서 바라본 발끝 rotation -> R7 
  // float_rleg_transform.translation() : float에서 바라본 발끝 좌표 -> P7
  // float_trunk_transform.translation() + float_trunk_transform.rotation()*D  : float에서 바라본 pelvis 좌표 + float 좌표계로 변환 * pelvis 좌표계에서 바라본 pelvis ~ hip까지 거리-> P2
  
  // R7.transpose * (P2 - P7) , P2 = P1 + R1*D
  double offset_hip_pitch = 24.0799945102*DEG2RAD;
  double offset_knee_pitch = 14.8197729791*DEG2RAD;
  double offset_ankle_pitch = 9.2602215311*DEG2RAD;

  Eigen::Vector3d R_r, R_D, L_r, L_D ;

  L_D << 0 , +0.105, -0.1119;
  R_D << 0 , -0.105, -0.1119;
  
  L_r = float_lleg_transform.rotation().transpose() * (float_trunk_transform.translation() + float_trunk_transform.rotation()*L_D - float_lleg_transform.translation());
  R_r = float_rleg_transform.rotation().transpose() * (float_trunk_transform.translation() + float_trunk_transform.rotation()*R_D - float_rleg_transform.translation());
  
  double R_C = 0, L_C = 0, L_upper = 0.3713, L_lower = 0.3728 , R_alpha = 0, L_alpha = 0;

  L_C = sqrt( pow(L_r(0),2) + pow(L_r(1),2) + pow(L_r(2),2) );
  R_C = sqrt( pow(R_r(0),2) + pow(R_r(1),2) + pow(R_r(2),2) );
  
  q_des(3) = (-acos((pow(L_upper,2) + pow(L_lower,2) - pow(L_C,2)) / (2*L_upper*L_lower))+ M_PI) ;
  q_des(9) = (-acos((pow(L_upper,2) + pow(L_lower,2) - pow(R_C,2)) / (2*L_upper*L_lower))+ M_PI) ;
  L_alpha = asin(L_upper / L_C * sin(M_PI - q_des(3)));
  R_alpha = asin(L_upper / R_C * sin(M_PI - q_des(9)));

  q_des(4)  = -atan2(L_r(0), sqrt(pow(L_r(1),2) + pow(L_r(2),2))) - L_alpha ;
  q_des(10) = -atan2(R_r(0), sqrt(pow(R_r(1),2) + pow(R_r(2),2))) - R_alpha ;
  
  // trunk_lleg_rotation -> R1.transpose * R7 
  // Ryaw * Rroll * Rpitch = R1.transpose * R7 * ~ 
  Eigen::Matrix3d R_Knee_Ankle_Y_rot_mat, L_Knee_Ankle_Y_rot_mat;
  Eigen::Matrix3d R_Ankle_X_rot_mat, L_Ankle_X_rot_mat;
  Eigen::Matrix3d R_Hip_rot_mat, L_Hip_rot_mat;

  L_Knee_Ankle_Y_rot_mat = DyrosMath::rotateWithY(-q_des(3)-q_des(4));
  L_Ankle_X_rot_mat = DyrosMath::rotateWithX(-q_des(5));
  R_Knee_Ankle_Y_rot_mat = DyrosMath::rotateWithY(-q_des(9)-q_des(10));
  R_Ankle_X_rot_mat = DyrosMath::rotateWithX(-q_des(11)); 
  
  L_Hip_rot_mat.setZero(); R_Hip_rot_mat.setZero();

  L_Hip_rot_mat = float_trunk_transform.rotation().transpose() * float_lleg_transform.rotation() * L_Ankle_X_rot_mat * L_Knee_Ankle_Y_rot_mat; 
  R_Hip_rot_mat = float_trunk_transform.rotation().transpose() * float_rleg_transform.rotation() * R_Ankle_X_rot_mat * R_Knee_Ankle_Y_rot_mat;

  q_des(0) = -atan2(-L_Hip_rot_mat(0,1),L_Hip_rot_mat(1,1)); // Hip yaw
  q_des(1) =  atan2(L_Hip_rot_mat(2,1), -L_Hip_rot_mat(0,1) * sin(q_des(0)) + L_Hip_rot_mat(1,1)*cos(q_des(0))); // Hip roll
  q_des(2) =  atan2(-L_Hip_rot_mat(2,0), L_Hip_rot_mat(2,2)) + offset_hip_pitch; // Hip pitch
  q_des(3) =  q_des(3) - 14.8197729791*DEG2RAD; // Knee pitch
  q_des(4) =  q_des(4) - 9.2602215311*DEG2RAD; // Ankle pitch
  q_des(5) =  atan2( L_r(1), L_r(2) ); // Ankle roll

  q_des(6) = -atan2(-R_Hip_rot_mat(0,1),R_Hip_rot_mat(1,1));
  q_des(7) =  atan2(R_Hip_rot_mat(2,1), -R_Hip_rot_mat(0,1) * sin(q_des(6)) + R_Hip_rot_mat(1,1)*cos(q_des(6)));
  q_des(8) = -atan2(-R_Hip_rot_mat(2,0), R_Hip_rot_mat(2,2)) - offset_hip_pitch;
  q_des(9) = -q_des(9) + 14.8197729791*DEG2RAD;
  q_des(10) = -q_des(10) + 9.2602215311*DEG2RAD; 
  q_des(11) =  atan2( R_r(1), R_r(2) );
  
  Des_q_tra_graph << q_des(0) << "," << q_des(1) << "," << q_des(2) << "," << q_des(3) << "," << q_des(4) << "," << q_des(5) << "," << DOB_IK_output_b_(0) << "," << DOB_IK_output_b_(1) << "," << DOB_IK_output_b_(2) << "," << DOB_IK_output_b_(3) << "," << DOB_IK_output_b_(4) << "," << DOB_IK_output_b_(5) << std::endl;
}

void WalkingController::computeJacobianControl(Eigen::Isometry3d float_lleg_transform, Eigen::Isometry3d float_rleg_transform, Eigen::Vector3d float_lleg_transform_euler, Eigen::Vector3d float_rleg_transform_euler, Eigen::Vector12d& desired_leg_q_dot)
{

  Eigen::Matrix6d jacobian_temp_l, jacobian_temp_r, current_leg_jacobian_l_inv, current_leg_jacobian_r_inv,
      J_damped, I_matrix;
  double wl, wr, w0, lambda, a;
  w0 = 0.001;
  lambda = 0.05;
  jacobian_temp_l=current_leg_jacobian_l_*current_leg_jacobian_l_.transpose();
  jacobian_temp_r=current_leg_jacobian_r_*current_leg_jacobian_r_.transpose();
  wr = sqrt(jacobian_temp_r.determinant());
  wl = sqrt(jacobian_temp_l.determinant());

  if (wr<=w0)
  { //Right Jacobi
    a = lambda * pow(1-wr/w0,2);
    J_damped = current_leg_jacobian_r_.transpose()*current_leg_jacobian_r_+a*Eigen::Matrix6d::Identity();
    J_damped = J_damped.inverse();

    cout << "Singularity Region of right leg: " << wr << endl;
    current_leg_jacobian_r_inv = J_damped*current_leg_jacobian_r_.transpose();
  }
  else
  {
    current_leg_jacobian_r_inv = (current_leg_jacobian_r_.transpose()*current_leg_jacobian_r_).inverse()*current_leg_jacobian_r_.transpose();
  }

  if (wl<=w0)
  {
    a = lambda*pow(1-wl/w0,2);
    J_damped = current_leg_jacobian_l_.transpose()*current_leg_jacobian_l_+a*Eigen::Matrix6d::Identity();
    J_damped = J_damped.inverse();

    cout << "Singularity Region of right leg: " << wr << endl;
    current_leg_jacobian_l_inv = J_damped*current_leg_jacobian_l_.transpose();
  }
  else
  {
    current_leg_jacobian_l_inv = (current_leg_jacobian_l_.transpose()*current_leg_jacobian_l_).inverse()*current_leg_jacobian_l_.transpose();
  }

  Eigen::Matrix6d kp; // for setting CLIK gains
  kp.setZero();
  kp(0,0) = 160;
  kp(1,1) = 160;
  kp(2,2) = 160;
  kp(3,3) = 100;
  kp(4,4) = 100;
  kp(5,5) = 100;

  Eigen::Vector6d lp, rp, cubic_xr, cubic_xl;
  lp.setZero(); rp.setZero(), cubic_xr.setZero(), cubic_xl.setZero();
  lp.topRows<3>() = (-lfoot_float_current_.translation() + float_lleg_transform.translation());
  rp.topRows<3>() = (-rfoot_float_current_.translation() + float_rleg_transform.translation());

  for(int i=0;i<3;i++)
  {
    cubic_xl(i) = float_lleg_transform.translation()(i);
    cubic_xl(i+3) = float_lleg_transform_euler(i);
  }

  for(int i=0;i<3;i++)
  {
    cubic_xr(i) = float_rleg_transform.translation()(i);
    cubic_xr(i+3) = float_rleg_transform_euler(i);
  }
  Eigen::Vector3d r_leg_phi, l_leg_phi;
  l_leg_phi = DyrosMath::legGetPhi(lfoot_float_current_, lfoot_float_init_, cubic_xl);
  r_leg_phi = DyrosMath::legGetPhi(rfoot_float_current_, rfoot_float_init_, cubic_xr);

  lp.bottomRows<3>() = - l_leg_phi;
  rp.bottomRows<3>() = - r_leg_phi;

  Eigen::Vector6d q_lfoot_dot,q_rfoot_dot;
  q_lfoot_dot=current_leg_jacobian_l_inv*kp*lp;
  q_rfoot_dot=current_leg_jacobian_r_inv*kp*rp;

  for (int i=0; i<6; i++)
  {
    desired_leg_q_dot(i+6) = q_rfoot_dot(i);
    desired_leg_q_dot(i) = q_lfoot_dot(i);
  }
}

void WalkingController::modifiedPreviewControl_MJ()
{
  /////reference: http://www.tandfonline.com/doi/pdf/10.1080/0020718508961156?needAccess=true/////////////

  if(walking_tick_ == 0) // 보행 시작할 때 딱 한번만 구하면 됨.
  {
   //previewParam_MJ(1.0/hz_, 16*hz_/10, k_ ,com_support_init_, gi_, gp_l_, gx_, a_, b_, c_);
   previewParam_MJ_Act(1.0/hz_, 16*hz_/10, K_act_ ,com_support_init_, Gi_, Gd_, Gx_, A_, B_, C_, D_, A_bar_, B_bar_);
   //previewParam_MJ_CPM(1.0/hz_, 16*hz_/10, K_ ,com_support_init_, Gi_, Gd_, Gx_, A_, B_, C_, D_, A_bar_, B_bar_);

   UX_ = com_support_init_(0);
   UY_ = com_support_init_(1);
   xs_(0) = xi_; xs_(1) = 0; xs_(2) = 0;
   ys_(0) = yi_; ys_(1) = 0; xs_(2) = 0;
  }

  if(current_step_num_ == 0)
    zmp_start_time_ = 0.0;
  else
    zmp_start_time_ = t_start_;
  
  //preview_MJ(1.0/hz_, 16*hz_/10, walking_tick_-zmp_start_time_, xi_, yi_, xs_, ys_, ux_, uy_, gi_, gp_l_, gx_, a_, b_, c_, xd_, yd_);  
  preview_MJ_Act(1.0/hz_, 16*hz_/10, walking_tick_-zmp_start_time_, xi_, yi_, xs_, ys_, UX_, UY_, Gi_, Gd_, Gx_, A_, B_, xd_, yd_);
  //preview_MJ_CPM(1.0/hz_, 16*hz_/10, walking_tick_-zmp_start_time_, xi_, yi_, xs_, ys_, UX_, UY_, Gi_, Gd_, Gx_, A_, B_, A_bar_, B_bar_, XD_, YD_, X_bar_p_, Y_bar_p_); 
}                                                                   

void WalkingController::previewParam_MJ_CPM(double dt, int NL, Eigen::Matrix3d& K, Eigen::Vector3d com_support_init_, Eigen::MatrixXd& Gi, Eigen::VectorXd& Gd, Eigen::MatrixXd& Gx, 
  Eigen::MatrixXd& A, Eigen::VectorXd& B, Eigen::MatrixXd& C, Eigen::MatrixXd& D, Eigen::MatrixXd& A_bar, Eigen::VectorXd& B_bar)
  { 
    //double Kp = 2098; double Kv = 32.2;
    double Kp = 10.35;
    double Kv = 2.2012;
    
    A.resize(2,2);
    A(0,0) = 1.0 - Kp*dt*dt*0.5;
    A(0,1) = dt - 0.5*Kv*dt*dt;
    A(1,0) = -Kp*dt;
    A(1,1) = 1.0 - Kv*dt;
    
    B.resize(2);
    B(0) = 0.5*Kp*dt*dt;
    B(1) = Kp*dt;
    
    C.resize(1,2);
    C(0,0) = 1 + zc_/GRAVITY*Kp;
    C(0,1) = zc_/GRAVITY*Kv;
    //std::cout << GRAVITY << "," << zc_/GRAVITY*Kp << std::endl;
    D.resize(1,1);
    D(0,0) = -zc_/GRAVITY*Kp;
    // A, B, C, D discrete matrix 정의

    B_bar.resize(3,1);
    B_bar(0,0) = D(0,0);
    B_bar(1,0) = B(0);
    B_bar(2,0) = B(1);
    
    Eigen::Matrix1x3d B_bar_tran;
    B_bar_tran = B_bar.transpose();
    
    Eigen::MatrixXd I_bar;
    Eigen::MatrixXd F_bar;
    A_bar.resize(3,3);
    I_bar.resize(3,1);
    F_bar.resize(3,2);
    F_bar.setZero();

    F_bar(0,0) = C(0,0);
    F_bar(0,1) = C(0,1);

    F_bar(1,0) = A(0,0);
    F_bar(1,1) = A(0,1);
    F_bar(2,0) = A(1,0);
    F_bar(2,1) = A(1,1);
    
    I_bar.setZero();
    I_bar(0,0) = 1.0;

    A_bar(0,0) = I_bar(0,0);
    A_bar(1,0) = I_bar(1,0);
    A_bar(2,0) = I_bar(2,0);

    A_bar(0,1) = F_bar(0,0);
    A_bar(0,2) = F_bar(0,1);
    A_bar(1,1) = F_bar(1,0);
    A_bar(1,2) = F_bar(1,1);
    A_bar(2,1) = F_bar(2,0);
    A_bar(2,2) = F_bar(2,1);
   
    Eigen::MatrixXd Qe;
    Qe.resize(1,1);
    Qe(0,0) = 1.0;

    Eigen::MatrixXd R;
    R.resize(1,1);
    R(0,0) = 0.000001;

    Eigen::MatrixXd Qx;
    Qx.resize(3,3);
    Qx.setZero();

    Eigen::MatrixXd Q_bar;
    Q_bar.resize(3,3);
    Q_bar.setZero();
    Q_bar(0,0) = Qe(0,0);
    
    //K = discreteRiccatiEquationLQR(A_bar, B_bar, R, Q_bar);
    /*K(0,0) = 110.006257148045; 
    K(0,1) = 5995.685176876064;  
    K(0,2) = 1618.928784700213; 
    K(1,1) = 329781.443215973151; 
    K(1,2) = 89046.148242507566;  
    K(2,2) = 24043.852923290622;
    
    K(1, 0) = K(0, 1);
    K(2, 0) = K(0, 2);
    K(2, 1) = K(1, 2);*/
    K(0,0) = 110.073395244516;
    K(0,1) = 6002.539472194625;
    K(0,2) = 1619.933800848508;
    K(1,1) = 330730.600140600640;
    K(1,2) = 89149.373787370059;
    K(2,2) = 24059.898910228723;
    K(1, 0) = K(0, 1);
    K(2, 0) = K(0, 2);
    K(2, 1) = K(1, 2);
    
    Eigen::MatrixXd Temp_mat;
    Eigen::MatrixXd Temp_mat_inv;
    Eigen::MatrixXd Ac_bar;
    Temp_mat.resize(1,1);
    Temp_mat.setZero();
    Temp_mat_inv.resize(1,1);
    Temp_mat_inv.setZero();
    Ac_bar.setZero();
    Ac_bar.resize(3,3);

    Temp_mat = R + B_bar_tran * K * B_bar;
    Temp_mat_inv = Temp_mat.inverse();
    //std::cout << R << ","<< B_bar_tran << "," << Temp_mat << "," << Temp_mat_inv << std::endl;
    Ac_bar = A_bar - B_bar * Temp_mat_inv * B_bar_tran * K * A_bar;
    
    Eigen::MatrixXd Ac_bar_tran(3,3);
    Ac_bar_tran = Ac_bar.transpose();

    //std::cout << Temp_mat <<"," << Temp_mat_inv << std::endl;
    Gi.resize(1,1); Gx.resize(1,2);
    Gi = Temp_mat_inv * B_bar_tran * K * I_bar ;
    Gx = Temp_mat_inv * B_bar_tran * K * F_bar ;   
    
    Eigen::MatrixXd X_bar;
    Eigen::Vector3d X_bar_col;
    X_bar.resize(3, NL); 
    X_bar.setZero();
    X_bar_col.setZero();
    X_bar_col = - Ac_bar_tran * K * I_bar;

    for(int i = 0; i < NL; i++)
    {
      X_bar.block<3,1>(0,i) = X_bar_col;
      X_bar_col = Ac_bar_tran*X_bar_col;
    }           

    Gd.resize(NL);
    Eigen::VectorXd Gd_col(1);
    Gd_col(0) = -Gi(0,0);
    
    for(int i = 0; i < NL; i++)
    {
      Gd.segment(i,1) = Gd_col;
      Gd_col = Temp_mat_inv * B_bar_tran * X_bar.col(i) ;
    }
     
  }
  
void WalkingController::preview_MJ_CPM(double dt, int NL, int tick, double x_i, double y_i, Eigen::Vector3d xs, Eigen::Vector3d ys, double& UX, double& UY, 
       Eigen::MatrixXd Gi, Eigen::VectorXd Gd, Eigen::MatrixXd Gx, Eigen::MatrixXd A, Eigen::VectorXd B, Eigen::MatrixXd A_bar, Eigen::VectorXd B_bar, Eigen::Vector2d &XD, Eigen::Vector2d &YD, Eigen::VectorXd& X_bar_p, Eigen::VectorXd& Y_bar_p)
{
  int zmp_size;
  zmp_size = ref_zmp_.col(1).size(); // 보행 중 720개 (240 * 3)
  Eigen::VectorXd px_ref, py_ref;
  px_ref.resize(zmp_size);
  py_ref.resize(zmp_size);
  
  for(int i = 0; i < zmp_size; i++)
  {
    px_ref(i) = ref_zmp_(i,0);
    py_ref(i) = ref_zmp_(i,1);
  }
  
  Eigen::Matrix1x3d C;
  C(0,0) = 1; C(0,1) = 0; C(0,2) = -zc_/GRAVITY;
  
  Eigen::VectorXd px, py;
  px.resize(1); py.resize(1);

  X_bar_p.resize(3); Y_bar_p.resize(3);

  if(tick == 0 && current_step_num_ == 0)
  { 
    Preview_X_b(0) = x_i; // 이때 before 값이 없으면 첫 tick에서 delta x의 에러가 갑작스럽게 생겨서 발산
    Preview_Y_b(0) = y_i;
    Preview_X(0) = x_i;
    Preview_Y(0) = y_i;
    Preview_X(1) = 0;
    Preview_Y(1) = 0;
    X_bar_p.setZero(); Y_bar_p.setZero();
  }
  else
  {  
    Preview_X(0) = xs(0); Preview_Y(0) = ys(0);    
    Preview_X(1) = xs(1); Preview_Y(1) = ys(1);
  }  
  
  Eigen::VectorXd Temp_mat_X, Temp_mat_Y;
  Temp_mat_X.resize(3); Temp_mat_Y.resize(3);
  Temp_mat_X.setZero(); Temp_mat_Y.setZero();  
    
  Temp_mat_X(0) = Preview_X(0); Temp_mat_Y(0) = Preview_Y(0); // preview_x(0)이 x_i를 안넣고 xs를??
  Temp_mat_X(2) = X_bar_p(2)/dt; Temp_mat_Y(2) = Y_bar_p(2)/dt;
  //Temp_mat_X(2) = com_support_ddot_current_(0); Temp_mat_Y(2) = com_support_ddot_current_(1); // 값 2개 완전히 같음.  
  
  px = C*Temp_mat_X;
  py = C*Temp_mat_Y;
   
  X_bar_p(0) = px(0) - px_ref(tick); //e(i) 정의
  Y_bar_p(0) = py(0) - py_ref(tick);
   
  double sum_Gd_px_ref = 0, sum_Gd_py_ref = 0;

  for(int i = 0; i < NL; i++) // Preview Step 개수만큼 더함.
  {
    sum_Gd_px_ref = sum_Gd_px_ref + Gd(i)*(px_ref(tick + 1 + i) - px_ref(tick + i));
    sum_Gd_py_ref = sum_Gd_py_ref + Gd(i)*(py_ref(tick + 1 + i) - py_ref(tick + i));
  }

  Eigen::MatrixXd temp; temp.resize(2, 1);
  Eigen::VectorXd GxX(1); Eigen::VectorXd GxY(1);
  temp(0, 0) = X_bar_p(1);
  temp(1, 0) = X_bar_p(2);
  GxX = Gx*temp;

  temp(0, 0) = Y_bar_p(1);
  temp(1, 0) = Y_bar_p(2);
  GxY = Gx*temp;

  Eigen::MatrixXd del_ux(1,1);
  Eigen::MatrixXd del_uy(1,1);
  del_ux.setZero();
  del_uy.setZero();

  del_ux(0,0) = -(X_bar_p(0) * Gi(0,0)) - GxX(0) - sum_Gd_px_ref;
  del_uy(0,0) = -(Y_bar_p(0) * Gi(0,0)) - GxY(0) - sum_Gd_py_ref;

  X_bar_p = A_bar*X_bar_p + B_bar*del_ux;
  Y_bar_p = A_bar*Y_bar_p + B_bar*del_uy;  

  UX = UX + del_ux(0,0);
  UY = UY + del_uy(0,0);
  
  Preview_X_b(1) = Preview_X(1); Preview_Y_b(1) = Preview_Y(1);  
 
  XD = A*Preview_X + B*UX;
  YD = A*Preview_Y + B*UY;
  
  Preview_X(1) = XD(1); Preview_Y(1) = YD(1);

  Com_tra_graph << XD(0) << "," << YD(0) << "," << px(0) << "," << py(0) << "," << px_ref(tick) << "," << py_ref(tick) << "," << X_bar_p(2) << "," << Y_bar_p(2) << std::endl;
       
}

void WalkingController::previewParam_MJ_Act(double dt, int NL, Eigen::Matrix4d& K, Eigen::Vector3d com_support_init_, Eigen::MatrixXd& Gi, Eigen::VectorXd& Gd, Eigen::MatrixXd& Gx, 
  Eigen::MatrixXd& A, Eigen::VectorXd& B, Eigen::MatrixXd& C, Eigen::MatrixXd& D, Eigen::MatrixXd& A_bar, Eigen::VectorXd& B_bar)
  { 
    A.resize(3,3);
    A(0,0) = 1.0;
    A(0,1) = dt;
    A(0,2) = dt*dt*0.5;
    A(1,0) = 0;
    A(1,1) = 1.0;
    A(1,2) = dt;
    A(2,0) = 0;
    A(2,1) = 0;
    A(2,2) = 1;
    
    B.resize(3);
    B(0) = dt*dt*dt/6;
    B(1) = dt*dt/2;
    B(2) = dt;
    
    C.resize(1,3);
    C(0,0) = 1;
    C(0,1) = 0;
    C(0,2) = -zc_/GRAVITY;

    B_bar.resize(4);    
    B_bar.segment(0,1) = C*B; 
    B_bar.segment(1,3) = B;
    
    Eigen::Matrix1x4d B_bar_tran;
    B_bar_tran = B_bar.transpose();
    
    Eigen::MatrixXd I_bar;
    Eigen::MatrixXd F_bar;
    A_bar.resize(4,4);
    I_bar.resize(4,1);
    F_bar.resize(4,3);
    F_bar.setZero();

    F_bar.block<1,3>(0,0) = C*A;
    F_bar.block<3,3>(1,0) = A;
    
    I_bar.setZero();
    I_bar(0,0) = 1.0;

    A_bar.block<4,1>(0,0) = I_bar;
    A_bar.block<4,3>(0,1) = F_bar;
   
    Eigen::MatrixXd Qe;
    Qe.resize(1,1);
    Qe(0,0) = 1.0;

    Eigen::MatrixXd R;
    R.resize(1,1);
    R(0,0) = 0.000001;

    Eigen::MatrixXd Qx;
    Qx.resize(3,3);
    Qx.setZero();

    Eigen::MatrixXd Q_bar;
    Q_bar.resize(3,3);
    Q_bar.setZero();
    Q_bar(0,0) = Qe(0,0);
    
    K(0,0) = 110.946733178638; 
    K(0,1) = 6099.115434920124;  
    K(0,2) = 1670.206808355727; 
    K(0,3) = 4.277193660725; 
    K(1,0) = K(0,1); 
    K(1,1) = 342635.571113037353; 
    K(1,2) = 93854.611649038387; 
    K(1,3) = 247.309059159967; 
    K(2,0) = K(0,2); 
    K(2,1) = K(1,2); 
    K(2,2) = 25708.928259919834; 
    K(2,3) = 67.827431998836; 
    K(3,0) = K(0,3); 
    K(3,1) = K(1,3); 
    K(3,2) = K(2,3); 
    K(3,3) = 0.202193033700;
  
    Eigen::MatrixXd Temp_mat;
    Eigen::MatrixXd Temp_mat_inv;
    Eigen::MatrixXd Ac_bar;
    Temp_mat.resize(1,1);
    Temp_mat.setZero();
    Temp_mat_inv.resize(1,1);
    Temp_mat_inv.setZero();
    Ac_bar.setZero();
    Ac_bar.resize(4,4);

    Temp_mat = R + B_bar_tran * K * B_bar;
    Temp_mat_inv = Temp_mat.inverse();
    
    Ac_bar = A_bar - B_bar * Temp_mat_inv * B_bar_tran * K * A_bar;
    
    Eigen::MatrixXd Ac_bar_tran(4,4);
    Ac_bar_tran = Ac_bar.transpose();
    
    Gi.resize(1,1); Gx.resize(1,3);
    Gi = Temp_mat_inv * B_bar_tran * K * I_bar ;
    Gx = Temp_mat_inv * B_bar_tran * K * F_bar ;   
    
    Eigen::MatrixXd X_bar;
    Eigen::Vector4d X_bar_col;
    X_bar.resize(4, NL); 
    X_bar.setZero();
    X_bar_col.setZero();
    X_bar_col = - Ac_bar_tran * K * I_bar;

    for(int i = 0; i < NL; i++)
    {
      X_bar.block<4,1>(0,i) = X_bar_col;
      X_bar_col = Ac_bar_tran*X_bar_col;
    }           

    Gd.resize(NL);
    Eigen::VectorXd Gd_col(1);
    Gd_col(0) = -Gi(0,0);
    
    for(int i = 0; i < NL; i++)
    {
      Gd.segment(i,1) = Gd_col;
      Gd_col = Temp_mat_inv * B_bar_tran * X_bar.col(i) ;
    }
     
  }
  
void WalkingController::preview_MJ_Act(double dt, int NL, int tick, double x_i, double y_i, Eigen::Vector3d xs, Eigen::Vector3d ys, double& UX, double& UY, 
       Eigen::MatrixXd Gi, Eigen::VectorXd Gd, Eigen::MatrixXd Gx, Eigen::MatrixXd A, Eigen::VectorXd B, Eigen::Vector3d &XD, Eigen::Vector3d &YD)
{
  int zmp_size;
  zmp_size = ref_zmp_.col(1).size(); // 보행 중 720개 (240 * 3)
  Eigen::VectorXd px_ref, py_ref;
  px_ref.resize(zmp_size);
  py_ref.resize(zmp_size);
  
  for(int i = 0; i < zmp_size; i++)
  {
    px_ref(i) = ref_zmp_(i,0);
    py_ref(i) = ref_zmp_(i,1);
  }
  
  Eigen::Matrix1x3d C;
  C(0,0) = 1; C(0,1) = 0; C(0,2) = -zc_/GRAVITY;
  
  Eigen::VectorXd px, py;
  px.resize(1); py.resize(1);

  if(tick == 0 && current_step_num_ == 0)
  { 
    preview_x_b(0) = x_i; // 이때 before 값이 없으면 첫 tick에서 delta x의 에러가 갑작스럽게 생겨서 발산
    preview_y_b(0) = y_i;
    preview_x(0) = x_i;
    preview_y(0) = y_i;
  }
  else
  {   
    preview_x_b = preview_x;
    preview_y_b = preview_y;   
    preview_x = xs; preview_y = ys;
  }     
     
  px = C*preview_x;
  py = C*preview_y;
  
  Com_tra_graph << xs_(0) << "," << ys_(0) << "," << px << "," << py << "," << px_ref(tick) << "," << py_ref(tick)<< std::endl;
  
  double sum_Gd_px_ref = 0, sum_Gd_py_ref = 0;

  for(int i = 0; i < NL; i++) // Preview Step 개수만큼 더함.
  {
    sum_Gd_px_ref = sum_Gd_px_ref + Gd(i)*(px_ref(tick + 1 + i) - px_ref(tick + i));
    sum_Gd_py_ref = sum_Gd_py_ref + Gd(i)*(py_ref(tick + 1 + i) - py_ref(tick + i));
  }
 
  Eigen::MatrixXd del_ux(1,1);
  Eigen::MatrixXd del_uy(1,1);
  del_ux.setZero();
  del_uy.setZero();
   
  Eigen::VectorXd GX_X(1); 
  GX_X = Gx * (preview_x - preview_x_b);
  Eigen::VectorXd GX_Y(1); 
  GX_Y = Gx * (preview_y - preview_y_b);
  
  del_ux(0,0) = -(px(0) - px_ref(tick))*Gi(0,0) - GX_X(0) - sum_Gd_px_ref;
  del_uy(0,0) = -(py(0) - py_ref(tick))*Gi(0,0) - GX_Y(0) - sum_Gd_py_ref;
   
  UX = UX + del_ux(0,0);
  UY = UY + del_uy(0,0);

  XD = A*preview_x + B*UX;
  YD = A*preview_y + B*UY;  
  
}
/*
 // Control Performance model 미반영
void WalkingController::preview_MJ(double dt, int NL, int tick, double x_i, double y_i, Eigen::Vector3d xs,Eigen::Vector3d ys, double& ux, double& uy, 
       double gi, Eigen::VectorXd gp_l, Eigen::Matrix1x3d gx, Eigen::Matrix3d a, Eigen::Vector3d b, Eigen::Matrix1x3d c, Eigen::Vector3d &xd, Eigen::Vector3d &yd)
{
  // dt = 0.005, NL = 320 (Preview Step 갯수), tick = 0~239, x_i,y_i = 지지발에서 본 X,Y CoM의 초기값, xs,ys = 이전 tick의 CoM desired 
  int zmp_size;
  zmp_size = ref_zmp_.col(1).size(); // 보행 중 720개 (240 * 3)
  Eigen::VectorXd px_ref, py_ref;
  px_ref.resize(zmp_size);
  py_ref.resize(zmp_size);

  for(int i = 0; i < zmp_size; i++)
  {
    px_ref(i) = ref_zmp_(i,0);
    py_ref(i) = ref_zmp_(i,1);
  }

  preview_x_b = preview_x;
  preview_y_b = preview_y;

  if(tick == 0 && current_step_num_ == 0)
  { 
    preview_x_b(0) = x_i; // 이때 before 값이 없으면 첫 tick에서 delta x의 에러가 갑작스럽게 생겨서 발산
    preview_y_b(0) = y_i;
    preview_x(0) = x_i;
    preview_y(0) = y_i;
  }
  else
  {  
    preview_x = xs;
    preview_y = ys;
  }
   
  preview_x_b(2) = preview_x(2);
  preview_y_b(2) = preview_y(2);
  
  double xzmp_err = 0, yzmp_err = 0;
  
  Eigen::Matrix<double, 1, 1> px, py;
  px = c * preview_x;
  py = c * preview_y;
  // e(i) 정의
  xzmp_err = px(0) - px_ref(tick);
  yzmp_err = py(0) - py_ref(tick);

  double sum_gp_px_ref = 0, sum_gp_py_ref = 0;

  for(int i = 0; i < NL; i++) // Preview Step 개수만큼 더함.
  {
    sum_gp_px_ref += gp_l(i) * (px_ref(tick+1 + i)- px_ref(tick + i));
    sum_gp_py_ref += gp_l(i) * (py_ref(tick+1 + i)- py_ref(tick + i));
  }
  
  double gx_x, gy_y, del_ux, del_uy;
 
  del_ux = -(xzmp_err*gi) - (gx*(preview_x - preview_x_b)) - sum_gp_px_ref;
  del_uy = -(yzmp_err*gi) - (gx*(preview_y - preview_y_b)) - sum_gp_py_ref;

  xd = a*preview_x + b*del_ux;
  yd = a*preview_y + b*del_uy;

}

void WalkingController::previewParam_MJ(double dt, int NL, Eigen::Matrix4d& k, Eigen::Vector3d com_support_init_, double& gi, Eigen::VectorXd& gp_l, Eigen::Matrix1x3d& gx, 
  Eigen::Matrix3d& a, Eigen::Vector3d& b, Eigen::Matrix1x3d& c)
  { 
    a.setIdentity();
    a(0,1) = dt;
    a(0,2) = dt*dt/2;
    a(1,2) = dt;
    
    b.setZero();
    b(0) = dt*dt*dt/6;
    b(1) = dt*dt/2;
    b(2) = dt;

    c(0,0) = 1;
    c(0,1) = 0;
    c(0,2) = -zc_/GRAVITY;
    // A, B, C discrete Matrix 정의

    Eigen::Vector4d b_bar;
    b_bar(0) = c*b; // (1x3)*(3x1) = 1x1
    b_bar.segment(1,3) = b;

    Eigen::Matrix1x4d b_bar_tran;
    b_bar_tran = b_bar.transpose();

    Eigen::Vector4d i_bar;
    i_bar.setZero();
    i_bar(0) = 1;

    Eigen::Matrix4x3d f_bar;
    f_bar.setZero();
    f_bar.block<1,3>(0,0) = c*a;
    f_bar.block<3,3>(1,0) = a;

    Eigen::Matrix4d a_bar;
    a_bar.block<4,1>(0,0) = i_bar;
    a_bar.block<4,3>(0,1) = f_bar;

    double qe;
    qe = 1;
    Eigen::Matrix<double, 1,1> r;
    r(0,0) = 0.000001;

    Eigen::Matrix3d qx;
    qx.setZero();
    Eigen::Matrix4d q_bar;
    q_bar.setZero();
    q_bar(0,0) = qe;
    q_bar.block<3,3>(1,1) = qx;

    k = discreteRiccatiEquationPrev(a_bar, b_bar, r, q_bar); 
    std::cout << k << std::endl;
    double temp_mat;
    temp_mat = r(0) + b_bar_tran * k * b_bar;
     
    Eigen::Matrix4d ac_bar;
    ac_bar.setZero();
    ac_bar = a_bar - b_bar/temp_mat*b_bar_tran*k*a_bar;

    gi = (b_bar_tran*k*i_bar);
    gi *= 1/temp_mat;
    gx = (b_bar_tran*k*f_bar)/temp_mat;

    Eigen::MatrixXd x_l(4, NL);
    Eigen::Vector4d x_l_column;
    x_l.setZero();
    x_l_column.setZero();
    x_l_column = -ac_bar.transpose()*k*i_bar;

    for(int i=0; i<NL; i++)
    {
      x_l.col(i) = x_l_column;
      x_l_column = ac_bar.transpose()*x_l_column;
    }

    gp_l.resize(NL);
    double gp_l_column;
    gp_l_column = -gi;
    
    for(int i=0; i<NL; i++)
    {
      gp_l(i) = gp_l_column;
      gp_l_column = b_bar_tran*x_l.col(i);
      gp_l_column = gp_l_column/temp_mat;
    }
  }
*/

void WalkingController::compensator()
{
  if(hip_compensator_mode_ == true)
  {
    hipCompensation(); 
    if(lqr_compensator_mode_ == false)
    { 
      hipCompensator();
    }
  }

  if(lqr_compensator_mode_ == true)
  {
    Eigen::Vector12d d_q;

    for (int i=0; i<12; i++)
      d_q(i) = desired_q_(i); // IK 풀어서 나온 Desired Joint angle

    Eigen::Vector12d lqr_joint_input;

    slowcalc_mutex_.lock();

    ad_copy_ = ad_right_;
    bd_copy_ = bd_right_;
    ad_total_copy_ = ad_total_right_;
    bd_total_copy_ = bd_total_right_;
    kkk_copy_ = kkk_motor_right_;
    slowcalc_mutex_.unlock();
    
    vibrationControl_MJ(d_q,lqr_joint_input);
    
    //LQR_q_tra_graph << desired_q_(0) << "," << desired_q_(1) << "," << desired_q_(2) << "," << desired_q_(3) << "," << desired_q_(4) << "," << desired_q_(5) << "," << DOB_IK_output_(0) << "," << DOB_IK_output_(1) << "," << DOB_IK_output_(2) << "," << DOB_IK_output_(3) << "," << DOB_IK_output_(4) << "," << DOB_IK_output_(5) << std::endl;
    
    double grav_gain_timing = 1.0;
    if(foot_step_(current_step_num_,6) == 1) // left foot support (right foot gain)
    {
      if(walking_tick_ > t_start_+t_rest_init_+t_double1_ && walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_)
        grav_gain_timing = DyrosMath::cubic(walking_tick_,t_start_+t_rest_init_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_,1.0,0.3,0.0,0.0);
      else
        grav_gain_timing = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_,t_start_+t_total_-t_rest_last_,0.3,1.0,0.0,0.0);
    }


    for (int n = 7 ; n < 12; n++) // right foot
    {
      if(abs(lqr_joint_input(n)-desired_q_(n)) > 20.0*DEG2RAD )
      {
      }
      else
      {
        if(n == 7 || n == 8)
          desired_q_(n) = desired_q_(n) - 0.0022*grav_gain_timing*grav_ground_torque_[n];
        else if (n == 9)
          desired_q_(n) = desired_q_(n) - 0.0010*grav_gain_timing*grav_ground_torque_[n];
        else
          desired_q_(n) = desired_q_(n);
      }
    }

    grav_gain_timing = 1.0;
    if(foot_step_(current_step_num_,6) == 0) // left foot (right foot gain)
    {
      if(walking_tick_ > t_start_+t_rest_init_+t_double1_ && walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_)
        grav_gain_timing = DyrosMath::cubic(walking_tick_,t_start_+t_rest_init_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_,1.0,0.3,0.0,0.0);
      else
        grav_gain_timing = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_,t_start_+t_total_-t_rest_last_,0.3,1.0,0.0,0.0);
    }

    for (int n = 1 ; n < 6; n++)
    {
      if(abs(lqr_joint_input(n)-desired_q_(n)) > 20.0*DEG2RAD )
      {
      }
      else
      {
        if ( n == 1 || n == 2)
          desired_q_(n)  = desired_q_(n) - 0.0022*grav_gain_timing*grav_ground_torque_[n];
        else if (n == 3)
          desired_q_(n) = desired_q_(n) - 0.0010*grav_gain_timing*grav_ground_torque_[n];
        else
          desired_q_(n) = desired_q_(n);
      }
    }
  }
}

void WalkingController::hipCompensator()
{
  double left_hip_angle = 4.0*DEG2RAD, right_hip_angle = 8.5*DEG2RAD, left_hip_angle_first_step = 4.0*DEG2RAD, right_hip_angle_first_step = 8.5*DEG2RAD,

      left_hip_angle_temp = 0.0, right_hip_angle_temp = 0.0, temp_time = 0.1*hz_, left_pitch_angle = 0.0*DEG2RAD;

  if (current_step_num_ == 0)
  {
    if(foot_step_(current_step_num_, 6) == 1) //left support foot
    {
      if(walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_real_+t_double1_+temp_time,0.0*DEG2RAD, left_hip_angle_first_step, 0.0, 0.0);
      else if(walking_tick_ >= t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-temp_time,t_start_+t_total_-t_rest_last_,left_hip_angle_first_step, 0.0, 0.0, 0.0);
      else
        left_hip_angle_temp = 0.0*DEG2RAD;
    }
    else if (foot_step_(current_step_num_, 6) == 0) // right support foot
    {
      if(walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_real_+t_double1_+temp_time,0.0*DEG2RAD, right_hip_angle_first_step, 0.0, 0.0);
      else if(walking_tick_ >= t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-temp_time,t_start_+t_total_-t_rest_last_,right_hip_angle_first_step, 0.0, 0.0, 0.0);
      else
        right_hip_angle_temp = 0.0*DEG2RAD;
    }
    else
    {
      left_hip_angle_temp = 0.0*DEG2RAD;
      right_hip_angle_temp = 0.0*DEG2RAD;
    }
  }
  else
  {
    if(foot_step_(current_step_num_, 6) == 1)
    {
      if(walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_real_+t_double1_+temp_time,0.0*DEG2RAD,left_hip_angle,0.0,0.0);
      else if (walking_tick_ >= t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-temp_time,t_start_+t_total_-t_rest_last_,left_hip_angle,0.0,0.0,0.0);
      else
        left_hip_angle_temp = 0.0*DEG2RAD;

    }
    else if(foot_step_(current_step_num_,6) == 0)
    {
      if(walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_real_+t_double1_+temp_time,0.0*DEG2RAD,right_hip_angle,0.0,0.0);
      else if(walking_tick_ >= t_start_+t_total_-t_rest_last_-t_double2_-temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-temp_time,t_start_+t_total_-t_rest_last_,left_hip_angle,0.0,0.0,0.0);
      else
        right_hip_angle_temp = 0.0*DEG2RAD;
    }
    else
    {
      left_hip_angle_temp = 0.0*DEG2RAD;
      right_hip_angle_temp = 0.0*DEG2RAD;
    }
  }
  desired_q_(1) = desired_q_(1) + left_hip_angle_temp;
  desired_q_(7) = desired_q_(7) - right_hip_angle_temp;
  joint_offset_angle_(1) = left_hip_angle_temp;
  joint_offset_angle_(7) = -right_hip_angle_temp;
}
void WalkingController::hipCompensation()
{
  double a_total, b_total, alpha_l, alpha_r, rq0, rq1, rq2, rq3, rq4, rq5, lq0, lq1, lq2, lq3, lq4, lq5, robotweight, fromright, fromleft, alpha, alpha_1, f_r, f_l; //alpha is weighting factor

  a_total = -0.0012;
  b_total = 0.00087420;
  //robotweight = 50.3082*9.81;

  robotweight = 46.892*9.81;

  lq0 = desired_q_(0);
  lq1 = desired_q_(1);
  lq2 = desired_q_(2);
  lq3 = desired_q_(3);
  lq4 = desired_q_(4);
  lq5 = desired_q_(5);
  
  rq0 = desired_q_(6);
  rq1 = desired_q_(7);
  rq2 = desired_q_(8);
  rq3 = desired_q_(9);
  rq4 = desired_q_(10);
  rq5 = desired_q_(11);

  fromright = cos(rq2)*sin(rq0)*(-1.454E-1)-sin(rq0)*sin(rq2)*(3.39E2/1.0E3)-cos(rq3)*(cos(rq2)*sin(rq0)+cos(rq0)*sin(rq1)*sin(rq2))*(3.0/5.0E1)-cos(rq3)*(sin(rq0)*sin(rq2)-cos(rq0)*cos(rq2)*sin(rq1))*(4.6E1/1.25E2)-sin(rq3)*(cos(rq2)*sin(rq0)+cos(rq0)*sin(rq1)*sin(rq2))*(4.6E1/1.25E2)+sin(rq3)*(sin(rq0)*sin(rq2)-cos(rq0)*cos(rq2)*sin(rq1))*(3.0/5.0E1)+cos(rq0)*cos(rq2)*sin(rq1)*(3.39E2/1.0E3)-cos(rq0)*sin(rq1)*sin(rq2)*1.454E-1-2.1E1/2.0E2;
  fromleft = cos(lq2)*sin(lq0)*(-1.454E-1)+sin(lq0)*sin(lq2)*(3.39E2/1.0E3)+cos(lq3)*(sin(lq0)*sin(lq2)+cos(lq0)*cos(lq2)*sin(lq1))*(4.6E1/1.25E2)-cos(lq3)*(cos(lq2)*sin(lq0)-cos(lq0)*sin(lq1)*sin(lq2))*(3.0/5.0E1)+sin(lq3)*(sin(lq0)*sin(lq2)+cos(lq0)*cos(lq2)*sin(lq1))*(3.0/5.0E1)+sin(lq3)*(cos(lq2)*sin(lq0)-cos(lq0)*sin(lq1)*sin(lq2))*(4.6E1/1.25E2)+cos(lq0)*cos(lq2)*sin(lq1)*(3.39E2/1.0E3)+cos(lq0)*sin(lq1)*sin(lq2)*1.454E-1+2.1E1/2.0E2;

  alpha = -fromleft / (fromright - fromleft);
  
  if(fromright >= 0)
  {
    alpha = 1;
  }
  if(fromleft <= 0)
  {
    alpha = 0;
  }
  std::cout << fromleft << "," << fromright << "," << alpha << std::endl;
  alpha_1 = 1-alpha;

  f_r = robotweight*alpha;
  f_l = robotweight*alpha_1;

  Eigen::Vector6d lTau, rTau, Jc2, Jc8;
  lTau.setZero();
  rTau.setZero();
  Jc2.setZero();
  Jc8.setZero();

  Jc2(0) = 0;
  Jc2(1) = cos(rq2)*sin(rq1)*(3.39E2/1.0E3)-sin(rq1)*sin(rq2)*1.454E-1-sin(rq1)*sin(rq2)*sin(rq3)*(4.6E1/1.25E2)+cos(rq2)*cos(rq3)*sin(rq1)*(4.6E1/1.25E2)-cos(rq2)*sin(rq1)*sin(rq3)*(3.0/5.0E1)-cos(rq3)*sin(rq1)*sin(rq2)*(3.0/5.0E1);
  Jc2(2) = cos(rq1)*cos(rq2)*1.454E-1+cos(rq1)*sin(rq2)*(3.39E2/1.0E3)+cos(rq1)*cos(rq2)*cos(rq3)*(3.0/5.0E1)+cos(rq1)*cos(rq2)*sin(rq3)*(4.6E1/1.25E2)+cos(rq1)*cos(rq3)*sin(rq2)*(4.6E1/1.25E2)-cos(rq1)*sin(rq2)*sin(rq3)*(3.0/5.0E1);
  Jc2(3) = cos(rq1)*cos(rq2)*cos(rq3)*(3.0/5.0E1)+cos(rq1)*cos(rq2)*sin(rq3)*(4.6E1/1.25E2)+cos(rq1)*cos(rq3)*sin(rq2)*(4.6E1/1.25E2)-cos(rq1)*sin(rq2)*sin(rq3)*(3.0/5.0E1);
  Jc2(4) = 0;
  Jc2(5) = 0;

  Jc8(0) = 0;
  Jc8(1) = cos(lq2)*sin(lq1)*(3.39E2/1.0E3)+sin(lq1)*sin(lq2)*1.454E-1-sin(lq1)*sin(lq2)*sin(lq3)*(4.6E1/1.25E2)+cos(lq2)*cos(lq3)*sin(lq1)*(4.6E1/1.25E2)+cos(lq2)*sin(lq1)*sin(lq3)*(3.0/5.0E1)+cos(lq3)*sin(lq1)*sin(lq2)*(3.0/5.0E1);
  Jc8(2) = cos(lq1)*cos(lq2)*(-1.454E-1)+cos(lq1)*sin(lq2)*(3.39E2/1.0E3)-cos(lq1)*cos(lq2)*cos(lq3)*(3.0/5.0E1)+cos(lq1)*cos(lq2)*sin(lq3)*(4.6E1/1.25E2)+cos(lq1)*cos(lq3)*sin(lq2)*(4.6E1/1.25E2)+cos(lq1)*sin(lq2)*sin(lq3)*(3.0/5.0E1);
  Jc8(3) = cos(lq1)*cos(lq2)*cos(lq3)*(-3.0/5.0E1)+cos(lq1)*cos(lq2)*sin(lq3)*(4.6E1/1.25E2)+cos(lq1)*cos(lq3)*sin(lq2)*(4.6E1/1.25E2)+cos(lq1)*sin(lq2)*sin(lq3)*(3.0/5.0E1);
  Jc8(4) = 0;
  Jc8(5) = 0;

  for(int i=0; i<6; i++)
  {
    rTau(i)=Jc2(i)*f_r;
    lTau(i)=Jc8(i)*f_l;
  }

  double rising = 1.0, timingtiming = 1.0, k = 0.2, k1 = 0.2;

  if (lqr_compensator_mode_ == false)
  {
    desired_q_(8)=desired_q_(8)+(a_total*rTau(2)+b_total)*rising*k1;//offwhenslow
    desired_q_(9)=desired_q_(9)+(a_total*rTau(3)+b_total)*rising*0.3;//offwhenslow
    desired_q_(10)=desired_q_(10)+(a_total*rTau(4)+b_total)*rising*k1;//offwhenslow
  }

  joint_offset_angle_(8) = (a_total*rTau(2)+b_total)*rising*k;
  joint_offset_angle_(9) = (a_total*rTau(3)+b_total)*rising*0.2;
  joint_offset_angle_(10) = (a_total*rTau(4)+b_total)*rising*k;

  if (lqr_compensator_mode_  == false)
  {
    desired_q_(2)=desired_q_(2)+(a_total*lTau(2)+b_total)*rising*k;//offwhenslow
    desired_q_(3)=desired_q_(3)+(a_total*lTau(3)+b_total)*rising*0.3;//offwhenslow
    desired_q_(4)=desired_q_(4)+(a_total*lTau(4)+b_total)*rising*k;//offwhenslow
  }

  joint_offset_angle_(2) = (a_total*lTau(2)+b_total)*rising*k1;
  joint_offset_angle_(3) = (a_total*lTau(3)+b_total)*rising*0.2;
  joint_offset_angle_(4) = (a_total*lTau(4)+b_total)*rising*k1;

  for (int i = 0; i < 6; i ++)
  {
    grav_ground_torque_(i) = lTau[i];
    grav_ground_torque_(i+6) = rTau[i];
  }

}

void WalkingController::vibrationControl_MJ(const Eigen::Vector12d desired_leg_q, Eigen::Vector12d &output)
{ 
  if(walking_tick_ == 0) // 모터 엔코더 각도, 링크에 부착 된 외부 엔코더의 각도 초기화.
  {
    pre_motor_q_leg_ = current_motor_q_leg_;
    pre_link_q_leg_ = current_link_q_leg_;
    //lqr_output_pre_ = current_motor_q_leg_;
    DOB_IK_output_b_ = current_motor_q_leg_;
  }
  x_bar_right_.setZero(); // X_bar는 e(t)와 theta_m_dot, j_dot, j_ddot 으로 구성된 48 x 1 상태 변수.
  
  for (int i = 0; i < 12; i++)
  { // Simulation이라서 그런지 getRobotState 에서 link_q와 motor_q에 같은 값이 들어가고 있다.
    x_bar_right_(i, 0) = current_link_q_leg_(i) - desired_leg_q(i); // e(t) // current_link_q_leg는 현재 외부 엔코더값이 아니라, Motor angle이 들어가고 있음. desired_leg_q는 IK를 풀어서 나온 Desired 관절 각
    x_bar_right_(i+12,0) = current_motor_q_leg_(i) - pre_motor_q_leg_(i); // theta_m_dot -> theta_m 의 변화량
    x_bar_right_(i+24,0) = 0.0; // current_link_q_leg_(i) - pre_link_q_leg_(i); // theta_j_dot -> theta_j 의 변화량, 민곤이형 실수로 계속 0들어가고 있었음.
    //x_bar_right_(i+36,0)은 의도적으로 노이즈때문에 뺐음.
  }

  Eigen::Vector12d del_u_right;
  del_u_right = -kkk_copy_*x_bar_right_; // (12 x 48) X (48 x 1) -> u_c dot이 discrete니까 del_u_c
  x_bar_right_ = ad_total_copy_*x_bar_right_ + bd_total_copy_*del_u_right;
  /////////////////////////////////////////////////////////////////////////////////////////////////////LQR 기반 위치제어 입력.
  
  // ad,bd_total 은 LQR Dynamics, ad,bd는 모터 제어 Dynamics
  Eigen::Vector12d current_u; // left right order
  
  for (int i = 0; i < 12; i++)
  { // (theta_m(k+1) - theta_m(k)) / dt = Kp (u - theta_m(k)) 
    current_u(i) = (current_motor_q_leg_(i) - ad_copy_(i, i)*pre_motor_q_leg_(i)) / bd_copy_(i, i); // ad_copy, bd_copy에 들어가는 Kp를 구하는게 애매함.
    // ad_copy는 A (36x36), bd_copy는 B (36x12)를 Discretization 한 것이고, 거기서 첫번째 상태 변수인 theta_m (12개)만 뽑아서 사용한 것.
  }
  
  Eigen::Vector12d d_hat; //left right order
  // d_hat = lqr_output_pre_ - current_u ; //lqr_output_pre_ : u' , current_u : (Pn^-1(s))*theta_m(motor angles) -> u' + d
  d_hat = DOB_IK_output_b_ - current_u ; // 부호가 반대
  // d_hat 을 이런식으로 관측함.

  if(walking_tick_ == 0)
    d_hat_b = d_hat;

  d_hat = 0.7*d_hat_b + 0.3*d_hat; // 필터링

  // Mingon's LQR contorller gain (using external encoder)
  double default_gain = 0.2; // Kp가 정확하다면 시뮬레이션이나 실제 로봇이나 0.2~1의 의미는 같다.
  double compliant_gain = 1;
  double compliant_tick = 0.1*hz_;
  double gain_temp;
  for (int i = 0; i < 12; i ++)
  {
    if(i < 6) //왼쪽 다리 관절
    {
      //double gain_temp = default_gain;
      gain_temp = default_gain;
      if(walking_enable_ == true)
      {
        if (foot_step_(current_step_num_,6) == 0) // 오른발 지지 상태
        { 
          if(walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_ - compliant_tick) // gain_temp -> 0.2
          { // t_total_: 240tick, t_rest_init_,last: 10tick (0.05초), t_double1,2_: 20tick (0.1초), compliant_tick : 20tick (0.1초)
            gain_temp = default_gain; // 1step 240tick 동안 0~190 tick까지 계속 기본 gain.
          }
          else if(walking_tick_ >= t_start_ + t_total_ - t_rest_last_ - t_double2_ - compliant_tick && walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_)
          { 
            gain_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_ - t_double2_ - compliant_tick, t_start_ + t_total_ - t_rest_last_ - t_double2_, default_gain, compliant_gain, 0.0, 0.0);
          } // 1step 240tick 동안 190~210 tick까지 기본 gain 0.2에서 1까지 올림.
          else
          {
            gain_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_, t_start_ + t_total_, compliant_gain, default_gain, 0.0, 0.0);
          } // 1step 240tick 동안 210~240 tick까지 기본 gain 1에서 0.2까지 낮춤.
        } 
        else // 왼발 지지 상태
        {
          gain_temp = default_gain;
        }
      }
      //lqr_output_(i) = lqr_output_pre_(i) + del_u_right(i, 0) - gain_temp*d_hat(i); // u(tk) = uc(tk - del t) + del_u_c(discrete) + d_hat
      DOB_IK_output_(i) = desired_leg_q(i) - gain_temp*d_hat(i);  // LQR 위치 대신 단순 IK 기반 위치 ( u_c + d_hat = u' (논문))
    }
    else //오른 다리 관절
    {
      //double gain_temp = default_gain;
      gain_temp = default_gain;

      if(walking_enable_ == true)
      {
        if (foot_step_(current_step_num_,6) == 1) // 왼발 지지 상태
        {
          if(walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_ - compliant_tick) // gain_temp -> 0.2
          {
            gain_temp = default_gain;
          }
          else if(walking_tick_ >= t_start_ + t_total_ - t_rest_last_ - t_double2_ - compliant_tick && walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_)
          {
            gain_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_ - t_double2_ - compliant_tick, t_start_ + t_total_ - t_rest_last_ - t_double2_, default_gain, compliant_gain, 0.0, 0.0);
          }
          else
          {
            gain_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_, t_start_ + t_total_, compliant_gain, default_gain, 0.0, 0.0);
          }
        }
        else // 오른발 지지 상태
        {
          gain_temp = default_gain;
        }
      }
      //lqr_output_(i) = lqr_output_pre_(i) + del_u_right(i, 0) - gain_temp*d_hat(i);
      DOB_IK_output_(i) = desired_leg_q(i) - gain_temp*d_hat(i);  // LQR 위치 대신 단순 IK 기반 위치
    }
  }
  
  lqr_output_pre_ = lqr_output_; // 시뮬에선 안씀.
  pre_motor_q_leg_ = current_motor_q_leg_;
  pre_link_q_leg_ = current_link_q_leg_;
  d_hat_b = d_hat;
  DOB_IK_output_b_ = DOB_IK_output_;

  for (int i=0; i<12; i++)
  {
    output(i) = DOB_IK_output_(i); 
  }
}
/*
void WalkingController::vibrationControl(const Eigen::Vector12d desired_leg_q, Eigen::Vector12d &output)
{
  if(walking_tick_ == 0)
  {
    pre_motor_q_leg_ = current_motor_q_leg_;
    pre_link_q_leg_ = current_link_q_leg_;
    lqr_output_pre_ = current_motor_q_leg_;
  }
  x_bar_right_.setZero();

  // right foot X_bar
  for (int i = 0; i < 12; i++)
  {
    x_bar_right_(i, 0) = current_link_q_leg_(i) - desired_leg_q(i); //left

    x_bar_right_(i+12,0) = current_motor_q_leg_(i) - pre_motor_q_leg_(i); //left

    x_bar_right_(i+24,0) = 0.0;//current_link_q_leg_(i) - pre_link_q_leg_(i); //left

  }

  Eigen::Vector12d del_u_right;
  del_u_right = -kkk_copy_*x_bar_right_;
  x_bar_right_ = ad_total_copy_*x_bar_right_ + bd_total_copy_*del_u_right;

  Eigen::Vector12d current_u;// left right order

  for (int i = 0; i < 12; i++)
  {
    //left
    current_u(i) = (current_motor_q_leg_(i) - ad_copy_(i, i)*pre_motor_q_leg_(i)) / bd_copy_(i, i);

    //right
    //current_u(i+6) = (current_motor_q_leg_(i+6) - ad_copy_(i, i)*pre_motor_q_leg_(i+6)) / bd_copy_(i, i);
  }


  Eigen::Vector12d dist;//left right order


  dist = lqr_output_pre_ - current_u; //lqr_output_pre_ : u' , current_u : u_c
  // d_hat 을 이런식으로 관측함.

  if(walking_tick_ == 0)
    dist_prev_ = dist;

  dist = 0.7*dist_prev_+0.3*dist;

  // Mingon's LQR contorller gain(using external encoder)
  double default_gain = 0.200;

  double compliant_gain = 1.0;

  double compliant_tick = 0.1*hz_;


  for (int i = 0; i < 12; i++)
  {
    if(i < 6) //left foot
    {
      double gain_temp = default_gain;


      if(walking_enable_ == true)
      {
        if (foot_step_(current_step_num_,6) == 0) // right support foot
        {
          if(walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_-compliant_tick) // gain_temp -> 0.2
          {
            gain_temp = default_gain;
          }
          else if(walking_tick_ >= t_start_+t_total_-t_rest_last_-t_double2_-compliant_tick && walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_)
          {
            gain_temp = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-compliant_tick,t_start_+t_total_-t_rest_last_-t_double2_,default_gain,compliant_gain,0.0,0.0);
          }
          else
          {
            gain_temp = DyrosMath::cubic(walking_tick_, t_start_+t_total_-t_rest_last_,t_start_+t_total_,compliant_gain,default_gain,0.0,0.0);
          }
        }
        else //left support foot
        {
          gain_temp = default_gain;

        }
      }
      lqr_output_(i) = lqr_output_pre_(i) + del_u_right(i, 0) - gain_temp*dist(i);
    }
    else // right foot
    {
      double gain_temp = default_gain;


      if(walking_enable_ == true)
      {
        if (foot_step_(current_step_num_,6) == 1) // left suppor foot
        {
          if(walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_-compliant_tick) // gain_temp -> 0.2
          {
            gain_temp = default_gain;
          }
          else if(walking_tick_ >= t_start_+t_total_-t_rest_last_-t_double2_-compliant_tick && walking_tick_ < t_start_+t_total_-t_rest_last_-t_double2_)
          {
            gain_temp = DyrosMath::cubic(walking_tick_,t_start_+t_total_-t_rest_last_-t_double2_-compliant_tick,t_start_+t_total_-t_rest_last_-t_double2_,default_gain,compliant_gain,0.0,0.0);
          }
          else
          {
            gain_temp = DyrosMath::cubic(walking_tick_, t_start_+t_total_-t_rest_last_,t_start_+t_total_,compliant_gain,default_gain,0.0,0.0);
          }
        }
        else //left foot support
        {
          gain_temp = default_gain;
        }
      }
      lqr_output_(i) = lqr_output_pre_(i) + del_u_right(i, 0) - gain_temp*dist(i);
    }
  }

  lqr_output_pre_ = lqr_output_;
  pre_motor_q_leg_ = current_motor_q_leg_;
  pre_link_q_leg_ = current_link_q_leg_;
  dist_prev_ = dist;

  for (int i=0; i<12; i++)
  {
    output(i) = lqr_output_(i);
  }
}*/

void WalkingController::massSpringMotorModel(double spring_k, double damping_d, double motor_k, Eigen::Matrix12d & mass, Eigen::Matrix<double, 36, 36>& a, Eigen::Matrix<double, 36, 12>& b, Eigen::Matrix<double, 12, 36>& c)
{
  int dof = 12;
  Eigen::Matrix12d spring_k_mat;
  Eigen::Matrix12d damping_d_mat;
  Eigen::Matrix12d motor_k_mat;

  spring_k_mat.setZero();
  damping_d_mat.setZero();
  motor_k_mat.setZero();

  for (int i = 0; i < dof; i++) // 전부 diagonal matrix니까
  {
    spring_k_mat(i, i) = spring_k; // K = 2000
    damping_d_mat(i, i) = damping_d; // D = 100
    motor_k_mat(i, i) = motor_k; // Kp = 10, 시뮬레이션에서는 30. 20이하 발산.
  }
  // knee joint
  motor_k_mat(3,3) = 18.0 * 3; // 기존 코드 18, 시뮬레이션 54
  motor_k_mat(9,9) = 18.0 * 3;

  Eigen::Matrix12d inv_mass;
  inv_mass = mass;

  Eigen::Matrix12d zero_12;
  zero_12.setZero();

  Eigen::Matrix12d eye_12;
  eye_12.setIdentity();

  Eigen::Matrix12d mass_temp;
  mass_temp = inv_mass;

  Eigen::Matrix12d a1;
  a1 = mass_temp*(spring_k_mat - damping_d_mat*motor_k_mat); // S_m*(K - D * K_p)

  Eigen::Matrix12d a2;
  a2 = -mass_temp*(spring_k_mat); // -S_m * K

  Eigen::Matrix12d a3;
  a3 = -mass_temp*(damping_d_mat); // -S_m * D

  a.setZero();
  a << -motor_k_mat, zero_12, zero_12, zero_12, zero_12, eye_12, a1, a2, a3; // a_right_mat

  Eigen::Matrix12d b1;
  b1 = mass_temp*damping_d_mat*motor_k_mat;


  b.setZero();
  b << motor_k_mat, zero_12, b1; // b_right_mat -> Kp, 0(nxn), S_m * D * K_p

  c.setZero();
  c << zero_12, eye_12, zero_12; // c_right_mat

}

void WalkingController::discreteModel(Eigen::Matrix<double, 36, 36>& a, Eigen::Matrix<double, 36, 12>& b, Eigen::Matrix<double, 12, 36>& c, int np, double dt,
                                      Eigen::Matrix<double, 36, 36>& ad, Eigen::Matrix<double, 36, 12>& bd, Eigen::Matrix<double, 12, 36>& cd,
                                      Eigen::Matrix<double, 48, 48>& ad_total, Eigen::Matrix<double, 48, 12>& bd_total)
{
  int n = a.rows(); //state \B0\B9\BC\F6  (36)
  int r = b.cols(); // Input \B0\B9\BC\F6  (12)
  int p = c.rows(); // output \B0\B9\BC\F6  (12)
  
  Eigen::Matrix<double, 36, 36> inv_a;
  inv_a = a.inverse();

  Eigen::Matrix<double, 36, 36> eye6;
  eye6.setIdentity();
  // X(k+1) = (I + A*dt)*X(k) + B*dt*U(k)
  ad = eye6 +a*dt; 
  bd = b*dt;
  cd = c;
  
  Eigen::Matrix<double, 12, 36> ca;
  ca = cd*ad;

  Eigen::Matrix<double, 12, 12> eye_p;
  eye_p.setIdentity();

  Eigen::Matrix<double, 36, 12> zero_n_p;
  zero_n_p.setZero();

  Eigen::Matrix<double, 12, 12> cb;
  cb = cd*bd;
  
  if (np < 1) // discretization (Preview 할때 Katayama 논문과 같은 방식, State에 Error가 추가 됐을 때 48 x 48 matrix는 다음과 같음.)
  {
    ad_total << eye_p, ca, zero_n_p, ad;
    bd_total << cb, bd;
  }
  else
  {
    ad_total.resize(n + p + np, n + p + np);
    bd_total.resize(n + p + np, r);

    Eigen::MatrixXd zero_temp1;
    zero_temp1.resize(p, (np - 1)*p); zero_temp1.setZero();

    Eigen::MatrixXd zero_temp2;
    zero_temp2.resize(n, np*p); zero_temp2.setZero();

    Eigen::MatrixXd zero_temp3;
    zero_temp3.resize(np*p, p + n); zero_temp3.setZero();

    Eigen::MatrixXd shift_register;
    shift_register.resize(np*p, np*p); shift_register.setZero();


    Eigen::MatrixXd zero_temp4;
    zero_temp4.resize(np*p, r); zero_temp4.setZero();
    for (int i = 0; i<p*(np - 1); i++)
      shift_register(i, i + p) = 1;

    ad_total << eye_p, ca, -eye_p, zero_temp1, zero_n_p, ad, zero_temp2, zero_temp3, shift_register;
    bd_total << cb, bd, zero_temp4;
  }

}

void WalkingController::riccatiGain(Eigen::Matrix<double, 48, 48>& ad_total, Eigen::Matrix<double, 48, 12>& bd_total, Eigen::Matrix<double, 48, 48>& q, Eigen::Matrix12d& r, Eigen::Matrix<double, 12, 48>& k)
{

  const int n = ad_total.rows(); //state \B0\B9\BC\F6 (48)
  const int m = bd_total.cols(); // Input \B0\B9\BC\F6 (12)

  static bool InitDare;

  if (!InitDare)
  {
    std::cout << "insiadfnisdaifnadisfasd" <<std::endl;
    discreteRiccatiEquationInitialize(ad_total, bd_total);
    InitDare = true;
  }

  kkk_ = discreteRiccatiEquationLQR(ad_total, bd_total, r, q);

  Eigen::Matrix<double, 12, 48> trans_bd;
  trans_bd = bd_total.transpose();
  Eigen::Matrix<double, 12, 12> temp_r_inv;
  temp_r_inv = (r + trans_bd*kkk_*bd_total);
  temp_r_inv = temp_r_inv.inverse();
  k = temp_r_inv *trans_bd*kkk_*ad_total;
}

void WalkingController::slowCalc()
{
  while(true)
  {
    if(ready_for_thread_flag_)
    {
      slowCalcContent();
      if (ready_for_compute_flag_ == false)
      {
        ready_for_compute_flag_ = true;
      }
    }

    this_thread::sleep_for(chrono::milliseconds(100));

  }
}

void WalkingController::slowCalcContent()
{
  Eigen::Vector28d qqq;
  slowcalc_mutex_.lock();

  for(int i =0; i<28; i++)
  {
    qqq(i) = thread_q_(i); // Motor의 Current angle (theta_m)
  }
  slowcalc_mutex_.unlock();

  Eigen::Matrix<double, 6, 18> contact_j;
  Eigen::Matrix6d lamda;
  Eigen::Matrix<double, 6, 18> j_c_bar;
  Eigen::Matrix<double, 18, 18> pc;
  Eigen::Matrix<double, 18, 18> eye_18;
  Eigen::Matrix<double, 18, 18> mass_inverse;
  Eigen::Matrix<double, 12, 12> temp22;
  Eigen::Matrix<double, 48, 48> q_mat;
  Eigen::Matrix12d r_mat;

  Eigen::Matrix<double, 12, 48> lqr_k; //lqr_k.resize(12, 12*4);

  if(thread_tick_ == 0)
  {
    mass_matrix_.setZero();
    mass_matrix_pc_.setZero();
    mass_matrix_sel_.setZero();
  }

  mass_matrix_ = model_.getLegInertia();

  if((walking_enable_) == true)
  {
    if(foot_step_(current_step_num_,6) == 0) // 오른발
    {
      contact_j = model_.getLegWithVLinkJacobian((DyrosJetModel::EndEffector)(1)); // 발끝(E.E) 자코비언, 6 X 18 (Virtual Joint 6 + 다리 관절 12)
    }
    else if(foot_step_(current_step_num_,6) == 1) // 왼발
    {
      contact_j = model_.getLegWithVLinkJacobian((DyrosJetModel::EndEffector)(0));
    }

    lamda = (contact_j*mass_matrix_.inverse()*contact_j.transpose());
    lamda = lamda.inverse();
    j_c_bar = lamda*contact_j*mass_matrix_.inverse();

    pc = contact_j.transpose()*j_c_bar;

    eye_18.setIdentity();


    mass_inverse = mass_matrix_.inverse();

    mass_matrix_pc_ = mass_inverse*(eye_18-pc);

    for (int i=0; i<12; i++)
    {
      for (int j=0; j<12; j++)
        mass_matrix_sel_(i,j) = mass_matrix_pc_(i+6,j+6);
    }

  }
  else // 공중에 떠있을 때, Contact Jacobian 없을 때 사용했음.
  {
    mass_inverse = mass_matrix_.inverse();

    for (int i=0; i<12; i++)
    {
      for (int j=0; j<12; j++)
        mass_matrix_sel_(i,j) = mass_inverse(i+6,j+6);
    }
  }

  temp22.setZero();
  temp22 = mass_matrix_sel_;

  mass_matrix_sel_.setZero();
  for (int i=0; i<12; i++)
  {
    mass_matrix_sel_(i,i) = temp22(i,i); // Diagonal 부분만 빼옴.
  }

  double spring_k = 2000.0;
  double damping_d = 100.0;
  double motor_k = 10.0*3; // 기존 코드 10, Simulation 에서 30
  massSpringMotorModel(spring_k, damping_d, motor_k, mass_matrix_sel_, a_right_mat_, b_right_mat_, c_right_mat_);
  discreteModel(a_right_mat_, b_right_mat_, c_right_mat_, 0, 1.0/hz_, a_disc_, b_disc_, c_right_mat_, a_disc_total_, b_disc_total_);

  q_mat.setIdentity();
  q_mat = q_mat*1.0;

  for (int i=0; i<12; i++)
  {
    q_mat(i,i) = 10.0;
  }

  r_mat.setIdentity();
  r_mat = r_mat * 1.0;

  riccatiGain(a_disc_total_, b_disc_total_, q_mat, r_mat, lqr_k); // q_mat (48 x 48) r_mat (12 x 12) lqr_k (12 x 48)

  if(ready_for_compute_flag_==false)
  {
  }

  calc_update_flag_ = true;
  thread_tick_++;

  slowcalc_mutex_.lock();

  ad_right_ = a_disc_; // Ad, Bd -> 36x36 Matrix, State가 theta_m, theta_j, theta_j_dot 을 Discretization 한 Matrix
  bd_right_ = b_disc_;
  ad_total_right_ = a_disc_total_; // e(t)를 상태변수로 추가 한 Discrete 모델.
  bd_total_right_ = b_disc_total_;
  kkk_motor_right_ = lqr_k;
  slowcalc_mutex_.unlock();
}

void WalkingController::discreteRiccatiEquationInitialize(Eigen::MatrixXd a, Eigen::MatrixXd b)
{
  int n=a.rows(); //number of rows
  int m=b.cols(); //number of columns

  Z11.resize(n, n);
  Z12.resize(n, n);
  Z21.resize(n, n);
  Z22.resize(n, n);
  temp1.resize(m, m);
  temp2.resize(n, n);
  temp3.resize(m, n);

  eigVal_real.resize(2 * n); //eigen value\C0\C7 real\B0\AA
  eigVal_img.resize(2 * n); //eigen value\C0\C7 img\B0\AA
  eigVec_real.resize(2 * n); //eigen vector\C0\C7 real\B0\AA
  eigVec_img.resize(2 * n); //eigen vector\C0\C7 img\B0\AA

  Z.resize(2 * n, 2 * n);

  deigVal_real.resize(2 * n);
  deigVal_img.resize(2 * n);

  deigVec_real.resize(2 * n, 2 * n);
  deigVec_img.resize(2 * n, 2 * n);

  tempZ_real.resize(2 * n, n);
  tempZ_img.resize(2 * n, n);

  U11_inv.resize(n, n);
  X.resize(n, n);
  X_sol.resize(n, n);

  tempZ_comp.resize(n * 2, n);
  U11.resize(n, n);
  U21.resize(n, n);

}

Eigen::MatrixXd WalkingController::discreteRiccatiEquationLQR(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd R, Eigen::MatrixXd Q)
{
  int n=A.rows(); //number of rows
  int n2 = n * 2;
  int m=B.cols(); //number of columns

  for (int i = 0; i<2 * n; i++) //\B0\AA\C0\C7 \C3ʱ\E2ȭ
  {
    eigVec_real[i].resize(n2);
    eigVec_real[i].setZero();
    eigVec_img[i].resize(n2);
    eigVec_img[i].setZero();
  }

  deigVal_real.setZero();
  deigVal_img.setZero();
  deigVec_real.setZero();
  deigVec_img.setZero();
  tempZ_real.setZero();
  tempZ_img.setZero();
  U11_inv.setZero();
  X.setZero();
  X_sol.setZero();


  Z11 = A.inverse();
  temp1 = R.inverse();
  Z21 = Q*Z11;

  Eigen::MatrixXd B_T;
  B_T = B.transpose();
  temp2 = B * R.inverse() * B.transpose();     //B*inv(R)*B'
  Z12 = Z11*temp2;
  Eigen::MatrixXd A_T;
  A_T = A.transpose();
  Z22 = A.transpose() + Z21*temp2;
  Z.setZero();

  //ggory15
  Z.topLeftCorner(n, n) = Z11;
  Z.topRightCorner(n, n) = Z12;
  Z.bottomLeftCorner(n, n) = Z21;
  Z.bottomRightCorner(n, n) = Z22;

  using namespace std;

  Eigen::MatrixXd Z_evr = Z; // \C0ӽ\C3 \C0\FA\C0\E5, rmath\C0\C7 evr\C0\BB \C7ϸ\E9 \BF\F8\BA\BB Z matrix\B0\A1 \BA\AF\C7\FC\B5\CA

  /////////////////////////
  Eigen::EigenSolver<Eigen::MatrixXd> es(Z_evr);    //8ms


  Z_eig = Z.eigenvalues();    //5ms
  es_eig = es.eigenvectors();

  deigVal_real = Z_eig.real();
  deigVal_img = Z_eig.imag();


  for (int i = 0; i<n2; i++)
  {
    for (int ii = 0; ii<n2; ii++)
    {
      deigVec_real(ii, i) = es_eig.col(i)(ii).real();
      deigVec_img(ii, i) = es_eig.col(i)(ii).imag();
    }
  }

  int c1 = 0;

  for (int i = 0; i<n2; i++)
  {
    if ((deigVal_real(i)*deigVal_real(i) + deigVal_img(i)*deigVal_img(i))>1.0) //outside the unit cycle
    {
      for (int j = 0; j<n2; j++)
      {
        tempZ_real(j, c1) = deigVec_real(j, i);
        tempZ_img(j, c1) = deigVec_img(j, i);
      }
      c1++;
    }
  }

  using namespace Eigen;

  for (int i = 0; i<n2; i++)
  {
    for (int j = 0; j<n; j++)
    {
      tempZ_comp.real()(i, j) = tempZ_real(i, j);
      tempZ_comp.imag()(i, j) = tempZ_img(i, j);
    }
  }

  for (int i = 0; i<n; i++)
  {
    for (int j = 0; j<n; j++)
    {
      U11(i, j) = tempZ_comp(i, j);
      U21(i, j) = tempZ_comp(i + n, j);
    }
  }

  U11_inv = U11.inverse();
  X = U21*(U11_inv);

  for (int i = 0; i<n; i++)
  {
    for (int j = 0; j<n; j++)
    {
      X_sol(i, j) = X.real()(i, j);
    }
  }

  return X_sol;
}

Eigen::MatrixXd WalkingController::discreteRiccatiEquationPrev(Eigen::MatrixXd a, Eigen::MatrixXd b, Eigen::MatrixXd r, Eigen::MatrixXd q)
{
  int n=a.rows(); //number of rows
  int m=b.cols(); //number of columns
  
  Eigen::MatrixXd z11(n, n), z12(n, n), z21(n, n), z22(n, n);

  z11 = a.inverse();
  z12 = a.inverse()*b*r.inverse()*b.transpose();
  z21 = q*a.inverse();
  z22 = a.transpose() + q*a.inverse()*b*r.inverse()*b.transpose();

  Eigen::MatrixXd z; z.resize(2*n, 2*n);
  z.setZero();
  z.topLeftCorner(n,n) = z11;
  z.topRightCorner(n,n) = z12;
  z.bottomLeftCorner(n,n) = z21;
  z.bottomRightCorner(n,n) = z22;
    
  std::vector<Eigen::VectorXd> eigVec_real(2*n);
  std::vector<Eigen::VectorXd> eigVec_img(2*n);

  for(int i=0; i<8; i++)
  {
    eigVec_real[i].resize(2*n);
    eigVec_real[i].setZero();
    eigVec_img[i].resize(2*n);
    eigVec_img[i].setZero();
  }

  Eigen::VectorXd deigVal_real(2*n);
  Eigen::VectorXd deigVal_img(2*n);
  deigVal_real.setZero();
  deigVal_img.setZero();
  Eigen::MatrixXd deigVec_real(2*n,2*n);
  Eigen::MatrixXd deigVec_img(2*n,2*n);
  deigVec_real.setZero();
  deigVec_img.setZero();

  deigVal_real = z.eigenvalues().real();
  deigVal_img = z.eigenvalues().imag();
  
  Eigen::EigenSolver<Eigen::MatrixXd> ev(z);

  for(int i=0;i<2*n; i++)
  {
    for(int j=0; j<2*n; j++)
    {
      deigVec_real(j,i) = ev.eigenvectors().col(i)(j).real();
      deigVec_img(j,i) = ev.eigenvectors().col(i)(j).imag();
    }
  }
  
  //Order the eigenvectors
  //move e-vectors correspnding to e-value outside the unite circle to the left

  Eigen::MatrixXd tempZ_real(2*n, n), tempZ_img(2*n, n);
  tempZ_real.setZero();
  tempZ_img.setZero();
  int c=0;

  for (int i=0;i<2*n;i++)
  {
    if ((deigVal_real(i)*deigVal_real(i)+deigVal_img(i)*deigVal_img(i))>1.0) //outside the unit cycle
    {
      for(int j=0; j<2*n; j++)
      {
        tempZ_real(j,c) = deigVec_real(j,i);
        tempZ_img(j,c) = deigVec_img(j,i);
      }
      c++;
    }
  }

  Eigen::MatrixXcd tempZ_comp(2*n, n);
  for(int i=0;i<2*n;i++)
  {
    for(int j=0;j<n;j++)
    {
      tempZ_comp.real()(i,j) = tempZ_real(i,j);
      tempZ_comp.imag()(i,j) = tempZ_img(i,j);
    }
  }

  Eigen::MatrixXcd U11(n, n), U21(n, n), X(n, n);
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      U11(i,j) = tempZ_comp(i,j);
      U21(i,j) = tempZ_comp(i+n,j);
    }
  }
  X = U21*(U11.inverse());

  Eigen::MatrixXd X_sol(n, n);
  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      X_sol(i,j) = X.real()(i,j);
    }
  }

  return X_sol;
}

}
 