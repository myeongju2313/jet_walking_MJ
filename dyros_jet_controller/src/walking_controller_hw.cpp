#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/dyros_jet_model.h"
#include "dyros_jet_controller/walking_controller_hw.h"
#include "cvxgen_6_8_0/cvxgen/solver.h"
#include <fstream>
#include <tf/tf.h>

Vars vars;
Params params;
Workspace work;
Settings settings;

namespace dyros_jet_controller
{ 
  ofstream MJ_graph("/home/myeongju/MJ_graph.txt");

void WalkingController::compute()
{   
  if(walking_enable_ == true)
  {
    updateInitialState();  
    getRobotState();  
    floatToSupportFootstep(); 

    if(ready_for_thread_flag_ == false)
    { ready_for_thread_flag_ = true; }
    
    if(ready_for_compute_flag_ == true)
    {
      if(current_step_num_< total_step_num_)
      {                  
        getZmpTrajectory();
        getComTrajectory();
        getPelvTrajectory();     
        getFootTrajectory();   
        supportToFloatPattern();
        computeIkControl_MJ(pelv_trajectory_float_, lfoot_trajectory_float_, rfoot_trajectory_float_, q_des);

        Eigen::Vector12d d_q;  
        for(int i=0; i<12; i++)
        { desired_q_(i) = q_des(i); } 
        desired_q_not_compensated_ = desired_q_ ;  
        
        if(hip_compensator_mode_ == true)
        { 
          hip_compensator(); 
          for(int i = 0; i < 12; i++)
          { d_q(i) = desired_q_(i); }
          Compliant_control(d_q); 
        }
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
  is_right_foot_swing_ = is_right_foot_swing;  
  walkingPatternDCM_ = walking_pattern; 

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
    if( mask[i] >= PRIORITY && mask[i] < PRIORITY * 2 )
    {
      if(hip_compensator_mode_ == true)
      {
        if(walking_tick_ == 0)
        { desired_q(i) = desired_q_(i); }
        else
        { desired_q(i) = DOB_IK_output_b_(i); }
      }
      else
      { desired_q(i) = desired_q_(i); }            
    }         
  }  
}

void WalkingController::parameterSetting()
{
  t_double1_ = 0.10*hz_; 
  t_double2_ = 0.10*hz_;
  t_rest_init_ = 0.05*hz_;
  t_rest_last_ = 0.05*hz_;
  t_total_= 1.1*hz_;
  t_temp_ = 3.0*hz_;
  t_last_ = t_total_ + t_temp_ ; //4.2*hz_;
  t_start_ = t_temp_ + 1 ;
 
  t_start_real_ = t_start_ + t_rest_init_;

  current_step_num_ = 0;
  walking_tick_ = 0; 
  foot_height_ = 0.05; 
  gyro_frame_flag_ = false;
  com_control_mode_ = true;
  estimator_flag_ = false; 
}

void WalkingController::getRobotState()
{
  Eigen::Matrix<double, DyrosJetModel::MODEL_WITH_VIRTUAL_DOF, 1> q_temp, qdot_temp;
  q_temp.setZero();
  qdot_temp;  
  q_temp.segment<28>(6) = current_q_.segment<28>(0);   
  qdot_temp.segment<28>(6)= current_qdot_.segment<28>(0);

  //if(walking_tick_ > 0) 
  //{ q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);}

  model_.updateKinematics(q_temp, qdot_temp);
  com_float_current_ = model_.getCurrentCom();
  com_float_current_dot_= model_.getCurrentComDot();
  lfoot_float_current_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)0); 
  rfoot_float_current_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)1);
    
  if(foot_step_(current_step_num_, 6) == 0) 
  { supportfoot_float_current_ = rfoot_float_current_; }
  else if(foot_step_(current_step_num_, 6) == 1)
  { supportfoot_float_current_ = lfoot_float_current_; }

  pelv_support_current_ = DyrosMath::inverseIsometry3d(supportfoot_float_current_);
  pelv_float_current_.setIdentity();
  lfoot_support_current_ = DyrosMath::multiplyIsometry3d(pelv_support_current_,lfoot_float_current_);
  rfoot_support_current_ = DyrosMath::multiplyIsometry3d(pelv_support_current_,rfoot_float_current_);     
  com_support_current_ =  DyrosMath::multiplyIsometry3dVector3d(pelv_support_current_, com_float_current_);
  
  current_motor_q_leg_ = current_q_.segment<12>(0);
}

void WalkingController::calculateFootStepTotal()
{
  /*
   * this function calculate foot steps which the robot should put on
   * algorithm: set robot orientation to the destination -> go straight -> set target orientation on the destination
   *
   * foot_step_(current_step_num_, i) is the foot step where the robot will step right after
   * foot_step_(current_step_num_, 6) = 0 means swingfoot is left(support foot is right)
   */ 
   
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
  double length_to_target = sqrt(target_x_*target_x_ + target_y_*target_y_); // 목표까지의 최종 거리.
  double dlength = step_length_x_; 
  unsigned int middle_total_step_number = length_to_target/dlength; // 목표 최종 거리를 Step 보폭으로 나눈 정수로 표현한 스텝 갯수.
  double middle_residual_length = length_to_target - middle_total_step_number*dlength; // 나머지 거리.
  
  if(length_to_target == 0) // 제자리 걸음
  {
    middle_total_step_number = 10; //walking on the spot 10 times
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
      number_of_foot_step = number_of_foot_step + 2*del_size; // 예를 들어, 45도 회전일 경우, initial_total_step_number = 4, 5번째 발이 목표에 오고, 6번째 발이 5번째 발과 맞춰짐 
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

void WalkingController::floatToSupportFootstep()
{
  Eigen::Isometry3d reference;

  if(current_step_num_ == 0) 
  {
    if(foot_step_(0,6) == 0) //right support
    {
      reference.translation() = rfoot_float_init_.translation(); // (1) -> 정확히 -0.127803은 아님. (실제 오차가 있기 때문에.)  
      reference.translation()(2) = 0.0;
      reference.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(rfoot_float_init_.linear())(2));
      reference.translation()(0) = 0.0;
    }
    else //left support
    {
      reference.translation() = lfoot_float_init_.translation(); // (1) -> 0.127803  
      reference.translation()(2) = 0.0; // Z는 Pelvis 기준으로 발까지 - 변위 였고,
      reference.linear() = DyrosMath::rotateWithZ(DyrosMath::rot2Euler(lfoot_float_init_.linear())(2)); // Pelvis부터 왼쪽 발의 Rotation Matirx를 alpha, beta, gamma로 바꾸고, gamma(Z 축)만 빼서 Rotation matrix로 다시 바꿈.
      reference.translation()(0) = 0.0; // X는 Pevlis 기준으로 Offset 되어 있는 변위 였다.
    } 
  }
  else
  {
    reference.linear() = DyrosMath::rotateWithZ(foot_step_(current_step_num_-1,5)); // 현재 current_step_num에 맞는 계획되어 있는 Global에서 Foot step의 Z방향 회전 각도를 Z회전 Rotation matrix의 변수에 넣는다.
    for(int i=0 ;i<3; i++)
      reference.translation()(i) = foot_step_(current_step_num_-1,i); // 현재 current_step_num에 맞는 계획되어 있는 Global에서 Foot step의 X,Y,Z 좌표를 넣음.
  } 

  Eigen::Vector3d temp_local_position;
  Eigen::Vector3d temp_global_position;
  // temp_local,global_position, 여기서는 Global에서 본 Foot step 좌표들을 Current_step_num_에 따른 현재 지지발 좌표계에서 보기위해 사용됨. 
 
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
    foot_step_support_frame_(i,5) = foot_step_(i,5) - foot_step_(current_step_num_-1,5); // 지지발 기준에서 본 거니까 계속 빼줘야함.      

    if(current_step_num_ == 0)
    { foot_step_support_frame_(i,5) = foot_step_(i,5) - supportfoot_float_init_(5); }
  } 

  for(int j=0;j<3;j++)
    temp_global_position(j) = supportfoot_float_init_(j);

  temp_local_position = reference.linear().transpose()*(temp_global_position - reference.translation());
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
 
    com_float_init_ = model_.getCurrentCom();
    pelv_float_init_.setIdentity();
    lfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(0));
    rfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(1));    
    
    Eigen::Isometry3d ref_frame;

    if(foot_step_(0, 6) == 0)  //right foot support
    { ref_frame = rfoot_float_init_; }    
    else if(foot_step_(0, 6) == 1)
    { ref_frame = lfoot_float_init_; }
    
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
    pelv_support_start_ = pelv_support_init_;
    total_step_num_ = foot_step_.col(1).size();
    xi_ = com_support_init_(0); // preview parameter
    yi_ = com_support_init_(1);
    zc_ = com_support_init_(2);     
    
  }
  else if(current_step_num_ != 0 && walking_tick_ == t_start_) // step change 
  {  
    Eigen::Matrix<double, DyrosJetModel::MODEL_WITH_VIRTUAL_DOF, 1> q_temp, qdot_temp;
    q_temp.setZero();
    qdot_temp;
    q_temp.segment<28>(6) = current_q_.segment<28>(0);
    qdot_temp.segment<28>(6)= current_qdot_.segment<28>(0);  
    //q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);
     
    model_.updateKinematics(q_temp, qdot_temp);

    lfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(0));
    rfoot_float_init_ = model_.getCurrentTransform((DyrosJetModel::EndEffector)(1));
    com_float_init_ = model_.getCurrentCom();
    pelv_float_init_.setIdentity();  

    Eigen::Isometry3d ref_frame;

    if(foot_step_(current_step_num_, 6) == 0)  //right foot support
    { ref_frame = rfoot_float_init_; }
    else if(foot_step_(current_step_num_, 6) == 1)
    { ref_frame = lfoot_float_init_; }

    pelv_support_init_ = DyrosMath::inverseIsometry3d(ref_frame);
    com_support_init_ = pelv_support_init_.linear()*com_float_init_ + pelv_support_init_.translation(); 
    pelv_support_euler_init_ = DyrosMath::rot2Euler(pelv_support_init_.linear()); 

    lfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),lfoot_float_init_);
    rfoot_support_init_ = DyrosMath::multiplyIsometry3d(DyrosMath::inverseIsometry3d(ref_frame),rfoot_float_init_);    
    rfoot_support_euler_init_ = DyrosMath::rot2Euler(rfoot_support_init_.linear());
    lfoot_support_euler_init_ = DyrosMath::rot2Euler(lfoot_support_init_.linear()); 
  }  
} 

void WalkingController::getZmpTrajectory()
{
  unsigned int planning_step_number = 3;
  unsigned int norm_size = 0;
  
  if(current_step_num_ >= total_step_num_ - planning_step_number)
    norm_size = (t_last_ - t_start_ + 1)*(total_step_num_ - current_step_num_) + 20*hz_;
  else
    norm_size = (t_last_ - t_start_ + 1)*(planning_step_number); 
  if(current_step_num_ == 0)
    norm_size = norm_size + t_temp_ + 1;
 
  zmpGenerator(norm_size, planning_step_number);   
}

void WalkingController::zmpGenerator(const unsigned int norm_size, const unsigned planning_step_num)
{ 
  ref_zmp_.resize(norm_size, 2); 
  Eigen::VectorXd temp_px;
  Eigen::VectorXd temp_py;
  
  unsigned int index = 0;
  // 매 tick 마다 zmp가 3발 앞까지 계산 된다. 

  if(current_step_num_ == 0) // Walking을 수행 할 때, 정지 상태 일때 3초 동안 Ref X ZMP를 0으로 보냄. Y ZMP는 제자리 유지.  
  {
    for (int i = 0; i <= t_temp_; i++) //600 tick
    {
      if(i < 0.5*hz_) 
      {
        ref_zmp_(i,0) = com_support_init_(0) ;
        ref_zmp_(i,1) = com_support_init_(1) ;
      }
      else if(i < 1.5*hz_) 
      {
        double del_x = i - 0.5*hz_;
        ref_zmp_(i,0) = com_support_init_(0) - del_x * com_support_init_(0)/(1.0*hz_);
        ref_zmp_(i,1) = com_support_init_(1) ;
      }
      else 
      {
        ref_zmp_(i,0) = 0.0;
        ref_zmp_(i,1) = com_support_init_(1) ;
      }      
      index++;
    }    
  }
  if(current_step_num_ >= total_step_num_ - planning_step_num)
  {  
    for(unsigned int i = current_step_num_; i < total_step_num_; i++)
    {
      onestepZmp(i, temp_px, temp_py);
     
      for(unsigned int j = 0; j < t_total_; j++)
      {
        ref_zmp_(index + j, 0) = temp_px(j);
        ref_zmp_(index + j, 1) = temp_py(j);    
      }
      index = index + t_total_;
    }
    
    for(unsigned int j = 0; j < 20*hz_; j++)
    {
      ref_zmp_(index + j, 0) = ref_zmp_(index -1, 0);
      ref_zmp_(index + j, 1) = ref_zmp_(index -1, 1);
    }
    index = index + 20*hz_;      
  }
  else // 보행 중 사용 하는 Ref ZMP
  { 
    for(unsigned int i = current_step_num_; i < current_step_num_ + planning_step_num; i++)  
    {
      onestepZmp(i, temp_px, temp_py);
      for (unsigned int j = 0; j < t_total_; j++) // 1 step 보행은 1.2초, 240 tick
      {
        ref_zmp_(index+j,0) = temp_px(j);
        ref_zmp_(index+j,1) = temp_py(j);
      }      
      index = index + t_total_; // 참조 zmp가 이만큼 쌓였다.      
      // 결국 실제 로봇 1Hz마다 720개의 ref_zmp를 생성함. 3.6초
    }   
  }   
}

void WalkingController::onestepZmp(unsigned int current_step_number, Eigen::VectorXd& temp_px, Eigen::VectorXd& temp_py)
{
  temp_px.resize(t_total_); // 함수가 실행 될 때 마다, 240 tick의 참조 ZMP를 담는다. Realtime = 1.2초
  temp_py.resize(t_total_);
  temp_px.setZero();
  temp_py.setZero();

  double Kx = 0, Ky = 0, A = 0, B = 0, wn = 0;
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
      else if(i >= t_total_ - t_rest_last_ - t_double2_  && i < t_total_) //1.05 ~ 1.15초 , 210 ~ 230 tick 
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

void WalkingController::preview_Parameter_CPM(double dt, int NL, Eigen::Matrix3d& K, Eigen::Vector3d com_support_init_, Eigen::MatrixXd& Gi, Eigen::VectorXd& Gd, Eigen::MatrixXd& Gx, 
  Eigen::MatrixXd& A, Eigen::VectorXd& B, Eigen::MatrixXd& C, Eigen::MatrixXd& D, Eigen::MatrixXd& A_bar, Eigen::VectorXd& B_bar)
  {      
    double Kp = 100;
    double Kv = 10;
    
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
    D.resize(1,1);
    D(0,0) = -zc_/GRAVITY*Kp; 

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
    R(0,0) = 1.0;

    Eigen::MatrixXd Qx;
    Qx.resize(3,3);
    Qx.setZero();

    Eigen::MatrixXd Q_bar;
    Q_bar.resize(3,3);
    Q_bar.setZero();
    Q_bar(0,0) = Qe(0,0);
        
    K(0,0) = 110.075528194525;
    K(0,1) = 6002.773189475650;
    K(0,2) = 1620.941388698153;
    K(1,1) = 330547.378525671258;
    K(1,2) = 89255.846463209440;
    K(2,2) = 24102.882783488018;
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
    Ac_bar = A_bar - B_bar * Temp_mat_inv * B_bar_tran * K * A_bar;
    
    Eigen::MatrixXd Ac_bar_tran(3,3);
    Ac_bar_tran = Ac_bar.transpose();
 
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

void WalkingController::previewcontroller_CPM(double dt, int NL, int tick, double x_i, double y_i, Eigen::Vector3d xs, Eigen::Vector3d ys, double& UX, double& UY, 
       Eigen::MatrixXd Gi, Eigen::VectorXd Gd, Eigen::MatrixXd Gx, Eigen::MatrixXd A, Eigen::VectorXd B, Eigen::MatrixXd A_bar, Eigen::VectorXd B_bar, Eigen::Vector2d &XD, Eigen::Vector2d &YD, Eigen::VectorXd& X_bar_p, Eigen::VectorXd& Y_bar_p)
{
  int zmp_size;
  zmp_size = ref_zmp_.col(1).size();  
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
    Preview_X_b(0) = x_i;  
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
    Preview_X(1) = X_bar_p(1)/dt; Preview_Y(1) = Y_bar_p(1)/dt;     
  }  
  
  Eigen::VectorXd Temp_mat_X, Temp_mat_Y;
  Temp_mat_X.resize(3); Temp_mat_Y.resize(3);
  Temp_mat_X.setZero(); Temp_mat_Y.setZero();  
    
  Temp_mat_X(0) = Preview_X(0); Temp_mat_Y(0) = Preview_Y(0); 
  Temp_mat_X(2) = X_bar_p(2)/dt; Temp_mat_Y(2) = Y_bar_p(2)/dt;  
  
  px = C*Temp_mat_X;
  py = C*Temp_mat_Y; 
   
  X_bar_p(0) = px(0) - px_ref(tick); 
  Y_bar_p(0) = py(0) - py_ref(tick);
   
  double sum_Gd_px_ref = 0, sum_Gd_py_ref = 0;

  for(int i = 0; i < NL; i++) 
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
  
  XD = A*Preview_X + B*UX;
  YD = A*Preview_Y + B*UY;

         
}

void WalkingController::getComTrajectory()
{
  if(walking_tick_ == 0)  
  { 
    preview_Parameter_CPM(1.0/hz_, 16*hz_/10, K_ ,com_support_init_, Gi_, Gd_, Gx_, A_, B_, C_, D_, A_bar_, B_bar_); 
    UX_ = com_support_init_(0);
    UY_ = com_support_init_(1);
    xs_(0) = xi_; xs_(1) = 0; xs_(2) = 0;
    ys_(0) = yi_; ys_(1) = 0; xs_(2) = 0;
  }

  if(current_step_num_ == 0)
  { zmp_start_time_ = 0.0; }
  else
  { zmp_start_time_ = t_start_; }
          
  previewcontroller_CPM(1.0/hz_, 16*hz_/10, walking_tick_-zmp_start_time_, xi_, yi_, xs_, ys_, UX_, UY_, Gi_, Gd_, Gx_, A_, B_, A_bar_, B_bar_, XD_, YD_, X_bar_p_, Y_bar_p_);
  
  xs_(0) = com_support_current_(0); ys_(0) = com_support_current_(1); 
  xs_(1) = XD_(1); ys_(1) = YD_(1); 

  com_desired_(0) = UX_;
  com_desired_(1) = UY_;
  com_desired_(2) = pelv_support_start_.translation()(2);

  if (walking_tick_ == t_start_ + t_total_-1 && current_step_num_ != total_step_num_-1)  
  { 
    Eigen::Vector3d com_pos_prev;
    Eigen::Vector3d com_pos;
    Eigen::Vector3d com_vel_prev;
    Eigen::Vector3d com_vel;
    Eigen::Vector3d com_acc_prev;
    Eigen::Vector3d com_acc;
    //
    Eigen::Vector3d com_u_prev;
    Eigen::Vector3d com_u;
    // 
    Eigen::Matrix3d temp_rot;
    Eigen::Vector3d temp_pos;
    
    temp_rot = DyrosMath::rotateWithZ(-foot_step_support_frame_(current_step_num_,5)); 
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

    com_u_prev(0) = UX_;
    com_u_prev(1) = UY_;
    com_u = temp_rot*(com_u_prev - temp_pos);

    xs_(0) = com_pos(0);
    ys_(0) = com_pos(1);
    xs_(1) = com_vel(0);
    ys_(1) = com_vel(1);
    xs_(2) = com_acc(0);
    ys_(2) = com_acc(1); 

    UX_ = com_u(0);
    UY_ = com_u(1);  
  } 
}

void WalkingController::getPelvTrajectory()
{
  double z_rot = foot_step_support_frame_(current_step_num_,5);
   
  pelv_trajectory_support_.translation()(0) = pelv_support_current_.translation()(0) + 2.0*(com_desired_(0) - com_support_current_(0));
  pelv_trajectory_support_.translation()(1) = pelv_support_current_.translation()(1) + 2.0*(com_desired_(1) - com_support_current_(1));
  pelv_trajectory_support_.translation()(2) = com_desired_(2);          
       
  Eigen::Vector3d Trunk_trajectory_euler;
  Trunk_trajectory_euler.setZero();

  if(walking_tick_ < t_start_ + t_rest_init_ + t_double1_)
  { Trunk_trajectory_euler(2) = pelv_support_euler_init_(2); }
  else if(walking_tick_ >= t_start_ + t_rest_init_ + t_double1_ && walking_tick_ < t_start_ + t_total_ - t_double2_ - t_rest_last_)
  { Trunk_trajectory_euler(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_double2_-t_rest_last_, pelv_support_euler_init_(2),z_rot/2.0,0.0,0.0); }
  else
  { Trunk_trajectory_euler(2) = z_rot/2.0; } 
  
  pelv_trajectory_support_.linear() = DyrosMath::rotateWithZ(Trunk_trajectory_euler(2))*DyrosMath::rotateWithY(Trunk_trajectory_euler(1))*DyrosMath::rotateWithX(Trunk_trajectory_euler(0));
}

void WalkingController::supportToFloatPattern()
{  
  pelv_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*pelv_trajectory_support_;
  lfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*lfoot_trajectory_support_;
  rfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*rfoot_trajectory_support_;
}

void WalkingController::getFootTrajectory()
{
  Eigen::Vector6d target_swing_foot;

  for(int i=0; i<6; i++)
  { target_swing_foot(i) = foot_step_support_frame_(current_step_num_,i); }
 
  if(walking_tick_ < t_start_ + t_rest_init_ + t_double1_) 
  { 
    lfoot_trajectory_support_.translation() = lfoot_support_init_.translation();  
    rfoot_trajectory_support_.translation() = rfoot_support_init_.translation();
    
    lfoot_trajectory_euler_support_ = lfoot_support_euler_init_;
    rfoot_trajectory_euler_support_ = rfoot_support_euler_init_;

    if(foot_step_(current_step_num_,6) == 1)  
    { 
      lfoot_trajectory_support_.translation().setZero();
      lfoot_trajectory_euler_support_.setZero();

      rfoot_trajectory_support_.translation() = rfoot_support_init_.translation();
      rfoot_trajectory_support_.translation()(2) = 0;
      rfoot_trajectory_euler_support_ = rfoot_support_euler_init_;
    }     
    else if(foot_step_(current_step_num_,6) == 0)  
    { 
      rfoot_trajectory_support_.translation().setZero();
      rfoot_trajectory_euler_support_.setZero();

      lfoot_trajectory_support_.translation() = lfoot_support_init_.translation();
      lfoot_trajectory_support_.translation()(2) = 0;
      lfoot_trajectory_euler_support_ = lfoot_support_euler_init_; 
    }     

    lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
    rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
  }
  
  else if(walking_tick_ >= t_start_ + t_rest_init_ + t_double1_ && walking_tick_ < t_start_ + t_total_ - t_double2_ - t_rest_last_)  
  {   
    if(foot_step_(current_step_num_,6) == 1) 
    {
      lfoot_trajectory_support_.translation() = lfoot_support_init_.translation();             
      lfoot_trajectory_euler_support_.setZero(); 

      lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
      
      if(walking_tick_ < t_start_ + t_rest_init_ + t_double1_ + (t_total_ - t_rest_init_ - t_rest_last_ - t_double1_ - t_double2_)/2.0) 
      { rfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_+ t_rest_init_ + t_double1_,t_start_real_ + t_double1_ + (t_total_ - t_rest_init_ - t_rest_last_ - t_double1_ - t_double2_)/2.0,0,foot_height_,0.0,0.0); }  
      else
      { rfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_)/2.0,t_start_+t_total_-t_rest_last_-t_double2_,foot_height_,target_swing_foot(2),0.0,0.0); }
      
      for(int i=0; i<2; i++) // X, Y 방향 궤적 생성, 여기서 불연속 문제를 회복해야함. 스윙발이 완벽히 따라가지 못하고 지지발이 되었을때, 현재 지지발 기준으로 봤을때 스윙발이 오차가 생긴것이기 때문에 3차곡선으로 이어줘야함. 
      { rfoot_trajectory_support_.translation()(i) = DyrosMath::cubic(walking_tick_,t_start_real_ + t_double1_ , t_start_+t_total_-t_rest_last_-t_double2_, rfoot_support_init_.translation()(i),target_swing_foot(i),0.0,0.0); } 
      
      rfoot_trajectory_euler_support_(0) = 0;
      rfoot_trajectory_euler_support_(1) = 0;
      rfoot_trajectory_euler_support_(2) = DyrosMath::cubic(walking_tick_,t_start_ + t_rest_init_ + t_double1_,t_start_ + t_total_ - t_rest_last_ - t_double2_,rfoot_support_euler_init_(2),target_swing_foot(5),0.0,0.0);
      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
    }
    else if(foot_step_(current_step_num_,6) == 0) // 오른발이 지지발일때, 지지발은 고정, 왼발은 목표 위치로 스윙
    { 
      rfoot_trajectory_support_.translation() = rfoot_support_init_.translation(); 
      rfoot_trajectory_euler_support_.setZero(); 

      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
 
      if(walking_tick_ < t_start_ + t_rest_init_ + t_double1_ + (t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_)/2.0)
      { lfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_)/2.0,0,foot_height_,0.0,0.0); }
      else
      { lfoot_trajectory_support_.translation()(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_+(t_total_-t_rest_init_-t_rest_last_-t_double1_-t_double2_)/2.0,t_start_+t_total_-t_rest_last_-t_double2_,foot_height_,target_swing_foot(2),0.0,0.0); }
         
      for(int i=0; i<2; i++)
      { lfoot_trajectory_support_.translation()(i) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_,lfoot_support_init_.translation()(i),target_swing_foot(i),0.0,0.0); }

      lfoot_trajectory_euler_support_(0) = 0;
      lfoot_trajectory_euler_support_(1) = 0;  
      lfoot_trajectory_euler_support_(2) = DyrosMath::cubic(walking_tick_,t_start_real_+t_double1_,t_start_+t_total_-t_rest_last_-t_double2_,lfoot_support_euler_init_(2),target_swing_foot(5),0.0,0.0);
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
      
      for(int i=0; i<3; i++)
      {
        rfoot_trajectory_support_.translation()(i) = target_swing_foot(i);
        rfoot_trajectory_euler_support_(i) = target_swing_foot(i+3);
      }
      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));
    }
    else if (foot_step_(current_step_num_,6) == 0) 
    {
      //rfoot_trajectory_support_.translation() = rfoot_support_init_.translation();
      rfoot_trajectory_euler_support_.setZero();
      rfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(rfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(rfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(rfoot_trajectory_euler_support_(0));

      for(int i=0; i<3; i++)
      {
        lfoot_trajectory_support_.translation()(i) = target_swing_foot(i);
        lfoot_trajectory_euler_support_(i) = target_swing_foot(i+3);
      }
      lfoot_trajectory_support_.linear() = DyrosMath::rotateWithZ(lfoot_trajectory_euler_support_(2))*DyrosMath::rotateWithY(lfoot_trajectory_euler_support_(1))*DyrosMath::rotateWithX(lfoot_trajectory_euler_support_(0));
    }
  }
}

void WalkingController::computeIkControl_MJ(Eigen::Isometry3d float_trunk_transform, Eigen::Isometry3d float_lleg_transform, Eigen::Isometry3d float_rleg_transform, Eigen::Vector12d& q_des)
{
  //float = World/ trunk = pelvis 
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
  
}

void WalkingController::updateNextStepTime()
{
  if(walking_tick_ == t_last_)  
  { 
    if(current_step_num_ != total_step_num_-1)
    { 
      t_start_ = t_last_ + 1 ;  
      t_start_real_ = t_start_ + t_rest_init_;  
      t_last_ = t_start_ + t_total_ -1;    
      current_step_num_ ++; 
    }    
  }
   if(current_step_num_ == total_step_num_-1 && walking_tick_ >= t_last_ + t_total_) 
   {
     walking_enable_ = false;
   }
   walking_tick_ ++;
}

void WalkingController::slowCalc()
{
  while(true)
  {
    if(ready_for_thread_flag_)
    {
       
      if (ready_for_compute_flag_ == false)
      {
        ready_for_compute_flag_ = true;
      }
    }
    this_thread::sleep_for(chrono::milliseconds(100));
  }
}

void WalkingController::hip_compensator()
{
  double left_hip_angle = 3.0*DEG2RAD, right_hip_angle = 3.0*DEG2RAD, left_hip_angle_first_step = 3.0*DEG2RAD, right_hip_angle_first_step = 3.0*DEG2RAD,
      left_hip_angle_temp = 0.0, right_hip_angle_temp = 0.0, temp_time = 0.1*hz_;

  if(current_step_num_ == 0)
  {
    if(foot_step_(current_step_num_, 6) == 1) //left support foot
    {
      if(walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_rest_init_ + t_double1_, t_start_ + t_rest_init_ + t_double1_ + temp_time, 0, left_hip_angle_first_step, 0.0, 0.0);
      else if(walking_tick_ >= t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time, t_start_ + t_total_ - t_rest_last_,left_hip_angle_first_step, 0.0, 0.0, 0.0);
      else
        left_hip_angle_temp = 0;
    }
    else if(foot_step_(current_step_num_, 6) == 0) // right support foot
    {
      if(walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_rest_init_ + t_double1_, t_start_ + t_rest_init_ + t_double1_ + temp_time, 0, right_hip_angle_first_step, 0.0, 0.0);
      else if(walking_tick_ >= t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time, t_start_ + t_total_ - t_rest_last_,right_hip_angle_first_step, 0.0, 0.0, 0.0);
      else
        right_hip_angle_temp = 0;
    }
  }
  else
  {
    if(foot_step_(current_step_num_, 6) == 1) //left support foot
    {
      if(walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_rest_init_ + t_double1_, t_start_ + t_rest_init_ + t_double1_ + temp_time, 0, left_hip_angle, 0.0, 0.0);
      else if(walking_tick_ >= t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        left_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time, t_start_ + t_total_ - t_rest_last_,left_hip_angle, 0.0, 0.0, 0.0);
      else
        left_hip_angle_temp = 0;
    }
    else if(foot_step_(current_step_num_, 6) == 0) // right support foot
    {
      if(walking_tick_ < t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_rest_init_ + t_double1_, t_start_ + t_rest_init_ + t_double1_ + temp_time, 0, right_hip_angle, 0.0, 0.0);
      else if(walking_tick_ >= t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time)
        right_hip_angle_temp = DyrosMath::cubic(walking_tick_, t_start_ + t_total_ - t_rest_last_ - t_double2_ - temp_time, t_start_ + t_total_ - t_rest_last_,right_hip_angle, 0.0, 0.0, 0.0);
      else
        right_hip_angle_temp = 0;
    }
  }
  desired_q_(1) = desired_q_(1) + left_hip_angle_temp;
  desired_q_(7) = desired_q_(7) - right_hip_angle_temp;
}

void WalkingController::Compliant_control(Eigen::Vector12d desired_leg_q)
{
  if(walking_tick_ == 0)  
  {
    pre_motor_q_leg_ = current_motor_q_leg_; 
    DOB_IK_output_b_ = current_motor_q_leg_;
  }
      
  Eigen::Vector12d current_u; // left right order
  double del_t = 0, Kp = 0;
  del_t = 1/hz_; Kp = 30;

  for (int i = 0; i < 12; i++)
  { // (theta_m(k+1) - theta_m(k)) / dt = Kp (u - theta_m(k)) 
    current_u(i) = (current_motor_q_leg_(i) - (1 - Kp*del_t)*pre_motor_q_leg_(i)) / (Kp*del_t);
  }
  
  Eigen::Vector12d d_hat;    
  d_hat = current_u - DOB_IK_output_b_ ; //current_u -> u' + d , DOB_IK_output_b_ -> u'= IK output + d hat
  
  if(walking_tick_ == 0)
    d_hat_b = d_hat;

  d_hat = 0.7*d_hat_b + 0.3*d_hat; // 필터링

  // Mingon's LQR contorller gain (using external encoder)
  double default_gain = 0.2; // Kp가 정확하다면 시뮬레이션이나 실제 로봇이나 0.2~1의 의미는 같다.
  double compliant_gain = 1.5;
  double compliant_tick = 0.1*hz_;
  double gain_temp;
  for (int i = 0; i < 12; i ++)
  {
    if(i < 6) //왼쪽 다리 관절
    { 
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
      DOB_IK_output_(i) = desired_leg_q(i) + gain_temp*d_hat(i);  // LQR 위치 대신 단순 IK 기반 위치 ( u_c + d_hat = u' (논문))
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
      DOB_IK_output_(i) = desired_leg_q(i) + gain_temp*d_hat(i);  // LQR 위치 대신 단순 IK 기반 위치
    }
  }  
  pre_motor_q_leg_ = current_motor_q_leg_; 
  d_hat_b = d_hat; 
  DOB_IK_output_b_ = DOB_IK_output_; 
}

}
 