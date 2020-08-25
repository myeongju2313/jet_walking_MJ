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
  if((walking_enable_ == true))
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
        circling_motion();
        //supportToFloatPattern();
        computeIkControl_MJ(pelv_trajectory_float_, lfoot_trajectory_float_, rfoot_trajectory_float_, q_des);
          
        for(int i=0; i<12; i++)
        { desired_q_(i) = q_des(i); }
        desired_q_not_compensated_ = desired_q_ ;  
        
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
      desired_q(i) = desired_q_(i);
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
  // 실시간 실제 CoM Update
  q_temp.segment<28>(6) = current_q_.segment<28>(0);   
  qdot_temp.segment<28>(6)= current_qdot_.segment<28>(0);

  if(walking_tick_ > 0) 
  { q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);}

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
    middle_total_step_number = 5; //walking on the spot 10 times
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
    q_temp.segment<12>(6) = desired_q_not_compensated_.segment<12>(0);
     
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

void WalkingController::supportToFloatPattern()
{
  pelv_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*pelv_trajectory_support_;
  lfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*lfoot_trajectory_support_;
  rfoot_trajectory_float_ = DyrosMath::inverseIsometry3d(pelv_trajectory_support_)*rfoot_trajectory_support_;
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
 
void WalkingController::circling_motion()
{
  pelv_trajectory_float_.translation().setZero();
  pelv_trajectory_float_.linear().setIdentity(); 
  
  rfoot_trajectory_float_.translation()(0) = -0.05 * cos(0.5*M_PI*walking_tick_/hz_);
  rfoot_trajectory_float_.translation()(1) = -0.12782;
  rfoot_trajectory_float_.translation()(2) = -0.76548 + 0.05 * sin(0.5*M_PI*walking_tick_/hz_); 
  
  lfoot_trajectory_float_.translation()(0) = -0.05 * cos(0.5*M_PI*(walking_tick_/hz_ - 1.0));
  lfoot_trajectory_float_.translation()(1) =  0.12782;
  lfoot_trajectory_float_.translation()(2) = -0.76548 + 0.05 * sin(0.5*M_PI*(walking_tick_/hz_ - 1.0));    
  
  lfoot_trajectory_float_.linear().setIdentity();
  rfoot_trajectory_float_.linear().setIdentity(); 
}

void WalkingController::updateNextStepTime()
{
  /*if(walking_tick_ == t_last_) // 840 tick 부터 시작해서 240 tick 주기
  { 
    if(current_step_num_ != total_step_num_-1)
    { 
      t_start_ = t_last_ + 1 ; //841 tick
      t_start_real_ = t_start_ + t_rest_init_; // 841 + 10 = 851 tick
      t_last_ = t_start_ + t_total_ -1; // 840 + 240 = 1080 tick  
      current_step_num_ ++; // 시작 단계에서, walking tick이 600, 한 step 걸었을 때, 즉 t_last의 시간과 같게 되면 current_step_num_ 이 증가, 매 step은 240tick, 1.2초가 걸림.
 
    }    
  }
   if(current_step_num_ == total_step_num_-1 && walking_tick_ >= t_last_ + t_total_)// + 160) // 여기에 tick을 늘려주는 만큼 보행 완료 대기 시간을 늘릴 수 있음.
   {
     walking_enable_ = false;
   }*/
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

}
 
