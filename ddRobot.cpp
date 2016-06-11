/**
 * Robô de direção diferencial
 * Disciplina de Robótica CIn/UFPE
 * 
 * @autor Prof. Hansenclever Bassani
 * 
 * Este código é proporcionado para facilitar os passos iniciais da programação.
 * Porém, não há garantia de seu correto funcionamento.
 * 
 * Testado em: Ubuntu 14.04 + Netbeans
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <unistd.h>

#define V_REP_IP_ADDRESS "127.0.0.1"
#define V_REP_PORT 19997//1999;
#define PI 3.14
extern "C" {
#include "extApi.h"
    /*	#include "extApiCustom.h" if you wanna use custom remote API functions! */
}

simxInt ddRobotHandle;
simxInt leftMotorHandle;
simxInt rightMotorHandle;
simxInt sensorHandle;
simxInt graphOdometryHandle;
simxFloat pos[3], orientation[3];
int g = 1;

simxFloat goals[8][3] = {{-1.0, 1.5, 3.14/2}, {-0.25, 0.55, 0}, {1, -0.3, -3.14/2}, {1.0, -1.5, -3.14/2}, {1.0, 1.5, 3.14/2}, {0.25, -0.5, -3.14}, {-1.0, -1.5, -3.14/2}, {0.35, 0, 3.14/2}};

simxInt goal= 0;
int state = 0;
simxFloat k_rho = 0.1, k_alpha = 2, k_beta = -0.5;
simxFloat r = 0.02, l = 0.1;

void getPosition(int clientID, simxFloat pos[]) { //[x,y,theta]

    simxInt ret = simxGetObjectPosition(clientID, ddRobotHandle, -1, pos, simx_opmode_oneshot_wait);
    if (ret > 0) {
        printf("Error reading robot position\n");
        return;
    }

    simxFloat orientation[3];
    ret = simxGetObjectOrientation(clientID, ddRobotHandle, -1, orientation, simx_opmode_oneshot_wait);
    if (ret > 0) {
        printf("Error reading robot orientation\n");
        return;
    }

    simxFloat theta = orientation[2];
    pos[2] = theta;
}

simxInt getSimTimeMs(int clientID) { //In Miliseconds
    return simxGetLastCmdTime(clientID);
}

simxFloat to180range(simxFloat angulo){
	angulo = angulo*180.0/PI;
	if(angulo <= -180){
		angulo = angulo + 360.0;
    }if(angulo >180){
         angulo = angulo - 360.0;
    }
    angulo = (PI*angulo)/180.0;
	return angulo;
}

float to_positive_angle(float angle) {

    angle = fmod(angle, 2 * M_PI);
    while (angle < 0) {
        angle = angle + 2 * M_PI;
    }
    return angle;
}

float smallestAngleDiff(float target, float source) {
    float a;
    a = to_positive_angle(target) - to_positive_angle(source);

    if (a > M_PI) {
        a = a - 2 * M_PI;
    } else if (a < -M_PI) {
        a = a + 2 * M_PI;
    }
    return a;
}

void readOdometers(int clientID, simxFloat &dPhiL, simxFloat &dPhiR) {
    //old joint angle position
    static simxFloat lwprev=0; 
    static simxFloat rwprev=0;
    
    //current joint angle position
    simxFloat lwcur=0;
    simxFloat rwcur=0;

    simxGetJointPosition(clientID, leftMotorHandle, &lwcur, simx_opmode_oneshot);
    simxGetJointPosition(clientID, rightMotorHandle, &rwcur, simx_opmode_oneshot);

    dPhiL = smallestAngleDiff(lwcur, lwprev);
    dPhiR = smallestAngleDiff(rwcur, rwprev);
    lwprev = lwcur;
    rwprev = rwcur;
}

void setTargetSpeed(int clientID, simxFloat phiL, simxFloat phiR) {
    simxSetJointTargetVelocity(clientID, leftMotorHandle, phiL, simx_opmode_oneshot);
    simxSetJointTargetVelocity(clientID, rightMotorHandle, phiR, simx_opmode_oneshot);   
}

inline double to_deg(double radians) {
    return radians * (180.0 / M_PI);
}


void nextGoal(){
	goal = goal+1;
    if (goal>7){
        goal =1;
    }
    state = 0;
}

void motionControl(int clientID, simxFloat &phiL, simxFloat &phiR)
{
	
	simxFloat v, omega, wR, wL;
	simxFloat dx, dy, dtheta, atg, rho;
	simxFloat theta, beta, alpha;
	
    //Localizacao atual do robo [x, y, z]
    getPosition(clientID, pos);
	
    //Angulo atual do robo
    if(simxGetObjectOrientation(clientID, ddRobotHandle, -1, orientation, simx_opmode_oneshot_wait) > 0) {
        printf("Error reading robot orientation\n");
        return;
    }
	
    theta = orientation[2];
    
    //Computar dx, dy, rho, alpha e beta
    dx = goals[goal][0] - pos[0];
    dy = goals[goal][1] - pos[1];
    dtheta = smallestAngleDiff(goals[goal][2], theta);
    rho = sqrt(dx*dx + dy*dy);
    
    atg = atan2(dy, dx);
    atg = to180range(atg);
    alpha = smallestAngleDiff(atg,theta);
    alpha = to180range(alpha);
    beta = (goals[goal][2])-theta - alpha;
    beta = to180range(beta);

	v = k_rho*rho;

    if (v<0.2) v = 0.2;
	
    omega = k_alpha*alpha + k_beta*(beta);

    wR = v + l*omega;
    wL = v - l*omega;

    phiR = wR/r;
    phiL = wL/r;

    // o controle eh instavel com rho ~= 0
    if (rho < 0.05 || state > 0){
        state = 1;
        // faz o robo apenas girar para corrigir theta
        if (abs(dtheta)>0.5){
            phiR = 2*dtheta; // proporcionalmente ao que falta
            phiL = -2*dtheta;
		}
        else{
            phiR = 0;
            phiL = 0;
            state = 2;
		}
	}
}
            
int main(int argc, char* argv[]) {

    std::string ipAddr = V_REP_IP_ADDRESS;
    int portNb = V_REP_PORT;

    if (argc > 1) {
        ipAddr = argv[1];
    }

    printf("Iniciando conexao com: %s...\n", ipAddr.c_str());

    int clientID = simxStart((simxChar*) (simxChar*) ipAddr.c_str(), portNb, true, true, 2000, 5);
	
    if (clientID != -1) {
        printf("Conexao efetuada\n");
        
		simxFloat phiL = 0; //rad/s
        simxFloat phiR = 0; //rad/s
		
        //Get handles for robot parts, actuators and sensores:
        simxGetObjectHandle(clientID, "RobotFrame#", &ddRobotHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "LeftMotor#", &leftMotorHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "RightMotor#", &rightMotorHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "GraphOdometry#", &graphOdometryHandle, simx_opmode_oneshot_wait);

        printf("RobotFrame: %d\n", ddRobotHandle);
        printf("LeftMotor: %d\n", leftMotorHandle);
        printf("RightMotor: %d\n", rightMotorHandle);
        printf("GraphOdometry: %d\n", graphOdometryHandle);

        //start simulation
        int ret = simxStartSimulation(clientID, simx_opmode_oneshot_wait);
        
        if (ret==-1) {
            printf("Não foi possível iniciar a simulação.\n");
            return -1;
        }
        
        printf("Simulação iniciada.\n");

        //While is connected:
        while (simxGetConnectionId(clientID) != -1) {
			if(state == 2){
				nextGoal();
			}

			motionControl(clientID, phiL, phiR);
			setTargetSpeed(clientID, phiL, phiR);
			/*
            //Read current position:
            simxFloat pos[3]; //[x,y,theta] in [cm cm rad]
            getPosition(clientID, pos);

            //Read simulation time of the last command:
            simxInt time = getSimTimeMs(clientID); //Simulation time in ms or 0 if sim is not running
            //stop the loop if simulation is has been stopped:
            if (time == 0) break;             
            printf("Posicao: [%.2f %.2f %.2fº], time: %dms\n", pos[0], pos[1], to_deg(pos[2]), time);
            
            //Read current wheels angle variation:
            simxFloat dPhiL, dPhiR; //rad
            readOdometers(clientID, dPhiL, dPhiR);
            printf("dPhiL: %.2f dPhiR: %.2f\n", dPhiL, dPhiR);
            
            //Set new target speeds: robot going in a circle:
            simxFloat phiL = 5; //rad/s
            simxFloat phiR = 20; //rad/s
            setTargetSpeed(clientID, phiL, phiR);
			motionControl(clientID, phiL, phiR);

            //Let some time for V-REP do its work:
            extApi_sleepMs(2);
*/
        }
        
        //Stop the robot and disconnect from V-Rep;
        setTargetSpeed(clientID, 0, 0);
        simxPauseSimulation(clientID, simx_opmode_oneshot_wait);
        simxFinish(clientID);
        
    } else {
        printf("Nao foi possivel conectar.\n");
        return -2;
    }
    
    return 0;
}


